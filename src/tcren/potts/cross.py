r"""Every TCR read against every pMHC: the all-vs-all Cα distance table.

:func:`tcren.potts.geom.fit_geom_prior` sees one geometry per crystal --- 359 of them, covering
whichever ``(loop length, peptide length)`` combinations happen to have been solved. A CDR3 of 16
against an 11-mer may have no template at all, and the estimator's answer is to normalise length
away into graphon coordinates, which is where its placement goes (Jensen--Shannon 0.4200 against
the crystal's own 0.7052).

This module removes the scarcity instead. ``Canonical2026`` is ``Native2026`` after canonical MHC
superposition, so **every complex already sits in one frame** --- which is exactly why
:func:`tcren.shuffle.graft_tcr` can splice a pMHC and a foreign TCR with no coordinate transform at
all. Reading the TCR of crystal *s* against the pMHC of crystal *t* is therefore a legitimate
pseudo-complex, and within MHC class there are :math:`265^2 + 94^2 = 79{,}061` of them against 359
crystals.

**The cross product over structures is the cross product over residues**, so there is no loop over
the 79,061 pairs: group the two residue tables by ``(class, chain, region, loop length, peptide
length)``, and each group is one ``(n_s, n_t, L, M)`` broadcast.

What that buys, and what it does not
------------------------------------

It buys **length coverage** --- a length-16 CDR3 against every peptide length any crystal carries,
which no single crystal set provides --- and **hundreds of samples per cell** where the crystals
give a handful. Those are the two things a conditional :math:`P(d_{C\alpha} \mid \text{chain},
\text{region}, i, L, j, M, \text{class})` needs.

It does **not** buy receptor--peptide dependence beyond the 359 solved pairings. Crossing every TCR
with every pMHC *is* the product measure, so the receptor x peptide cross-covariance it implies is
zero by construction. That block is carried by the ``is_native`` rows and estimated from the
crystals. This is a property of the design, not a defect to work around.
"""
from __future__ import annotations

import time
from collections.abc import Iterable, Iterator
from pathlib import Path

import numpy as np
import polars as pl

from ..structure.model import PEPTIDE_TYPE, Structure
from .geom import _CHAINS, REGIONS_MODELLED

__all__ = ["CROSS_COLUMNS", "MIN_CA", "PAR_COLUMNS", "REC_COLUMNS", "cross_from_structures",
           "cross_pairs", "partner_residues", "receptor_residues", "residue_tables",
           "write_cross_pairs"]

#: One row per resolved CDR Cα. ``len.rec`` is the *resolved* region length --- the ``L`` the
#: geometry model conditions on --- carried on the row so nothing downstream re-derives it.
REC_COLUMNS: tuple[str, ...] = (
    "pdb.id", "mhc_class", "chain.rec", "region.rec", "pos.rec", "len.rec", "aa.rec",
    "ca_x", "ca_y", "ca_z", "cb_x", "cb_y", "cb_z", "has_cb",
)

#: One row per resolved peptide Cα.
PAR_COLUMNS: tuple[str, ...] = (
    "pdb.id", "mhc_class", "pos.par", "par.len", "aa.par", "role.par", "ca_x", "ca_y", "ca_z",
)

#: The cross table. ``pdb.id.rec`` and ``pdb.id.par`` stay **separate columns** so a caller can hold
#: out by structure on either axis; a merged id would make that impossible, and every held-out
#: number this library reports is held out by structure.
CROSS_COLUMNS: tuple[str, ...] = (
    "pdb.id.rec", "pdb.id.par", "mhc_class", "chain.rec", "region.rec", "pos.rec", "len.rec",
    "aa.rec", "pos.par", "par.len", "aa.par", "role.par", "d_ca", "cos_cb", "has_cb", "is_native",
)

#: Closest a receptor Cα may sit to a peptide Cα, Å. **Measured, not assumed**: over the 175,358
#: receptor:peptide Cα pairs of the 359 solved complexes the minimum is **3.366 Å**, with one pair
#: under 3.5 and 23 under 4.0. This sits just below that, so filtering at it drops no row a real
#: complex produced and removes only pairings no chain could adopt --- a foreign TCR dropped onto a
#: groove clashes, and 0.062 % of cross rows sit under 3.0 Å where no native row does at all.
#:
#: The filter is per **row**, not per pairing: a pseudo-complex with one impossible pair still
#: carries hundreds of ordinary ones, and its geometry is worth fitting.
MIN_CA: float = 3.35

#: Columns the cross product is grouped by. Within a group the resample matrices are constant and
#: the coordinates are a dense tensor, which is what makes the whole thing a handful of broadcasts.
_BUCKET: tuple[str, ...] = ("mhc_class", "chain.rec", "region.rec", "len.rec", "par.len")


def _mhc_class(structure: Structure) -> str:
    """``"MHCII"`` / ``"MHCI"`` exactly as :func:`tcren.potts.geom.fit_geom_prior` reads it."""
    return "MHCII" if any(getattr(c, "mhc_class", None) == "MHCII"
                          for c in getattr(structure, "mhc_calls", []) or []) else "MHCI"


def receptor_residues(structure: Structure, regions: Iterable[str] = REGIONS_MODELLED,
                      chains: Iterable[str] = _CHAINS) -> pl.DataFrame:
    """One row per resolved CDR Cα of a chain-typed, MHC-annotated structure.

    Args:
        structure: A structure already through ``classify_chains`` and ``annotate_mhc``. In the
            canonical frame if the rows are destined for :func:`cross_pairs`.
        regions: Receptor regions to emit.
        chains: Receptor chain types to emit.

    Returns:
        A frame with :data:`REC_COLUMNS`, empty if the structure carries none of these regions.
    """
    cls = _mhc_class(structure)
    rows: list[dict] = []
    wanted, want_chains = set(regions), set(chains)
    for chain in structure.chains:
        if getattr(chain, "chain_type", None) not in want_chains:
            continue
        for region in getattr(chain, "regions", []) or []:
            if region.region_type not in wanted:
                continue
            resolved = [r for r in region.residues if r.ca is not None]
            for pos, residue in enumerate(resolved):
                # Glycine has no Cβ. Carry a flag and put the Cα there rather than dropping the row:
                # the residue still has a backbone position and still makes contacts.
                cb = residue.cb
                rows.append({
                    "pdb.id": structure.pdb_id, "mhc_class": cls,
                    "chain.rec": chain.chain_type, "region.rec": region.region_type,
                    "pos.rec": pos, "len.rec": len(resolved), "aa.rec": residue.aa,
                    "ca_x": residue.ca[0], "ca_y": residue.ca[1], "ca_z": residue.ca[2],
                    "cb_x": (cb if cb is not None else residue.ca)[0],
                    "cb_y": (cb if cb is not None else residue.ca)[1],
                    "cb_z": (cb if cb is not None else residue.ca)[2],
                    "has_cb": cb is not None,
                })
    return _frame(rows, REC_COLUMNS)


def partner_residues(structure: Structure) -> pl.DataFrame:
    """One row per resolved peptide Cα, with the anchor call the site table uses.

    ``role.par`` comes from :func:`tcren.refine.anchors.predict_anchors` with the structure passed,
    so the real MHC-class call is used rather than the peptide-length heuristic.

    Args:
        structure: A chain-typed, MHC-annotated structure.

    Returns:
        A frame with :data:`PAR_COLUMNS`, empty if there is no peptide chain.
    """
    from ..refine.anchors import predict_anchors

    chain = next((c for c in structure.chains
                  if getattr(c, "chain_type", None) == PEPTIDE_TYPE), None)
    if chain is None:
        return _frame([], PAR_COLUMNS)
    resolved = [r for r in chain.residues if r.ca is not None]
    if not resolved:
        return _frame([], PAR_COLUMNS)
    peptide = "".join(r.aa for r in resolved)
    anchors = set(predict_anchors(peptide, structure=structure).anchors)
    cls = _mhc_class(structure)
    rows = [{"pdb.id": structure.pdb_id, "mhc_class": cls, "pos.par": pos,
             "par.len": len(resolved), "aa.par": residue.aa,
             "role.par": "anchor" if pos in anchors else "tcr_facing",
             "ca_x": residue.ca[0], "ca_y": residue.ca[1], "ca_z": residue.ca[2]}
            for pos, residue in enumerate(resolved)]
    return _frame(rows, PAR_COLUMNS)


def _frame(rows: list[dict], columns: tuple[str, ...]) -> pl.DataFrame:
    """Build a frame with a stable schema, so an empty structure concatenates with a full one."""
    schema = {c: (pl.Boolean if c == "has_cb"
                  else pl.Int32 if c.startswith(("pos.", "len.", "par.len"))
                  else pl.Float32 if c[:2] in ("ca", "cb")
                  else pl.Utf8) for c in columns}
    return pl.DataFrame(rows, schema=schema) if rows else pl.DataFrame(schema=schema)


def residue_tables(structures: Iterable[Structure]) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Stack :func:`receptor_residues` and :func:`partner_residues` over a structure set.

    A structure contributing no CDR loop or no peptide is dropped from **both** tables, so the two
    always describe the same set of complexes.

    Args:
        structures: Chain-typed, MHC-annotated structures --- pass
            :func:`tcren.annotation.batch.iter_annotated_set` through
            :func:`tcren.mhc.regions.annotate_mhc_batch`, never annotating in a loop.

    Returns:
        ``(receptor, partner)``.

    Raises:
        ValueError: if no structure yielded both a loop and a peptide.
    """
    rec, par = [], []
    for structure in structures:
        r, p = receptor_residues(structure), partner_residues(structure)
        if r.is_empty() or p.is_empty():
            continue
        rec.append(r)
        par.append(p)
    if not rec:
        raise ValueError("no structure yielded both a CDR loop and a peptide — were the chains "
                         "typed (classify_chains) and the MHC annotated (annotate_mhc_batch) "
                         "before this was called?")
    return pl.concat(rec), pl.concat(par)


def _bucket_frames(rec: pl.DataFrame, par: pl.DataFrame, radius: float, include_self: bool,
                   min_dist: float) -> Iterator[pl.DataFrame]:
    """Yield one frame per ``(class, chain, region, L, M)`` bucket. One broadcast each."""
    # Buckets are walked in sorted key order and each is emitted with its rows already sorted, so
    # the concatenation is globally sorted without a second pass -- which is what lets
    # `write_cross_pairs` stream parts that still read back in the declared order.
    rec_groups = dict(sorted(
        (k, g.sort("pdb.id", "pos.rec"))
        for k, g in rec.group_by(["mhc_class", "chain.rec", "region.rec", "len.rec"])))
    par_groups = dict(sorted(
        (k, g.sort("pdb.id", "pos.par")) for k, g in par.group_by(["mhc_class", "par.len"])))
    for (cls, chain, region, loop_len), r in rec_groups.items():
        s_ids = r["pdb.id"].unique(maintain_order=True).to_list()
        n_s = len(s_ids)
        xyz = r.select("ca_x", "ca_y", "ca_z").to_numpy().reshape(n_s, loop_len, 3)
        cb = r.select("cb_x", "cb_y", "cb_z").to_numpy().reshape(n_s, loop_len, 3)
        for (par_cls, pep_len), p in par_groups.items():
            if par_cls != cls:                       # never cross an MHC class
                continue
            t_ids = p["pdb.id"].unique(maintain_order=True).to_list()
            n_t = len(t_ids)
            pep = p.select("ca_x", "ca_y", "ca_z").to_numpy().reshape(n_t, pep_len, 3)
            # (n_s, n_t, L, M): every receptor residue against every partner residue, for every
            # pairing of a TCR with a pMHC. This is the whole engine.
            delta = xyz[:, None, :, None, :] - pep[None, :, None, :, :]
            d_ca = np.linalg.norm(delta, axis=-1)
            # "Will the residue turn to the peptide": cos of the angle between the side-chain
            # direction Ca->Cb and the direction to the partner Ca. +1 points at it, -1 away.
            side = (cb - xyz)[:, None, :, None, :]
            side_norm = np.linalg.norm(side, axis=-1)
            with np.errstate(invalid="ignore", divide="ignore"):
                cos_cb = np.einsum("stlmk,stlmk->stlm", side, -delta) / (side_norm * d_ca)
            cos_cb = np.where(side_norm > 1e-6, cos_cb, 0.0)
            keep = (d_ca <= radius) & (d_ca >= min_dist)
            if not include_self:
                same = np.equal.outer(np.asarray(s_ids), np.asarray(t_ids))
                keep &= ~same[:, :, None, None]
            if not keep.any():
                continue
            si, ti, li, mi = np.nonzero(keep)
            rec_i, par_i = si * loop_len + li, ti * pep_len + mi
            yield pl.DataFrame({
                "pdb.id.rec": r["pdb.id"].to_numpy()[rec_i],
                "pdb.id.par": p["pdb.id"].to_numpy()[par_i],
                "mhc_class": np.repeat(cls, len(si)),
                "chain.rec": np.repeat(chain, len(si)),
                "region.rec": np.repeat(region, len(si)),
                "pos.rec": li.astype(np.int32),
                "len.rec": np.full(len(si), loop_len, dtype=np.int32),
                "aa.rec": r["aa.rec"].to_numpy()[rec_i],
                "pos.par": mi.astype(np.int32),
                "par.len": np.full(len(si), pep_len, dtype=np.int32),
                "aa.par": p["aa.par"].to_numpy()[par_i],
                "role.par": p["role.par"].to_numpy()[par_i],
                "d_ca": d_ca[keep].astype(np.float32),
                "cos_cb": cos_cb[keep].astype(np.float32),
                "has_cb": r["has_cb"].to_numpy()[rec_i],
            }).with_columns((pl.col("pdb.id.rec") == pl.col("pdb.id.par")).alias("is_native"))


def cross_pairs(rec: pl.DataFrame, par: pl.DataFrame, *, radius: float = float("inf"),
                include_self: bool = True, min_dist: float = MIN_CA) -> pl.DataFrame:
    """Every within-class ``(TCR, pMHC)`` pairing's Cα distances, in memory.

    Args:
        rec: Receptor residues, :data:`REC_COLUMNS`.
        par: Partner residues, :data:`PAR_COLUMNS`.
        radius: Keep pairs at or under this Cα distance, Å. ``inf`` keeps complete ``(L, M)``
            blocks, which is what a geometry fit needs; a finite radius gives the smaller analysis
            table and **cannot** be fitted from.
        include_self: Keep the ``is_native`` diagonal --- the 359 real complexes, which are the only
            rows carrying a receptor:peptide pairing that was actually solved.
        min_dist: Drop pairs closer than this, Å. See :data:`MIN_CA`; pass 0 to keep the clashes.

    Returns:
        A frame with :data:`CROSS_COLUMNS`, sorted by :data:`_BUCKET` then
        ``pdb.id.rec, pdb.id.par, pos.rec, pos.par`` --- so a bucket is a contiguous, dense
        ``(pairs, L, M)`` tensor that a fit reshapes rather than pivots.

    Raises:
        ValueError: if no pairing survived.
    """
    frames = list(_bucket_frames(rec, par, radius, include_self, min_dist))
    if not frames:
        raise ValueError(f"no receptor:peptide pair in [{min_dist}, {radius}] Å — are both tables "
                         "non-empty and do they share an MHC class?")
    return pl.concat(frames).select(CROSS_COLUMNS)


def write_cross_pairs(rec: pl.DataFrame, par: pl.DataFrame, out: str | Path, *,
                      radius: float = float("inf"), include_self: bool = True,
                      min_dist: float = MIN_CA, part_rows: int = 4_000_000) -> tuple[int, int]:
    """Write :func:`cross_pairs` as a partitioned parquet dataset.

    ``out`` is a **directory** of ``part-0000.parquet`` … , read back as one lazy frame with
    ``pl.scan_parquet(out / "*.parquet")``. One rule, with no branch on whether the table happens to
    fit in memory: the complete class-I product is about 31 M rows.

    Args:
        rec: Receptor residues.
        par: Partner residues.
        out: Output directory; created if absent, and existing ``part-*.parquet`` are replaced.
        radius: See :func:`cross_pairs`.
        include_self: See :func:`cross_pairs`.
        min_dist: See :func:`cross_pairs`.
        part_rows: Rows per part.

    Returns:
        ``(rows, parts)``.
    """
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    for stale in out.glob("part-*.parquet"):
        stale.unlink()
    rows = parts = 0
    buffer: list[pl.DataFrame] = []
    held = 0

    def flush() -> None:
        nonlocal parts, buffer, held
        if held:
            pl.concat(buffer).select(CROSS_COLUMNS).write_parquet(out / f"part-{parts:04d}.parquet")
            parts += 1
            buffer, held = [], 0

    for frame in _bucket_frames(rec, par, radius, include_self, min_dist):
        buffer.append(frame)
        held += len(frame)
        rows += len(frame)
        if held >= part_rows:
            flush()
    flush()
    if not rows:
        raise ValueError(f"no receptor:peptide pair in [{min_dist}, {radius}] Å — nothing written")
    return rows, parts


def cross_from_structures(structures: Iterable[Structure], out: str | Path, *,
                          radius: float = float("inf"), include_self: bool = True,
                          min_dist: float = MIN_CA, part_rows: int = 4_000_000,
                          verbose: bool = False) -> dict:
    """:func:`residue_tables` then :func:`write_cross_pairs`, with a summary of what was written.

    Args:
        structures: Chain-typed, MHC-annotated structures in one canonical frame.
        out: Output directory.
        radius: See :func:`cross_pairs`.
        include_self: See :func:`cross_pairs`.
        min_dist: See :func:`cross_pairs`.
        part_rows: Rows per part.
        verbose: Print timings as it goes.

    Returns:
        ``{"structures", "pairs", "rows", "parts", "bytes", "seconds", "by_class"}``, where
        ``pairs`` is the number of pseudo-complexes --- :math:`\\sum_c n_c^2` over MHC classes when
        ``include_self``, and :math:`\\sum_c n_c(n_c - 1)` without.
    """
    start = time.time()
    rec, par = residue_tables(structures)
    n_by_class = dict(rec.group_by("mhc_class").agg(pl.col("pdb.id").n_unique().alias("n"))
                      .iter_rows())
    pairs = sum(n * (n if include_self else n - 1) for n in n_by_class.values())
    if verbose:
        print(f"{len(rec):,} receptor residues, {len(par):,} peptide residues, "
              f"{pairs:,} pseudo-complexes ({time.time() - start:.0f} s)", flush=True)
    rows, parts = write_cross_pairs(rec, par, out, radius=radius, include_self=include_self,
                                    min_dist=min_dist, part_rows=part_rows)
    size = sum(f.stat().st_size for f in Path(out).glob("part-*.parquet"))
    return {"structures": rec["pdb.id"].n_unique(), "pairs": pairs, "rows": rows, "parts": parts,
            "bytes": size, "seconds": time.time() - start, "by_class": n_by_class}
