"""Footprint shape: how a receptor's contacts are *distributed*, not what they score.

Every other scorer in ``tcren`` reads the interface as a sum over contacts. The same contact map
also has a **shape** — which of the six CDR loops touched what, and whether the resulting footprint
is one connected patch — and that shape is a different observable. It carries no potential, no
fitted parameter and no reference structure.

Two families, both computed from one contact map:

**Coverage.** Partition the TCR:pMHC residue contacts into cells and measure how evenly they are
spread. With ``p_i`` the fraction of contacts in cell ``i`` over ``k`` cells,

.. math::
    H = -\\frac{1}{\\ln k}\\sum_i p_i \\ln p_i, \\qquad D_q = \\Big(\\sum_i p_i^q\\Big)^{1/(1-q)}

``H`` is the normalised Shannon entropy (1.0 = perfectly even) and ``D_q`` the Hill number of order
``q`` — the *effective number of engaged cells* (Hill 1973, `doi:10.2307/1934352
<https://doi.org/10.2307/1934352>`_; Jost 2006, `doi:10.1111/j.2006.0030-1299.14714.x
<https://doi.org/10.1111/j.2006.0030-1299.14714.x>`_). ``D_1 = exp H_raw`` is a monotone transform of ``H`` and
ranks identically; ``D_2 = 1/\\sum_i p_i^2`` discounts weakly populated cells and separates better.
Two partitions ship: the **12 cells** of the 6 CDR loops × {peptide, MHC}, and the **24 cells** that
additionally split the peptide into N-terminal, central and C-terminal bands. Refining the peptide
side helps; refining the MHC side into its helices does not, which is why it is not offered.

**Topology.** Join the contacted pMHC residues at a Cα threshold and build the flag (clique)
complex on them. ``b0`` counts disconnected footprint patches and ``b1`` its holes. The cyclomatic
number of the *bipartite contact graph* (``E - V + C``) is deliberately **not** the headline: with
of order thirty contacts among of order thirty residues it is dominated by ``E`` and simply tracks
interface size. The patch count is scale-free and is not redundant with the coverage entropy.

Everything here is invariant under rigid motion, so **no canonical orientation is required** — only
chain typing and CDR region markup (:func:`tcren.annotation.classify_chains`). MHC *region* markup
is not needed either, so the two-pass MHC annotation trap does not apply.

    >>> from tcren.footprint import footprint_features
    >>> row = footprint_features(structure)      # doctest: +SKIP
    >>> row["D2_pep24"], row["fp_b0_r7"]         # doctest: +SKIP
"""

from __future__ import annotations

import warnings
from collections.abc import Iterable, Sequence
from pathlib import Path

import numpy as np
import polars as pl

from .contactmap import ContactMap
from .structure.model import PEPTIDE_TYPE, Structure

__all__ = [
    "CELL_LOOPS",
    "FOOTPRINT_FEATURES",
    "FOOTPRINT_SIZE_FEATURES",
    "footprint_topology_features",
    "cell_counts",
    "footprint_batch",
    "footprint_features",
]

#: The six CDR loops, in the order the cell partition indexes them.
CELL_LOOPS: tuple[str, ...] = ("TRA:CDR1", "TRA:CDR2", "TRA:CDR3",
                               "TRB:CDR1", "TRB:CDR2", "TRB:CDR3")

#: Peptide band edges: residues 0-2 are N-terminal, 3-5 central, 6+ C-terminal (0-based).
_BANDS = (3, 6)

#: Haldane-Anscombe correction. A 9-mer leaves cells genuinely empty, so the canonical-preference
#: log-odds needs it to stay finite; the diversity measures use the raw counts.
_PSEUDO = 0.5

#: The three raw contact counts this module emits alongside the shape measures. They are interface
#: **size**, not shape, and are catalogued under ``interface`` in :data:`tcren.recognition.DESCRIPTORS`
#: for exactly that reason: a shape channel that carried the contact count would correlate with the
#: interface channel by construction, and the whole point of the coverage and topology measures is
#: that they are size-free.
#:
#: The total is ``n_loop_contacts``, not ``n_contacts``: through 2.19.0 it was written under the
#: latter name, which :mod:`tcren.potts` also emits for a *different* quantity — the available pairs
#: that engaged, 29 against this module's 66 on 1ao7. Whichever pass ran last won the column, so a
#: feature table built without ``potts`` handed the footprint tally to a read-out standardized on
#: the Potts population. The two now have two names.
FOOTPRINT_SIZE_FEATURES: tuple[str, ...] = ("n_loop_contacts", "n_pep_contacts", "n_mhc_contacts")

#: Every column :func:`footprint_features` guarantees, size columns included. The radius-tagged
#: Betti columns (``fp_b0_r7`` and friends) are named from the ``radii`` argument and so are not
#: listed here; :func:`footprint_topology_features` gives the shape-only subset.
FOOTPRINT_FEATURES: tuple[str, ...] = (
    *FOOTPRINT_SIZE_FEATURES,
    "H_cell", "D1_cell", "D2_cell", "S_cell", "J_cell",
    "H_loop", "D2_loop", "D2_pep24",
    "ab_imb", "ab_imb_pep", "ab_imb_mhc",
    "L_canon", "p_germ_mhc", "p_cdr3_pep",
    "pep_free_frac",
    "pep_cov_frac", "pep_cov_even", "pep_cov_d2n", "pep_cov_centre", "pep_cov_spread",
    "h0_pers_ent",
)

#: The peptide-coverage block, normalized by peptide length so class I and class II compare.
PEPTIDE_COVERAGE_FEATURES: tuple[str, ...] = (
    "pep_free_frac",
    "pep_cov_frac", "pep_cov_even", "pep_cov_d2n", "pep_cov_centre", "pep_cov_spread",
)


def footprint_topology_features(radii: Sequence[float] = (7.0, 8.0)) -> tuple[str, ...]:
    """The shape-only feature names: :data:`FOOTPRINT_FEATURES` without the size counts, plus the
    radius-tagged Betti columns that ``radii`` produces.
    """
    tags = [f"r{r:g}" for r in radii]
    return tuple([f for f in FOOTPRINT_FEATURES if f not in FOOTPRINT_SIZE_FEATURES]
                 + [f"fp_{k}_{t}" for t in tags for k in ("b0", "b1", "chi", "b0_frac")])


# --- the cell partition ------------------------------------------------------------------------

def cell_counts(structure: Structure, cutoff: float = 5.0) -> pl.DataFrame:
    """Long ``(loop, target, band, n)`` tally of TCR:pMHC residue contacts.

    One :class:`~tcren.contactmap.ContactMap` build covers both TCR interfaces. ``loop`` is
    ``"<chain>:<region>"`` restricted to :data:`CELL_LOOPS`; ``target`` is ``"pep"`` or ``"mhc"``;
    ``band`` is the peptide third (``"pN"``/``"pM"``/``"pC"``) or ``"mhc"``. Counting is a single
    polars ``group_by`` — no Python loop over contacts.

    Args:
        structure: a chain-typed, CDR-region-annotated TCR-pMHC structure.
        cutoff: heavy-atom contact threshold in Angstrom.

    Returns:
        A frame with columns ``loop``, ``target``, ``band``, ``n``. Empty if the structure makes no
        CDR-loop contact with the pMHC.
    """
    # `classify_chains` leaves an MHC chain typed generically as "MHC"; `interface("tcr_mhc")`
    # matches only the supertype `annotate_mhc` assigns. Without that pass the MHC half of the
    # partition is unreachable and every measure here is computed on peptide contacts alone -- a
    # silent, plausible-looking wrong answer, so say it out loud.
    if any(getattr(c, "chain_type", None) == "MHC" for c in structure.chains):
        warnings.warn(f"{structure.pdb_id}: MHC chains are not annotated (chain_type == 'MHC'); "
                      f"run tcren.mhc.annotate_mhc first or the 6 MHC cells stay empty",
                      RuntimeWarning, stacklevel=2)
    cm = ContactMap.from_structure(structure, cutoff=cutoff)
    frames = []
    for target, iface in (("pep", "tcr_peptide"), ("mhc", "tcr_mhc")):
        d = cm.interface(iface)
        if d.is_empty():
            continue
        # `pos.to` is the peptide position when the map carries it, but it is null on a structure
        # whose peptide chain has no markup; the residue's own chain index is the same 0-based
        # count and is always present, so coalesce rather than silently collapsing every peptide
        # contact into one band.
        pos = pl.coalesce(pl.col("pos.to"), pl.col("residue.index.to"))
        band = (pl.lit("mhc") if target == "mhc" else
                pl.when(pos < _BANDS[0]).then(pl.lit("pN"))
                .when(pos < _BANDS[1]).then(pl.lit("pM")).otherwise(pl.lit("pC")))
        frames.append(d.select(
            pl.concat_str([pl.col("chain.type.from"), pl.col("region.type.from")],
                          separator=":").alias("loop"),
            pl.lit(target).alias("target"), band.alias("band")))
    if not frames:
        return pl.DataFrame(schema={"loop": pl.String, "target": pl.String,
                                    "band": pl.String, "n": pl.Int64})
    return (pl.concat(frames, how="vertical")
            .filter(pl.col("loop").is_in(list(CELL_LOOPS)))
            .group_by("loop", "target", "band", maintain_order=True)
            .agg(pl.len().cast(pl.Int64).alias("n"))
            .sort("loop", "target", "band"))


def _diversity(n: np.ndarray, k: int) -> dict[str, float]:
    """Normalised Shannon entropy, Hill numbers of order 1 and 2, richness and evenness."""
    tot = float(n.sum())
    if tot <= 0:
        return dict.fromkeys(("H", "D1", "D2", "S", "J"), float("nan"))
    p = n[n > 0] / tot
    h = float(-(p * np.log(p)).sum())
    s = float(len(p))
    return {"H": h / np.log(k), "D1": float(np.exp(h)), "D2": float(1.0 / (p ** 2).sum()),
            "S": s, "J": h / np.log(s) if s > 1 else float("nan")}


def _peptide_coverage(structure: Structure, cutoff: float) -> dict[str, float]:
    """How the TCR's contacts spread over the peptide, normalized by peptide length.

    ``n_pep_contacted`` is a raw count and ``D2_pep24`` splits the peptide into three fixed bands,
    so both are length-confounded: three contacted positions mean something different on a class I
    9-mer and a class II 15-mer, and a fixed third of a 9-mer is not a fixed third of a 15-mer.
    Every column here divides by the peptide's own length ``L``, so a class I 8-mer, a class I
    13-mer and a class II 18-mer are on one scale.

    Let ``c_i`` be the TCR contacts made to peptide position ``i`` and ``p_i = c_i / sum(c)``.

    ``pep_cov_frac``
        Contacted positions over ``L`` -- the share of the peptide the receptor touches at all.
    ``pep_free_frac``
        The share of the peptide the groove leaves for the receptor, ``L_eff / L``, where
        ``L_eff = sum_i a_i`` and ``a_i = n_TCR_i / (n_TCR_i + n_MHC_i)`` is position ``i``'s
        **accessibility**: its share of contacts that face the receptor rather than the groove.
        This is the threshold-free reading of "peptide without its MHC anchors" -- a P2 or
        C-terminal anchor buried in the groove has ``a_i`` near 0 and a solvent-exposed bulge near
        1, with no position index, no band and no cutoff anywhere in the definition. A binary
        anchor test cannot do this job: at a 5 A heavy-atom criterion **every** residue of a class
        I nonamer contacts the MHC, so the non-anchor set comes out empty.
    ``pep_cov_even``
        Pielou evenness ``H / ln L`` of ``q``, the TCR contact distribution **discounted by
        accessibility**: ``q_i`` is proportional to ``n_TCR_i * a_i``, so an anchor contributes
        almost nothing however many of its atoms sit within the cutoff of a CDR loop. The base is
        the peptide's own length, not the occupied-position count -- a receptor spread evenly over
        three of nine positions must score below one spread evenly over all nine, and normalising
        by the occupied count would make them equal. In ``[0, 1]`` for every peptide length.
    ``pep_cov_d2n``
        Hill number of order 2 of ``q`` over ``L`` -- the effective *share* of the peptide engaged,
        discounting positions carrying a single contact.
    ``pep_cov_centre``
        Contact-weighted mean position on ``[0, 1]`` from N- to C-terminus: where along the
        peptide the receptor sits. ``0.5`` is centred, which is the canonical diagonal docking.
    ``pep_cov_spread``
        Contact-weighted standard deviation of that position, doubled so a receptor reaching both
        termini approaches 1.

    Returns ``nan`` for every column on a peptide of one residue or with no TCR contact.
    """
    nan = dict.fromkeys(PEPTIDE_COVERAGE_FEATURES, float("nan"))
    pep = next((c for c in structure.chains if c.chain_type == PEPTIDE_TYPE), None)
    if pep is None or len(pep.residues) < 2:
        return nan
    length = len(pep.residues)

    cm = ContactMap.from_structure(structure, cutoff=cutoff)

    def per_position(frame: pl.DataFrame, side: str) -> np.ndarray:
        """Contacts per peptide position, on whichever side of the frame the peptide sits."""
        out = np.zeros(length, dtype=float)
        if frame.height == 0:
            return out
        # `pos.*` is the peptide position when the map carries it and null on a structure whose
        # peptide chain has no markup; the residue's own chain index is the same 0-based count and
        # is always present, so coalesce rather than dropping every contact -- the same fallback
        # `cell_counts` takes.
        tally = (frame.select(pl.coalesce(pl.col(f"pos.{side}"),
                                          pl.col(f"residue.index.{side}")).alias("pos"))
                 .drop_nulls("pos").group_by("pos").len())
        for pos, n in zip(tally["pos"].to_list(), tally["len"].to_list()):
            if pos is not None and 0 <= pos < length:
                out[int(pos)] = float(n)
        return out

    # The peptide is the `to` side against the TCR and the `from` side against the MHC.
    n_tcr = per_position(cm.interface("tcr_peptide"), "to")
    n_mhc = per_position(cm.interface("peptide_mhc"), "from")
    total = n_tcr.sum()
    if total <= 0:
        return nan

    both = n_tcr + n_mhc
    access = np.divide(n_tcr, both, out=np.zeros(length), where=both > 0)
    length_eff = float(access.sum())

    # The receptor-facing contact distribution: each position's TCR contacts discounted by how
    # much of it the groove holds. An anchor contributes almost nothing however many atoms of it
    # sit within 5 A of a CDR loop.
    weighted = n_tcr * access
    axis = np.arange(length, dtype=float) / (length - 1)
    if weighted.sum() <= 0:                       # every contacted position is fully buried
        return nan
    q = weighted / weighted.sum()
    nz = q[q > 0]
    h = float(-(nz * np.log(nz)).sum())
    centre = float((q * axis).sum())
    return {
        "pep_free_frac": length_eff / length,
        "pep_cov_frac": float((n_tcr > 0).sum()) / length,
        "pep_cov_even": h / float(np.log(length)),
        "pep_cov_d2n": float(1.0 / (q ** 2).sum()) / length,
        "pep_cov_centre": centre,
        "pep_cov_spread": 2.0 * float(np.sqrt((q * (axis - centre) ** 2).sum())),
    }


def _imbalance(a: float, b: float) -> float:
    """Signed contact imbalance in [-1, 1]; positive is alpha-shifted. NaN on an empty interface."""
    return (a - b) / (a + b) if (a + b) > 0 else float("nan")


# --- topology ----------------------------------------------------------------------------------

def _gf2_rank(M: np.ndarray) -> int:
    """Rank over GF(2) by Gaussian elimination. M is (triangles x edges) and small."""
    A = np.ascontiguousarray(M, dtype=np.uint8)
    r = 0
    for c in range(A.shape[1]):
        nz = np.nonzero(A[r:, c])[0]
        if not len(nz):
            continue
        i = r + int(nz[0])
        if i != r:
            A[[r, i]] = A[[i, r]]
        hit = np.nonzero(A[:, c])[0]
        A[hit[hit != r]] ^= A[r]
        r += 1
        if r == A.shape[0]:
            break
    return r


def _flag_betti(X: np.ndarray, radius: float) -> tuple[float, float]:
    """``(b0, b1)`` of the flag complex on points ``X`` with edges at distance <= ``radius``."""
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components

    n = len(X)
    if n < 2:
        return float(n), 0.0
    D = np.linalg.norm(X[:, None, :] - X[None], axis=-1)
    A = (D <= radius) & ~np.eye(n, dtype=bool)
    ei, ej = np.nonzero(np.triu(A))
    if not len(ei):
        return float(n), 0.0
    comps = connected_components(
        coo_matrix((np.ones(len(ei)), (ei, ej)), shape=(n, n)), directed=False)[0]
    b1 = len(ei) - n + comps
    eidx = {(int(a), int(b)): k for k, (a, b) in enumerate(zip(ei, ej))}
    tri = [(i, j, k) for i in range(n) for j in range(i + 1, n) if A[i, j]
           for k in range(j + 1, n) if A[i, k] and A[j, k]]
    if tri:
        B = np.zeros((len(tri), len(ei)), np.uint8)
        for t, (i, j, k) in enumerate(tri):
            B[t, [eidx[(i, j)], eidx[(i, k)], eidx[(j, k)]]] = 1
        b1 -= _gf2_rank(B)
    return float(comps), float(max(b1, 0))


def _h0_persistence_entropy(X: np.ndarray) -> float:
    """Persistence entropy of H0. The barcode's death times ARE the MST edge lengths, so this
    needs no filtration library: build the minimum spanning tree and take the normalised entropy
    of its edge lengths. Scale-free by construction."""
    from scipy.sparse.csgraph import minimum_spanning_tree

    if len(X) < 3:
        return float("nan")
    D = np.linalg.norm(X[:, None, :] - X[None], axis=-1)
    lengths = minimum_spanning_tree(D).toarray()
    lengths = lengths[lengths > 0]
    if len(lengths) < 2:
        return float("nan")
    p = lengths / lengths.sum()
    return float(-(p * np.log(p)).sum() / np.log(len(p)))


def _pmhc_footprint_ca(structure: Structure, cutoff: float) -> np.ndarray:
    """Calpha coordinates of the pMHC residues the receptor contacts."""
    cm = ContactMap.from_structure(structure, cutoff=cutoff)
    keys: set[tuple[str, int]] = set()
    for iface in ("tcr_peptide", "tcr_mhc"):
        d = cm.interface(iface)
        if not d.is_empty():
            keys |= set(zip(d["chain.id.to"].to_list(), d["residue.index.to"].to_list()))
    ca = {(c.chain_id, r.seq_index): r.ca
          for c in structure.chains for r in c.residues if r.ca is not None}
    pts = [ca[k] for k in sorted(keys) if k in ca]
    return np.asarray(pts, float) if pts else np.empty((0, 3))


# --- the feature row ----------------------------------------------------------------------------

def footprint_features(structure: Structure, *, cutoff: float = 5.0,
                       radii: Sequence[float] = (7.0, 8.0)) -> dict[str, float]:
    """Every coverage and topology feature of one structure, as a flat row.

    Args:
        structure: a chain-typed, CDR-region-annotated TCR-pMHC structure. No canonical
            orientation is needed — every feature is invariant under rigid motion.
        cutoff: heavy-atom contact threshold in Angstrom.
        radii: Calpha thresholds at which the footprint's flag complex is built. The patch count
            ``b0`` is most informative at 7 A and the hole count ``b1`` at 8 A, so both ship.

    Returns:
        ``{feature: value}`` over :data:`FOOTPRINT_FEATURES` plus ``fp_b0_r<r>``, ``fp_b1_r<r>``,
        ``fp_chi_r<r>`` and ``fp_b0_frac_r<r>`` for each radius. Values are ``nan`` where the
        structure gives them no support (no contacts, a single contacted residue).

    Note:
        ``n_loop_contacts`` and its two components count the contacts **the partition sees** —
        those made by the six CDR loops. Framework contacts are outside :data:`CELL_LOOPS` and are
        excluded by construction, so this is smaller than the full interface contact count, and it
        is a different quantity again from :mod:`tcren.potts`'s ``n_contacts``, which counts the
        available pairs that engaged rather than the residue pairs in reach of a loop. The
        topology features are not restricted this way: they are built on every contacted pMHC
        residue, framework-driven ones included, because the footprint is a region on the pMHC
        and does not care which part of the receptor produced it.
    """
    t = cell_counts(structure, cutoff=cutoff)
    row: dict[str, float] = {}
    if t.is_empty():
        row.update(dict.fromkeys(FOOTPRINT_FEATURES, float("nan")))
    else:
        # sorted, not just grouped: an unordered group_by makes `-sum p log p` add its terms in a
        # different order each run, which moved every entropy column by ~1e-15 between runs
        by = lambda *k: t.group_by(list(k)).agg(pl.col("n").sum()).sort(list(k))  # noqa: E731
        cell12 = by("loop", "target")
        cell24 = by("loop", "band")
        loop6 = by("loop")
        n_pep = int(t.filter(pl.col("target") == "pep")["n"].sum())
        n_mhc = int(t.filter(pl.col("target") == "mhc")["n"].sum())

        for tag, frame, k in (("cell", cell12, 12), ("loop", loop6, 6)):
            for stat, v in _diversity(frame["n"].to_numpy(), k).items():
                row[f"{stat}_{tag}"] = v
        row["D2_pep24"] = _diversity(cell24["n"].to_numpy(), 24)["D2"]
        row.update(_peptide_coverage(structure, cutoff))

        chain = pl.col("loop").str.slice(0, 3)
        side = lambda expr: {  # noqa: E731
            c: float(t.filter(expr & (chain == c))["n"].sum()) for c in ("TRA", "TRB")}
        allc, pepc, mhcc = (side(pl.lit(True)), side(pl.col("target") == "pep"),
                            side(pl.col("target") == "mhc"))
        # CDR-loop contacts ONLY, like every other measure here -- the cell partition is the
        # definition. A tally over *all* TCR residues gives a different number (framework contacts
        # are ~11% of the interface on the VDJdb benchmark), and the two must not be compared.
        row["ab_imb"] = _imbalance(allc["TRA"], allc["TRB"])
        row["ab_imb_pep"] = _imbalance(pepc["TRA"], pepc["TRB"])
        row["ab_imb_mhc"] = _imbalance(mhcc["TRA"], mhcc["TRB"])

        germ = pl.col("loop").str.contains("CDR1|CDR2")
        cnt = lambda g, tg: float(t.filter((germ if g else ~germ)  # noqa: E731
                                           & (pl.col("target") == tg))["n"].sum()) + _PSEUDO
        gm, gp, cm_, cp = cnt(1, "mhc"), cnt(1, "pep"), cnt(0, "mhc"), cnt(0, "pep")
        row["L_canon"] = float(np.log((gm * cp) / (gp * cm_)))
        row["p_germ_mhc"] = gm / (gm + gp)
        row["p_cdr3_pep"] = cp / (cp + cm_)
        row["n_loop_contacts"] = float(n_pep + n_mhc)
        row["n_pep_contacts"] = float(n_pep)
        row["n_mhc_contacts"] = float(n_mhc)

    X = _pmhc_footprint_ca(structure, cutoff)
    row["h0_pers_ent"] = _h0_persistence_entropy(X)
    for r in radii:
        tag = f"r{int(r)}" if float(r).is_integer() else f"r{r}"
        b0, b1 = _flag_betti(X, float(r)) if len(X) else (float("nan"), float("nan"))
        row[f"fp_b0_{tag}"] = b0
        row[f"fp_b1_{tag}"] = b1
        row[f"fp_chi_{tag}"] = b0 - b1
        row[f"fp_b0_frac_{tag}"] = b0 / len(X) if len(X) else float("nan")
    return row


def footprint_batch(structures: str | Path | Iterable[Structure], *, cutoff: float = 5.0,
                    radii: Sequence[float] = (7.0, 8.0),
                    organism: str = "human") -> pl.DataFrame:
    """One row per structure, over a folder / glob / archive or an iterable of structures.

    A path is resolved through :func:`tcren.paper.helpers.iter_annotated_set`, which sends every
    chain of every structure to arda in **one mmseqs call per organism**. Nothing here annotates
    per structure and nothing here uses a process pool: mmseqs is the parallel layer.

    The MHC pass then runs **after** chain typing, in one batched call. It is not optional and its
    order is not free: ``classify_chains`` leaves an MHC chain typed generically as ``"MHC"``, and
    ``ContactMap.interface("tcr_mhc")`` matches on the supertype that :func:`tcren.mhc.annotate_mhc`
    assigns. Skip it and every TCR:MHC contact vanishes without an error -- six of the twelve cells
    empty, ``p_germ_mhc`` collapses from ~0.78 to ~0.06, and ``H_cell`` is computed over a partition
    half of which is structurally unreachable.

    Args:
        structures: a directory, glob, ``.tar.gz`` or manifest of structures, or an iterable of
            already chain-typed :class:`~tcren.structure.model.Structure` objects.
        cutoff: heavy-atom contact threshold in Angstrom.
        radii: Calpha thresholds for the footprint flag complex.
        organism: organism for the single-structure path; ignored when a set is batched.

    Returns:
        A frame with ``pdb.id`` plus every feature of :func:`footprint_features`.
    """
    if isinstance(structures, (str, Path)):
        import os

        from .cli import _iter_typed
        from .mhc import annotate_mhc_batch

        structs = list(_iter_typed(Path(structures), organism=organism))
        annotate_mhc_batch(structs, threads=os.cpu_count() or 1)     # must follow chain typing
        it: Iterable[Structure] = structs
    else:
        it = structures
    rows = [{"pdb.id": s.pdb_id, **footprint_features(s, cutoff=cutoff, radii=radii)} for s in it]
    return pl.DataFrame(rows) if rows else pl.DataFrame(schema={"pdb.id": pl.String})



def _selfcheck() -> None:  # pragma: no cover - exercised by tests/unit/test_footprint.py
    """Assert the invariants that make these numbers meaningful."""
    n = np.array([10, 10, 10, 10], float)
    d = _diversity(n, 4)
    assert abs(d["H"] - 1.0) < 1e-12, d          # uniform over k cells is exactly 1
    assert abs(d["D2"] - 4.0) < 1e-9, d          # ... and engages exactly 4 effective cells
    skew = _diversity(np.array([37, 1, 1, 1], float), 4)
    assert skew["H"] < d["H"] and skew["D2"] < d["D2"], (skew, d)
    assert abs(_imbalance(3, 1) - 0.5) < 1e-12
    assert np.isnan(_imbalance(0, 0))

    # a square of 4 points at side 2: one patch, one hole at a radius that spans the side but
    # not the diagonal (2 < r < 2*sqrt(2)); at a radius spanning the diagonal the hole fills in
    sq = np.array([[0, 0, 0], [2, 0, 0], [2, 2, 0], [0, 2, 0]], float)
    assert _flag_betti(sq, 2.5) == (1.0, 1.0), _flag_betti(sq, 2.5)
    assert _flag_betti(sq, 3.0) == (1.0, 0.0), _flag_betti(sq, 3.0)
    far = np.array([[0, 0, 0], [1, 0, 0], [50, 0, 0], [51, 0, 0]], float)
    assert _flag_betti(far, 2.0)[0] == 2.0       # two patches, which is what a ragged footprint is
    assert _gf2_rank(np.array([[1, 1, 0], [0, 1, 1], [1, 0, 1]], np.uint8)) == 2  # sums to 0 mod 2


if __name__ == "__main__":  # pragma: no cover
    _selfcheck()
    print("footprint: self-check passed")
