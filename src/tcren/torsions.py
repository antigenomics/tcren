"""Backbone torsion angles (φ, ψ, ω) per residue, and the CDR3 loops in particular.

The contact map says *which* residues touch; the torsions say whether the backbone that puts them
there is one a real protein would adopt. That distinction is what makes torsions useful against
generated structures: a predictor can seat a side chain in a plausible contact while placing its
backbone in a region of the Ramachandran map that crystals essentially never visit.

Angles follow the IUPAC convention and are returned in **degrees** in ``(-180, 180]``:

* ``φ`` = C(i−1) – N(i) – Cα(i) – C(i)
* ``ψ`` = N(i) – Cα(i) – C(i) – N(i+1)
* ``ω`` = Cα(i−1) – C(i−1) – N(i) – Cα(i)   (≈ 180° trans, ≈ 0° cis-proline)

φ is undefined for the first residue of a chain and ψ for the last, so those come back as ``nan``
rather than being silently dropped — a caller counting residues must see the gap.

Example:
    >>> import tcren
    >>> s = tcren.parse_structure("1ao7.pdb.gz", pdb_id="1ao7")
    >>> tcren.annotation.classify_chains(s)
    >>> df = cdr3_torsions(s)
    >>> sorted(df.columns)[:3]
    ['aa', 'chain.id', 'chain.type']
"""

from __future__ import annotations

import numpy as np
import polars as pl

from .structure.model import Chain, Residue, Structure

__all__ = ["dihedral", "residue_torsions", "chain_torsions", "cdr3_torsions"]

#: Chain types that carry a CDR3 loop, mapped to the region markup's loop name.
_TCR_LOOP = {"TRA": "cdr3a", "TRB": "cdr3b", "TRD": "cdr3d", "TRG": "cdr3g"}


def dihedral(p0, p1, p2, p3) -> float:
    """Signed dihedral about the ``p1``–``p2`` axis, in degrees in ``(-180, 180]``.

    Returns ``nan`` if any coordinate is missing, so a residue with an unresolved backbone atom
    produces a gap rather than a fabricated angle.
    """
    if p0 is None or p1 is None or p2 is None or p3 is None:
        return float("nan")
    b0 = np.asarray(p0, float) - np.asarray(p1, float)
    b1 = np.asarray(p2, float) - np.asarray(p1, float)
    b2 = np.asarray(p3, float) - np.asarray(p2, float)
    n = np.linalg.norm(b1)
    if n < 1e-9:
        return float("nan")
    b1 = b1 / n
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    return float(np.degrees(np.arctan2(np.dot(np.cross(b1, v), w), np.dot(v, w))))


def _atom(res: Residue | None, name: str):
    if res is None:
        return None
    for a in res.atoms:
        if a.name == name:
            return a.coord
    return None


def residue_torsions(prev: Residue | None, cur: Residue,
                     nxt: Residue | None) -> tuple[float, float, float]:
    """``(phi, psi, omega)`` in degrees for one residue given its sequence neighbours."""
    phi = dihedral(_atom(prev, "C"), _atom(cur, "N"), _atom(cur, "CA"), _atom(cur, "C"))
    psi = dihedral(_atom(cur, "N"), _atom(cur, "CA"), _atom(cur, "C"), _atom(nxt, "N"))
    omega = dihedral(_atom(prev, "CA"), _atom(prev, "C"), _atom(cur, "N"), _atom(cur, "CA"))
    return phi, psi, omega


def chain_torsions(chain: Chain, structure_id: str = "") -> pl.DataFrame:
    """Per-residue ``phi``/``psi``/``omega`` for one chain.

    Residues are taken in ``seq_index`` order, so a residue whose neighbour is unresolved gets a
    ``nan`` angle: the peptide bond across a chain break does not exist and must not be invented.
    """
    res = sorted(chain.residues, key=lambda r: r.seq_index)
    rows = []
    for k, cur in enumerate(res):
        prev = res[k - 1] if k > 0 and res[k - 1].seq_index == cur.seq_index - 1 else None
        nxt = (res[k + 1] if k + 1 < len(res)
               and res[k + 1].seq_index == cur.seq_index + 1 else None)
        phi, psi, omega = residue_torsions(prev, cur, nxt)
        rows.append({"pdb.id": structure_id, "chain.id": chain.chain_id,
                     "chain.type": chain.chain_type, "seq_index": cur.seq_index,
                     "aa": cur.aa, "phi": phi, "psi": psi, "omega": omega})
    return pl.DataFrame(rows, schema={
        "pdb.id": pl.String, "chain.id": pl.String, "chain.type": pl.String,
        "seq_index": pl.Int64, "aa": pl.String, "phi": pl.Float64, "psi": pl.Float64,
        "omega": pl.Float64})


def cdr3_torsions(structure: Structure, region: str = "CDR3",
                  drop_incomplete: bool = True) -> pl.DataFrame:
    """Torsions of the CDR3 residues of every TCR chain in ``structure``.

    The structure must already be chain-typed (:func:`tcren.annotation.classify_chains`), since
    the CDR3 span comes from the region markup.

    Args:
        structure: a chain-typed structure.
        region: region markup type to select (``"CDR3"``; ``"CDR1"``/``"CDR2"`` also work).
        drop_incomplete: drop residues whose φ or ψ is ``nan`` (chain termini, breaks). Set
            ``False`` to see the gaps.

    Returns:
        One row per selected residue with ``pdb.id, chain.id, chain.type, loop, seq_index, aa,
        phi, psi, omega``. Empty (with the right schema) if the structure has no typed TCR chain.
    """
    frames = []
    for chain in structure.chains:
        loop = _TCR_LOOP.get(chain.chain_type or "")
        if loop is None:
            continue
        spans = [r for r in chain.regions if r.region_type == region]
        if not spans:
            continue
        keep = {i for r in spans for i in range(r.start_seq_index, r.end_seq_index + 1)}
        t = chain_torsions(chain, structure.pdb_id)
        frames.append(t.filter(pl.col("seq_index").is_in(list(keep)))
                       .with_columns(pl.lit(loop).alias("loop")))
    if not frames:
        return pl.DataFrame(schema={
            "pdb.id": pl.String, "chain.id": pl.String, "chain.type": pl.String,
            "seq_index": pl.Int64, "aa": pl.String, "phi": pl.Float64, "psi": pl.Float64,
            "omega": pl.Float64, "loop": pl.String})
    out = pl.concat(frames)
    if drop_incomplete:
        out = out.filter(pl.col("phi").is_not_nan() & pl.col("psi").is_not_nan())
    return out
