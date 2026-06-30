"""Fast ΔΔG of peptide point mutations (virtual-matrix path).

Implements the paper's fast ΔΔG: no atoms move and no re-docking is performed.
A mutation's effect is read straight off the substitution potential by re-scoring
the mutant sequence on the *same* contact map. ``ddg = E(native) - E(mutant)`` so a
positive value flags a destabilising mutation (the mutant binds less favourably,
i.e. has higher energy).
"""

from __future__ import annotations

from collections.abc import Iterable

import polars as pl

from .contactmap import ContactMap, Interface
from .potential import Potential
from .scoring import score_peptides


def _score_one(
    contact_map: ContactMap,
    peptide: str,
    potential: Potential,
    interface: Interface,
    tcr_regions: str,
    contact_weight: str = "residue",
) -> float:
    """Score a single peptide and return its scalar energy."""
    res = score_peptides(
        contact_map, [peptide], potential, interface=interface,
        tcr_regions=tcr_regions, contact_weight=contact_weight,
    )
    if res.height == 0:
        raise ValueError(
            f"peptide {peptide!r} was not scored "
            "(length mismatch with the structure's peptide?)"
        )
    return float(res["score"][0])


def ddg(
    contact_map: ContactMap,
    native: str,
    mutant: str,
    potential: Potential,
    *,
    interface: Interface = "tcr_peptide",
    tcr_regions: str = "all",
    contact_weight: str = "residue",
) -> float:
    """ΔΔG of a peptide mutation as ``E(native) - E(mutant)``.

    Args:
        contact_map: The structure's contact map.
        native: Native peptide sequence.
        mutant: Mutant peptide sequence (same length as ``native``).
        potential: Pairwise potential to score with.
        interface: Which interface to score over (default ``"tcr_peptide"``).
        tcr_regions: Which TCR regions to keep on the TCR side (passed through to
            ``score_peptides``).
        contact_weight: ``"residue"`` (default) or ``"atomic"``; passed through to
            ``score_peptides``.

    Returns:
        ``E(native) - E(mutant)``; positive means the mutation is destabilising.
    """
    e_native = _score_one(
        contact_map, native, potential, interface, tcr_regions, contact_weight
    )
    e_mutant = _score_one(
        contact_map, mutant, potential, interface, tcr_regions, contact_weight
    )
    return e_native - e_mutant


def alanine_scan(
    contact_map: ContactMap,
    native: str,
    potential: Potential,
    *,
    interface: Interface = "tcr_peptide",
    tcr_regions: str = "all",
    contact_weight: str = "residue",
) -> pl.DataFrame:
    """Alanine scan of the native peptide.

    Mutates each position of ``native`` to alanine in turn and reports the ΔΔG of
    that single substitution. One row per peptide position.

    Args:
        contact_map: The structure's contact map.
        native: Native peptide sequence.
        potential: Pairwise potential to score with.
        interface: Which interface to score over (default ``"tcr_peptide"``).
        tcr_regions: Which TCR regions to keep on the TCR side.
        contact_weight: ``"residue"`` (default) or ``"atomic"``; passed through to
            ``score_peptides``.

    Returns:
        Columns ``pos`` (0-based), ``wt_aa`` (native residue at that position) and
        ``ddG`` (``E(native) - E(Ala@pos)``). Positions without TCR contacts yield
        ``ddG == 0.0``.
    """
    e_native = _score_one(
        contact_map, native, potential, interface, tcr_regions, contact_weight
    )
    rows = []
    for pos, wt in enumerate(native):
        mutant = native[:pos] + "A" + native[pos + 1 :]
        e_mut = _score_one(
            contact_map, mutant, potential, interface, tcr_regions, contact_weight
        )
        rows.append({"pos": pos, "wt_aa": wt, "ddG": e_native - e_mut})
    return pl.DataFrame(rows, schema={"pos": pl.Int64, "wt_aa": pl.Utf8, "ddG": pl.Float64})


def neoantigen_ddg(
    contact_map: ContactMap,
    native: str,
    mutants: Iterable[str],
    potential: Potential,
    **kw,
) -> pl.DataFrame:
    """ΔΔG of candidate neoantigen mutants relative to a native peptide.

    Args:
        contact_map: The structure's contact map.
        native: Native peptide sequence.
        mutants: Candidate mutant peptides (each the same length as ``native``).
        potential: Pairwise potential to score with.
        **kw: Forwarded to :func:`ddg` (``interface``, ``tcr_regions``).

    Returns:
        Columns ``native``, ``mutant`` and ``ddG`` (``E(native) - E(mutant)``;
        positive means the mutation is destabilising), one row per mutant.
    """
    rows = [
        {"native": native, "mutant": m, "ddG": ddg(contact_map, native, m, potential, **kw)}
        for m in mutants
    ]
    return pl.DataFrame(
        rows, schema={"native": pl.Utf8, "mutant": pl.Utf8, "ddG": pl.Float64}
    )
