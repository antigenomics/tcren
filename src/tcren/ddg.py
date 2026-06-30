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

#: Interfaces whose contacts include the peptide. A peptide point mutation can only
#: change the energy of an interface that the peptide is part of; for any other
#: interface (e.g. ``"tcr_mhc"``) every per-position ΔΔG is exactly 0.
_PEPTIDE_INTERFACES: frozenset[str] = frozenset({"tcr_peptide", "peptide_mhc"})


def _score_one(
    contact_map: ContactMap,
    peptide: str,
    potential: Potential,
    interface: Interface,
    tcr_regions: str,
) -> float:
    """Score a single peptide and return its scalar energy."""
    res = score_peptides(
        contact_map, [peptide], potential, interface=interface, tcr_regions=tcr_regions
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

    Returns:
        ``E(native) - E(mutant)``; positive means the mutation is destabilising.
        Always ``0.0`` for interfaces that do not contain the peptide (e.g.
        ``"tcr_mhc"``), since a peptide mutation cannot affect them.
    """
    if interface not in _PEPTIDE_INTERFACES:
        return 0.0
    e_native = _score_one(contact_map, native, potential, interface, tcr_regions)
    e_mutant = _score_one(contact_map, mutant, potential, interface, tcr_regions)
    return e_native - e_mutant


def alanine_scan(
    contact_map: ContactMap,
    native: str,
    potential: Potential,
    *,
    interface: Interface = "tcr_peptide",
    tcr_regions: str = "all",
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

    Returns:
        Columns ``pos`` (0-based), ``wt_aa`` (native residue at that position) and
        ``ddG`` (``E(native) - E(Ala@pos)``). Positions without TCR contacts yield
        ``ddG == 0.0``. For interfaces that do not contain the peptide (e.g.
        ``"tcr_mhc"``) every position yields ``ddG == 0.0``.
    """
    peptide_iface = interface in _PEPTIDE_INTERFACES
    e_native = (
        _score_one(contact_map, native, potential, interface, tcr_regions)
        if peptide_iface
        else 0.0
    )
    rows = []
    for pos, wt in enumerate(native):
        if peptide_iface:
            mutant = native[:pos] + "A" + native[pos + 1 :]
            e_mut = _score_one(contact_map, mutant, potential, interface, tcr_regions)
            ddg_val = e_native - e_mut
        else:
            ddg_val = 0.0
        rows.append({"pos": pos, "wt_aa": wt, "ddG": ddg_val})
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
