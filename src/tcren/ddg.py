"""ΔΔG of peptide point mutations, virtually or on rebuilt coordinates.

``ddg = E(native) - E(mutant)``, and lower energy is a more favourable interface throughout
``tcren``, so a **positive** value flags a **stabilising** mutation: the mutant scores below the
native and binds better. A negative value is the destabilising one.

Two ways to get ``E(mutant)``, and they are not the same measurement.

**Virtual** (no ``structure`` argument) is the paper's fast path: no atoms move, and the mutant
sequence is re-indexed against the potential over the *native* contact map. It is exact for the
energy bookkeeping and wrong about geometry — a contact that exists only because a long arginine
reaches across is still counted after that arginine is notionally an alanine, whose Cβ stops 4 Å
short. On the 374 reference crystals only 54 % of 5 Å TCR:peptide residue pairs have both side
chains in range at all, so this is not a rare corner.

**Structural** (pass ``structure=``) rebuilds the mutant's coordinates with
:func:`tcren.refine.substitute.substitute_peptide`, recomputes its contact map, and scores that.
For an **alanine** target it is exact and needs no relaxation, because alanine's heavy atoms are
exactly backbone + Cβ: truncating at Cβ *is* the alanine, and a position mutated from glycine gets
an ideal-geometry Cβ built. Contacts the wild-type side chain alone was reaching then disappear, as
they physically must. For any other target the substituted residue is left as a Cβ stub, so its
reach is under-stated -- see :func:`tcren.refine.substitute.substitute_peptide` for what would
have to be built.
"""

from __future__ import annotations

from collections.abc import Iterable

import numpy as np
import polars as pl

from .contactmap import ContactMap, Interface
from .potential import Potential
from .scoring import score_peptides
from .structure.model import Structure

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
    contact_weight: str = "residue",
    weights: "np.ndarray | None" = None,
) -> float:
    """Score a single peptide and return its scalar energy."""
    res = score_peptides(
        contact_map, [peptide], potential, interface=interface,
        tcr_regions=tcr_regions, contact_weight=contact_weight, weights=weights,
    )
    if res.height == 0:
        raise ValueError(
            f"peptide {peptide!r} was not scored "
            "(length mismatch with the structure's peptide?)"
        )
    return float(res["score"][0])


def _mutant_map(structure: Structure, mutant: str, cutoff: float, sidechain: bool) -> ContactMap:
    """The mutant's own contact map, from rebuilt coordinates."""
    from .refine.substitute import substitute_peptide

    return ContactMap.from_structure(
        substitute_peptide(structure, mutant), cutoff=cutoff, sidechain=sidechain
    )


def ddg(
    contact_map: ContactMap,
    native: str,
    mutant: str,
    potential: Potential,
    *,
    interface: Interface = "tcr_peptide",
    tcr_regions: str = "all",
    contact_weight: str = "residue",
    structure: Structure | None = None,
    cutoff: float = 5.0,
    sidechain: bool = False,
    weights: np.ndarray | None = None,
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
        structure: When given, the mutant is **built** — its side chains are replaced and its
            contact map recomputed — rather than re-indexed on the native map. Exact for an
            alanine target; a Cβ stub for anything longer (see the module docstring).
        cutoff: Contact distance threshold for the rebuilt map (Å). Ignored when ``structure``
            is ``None``, in which case ``contact_map``'s own cutoff applies.
        sidechain: Passed to the rebuilt contact map, so a caller filtering on side-chain
            participation filters the mutant by the *mutant's* reach and not the native's.
        weights: An explicit per-contact multiplier, one value per row of the selected interface
            and in its row order, forwarded to :func:`tcren.scoring.score_peptides`. Its use here
            is to replace the map's hard 0/1 contact indicator with a contact **probability** --
            :func:`tcren.potts.contact_probabilities`' ``p_model``, or a rotamer-averaged
            occupancy -- so a substitution is scored against how often each pair actually touches
            rather than against one frozen snapshot of whether it did. ``None`` (default) leaves
            the result byte-identical. Ignored on a rebuilt mutant map (``structure=`` given),
            whose rows are its own and no longer align with the native's.

    Returns:
        ``E(native) - E(mutant)``; positive means the mutant has the LOWER energy, i.e. the
        mutation is stabilising. Negative means destabilising.
        Always ``0.0`` for interfaces that do not contain the peptide (e.g.
        ``"tcr_mhc"``), since a peptide mutation cannot affect them.
    """
    if interface not in _PEPTIDE_INTERFACES:
        return 0.0
    e_native = _score_one(
        contact_map, native, potential, interface, tcr_regions, contact_weight, weights
    )
    mutant_map = (contact_map if structure is None
                  else _mutant_map(structure, mutant, cutoff, sidechain))
    # A rebuilt mutant has its own contact rows, so a weight vector indexed on the native map
    # would be silently misaligned; drop it rather than mis-apply it.
    e_mutant = _score_one(
        mutant_map, mutant, potential, interface, tcr_regions, contact_weight,
        weights if structure is None else None
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
    structure: Structure | None = None,
    cutoff: float = 5.0,
    sidechain: bool = False,
) -> pl.DataFrame:
    """Alanine scan of the native peptide.

    Mutates each position of ``native`` to alanine in turn and reports the ΔΔG of
    that single substitution. One row per peptide position.

    With ``structure`` the mutant at each position is **built** and its contact map recomputed, so
    a position whose side chain was the only thing reaching the TCR loses those contacts, as it
    physically must. This is the case the structural path gets exactly right (see the module
    docstring), and it costs one contact-map rebuild per position.

    Args:
        contact_map: The structure's contact map.
        native: Native peptide sequence.
        potential: Pairwise potential to score with.
        interface: Which interface to score over (default ``"tcr_peptide"``).
        tcr_regions: Which TCR regions to keep on the TCR side.
        contact_weight: ``"residue"`` (default) or ``"atomic"``; passed through to
            ``score_peptides``.
        structure: Build each alanine mutant and rescore it on its own contact map.
        cutoff: Contact threshold for the rebuilt maps (Å).
        sidechain: Passed to the rebuilt maps.

    Returns:
        Columns ``pos`` (0-based), ``wt_aa`` (native residue at that position) and
        ``ddG`` (``E(native) - E(Ala@pos)``). Positions without TCR contacts yield
        ``ddG == 0.0``. For interfaces that do not contain the peptide (e.g.
        ``"tcr_mhc"``) every position yields ``ddG == 0.0``.
    """
    peptide_iface = interface in _PEPTIDE_INTERFACES
    e_native = (
        _score_one(contact_map, native, potential, interface, tcr_regions, contact_weight)
        if peptide_iface
        else 0.0
    )
    rows = []
    for pos, wt in enumerate(native):
        if peptide_iface:
            mutant = native[:pos] + "A" + native[pos + 1 :]
            mutant_map = (contact_map if structure is None
                          else _mutant_map(structure, mutant, cutoff, sidechain))
            e_mut = _score_one(
                mutant_map, mutant, potential, interface, tcr_regions, contact_weight
            )
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
        Columns ``native``, ``mutant`` and ``ddG`` (``E(native) - E(mutant)``; positive means
        the mutant has the lower energy, i.e. the substitution is stabilising), one row per
        mutant -- so ranking candidates by descending ``ddG`` puts the best first.
    """
    rows = [
        {"native": native, "mutant": m, "ddG": ddg(contact_map, native, m, potential, **kw)}
        for m in mutants
    ]
    return pl.DataFrame(
        rows, schema={"native": pl.Utf8, "mutant": pl.Utf8, "ddG": pl.Float64}
    )


def reference_delta(
    contact_map: ContactMap,
    peptide: str,
    potential: Potential,
    *,
    interface: Interface = "tcr_peptide",
    reference_aa: str = "A",
    tcr_regions: str = "all",
    contact_weight: str = "residue",
    structure: Structure | None = None,
    cutoff: float = 5.0,
    sidechain: bool = False,
) -> float:
    """Poly-alanine reference difference ΔΦ = Φ(peptide) − Φ(reference) on this contact map.

    ΔΦ is the full-peptide alanine-scan difference — the sum of the per-position native→Ala ΔΔGs of
    :func:`alanine_scan`. It subtracts the interface's identity-independent baseline Φ(reference), i.e.
    what the *pose geometry* scores when every peptide residue is ``reference_aa``, leaving the
    sequence-specific part.

    On a **fixed** contact map this is Φ(peptide) minus a constant, so it does not change the ranking of
    candidates threaded onto one structure. It differs from raw Φ only across candidates that each have
    their **own** structure (e.g. AlphaFold peptide-swap models), where it normalises out the per-pose
    interface geometry. That normalisation rescues forced / wrong-register poses whose geometry corrupts
    the raw contact energy (the CPL ila1 case: TCR-ranking ROC 0.35 → 0.83), at a small cost on clones
    where the generated geometry is itself informative — so it is a scoring mode for generated poses,
    not a default. It is **not** an affinity ΔΔG: a dimensionless contact-preference difference, not a
    free energy (see :mod:`tcren.refine.register` for the geometry defect it corrects). Empirically both
    raw Φ and ΔΦ are within-receptor *ranking* scores, not binding constants — on the ATLAS SPR set they
    correlate with ΔG/Kd/koff/kon only at ρ ≤ 0.3 in magnitude (off-rate comes from :mod:`tcren.mechanics`).

    Args:
        contact_map: The candidate's own contact map.
        peptide: The candidate peptide sequence.
        potential: Pairwise potential to score with.
        interface: Which interface to score over (default ``"tcr_peptide"``).
        reference_aa: The amino acid the reference peptide is made of (default alanine).
        tcr_regions: Which TCR regions to keep on the TCR side.
        contact_weight: ``"residue"`` (default) or ``"atomic"``.
        structure: **Build** the reference peptide and score it on its own contact map instead of
            re-indexing it on the candidate's. This changes what ΔΦ means. Virtually, the poly-Ala
            baseline is charged for every contact the real side chains make, so ΔΦ measures only
            the substitution of identities on a fixed contact set. Structurally, the baseline is
            what the *backbone plus Cβ* alone can reach, so ΔΦ measures what the side chains
            contribute at all -- which is what the poly-alanine reference is meant to mean, and on
            1ao7 is the difference between 29 TCR:peptide contacts and 14.
        cutoff: Contact threshold for the rebuilt reference map (Å).
        sidechain: Passed to the rebuilt reference map, so a side-chain-filtered score filters the
            reference by the reference's own reach.

    Returns:
        ΔΦ = Φ(peptide) − Φ(reference); more negative = the sequence adds more favourable contacts than
        the reference baseline. ``0.0`` for interfaces without the peptide (e.g. ``"tcr_mhc"``).
    """
    reference = reference_aa * len(peptide)
    return ddg(contact_map, peptide, reference, potential,
               interface=interface, tcr_regions=tcr_regions, contact_weight=contact_weight,
               structure=structure, cutoff=cutoff, sidechain=sidechain)
