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
#: Both peptide-bearing interfaces at once, each with its own potential. The CPL response
#: matrix has always summed them (:func:`tcren.cpl.response_matrix`); a whole peptide could
#: not be scored the same way, so a library ranking silently saw the receptor term alone.
_COMPLEX = "complex"


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
    """The mutant's own contact map, from rebuilt coordinates.

    Threads the whole peptide, so every position is truncated to backbone + Cβ. That is right when
    every position is genuinely mutated -- the poly-alanine reference -- and wrong for a single
    substitution, which must leave its neighbours' side chains alone
    (:func:`_point_mutant_map`).
    """
    from .refine.substitute import substitute_peptide

    return ContactMap.from_structure(
        substitute_peptide(structure, mutant), cutoff=cutoff, sidechain=sidechain
    )


def _peptide_sites(structure: Structure) -> tuple[str, list[int]]:
    """The peptide chain's id and its residues' ``seq_index`` values, in sequence order."""
    from .structure.model import PEPTIDE_TYPE

    pep = next((c for c in structure.chains if c.chain_type == PEPTIDE_TYPE), None)
    if pep is None:
        raise ValueError(f"no peptide chain in structure {structure.pdb_id!r}")
    return pep.chain_id, [r.seq_index for r in pep.residues]


def _point_mutant_map(
    structure: Structure, chain_id: str, seq_index: int, aa: str, cutoff: float, sidechain: bool
) -> ContactMap:
    """Contact map with exactly one residue re-typed in 3D; every other side chain untouched."""
    from .refine.substitute import substitute_residues

    return ContactMap.from_structure(
        substitute_residues(structure, {(chain_id, seq_index): aa}),
        cutoff=cutoff, sidechain=sidechain,
    )


def ddg(
    contact_map: ContactMap,
    native: str,
    mutant: str,
    potential: Potential,
    *,
    interface: Interface | str = "tcr_peptide",
    tcr_regions: str = "all",
    contact_weight: str = "residue",
    structure: Structure | None = None,
    cutoff: float = 5.0,
    sidechain: bool = False,
    weights: np.ndarray | None = None,
    mhc_potential: Potential | None = None,
) -> float:
    """ΔΔG of a peptide mutation as ``E(native) - E(mutant)``.

    Args:
        contact_map: The structure's contact map.
        native: Native peptide sequence.
        mutant: Mutant peptide sequence (same length as ``native``).
        potential: Pairwise potential to score with.
        interface: Which interface to score over (default ``"tcr_peptide"``). ``"complex"``
            scores BOTH peptide-bearing interfaces and sums them -- ``potential`` over
            TCR:peptide plus ``mhc_potential`` over peptide:MHC -- which is the convention
            :func:`tcren.cpl.response_matrix` has always used for a response-matrix cell, and the
            one an activation read-out needs: the assay fires only if the peptide is presented AND
            the receptor engages. Scoring ``"tcr_peptide"`` alone answers a recognition question
            and is blind to presentation, so a peptide whose anchors are destroyed scores like any
            other. Note the two channels are NOT separable in a library that varies every position.
        mhc_potential: The peptide:MHC potential used by ``interface="complex"``. ``None``
            (default) is Miyazawa-Jernigan, matching :func:`tcren.cpl.response_matrix`. Ignored
            for any single-interface call.
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
    if interface == _COMPLEX:
        from .potential import mj
        # `weights` reweights the RECEPTOR channel only, exactly as `response_matrix`'s
        # `tcr_weights` does; the presentation channel keeps the map's own indicator.
        common = dict(tcr_regions=tcr_regions, contact_weight=contact_weight,
                      structure=structure, cutoff=cutoff, sidechain=sidechain)
        return (ddg(contact_map, native, mutant, potential, interface="tcr_peptide",
                    weights=weights, **common)
                + ddg(contact_map, native, mutant, mhc_potential or mj(),
                      interface="peptide_mhc", **common))
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

    With ``structure`` **only that position** is truncated to alanine in 3D and the contact map
    recomputed, so a position whose side chain was the only thing reaching the TCR loses those
    contacts as it physically must, while its neighbours keep theirs. This is the case the
    structural path gets exactly right (see the module docstring), and it costs one contact-map
    rebuild per position.

    Before 2.25.0 this path threaded the whole peptide through
    :func:`~tcren.refine.substitute.substitute_peptide`, which truncates *every* residue to
    backbone + Cβ. The scan therefore measured each position against a poly-stub baseline rather
    than the native: on 1ao7 the native sequence threaded back through it kept 14 of 29
    TCR:peptide contacts, and the resulting offset appeared in every position, including
    positions with no contacts at all.

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
    chain_id, seq_indices = _peptide_sites(structure) if structure is not None else ("", [])
    if structure is not None and len(seq_indices) != len(native):
        raise ValueError(
            f"peptide chain has {len(seq_indices)} residues, native sequence has {len(native)}"
        )
    e_native = (
        _score_one(contact_map, native, potential, interface, tcr_regions, contact_weight)
        if peptide_iface
        else 0.0
    )
    rows = []
    for pos, wt in enumerate(native):
        if peptide_iface:
            mutant = native[:pos] + "A" + native[pos + 1 :]
            mutant_map = (
                contact_map if structure is None
                else _point_mutant_map(
                    structure, chain_id, seq_indices[pos], "A", cutoff, sidechain
                )
            )
            e_mut = _score_one(
                mutant_map, mutant, potential, interface, tcr_regions, contact_weight
            )
            ddg_val = e_native - e_mut
        else:
            ddg_val = 0.0
        rows.append({"pos": pos, "wt_aa": wt, "ddG": ddg_val})
    return pl.DataFrame(rows, schema={"pos": pl.Int64, "wt_aa": pl.Utf8, "ddG": pl.Float64})


#: TCR loops the receptor-side scan walks, and the per-loop aggregates it reports.
_CDR_LOOPS: tuple[str, ...] = ("CDR1", "CDR2", "CDR3")


def _tcr_mutant_map(
    structure: Structure, chain_id: str, seq_index: int, cutoff: float, sidechain: bool
) -> ContactMap:
    """The contact map of the structure with one receptor residue truncated to alanine, in 3D."""
    from .refine.substitute import substitute_residues

    return ContactMap.from_structure(
        substitute_residues(structure, {(chain_id, seq_index): "A"}),
        cutoff=cutoff, sidechain=sidechain,
    )


def _peptide_sequence(structure: Structure) -> str:
    from .structure.model import PEPTIDE_TYPE

    pep = next((c for c in structure.chains if c.chain_type == PEPTIDE_TYPE), None)
    if pep is None:
        raise ValueError(f"no peptide chain in structure {structure.pdb_id!r}")
    return "".join(r.aa for r in pep.residues)


def tcr_alanine_scan(
    contact_map: ContactMap,
    structure: Structure,
    potential: Potential,
    *,
    peptide: str | None = None,
    tcr_regions: str = "cdr",
    contact_weight: str = "residue",
    cutoff: float = 5.0,
    sidechain: bool = False,
) -> pl.DataFrame:
    """Alanine scan of the **receptor** side, on rebuilt coordinates.

    The mirror of :func:`alanine_scan`. Each contacted TCR residue is truncated to alanine **in
    3D** by :func:`tcren.refine.substitute.substitute_residues`, the contact map is recomputed and
    the interface rescored, so a loop residue whose side chain was the only thing bridging to the
    peptide loses those contacts exactly as it physically must. One rebuild per contacted residue.

    Only residues that actually contact the peptide are walked, because a residue with no contact
    has ``ddG == 0`` by construction and rebuilding it would cost a contact map for nothing.

    Scored over ``"tcr_peptide"`` alone: a receptor substitution cannot change the peptide:MHC
    energy, so the complex sum would only add a constant.

    Args:
        contact_map: The native structure's contact map.
        structure: The annotated structure the map came from. Required — there is no virtual
            path here, because truncating a receptor side chain without moving atoms would leave
            every contact it made in place, which is the failure mode this function exists to fix.
        potential: Pairwise potential to score with.
        peptide: Peptide sequence; taken from the structure's peptide chain when omitted.
        tcr_regions: Which TCR regions to walk — ``"cdr"`` (default), ``"cdr+fr"`` or ``"all"``.
        contact_weight: ``"residue"`` (default) or ``"atomic"``.
        cutoff: Contact threshold for the rebuilt maps (Å).
        sidechain: Passed to the rebuilt maps.

    Returns:
        One row per contacted receptor residue, with ``chain.id``, ``chain.type`` (TRA/TRB),
        ``region.type``, ``residue.index``, ``pos`` (0-based within its region), ``wt_aa`` and
        ``ddG`` = ``E(native) - E(Ala@residue)``. A **positive** ``ddG`` marks a stabilising
        residue: removing it costs energy, so it was earning its place.
    """
    peptide = peptide if peptide is not None else _peptide_sequence(structure)
    iface = contact_map.interface("tcr_peptide", tcr_regions=tcr_regions)
    sites = (
        iface.select("chain.id.from", "chain.type.from", "region.type.from",
                     "residue.index.from", "residue.aa.from", "pos.from")
        .unique()
        .sort("chain.type.from", "residue.index.from")
    )
    schema = {"chain.id": pl.Utf8, "chain.type": pl.Utf8, "region.type": pl.Utf8,
              "residue.index": pl.Int64, "pos": pl.Int64, "wt_aa": pl.Utf8, "ddG": pl.Float64}
    if sites.height == 0:
        return pl.DataFrame(schema=schema)

    e_native = _score_one(
        contact_map, peptide, potential, "tcr_peptide", tcr_regions, contact_weight
    )
    rows = []
    for chain_id, chain_type, region, index, wt, pos in sites.iter_rows():
        mutant_map = _tcr_mutant_map(structure, chain_id, index, cutoff, sidechain)
        e_mut = _score_one(
            mutant_map, peptide, potential, "tcr_peptide", tcr_regions, contact_weight
        )
        rows.append({"chain.id": chain_id, "chain.type": chain_type, "region.type": region,
                     "residue.index": index, "pos": pos, "wt_aa": wt,
                     "ddG": e_native - e_mut})
    return pl.DataFrame(rows, schema=schema)


def tcr_alanine_reference(scan: pl.DataFrame) -> dict[str, float]:
    """Per-loop poly-alanine references, summed from a :func:`tcr_alanine_scan`.

    Four numbers per structure: the germline loops together, each CDR3 on its own, and their
    total. Each is the **sum of the per-residue 3D ΔΔGs** of that loop, which is the additive
    reading and is defined whether or not a loop is engaged (an unengaged loop contributes 0).

    It is deliberately *not* the energy of mutating a whole loop to poly-alanine in one pass.
    Those differ once atoms move: truncating every side chain at once loses contacts that each
    residue alone retains, so the one-pass value is not the sum of the parts. The additive form
    is the one that says how much each residue earns.

    Args:
        scan: The frame returned by :func:`tcr_alanine_scan`.

    Returns:
        ``dPhi_ala_cdr12``, ``dPhi_ala_cdr3a``, ``dPhi_ala_cdr3b`` and ``dPhi_ala_tcr``.
    """
    def total(pred: pl.Expr) -> float:
        sel = scan.filter(pred)
        return float(sel["ddG"].sum()) if sel.height else 0.0

    region, chain = pl.col("region.type"), pl.col("chain.type")
    return {
        "dPhi_ala_cdr12": total(region.is_in(["CDR1", "CDR2"])),
        "dPhi_ala_cdr3a": total((region == "CDR3") & (chain == "TRA")),
        "dPhi_ala_cdr3b": total((region == "CDR3") & (chain == "TRB")),
        "dPhi_ala_tcr": float(scan["ddG"].sum()) if scan.height else 0.0,
    }


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
    interface: Interface | str = "tcr_peptide",
    reference_aa: str = "A",
    tcr_regions: str = "all",
    contact_weight: str = "residue",
    structure: Structure | None = None,
    cutoff: float = 5.0,
    sidechain: bool = False,
    mhc_potential: Potential | None = None,
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
        interface: Which interface to score over (default ``"tcr_peptide"``). ``"complex"``
            sums both peptide-bearing interfaces, which is the whole-complex ΔΦ a combinatorial
            library ranking needs -- the receptor term alone cannot see a destroyed anchor.
        mhc_potential: peptide:MHC potential for ``interface="complex"`` (default
            Miyazawa-Jernigan).
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
               structure=structure, cutoff=cutoff, sidechain=sidechain,
               mhc_potential=mhc_potential)
