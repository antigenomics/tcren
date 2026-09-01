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


#: Which interfaces a smoothed reference varies, and which side of each carries the varying chain.
#:
#: The complex Hamiltonian is
#: :math:`\Phi = c_{\mathrm{TP}}\Phi_{\mathrm{TCR:pep}} + c_{\mathrm{TM}}\Phi_{\mathrm{TCR:MHC}}
#: + c_{\mathrm{PM}}\Phi_{\mathrm{pep:MHC}}`, so a substitution on one chain leaves one whole term
#: untouched and that term drops out of the difference:
#:
#: * varying the **peptide** kills :math:`\Delta\Phi_{\mathrm{TCR:MHC}}` -- no peptide residue is in it;
#: * varying the **TCR** kills :math:`\Delta\Phi_{\mathrm{pep:MHC}}` -- no TCR residue is in it.
#:
#: Each remaining interface is scored with its own potential and divided by its own Native2026
#: scale, so the two surviving terms are commensurate before they are added.
SMOOTH_INTERFACES: dict[str, tuple[tuple[str, str], ...]] = {
    "peptide": (("tcr_peptide", "to"), ("peptide_mhc", "from")),
    "tcr": (("tcr_peptide", "from"), ("tcr_mhc", "from")),
}


def smoothed_reference(
    contact_map: ContactMap,
    potential: Potential,
    *,
    side: str = "peptide",
    beta: float = 1.0,
    background=None,
    tcr_regions: str = "all",
    chain: str | None = None,
    mhc_potential: Potential | None = None,
    weights: dict | None = None,
) -> dict:
    r"""Boltzmann-smoothed reference difference :math:`\delta\Phi` and its curvature.

    The hard reference of :func:`reference_delta` subtracts the energy of one arbitrary sequence
    (poly-alanine). This subtracts the **free energy of the residue background** instead, so the
    baseline is a distribution rather than a choice of amino acid.

    Both interfaces that contain the varying chain are summed; the third drops out identically (see
    :data:`SMOOTH_INTERFACES`). Because no interface energy carries a within-chain term,
    :math:`\Phi` is a sum of independent local fields over the varying positions,

    .. math::  \varphi_i(a) \;=\; \sum_{\text{interfaces } I} c_I
               \sum_{j \,:\, (i,j) \in C_I} e_I(a, y_j)

    -- position :math:`i` of the varying chain, amino acid :math:`a`, summed over that position's
    contacts with the frozen partners :math:`y`, each interface weighted by its own
    :math:`c_I = 1/\mathrm{sd}_{\mathrm{Native2026}}(\Phi_I)`. The partition function therefore
    factorizes exactly, and the reference free energy is available in closed form:

    .. math::  \Phi_{\mathrm{ref}} \;=\; -\frac{1}{\beta}\sum_i
               \log \sum_a p(a)\, e^{-\beta \varphi_i(a)}

    with :math:`p` the background composition over the 20 amino acids. Then

    .. math::  \delta\Phi \;=\; \Phi(\text{observed}) - \Phi_{\mathrm{ref}}, \qquad
               \operatorname{Var}\Phi \;=\; \sum_i \operatorname{Var}_{\beta}\!\left[\varphi_i\right]

    where the variance is taken under the tilted weights
    :math:`p(a)e^{-\beta\varphi_i(a)}/\sum_b p(b)e^{-\beta\varphi_i(b)}`.

    :math:`\delta\Phi` is a first difference in sequence, against a smooth baseline;
    :math:`\operatorname{Var}\Phi` is the second cumulant of the same log partition function, i.e. how
    sharply that position's energy responds to residue identity at all. A position whose twenty
    fields are equal contributes nothing to either; one with a single strongly preferred residue
    contributes to both.

    :math:`\beta` sets how much of the background is averaged over. :math:`\beta \to 0` gives the
    arithmetic mean field :math:`\varphi_i(a_i) - \langle\varphi_i\rangle_p`, which is the reference
    state a combinatorial peptide library actually realises (every other position held at an
    equimolar mixture); :math:`\beta \to \infty` gives :math:`\varphi_i(a_i) - \min_a \varphi_i(a)`,
    the distance from the best residue available at that position. The default :math:`\beta = 1` is
    the potential's own scale, since a Boltzmann-inverted potential is already in units of
    :math:`k_{\mathrm B}T`.

    Args:
        contact_map: the structure's contact map.
        potential: the TCR:peptide potential (TCRen2).
        side: ``"peptide"`` (vary the peptide, receptor frozen -- the peptide scan) or ``"tcr"``
            (vary the receptor, peptide frozen -- the TCR scan). See :data:`SMOOTH_INTERFACES`.
        beta: inverse temperature in the potential's units (default 1.0).
        background: 20-vector of amino-acid frequencies in :data:`tcren.scoring.RecognitionMatrix`
            column order, or ``None`` (default) for the equimolar background.
        tcr_regions: TCR-region filter, applied on the TCR side.
        chain: restrict a ``side="tcr"`` scan to one chain (``"TRA"`` or ``"TRB"``), so the two
            chains can be read apart rather than pooled. ``None`` (default) keeps both.
        mhc_potential: the potential for the presentation interface (default Miyazawa-Jernigan,
            which is what :mod:`tcren.pipeline` assigns there).
        weights: ``{interface: coefficient}`` override; ``None`` (default) reads the Native2026
            scales through :func:`tcren.pipeline._phi_scale`.

    Returns:
        ``{"dPhi": float, "varPhi": float, "n_positions": int}``. Both sums are ``0.0`` over an empty
        position set, which is what an interface with no contacts on that side should score.
    """
    import numpy as np

    from .pipeline import _phi_scale
    from .potential import mj
    from .scoring import recognition_matrix

    if side not in SMOOTH_INTERFACES:
        raise ValueError(f"side must be one of {sorted(SMOOTH_INTERFACES)}, got {side!r}")
    if beta <= 0:
        raise ValueError(f"beta must be positive, got {beta}")

    mhc_pot = mj() if mhc_potential is None else mhc_potential
    fields: dict[tuple, "np.ndarray"] = {}
    aa: tuple = ()
    for iface, which in SMOOTH_INTERFACES[side]:
        pot = potential if iface == "tcr_peptide" else mhc_pot
        c = (weights or {}).get(iface, 1.0 / _phi_scale(iface, pot))
        rm = recognition_matrix(contact_map, pot, interface=iface, side=which,
                                tcr_regions=tcr_regions)
        aa = aa or rm.aa
        for i, key in enumerate(rm.positions):
            if chain is not None and key[0] != chain:
                continue
            v = c * np.asarray(rm.energy, float)[i]
            fields[key] = v if key not in fields else fields[key] + v

    if not fields:
        return {"dPhi": 0.0, "varPhi": 0.0, "n_positions": 0}

    keys = list(fields)
    phi = np.vstack([fields[k] for k in keys])                       # (n_positions, 20)
    native = np.array([aa.index(k[3]) if k[3] in aa else -1 for k in keys])

    p = (np.full(20, 1.0 / 20.0) if background is None
         else np.asarray(background, float) / np.sum(background))

    # An amino acid the potential leaves undefined (NaN) carries no weight, rather than poisoning
    # the whole position -- the same rule `score_peptides` applies when it drops those contacts.
    ok = np.isfinite(phi)
    w = np.where(ok, p[None, :], 0.0)
    e = np.where(ok, phi, 0.0)

    shift = np.min(np.where(ok, phi, np.inf), axis=1, keepdims=True)   # log-sum-exp, per position
    u = w * np.exp(-beta * (e - shift))
    z = u.sum(axis=1)
    live = z > 0
    if not live.any():
        return {"dPhi": 0.0, "varPhi": 0.0, "n_positions": 0}

    f_ref = np.where(live, shift[:, 0] - np.log(np.where(live, z, 1.0)) / beta, 0.0)
    q = np.where(live[:, None], u / np.where(live, z, 1.0)[:, None], 0.0)   # tilted weights
    m1 = (q * e).sum(axis=1)
    var = (q * e * e).sum(axis=1) - m1 * m1

    got = live & (native >= 0)
    phi_obs = np.where(got, e[np.arange(len(native)), np.clip(native, 0, 19)], 0.0)
    return {"dPhi": float(np.sum(phi_obs[got] - f_ref[got])),
            "varPhi": float(np.sum(var[live])),
            "n_positions": int(got.sum())}
