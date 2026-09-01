"""Per-structure interface descriptors: the catalogue, and the extractor that fills it.

:data:`DESCRIPTORS` is the single registry -- every emitted column, its family and whether the
receptor enters its definition. :data:`INVARIANCE` says what each quantity is invariant under, and
:data:`DETAIL` its units and one-line meaning; the docs table is generated from those three, so a
descriptor cannot reach a feature table undocumented.

This module emits **descriptors only**. The fitted composites that used to ride along here -- the
Gaussian Bayes-net and Bayesian-logistic real-vs-shuffled recognizers, the frozen forced-pose
logistic, the fitted binder score and the cohort posterior -- were removed in 2.26.0: their
coefficients were frozen against training sets that no longer exist, which made them the one part
of the package a reader could not reproduce. Scoring is :mod:`tcren.reliability` on this table.
"""
from __future__ import annotations

import math
import warnings
from collections.abc import Sequence

import numpy as np

from .footprint import (FOOTPRINT_SIZE_FEATURES, footprint_topology_features)

_EPS = 1e-9


# ===================================================================================================
# Structure -> the 35-descriptor recognition vector the frozen recognizers consume, and P(real).
#
# Reproduces the extractor the shipped models were trained on (the manuscript's compute_features.py):
# docking geometry + per-interface TCRen/MJ energetics (F, poly-Ala ΔF) + contact-type tallies +
# biopython ΔSASA burial + MHC-class indicator. Heavy imports are function-local so that a bare
# ``import tcren`` (and ``import tcren.recognition``) stays dependency-light.
# ===================================================================================================

#: The core descriptor block ``recognize`` emits, and the vector the frozen recognizers consume.
#:
#: Every statistical-potential energy is named ``Phi_*`` — there is one potential per interface (TCRen
#: on TCR:peptide, MJ on the two presentation interfaces), so the potential's name does not belong in
#: the column's. Two exact duplicates were dropped in the 2026-07-28 audit: ``e_tcr_mhc`` (the same
#: number as ``Phi_tcr_mhc``) and ``ct_tp_hydrogen_bond`` (the same number as ``n_hbond``, which is the
#: name Eq. Q uses). The frozen models still ask for the old names and get the same values through
#: :data:`_FROZEN_ALIASES`, so their predictions are unchanged.
RECOGNITION_FEATURES = (
    "extent", "chain_balance", "pitch", "crossing", "crossing_signed", "dock_d", "dock_torsion",
    "dock_tcr_uy", "dock_tcr_uz", "dock_mhc_uy", "dock_mhc_uz", "Phi_cdr12", "Phi_cdr3a", "Phi_cdr3b",
    "Phi_tcr_pep", "Phi_tcr_mhc", "Phi_pep_mhc", "dPhi_tcr_pep", "dPhi_pep_mhc",
    "dPhi_pep_soft", "varPhi_pep_soft", "dPhi_tcr_soft", "varPhi_tcr_soft",
    "dPhi_tra_soft", "dPhi_trb_soft", "n_contacts_tp",
    "n_pep_contacted", "n_contacts_tm", "ct_tp_salt_bridge", "ct_tm_salt_bridge",
    "ct_tm_hydrogen_bond", "ct_tp_aromatic", "ct_tm_aromatic",
    "ct_tp_hydrophobic", "ct_tm_hydrophobic", "ct_tp_other", "ct_tm_other", "n_hbond",
    "burial", "mhc_class_bin",
)

_CT_TYPES = ("salt_bridge", "hydrogen_bond", "aromatic", "hydrophobic", "other")
_TCR_TYPES = ("TRA", "TRB", "TRG", "TRD")

#: Interface-symmetry descriptors from per-loop TCR:peptide contact **counts** (not energies), emitted as
#: extra ``recognize`` output columns — **not** part of :data:`RECOGNITION_FEATURES` (the frozen models'
#: 35-vector is fixed). ``cdr3_dominance`` = CDR3(α+β) share of CDR contacts (higher = CDR3-dominated,
#: oriented positive); ``cdr3_ab_imbalance`` = ``|CDR3α−CDR3β|`` normalised (absolute); ``chain_cdr_imbalance``
#: = ``|α−β|`` whole-CDR normalised (absolute). See :func:`_interface_symmetry`.
INTERFACE_SYMMETRY_FEATURES = ("cdr3_dominance", "cdr3_ab_imbalance", "chain_cdr_imbalance")

#: Where the receptor body sits over the groove (:func:`tcren.orient.tcr_placement`), emitted as extra
#: output columns — **not** part of :data:`RECOGNITION_FEATURES` (the frozen models' vector is fixed).
#: ``height`` = elevation of the CDR centroid above the groove plane, ``shift_u``/``shift_w`` its
#: in-plane displacement from the peptide centroid along the groove long/short axes, ``offset`` the
#: in-plane distance. These are the translational degrees of freedom no docking *angle* can see, and
#: the mechanism behind the coverage entropy (uniform coverage = riding low).
TCR_PLACEMENT_FEATURES = ("height", "shift_u", "shift_w", "offset")

#: The intra-peptide term, emitted as extra ``recognize --full`` output columns — **not** part of
#: :data:`RECOGNITION_FEATURES` (the frozen models' 35-vector is fixed). ``Phi_pep_int`` = the peptide's
#: MJ contact energy with **itself** (:func:`tcren.intra_peptide_energy`), the term every interface
#: sum omits; ``n_pep_int`` = how many such contacts there are. Both are properties of the pMHC alone
#: — no receptor enters them — so they carry cohort identity; see :data:`DESCRIPTORS`.
PEPTIDE_INTERNAL_FEATURES = ("Phi_pep_int", "n_pep_int")

#: CDR3-local frame features (18), the FramePose layer the whole-TCR :data:`RECOGNITION_FEATURES` miss.
#: Per loop, relative to the pMHC groove frame (u, w, n; origin = peptide Cα centroid):
#: ``reach`` = |loop centroid − origin|; ``o{u,w,n}`` = unit(centroid−origin)·(u,w,n) (where over the
#: groove the loop sits); ``a{u,w,n}`` = unit(Cα_N→Cα_C)·(u,w,n) (loop orientation over the groove);
#: ``topep`` = min Cα-Cα distance loop→peptide (engagement depth); ``ext`` = |Cα_C − Cα_N| (extension).
_CDR3_FRAME_KEYS = ("reach", "ou", "ow", "on", "au", "aw", "an", "topep", "ext")
CDR3_FRAME_FEATURES = tuple(f"{loop}_{k}" for loop in ("cdr3a", "cdr3b") for k in _CDR3_FRAME_KEYS)

#: The full feature vector: the core recognition descriptors + the 18 CDR3-frame descriptors.
#:
#: The 12 "matrix-swap" columns (``tcren_{g}``/``mj_{g}``/``d_{g}`` for ``g`` ∈ {tp, cdr12, cdr3a,
#: cdr3b}) were removed in the 2026-07-28 audit. The ``tcren_*`` four were exact duplicates of the
#: ``Phi_*`` energies; the ``mj_*`` four scored TCR:peptide contacts under the generic MJ potential,
#: which is not the potential this method uses on that interface, and the ``d_*`` four were their
#: difference. Nothing consumed them.
FULL_FEATURES = RECOGNITION_FEATURES + CDR3_FRAME_FEATURES

# ===================================================================================================
# The descriptor catalogue: what every emitted column is, and whether the receptor enters it.
# ===================================================================================================
#: Family of each descriptor, and whether the TCR enters its definition.
#:
#: Five families, split by **what each quantity is invariant under** — which is also the axis along
#: which they carry independent evidence:
#:
#: * ``placement`` — where the receptor sits, expressed in the pMHC groove frame: docking angles,
#:   the TCRdock rigid-body parameters, the ride height/shift/offset of the receptor body, and the
#:   per-loop CDR3 frame descriptors. **Frame-dependent**: these change if the groove frame does.
#: * ``interface`` — how much contact there is and of what chemical kind: buried area, contact
#:   counts and types, hydrogen bonds, clashes, chain and loop balance. SE(3)-invariant. This is the
#:   channel Eq. Q is built from.
#: * ``topology`` — the *shape* of the contact set, independent of both its size and its chemistry:
#:   coverage entropy and Hill numbers over the CDR-loop x target cells, the footprint's Betti
#:   numbers and persistence entropy, the canonical germline/CDR3 preference. SE(3)-invariant, which
#:   is why these need no canonical orientation (:mod:`tcren.footprint`).
#: * ``energetics`` — statistical-potential interface energies ``F`` and their poly-alanine
#:   references ``dF``. Lower is more favourable. SE(3)-invariant.
#: * ``kinetics`` — the interface as a network of breakable springs: stiffness, anisotropy, strain,
#:   rupture, and the residues that couple the pre-formed scaffold to the interface.
#:
#: ``placement`` and ``interface`` were one ``geometry`` family until 2026-08-24. Splitting them is
#: what lets the three-channel claim be stated at all: the coverage entropy is coupled to the ride
#: height (Spearman -0.559 / -0.525) and so is *not* independent of ``placement``, while it is a
#: different question whether it is independent of ``interface``. :func:`descriptors` keeps
#: ``"geometry"`` and ``"physics"`` working as aliases.
#:
#: ``involves_tcr`` is ``False`` for a quantity computed from the peptide and the MHC alone. Such a
#: column is a property of the *cohort*, not of the receptor: two structures of the same epitope on
#: the same allele share its value whatever their TCR. A model handed one can reach a cohort-level
#: label through epitope or allele identity instead of through interface physics, so any analysis
#: whose question is about receptors must select with ``tcr_only=True`` (:func:`descriptors`).
#:
#: Fitted and cohort-relative composites (``p_real``, ``p_real_bn``, ``p_forced``, ``p_bind``,
#: ``q_bind``, ``s_strain``) are listed under ``score``. They are outputs built from the descriptors
#: above and must never be fed back in as inputs; :func:`descriptors` excludes them by default.
DESCRIPTORS: dict[str, tuple[str, bool]] = {
    # -- placement: rigid-body pose in the groove frame ----------------------------------------
    "pitch": ("placement", True),            # incident/tilt angle out of the groove plane
    "crossing": ("placement", True),         # crossing (scanning) angle against the groove long axis
    "crossing_signed": ("placement", True),  # the same, signed: carries the docking polarity
    "dock_d": ("placement", True),
    "dock_torsion": ("placement", True),
    "dock_tcr_uy": ("placement", True),
    "dock_tcr_uz": ("placement", True),
    "dock_mhc_uy": ("placement", True),
    "dock_mhc_uz": ("placement", True),
    # where the receptor body sits over the groove (`tcren.orient.tcr_placement`)
    "height": ("placement", True),
    "shift_u": ("placement", True),
    "shift_w": ("placement", True),
    "offset": ("placement", True),
    # the CDR3 loops' own placement in the same frame (the FramePose layer)
    **{f: ("placement", True) for f in CDR3_FRAME_FEATURES},
    # -- interface: contact size and chemistry --------------------------------------------------
    "burial": ("interface", True),
    "extent": ("interface", True),
    "chain_balance": ("interface", True),
    "n_contacts_tp": ("interface", True),
    "n_contacts_tm": ("interface", True),
    "n_pep_contacted": ("interface", True),
    "n_hbond": ("interface", True),
    "ct_tp_salt_bridge": ("interface", True),
    "ct_tp_aromatic": ("interface", True),
    "ct_tp_hydrophobic": ("interface", True),
    "ct_tp_other": ("interface", True),
    "ct_tm_salt_bridge": ("interface", True),
    "ct_tm_hydrogen_bond": ("interface", True),
    "ct_tm_aromatic": ("interface", True),
    "ct_tm_hydrophobic": ("interface", True),
    "ct_tm_other": ("interface", True),
    "cdr3_dominance": ("interface", True),
    "cdr3_ab_imbalance": ("interface", True),
    "chain_cdr_imbalance": ("interface", True),
    "n_clashes": ("interface", True),
    "clash_score": ("interface", True),
    # the MHC class indicator is a property of the presenting molecule alone
    "mhc_class_bin": ("interface", False),
    # the raw footprint contact counts are interface SIZE, not shape -- see FOOTPRINT_SIZE_FEATURES.
    # The total among them is `n_loop_contacts`; bare `n_contacts` belongs to `potts` below.
    **{f: ("interface", True) for f in FOOTPRINT_SIZE_FEATURES},
    # -- topology: the shape of the contact set (`tcren.footprint`) ------------------------------
    **{f: ("topology", True) for f in footprint_topology_features()},
    # -- energetics: interface energies ----------------------------------------------------------
    "Phi_tcr_pep": ("energetics", True),
    "Phi_tcr_mhc": ("energetics", True),
    "Phi_cdr12": ("energetics", True),
    "Phi_cdr3a": ("energetics", True),
    "Phi_cdr3b": ("energetics", True),
    "dPhi_tcr_pep": ("energetics", True),
    # the same first difference against a smoothed background, both directions, TCR side split by
    # chain so a linear model can form the TRB-TRA contrast rather than being handed it
    "dPhi_pep_soft": ("energetics", True),
    "varPhi_pep_soft": ("energetics", True),
    "dPhi_tcr_soft": ("energetics", True),
    "varPhi_tcr_soft": ("energetics", True),
    "dPhi_tra_soft": ("energetics", True),
    "dPhi_trb_soft": ("energetics", True),
    # peptide:MHC energy and its poly-alanine reference: presentation, no receptor
    "Phi_pep_mhc": ("energetics", False),
    "dPhi_pep_mhc": ("energetics", False),
    # the peptide's contacts with itself: a property of the epitope's bound conformation alone
    "Phi_pep_int": ("energetics", False),
    "n_pep_int": ("interface", False),
    # -- potts: the contact map read against the coupled model (`tcren.potts`) -------------------
    # The same energy, referenced against the partition function instead of a poly-alanine
    # interface: `neg_energy = log_z + log_lik`, so the three carry capacity, typicality and their
    # sum. `n_sites` is how many pairs the backbone put in reach, `n_contacts` how many engaged.
    #
    # `n_contacts` is catalogued HERE and nowhere else. Through 2.19.0 it was also the name of the
    # footprint's CDR-loop tally, which is a different count on the same structure (1ao7: 29 here,
    # 66 there), and since the topology pass runs before the potts pass the emitted column meant
    # whichever family the caller happened to ask for. `tcren.reliability` standardizes it against
    # the Potts population, so the footprint tally reached the correction as a ~6-sd outlier with
    # no warning. The footprint one is now `n_loop_contacts`.
    "neg_energy": ("potts", True),
    "log_z": ("potts", True),
    "log_lik": ("potts", True),
    "psi": ("potts", True),
    "n_contacts": ("potts", True),
    # -- kinetics: contact fragility (``recognize``) -------------------------------------------
    "exp_lost": ("kinetics", True),
    "mean_margin": ("kinetics", True),
    "frac_robust": ("kinetics", True),
    # -- kinetics: the spring network (``mechanics``) ------------------------------------------
    "K_tens": ("kinetics", True),
    "K_shear": ("kinetics", True),
    "aniso": ("kinetics", True),
    "lam_max": ("kinetics", True),
    "lam_min": ("kinetics", True),
    "n_spring": ("kinetics", True),
    "S_tot": ("kinetics", True),
    "rupture_force": ("kinetics", True),
    "rupture_work": ("kinetics", True),
    "couple_pep": ("kinetics", True),
    "couple_mhc": ("kinetics", True),
    "couple_tcr": ("kinetics", True),
    "couple_total": ("kinetics", True),
    "n_interface": ("kinetics", True),
}

#: What each descriptor is invariant under -- the axis along which ``geometry`` and ``topology``
#: are different questions rather than two names for the contact set.
#:
#: *Geometry is the study of properties preserved by distance-preserving transformations;
#: topology is the study of properties preserved by continuous deformation.* Applied here:
#:
#: ``"geometric"``
#:     A continuous quantity in physical units -- a length (A), an area (A^2), an angle, or a
#:     direction cosine. Preserved by isometry, destroyed by deformation. **This is the docking**:
#:     where the receptor sits on the groove and how it leans.
#: ``"topological"``
#:     An invariant of the contact complex under continuous deformation -- Betti numbers, the
#:     Euler characteristic, and their size-normalized forms. **This is the interface surface**:
#:     how many patches it falls into and how many holes it has, whatever its shape.
#: ``"compositional"``
#:     A count over the *labelled* contact set, or a ratio, share, entropy or Hill number built
#:     from such counts. Preserved by both, because it reads the labelling rather than the shape.
#: ``"energetic"``
#:     A statistical-potential or Potts energy.
#: ``"categorical"``
#:     Which class of MHC presents the peptide -- class I or class II.
#:
#: Two consequences worth knowing before building a block from a family. The **topology family is
#: mostly compositional**: 20 of its 29 columns are diversity or coverage measures over labelled
#: cells and positions, and only 8 are topological invariants. And ``h0_pers_ent`` is filed
#: ``"geometric"``, not ``"topological"``, because the H0 barcode's bar lengths *are* the minimum
#: spanning tree's edge lengths in angstroms -- persistent homology is a metric construction, and
#: the entropy of a length distribution is not a homeomorphism invariant.
INVARIANCE: dict[str, str] = {
    # placement is metric throughout: distances, angles and direction cosines in the groove frame
    **{d: "geometric" for d, (fam, _) in DESCRIPTORS.items() if fam == "placement"},
    # interface: two continuous quantities, the rest counts and count ratios
    "burial": "geometric",              # dSASA, A^2
    "clash_score": "geometric",         # summed heavy-atom overlap, A
    "mhc_class_bin": "categorical",
    **{d: "compositional" for d in (
        "extent", "n_contacts_tp", "n_contacts_tm", "n_pep_contacted", "n_hbond",
        "ct_tp_salt_bridge", "ct_tp_aromatic", "ct_tp_hydrophobic", "ct_tp_other",
        "ct_tm_salt_bridge", "ct_tm_hydrogen_bond", "ct_tm_aromatic", "ct_tm_hydrophobic",
        "ct_tm_other", "n_clashes", "n_loop_contacts", "n_pep_contacts", "n_mhc_contacts",
        "n_pep_int", "chain_balance", "cdr3_dominance", "cdr3_ab_imbalance",
        "chain_cdr_imbalance",
    )},
    # topology: the Betti/Euler block is topological; the diversity block reads the labelling
    **{d: "topological" for d in (
        "fp_b0_r7", "fp_b1_r7", "fp_chi_r7", "fp_b0_frac_r7",
        "fp_b0_r8", "fp_b1_r8", "fp_chi_r8", "fp_b0_frac_r8",
    )},
    "h0_pers_ent": "geometric",
    **{d: "compositional" for d in (
        "H_cell", "D1_cell", "D2_cell", "S_cell", "J_cell", "H_loop", "D2_loop", "D2_pep24",
        "ab_imb", "ab_imb_pep", "ab_imb_mhc", "L_canon", "p_germ_mhc", "p_cdr3_pep",
        "pep_free_frac", "pep_cov_frac", "pep_cov_even", "pep_cov_d2n", "pep_cov_centre",
        "pep_cov_spread",
    )},
    # energies
    **{d: "energetic" for d, (fam, _) in DESCRIPTORS.items() if fam in ("energetics", "potts")},
    # kinetics: stiffnesses, forces and work are continuous; the coupling tallies are counts
    **{d: "geometric" for d in (
        "K_tens", "K_shear", "S_tot", "aniso", "lam_max", "lam_min",
        "rupture_force", "rupture_work", "mean_margin", "exp_lost",
    )},
    **{d: "compositional" for d in (
        "n_spring", "n_interface", "couple_pep", "couple_mhc", "couple_tcr", "couple_total",
        "frac_robust",
    )},
}

#: Units and a one-line definition for every descriptor. The single source the docs table is
#: generated from, so a new descriptor cannot reach a feature table undocumented.
#:
#: ``units`` is what the number is measured in -- ``A``, ``A^2``, ``deg``, ``rad``, ``kT``, ``N/m``,
#: or one of the dimensionless kinds ``count``, ``fraction``, ``signed fraction``, ``ratio``,
#: ``cosine``, ``log-odds``, ``indicator``. It is what a transform has to respect: a count is
#: variance-stabilized by a square root, a fraction by the arcsine (the classical angular
#: transformation), and an unbounded continuous quantity by neither.
DETAIL: dict[str, tuple[str, str]] = {
 "pitch": ("deg", "Incident angle of the TCR out of the groove plane. **Banned as a feature**: it reproduces no clean geometric angle yet out-discriminates every one of them, which is AlphaFold-confidence contamination rather than geometry."),
 "crossing": ("deg", "Crossing (scanning) angle between the Valpha->Vbeta axis projected into the groove plane and the groove long axis."),
 "crossing_signed": ("deg", "The same angle on [-180, 180); its sign is the docking polarity, canonical or reversed."),
 "dock_d": ("A", "MHC-stub to TCR-stub rigid-body separation."),
 "dock_torsion": ("rad", "Rigid-body dihedral of the TCR about the MHC stub; the docking twist. Circular, wraps at +-pi."),
 "dock_tcr_uy": ("cosine", "y component of the TCR stub unit vector in the MHC frame."),
 "dock_tcr_uz": ("cosine", "z component of the TCR stub unit vector; how high the receptor body rides over the groove."),
 "dock_mhc_uy": ("cosine", "y component of the MHC stub unit vector."),
 "dock_mhc_uz": ("cosine", "z component of the MHC stub unit vector."),
 "height": ("A", "Elevation of the CDR Calpha centroid above the groove plane."),
 "shift_u": ("A", "In-plane displacement of that centroid from the peptide centroid along the groove long axis."),
 "shift_w": ("A", "The same along the groove short axis."),
 "offset": ("A", "Length of the in-plane displacement; lateral shift whatever its direction."),
}
_LOOP_DETAIL: dict[str, tuple[str, str]] = {
 "reach": ("A", "Distance from the loop's Calpha centroid to the peptide Calpha centroid; how far {loop} reaches."),
 "topep": ("A", "Minimum Calpha-Calpha distance from {loop} to the peptide; its engagement depth."),
 "ext": ("A", "End-to-end extension of {loop}, the Calpha_N to Calpha_C distance."),
 "ou": ("cosine", "Where {loop} sits over the groove, along the long axis u."),
 "ow": ("cosine", "Where {loop} sits over the groove, along the short axis w."),
 "on": ("cosine", "Where {loop} sits over the groove, along the groove normal n."),
 "au": ("cosine", "Orientation of {loop}'s N->C axis against the groove long axis u."),
 "aw": ("cosine", "Orientation of {loop}'s N->C axis against the groove short axis w."),
 "an": ("cosine", "Orientation of {loop}'s N->C axis against the groove normal n."),
}
for _loop, _name in (("cdr3a", "CDR3alpha"), ("cdr3b", "CDR3beta")):
    for _k, (_u, _d) in _LOOP_DETAIL.items():
        DETAIL[f"{_loop}_{_k}"] = (_u, _d.format(loop=_name))

DETAIL.update({
 # --- interface -------------------------------------------------------------------------------
 "burial": ("A^2", "Interface buried surface, SASA(TCR) + SASA(pMHC) - SASA(complex), by Shrake-Rupley."),
 "extent": ("count", "Distinct TCR residues contacting the pMHC over both receptor interfaces."),
 "chain_balance": ("fraction", "min(a,b)/(a+b) over TCR:peptide contacts by chain; 0.5 when both chains engage equally, 0 when only one does."),
 "n_contacts_tp": ("count", "TCR-peptide residue-residue contacts."),
 "n_contacts_tm": ("count", "TCR-MHC residue-residue contacts."),
 "n_pep_contacted": ("count", "Distinct peptide residues the TCR contacts."),
 "n_hbond": ("count", "Polar N/O atom pairs within 3.5 A across TCR:peptide."),
 "ct_tp_salt_bridge": ("count", "Cationic-N / anionic-O pairs within 4 A across TCR:peptide."),
 "ct_tp_aromatic": ("count", "Ring-atom pairs between aromatic residues across TCR:peptide."),
 "ct_tp_hydrophobic": ("count", "Apolar C-C pairs between apolar residues across TCR:peptide."),
 "ct_tp_other": ("count", "Remaining classified TCR:peptide contacts."),
 "ct_tm_salt_bridge": ("count", "Salt bridges across TCR:MHC."),
 "ct_tm_hydrogen_bond": ("count", "Hydrogen bonds across TCR:MHC."),
 "ct_tm_aromatic": ("count", "Aromatic contacts across TCR:MHC."),
 "ct_tm_hydrophobic": ("count", "Hydrophobic contacts across TCR:MHC."),
 "ct_tm_other": ("count", "Remaining classified TCR:MHC contacts."),
 "cdr3_dominance": ("fraction", "CDR3(alpha+beta) share of all CDR TCR:peptide contacts."),
 "cdr3_ab_imbalance": ("fraction", "abs(CDR3a - CDR3b) / (CDR3a + CDR3b); how one-sided the CDR3 engagement is."),
 "chain_cdr_imbalance": ("fraction", "abs(a - b) / (a + b) over all CDR contacts; the chain-level mirror of chain_balance."),
 "n_clashes": ("count", "Peptide-partner heavy-atom pairs overlapping by more than 0.4 A on Bondi radii."),
 "clash_score": ("A", "Summed overlap depth of those clashing pairs; the steric burden of a forced pose."),
 "mhc_class_bin": ("class I / II", "Which class of MHC presents the peptide: 0 for class I, 1 for class II. Class I and class II grooves differ in shape and in how they hold a peptide, so this conditions the other descriptors rather than being scored beside them; a coefficient fitted across both classes without it is fitted to a mixture."),
 "n_loop_contacts": ("count", "Contacts the six-CDR-loop partition sees; framework contacts are outside it by construction."),
 "n_pep_contacts": ("count", "Of those loop contacts, the ones reaching the peptide."),
 "n_mhc_contacts": ("count", "Of those loop contacts, the ones reaching the MHC."),
 "n_pep_int": ("count", "Intra-peptide residue contacts, at 5 A with a sequence separation of at least three."),
 # --- topology --------------------------------------------------------------------------------
 "H_cell": ("fraction", "Normalized Shannon entropy of the contact composition over the twelve cells (six CDR loops x {peptide, MHC})."),
 "D1_cell": ("count", "Hill number of order 1 over the same cells, exp(H); monotone in H_cell."),
 "D2_cell": ("count", "Hill number of order 2, 1/sum(p^2); the effective number of engaged cells, discounting the weakly populated ones."),
 "S_cell": ("count", "Richness: how many of the twelve cells are occupied."),
 "J_cell": ("fraction", "Pielou evenness over the occupied cells. NaN when one cell is occupied."),
 "H_loop": ("fraction", "Normalized entropy over the six CDR loops alone, ignoring which target each contact reaches."),
 "D2_loop": ("count", "Hill number of order 2 over the six loops."),
 "D2_pep24": ("count", "Hill number of order 2 over the twenty-four-cell partition, the peptide split into N-terminal, central and C-terminal bands."),
 "ab_imb": ("signed fraction", "Signed (TRA - TRB)/(TRA + TRB) over CDR-loop contacts; positive is alpha-shifted."),
 "ab_imb_pep": ("signed fraction", "The same restricted to peptide-side contacts."),
 "ab_imb_mhc": ("signed fraction", "The same restricted to MHC-side contacts."),
 "L_canon": ("log-odds", "Canonical-docking log odds-ratio of loop class (germline, CDR3) against target (MHC, peptide), Haldane-Anscombe corrected. High when CDR3 sits on the peptide and the germline loops on the helices."),
 "p_germ_mhc": ("fraction", "Share of germline (CDR1/CDR2) contacts that reach the MHC."),
 "p_cdr3_pep": ("fraction", "Share of CDR3 contacts that reach the peptide."),
 "pep_free_frac": ("fraction", "Share of the peptide the groove leaves for the receptor: mean over positions of n_TCR/(n_TCR + n_MHC). The threshold-free reading of 'peptide without its MHC anchors'."),
 "pep_cov_frac": ("fraction", "Peptide positions the TCR contacts, over peptide length."),
 "pep_cov_even": ("fraction", "Pielou evenness of the accessibility-discounted contact distribution, base ln(peptide length); how evenly the receptor uses the peptide it can reach."),
 "pep_cov_d2n": ("fraction", "Hill number of order 2 of that distribution over peptide length; the effective share of the peptide engaged."),
 "pep_cov_centre": ("fraction", "Contact-weighted mean position on [0, 1] from N- to C-terminus; 0.5 is centred."),
 "pep_cov_spread": ("fraction", "Contact-weighted standard deviation of that position, doubled; approaches 1 when the receptor reaches both termini."),
 "h0_pers_ent": ("fraction", "Normalized entropy of the H0 barcode of the contacted pMHC Calpha cloud. The bar lengths are the minimum spanning tree's edges, so no filtration is chosen."),
 # --- energetics ------------------------------------------------------------------------------
 "Phi_tcr_pep": ("log-odds", "Phi over TCR-peptide contacts under TCRen2, summed over all TCR regions. Lower is more favourable."),
 "Phi_tcr_mhc": ("log-odds", "Phi over TCR-MHC contacts under Miyazawa-Jernigan."),
 "Phi_pep_mhc": ("log-odds", "Phi over peptide-MHC contacts under Miyazawa-Jernigan. Computed without the receptor."),
 "Phi_cdr12": ("log-odds", "The CDR1 + CDR2 part of the TCR:peptide energy, both chains."),
 "Phi_cdr3a": ("log-odds", "The CDR3alpha part of the TCR:peptide energy."),
 "Phi_cdr3b": ("log-odds", "The CDR3beta part of the TCR:peptide energy."),
 "dPhi_tcr_pep": ("log-odds", "Poly-alanine reference delta of the TCR:peptide energy; the pose-geometry baseline removed."),
 "dPhi_pep_mhc": ("log-odds", "The same reference across peptide:MHC. Computed without the receptor."),
 "dPhi_pep_soft": ("log-odds", "Smoothed reference delta, peptide direction: the peptide's energy minus the free energy of the residue background at each peptide position, receptor frozen."),
 "varPhi_pep_soft": ("log-odds^2", "Variance of the local field under the background, peptide direction: how sharply each peptide position's energy responds to residue identity, summed over positions. NOT a ddG -- it is a second cumulant, not a difference of differences."),
 "dPhi_tcr_soft": ("log-odds", "Smoothed reference delta, receptor direction: the receptor's energy minus the free energy of the residue background at each contacted TCR position, peptide frozen."),
 "varPhi_tcr_soft": ("log-odds^2", "Variance of the local field under the background, receptor direction, summed over contacted TCR positions."),
 "dPhi_tra_soft": ("log-odds", "The alpha-chain part of the receptor-direction smoothed reference delta."),
 "dPhi_trb_soft": ("log-odds", "The beta-chain part of the receptor-direction smoothed reference delta."),
 "Phi_pep_int": ("log-odds", "The peptide's own intra-chain contact energy. Computed without the receptor."),
 # --- potts -----------------------------------------------------------------------------------
 "neg_energy": ("kT", "-E of the observed contact map under the coupled Potts model; higher is more native-like. Exactly log_z + log_lik."),
 "log_z": ("kT", "Log partition function over every contact map the geometry admits; the interface's capacity."),
 "log_lik": ("kT", "Log probability of the observed contact map; its typicality."),
 "psi": ("kT/site", "log_lik per available site, so interfaces of different size compare."),
 "n_contacts": ("count", "Available residue pairs that engaged. Distinct from the footprint's loop tally."),
 # --- kinetics --------------------------------------------------------------------------------
 "exp_lost": ("count", "Expected TCR:peptide contacts lost under a 1 A isotropic shift."),
 "mean_margin": ("A", "Mean contact margin, cutoff minus minimum heavy-atom distance."),
 "frac_robust": ("fraction", "Share of TCR:peptide contacts with at least 1 A of margin."),
 "n_spring": ("count", "Springs in the interface network; every other kinetics column is NaN below three."),
 "S_tot": ("N/m", "Trace of the stiffness tensor; total interface stiffness."),
 "K_tens": ("N/m", "Tensile stiffness along the docking axis."),
 "K_shear": ("N/m", "In-plane stiffness, S_tot minus K_tens."),
 "aniso": ("ratio", "K_shear / K_tens; how much stiffer the interface is along the pull than across it."),
 "lam_max": ("N/m", "Largest eigenvalue of the stiffness tensor."),
 "lam_min": ("N/m", "Smallest eigenvalue of the stiffness tensor."),
 "rupture_force": ("N", "Peak resisting force under steered separation along the weaker axis."),
 "rupture_work": ("J", "Force integrated to full separation; the off-rate proxy, and a geometry-only quantity no potential enters."),
 "couple_pep": ("count", "Peptide residues contacting both the MHC and the TCR."),
 "couple_mhc": ("count", "MHC residues contacting both the peptide and the TCR."),
 "couple_tcr": ("count", "TCR residues in the Valpha-Vbeta interface that also contact the pMHC."),
 "couple_total": ("count", "Sum of the three coupling counts."),
 "n_interface": ("count", "Interface residue count; the size denominator for the coupling counts."),
})

# The radius-tagged footprint columns are named from the `radii` argument, so their entries are
# generated the same way the columns are.
for _r in (7, 8):
    DETAIL.update({
        f"fp_b0_r{_r}": ("count",
            f"Betti-0 of the flag complex on the contacted pMHC Calpha atoms at {_r} A: how many "
            "disconnected patches the footprint falls into."),
        f"fp_b1_r{_r}": ("count",
            f"Betti-1 at {_r} A: how many holes the footprint encloses."),
        f"fp_chi_r{_r}": ("count", f"Euler characteristic b0 - b1 at {_r} A."),
        f"fp_b0_frac_r{_r}": ("fraction",
            f"Patches per contacted residue at {_r} A; the size-free form of Betti-0."),
    })


#: The invariance classes, in the order the catalogue reports them.
INVARIANCE_CLASSES: tuple[str, ...] = (
    "geometric", "topological", "compositional", "energetic", "categorical",
)

#: Descriptors that need a second look before they are used, and why. A name absent from here has
#: no known defect; presence is **not** a reason to drop the column, only to know what it is.
#:
#: Two flags:
#:
#: * ``"suspicious"`` -- the quantity is not measuring what its family name suggests. Either it
#:   reads the generator rather than the interface, or it is fixed by an exact identity over other
#:   columns, or it identifies the cohort rather than the complex.
#: * ``"stalled"`` -- the quantity is defined but does not move: near-zero spread, or undefined on
#:   most of the corpus, so nothing downstream can use it.
#:
#: The identities were each verified to float tolerance on both receptor benchmarks (max relative
#: difference 3.6e-15), so they are algebra rather than correlation. A determined column is exact
#: information the model already has -- harmless in a report, and a rank deficiency in a fit.
STATUS: dict[str, tuple[str, str]] = {
    "pitch": ("suspicious", "reads the generator's confidence rather than the interface: it is "
                            "docking_angles' incident_angle, and it out-discriminates every clean "
                            "docking angle for that reason. Never use it as a feature."),
    # determined by an exact identity over other emitted columns
    "fp_chi_r7": ("suspicious", "determined: chi = b0 - b1 at the same radius."),
    "fp_chi_r8": ("suspicious", "determined: chi = b0 - b1 at the same radius."),
    "D1_cell": ("suspicious", "determined: D1 = 12 ** H_cell."),
    "J_cell": ("suspicious", "determined: J = H_cell * ln 12 / ln S_cell."),
    "offset": ("suspicious", "determined: offset = hypot(shift_u, shift_w)."),
    "n_loop_contacts": ("suspicious", "determined: n_pep_contacts + n_mhc_contacts."),
    "neg_energy": ("suspicious", "determined: log_z + log_lik."),
    "S_tot": ("suspicious", "determined by K_tens and K_shear."),
    "aniso": ("suspicious", "determined by K_tens and K_shear."),
    "couple_total": ("suspicious", "determined: couple_pep + couple_mhc + couple_tcr."),
    "crossing": ("suspicious", "determined: abs(crossing_signed)."),
    # computed without the receptor: constant within an epitope-allele cohort
    "Phi_pep_mhc": ("suspicious", "no receptor: constant across every structure of one epitope on "
                                  "one allele, so a receptor-ranking model reading it reaches the "
                                  "cohort label without reading an interface."),
    "dPhi_pep_mhc": ("suspicious", "no receptor; see Phi_pep_mhc."),
    "Phi_pep_int": ("suspicious", "no receptor; see Phi_pep_mhc."),
    "n_pep_int": ("suspicious", "no receptor; see Phi_pep_mhc."),
    "mhc_class_bin": ("suspicious", "no receptor: it is the MHC class, I or II. It is also "
                                    "constant on any single-class cohort -- both receptor "
                                    "benchmarks are class I -- so it contributes nothing there "
                                    "and separates the classes everywhere else."),
    # measured on the modelled corpus rather than argued from the definition
    "ct_tp_salt_bridge": ("stalled", "only 3 distinct values over 1,707 modelled complexes: a "
                                     "salt bridge across the TCR:peptide interface is rare enough "
                                     "that the count is almost always 0. The TCR:MHC counterpart "
                                     "ct_tm_salt_bridge does move."),
}

FAMILIES = ("placement", "interface", "topology", "energetics", "potts", "kinetics")

#: Retired family names kept working for callers written before the 2026-08-24 split.
_FAMILY_ALIASES = {"geometry": ("placement", "interface"), "physics": ("energetics",)}


def descriptors(family: str | None = None, *, tcr_only: bool = False,
                invariance: str | None = None) -> tuple[str, ...]:
    """Descriptor names from :data:`DESCRIPTORS`, filtered by family and receptor involvement.

    Args:
        family: keep one of :data:`FAMILIES` (``"placement"``, ``"interface"``, ``"topology"``,
            ``"energetics"``, ``"kinetics"``), or all of them if ``None``. The retired names
            ``"geometry"`` (= ``placement`` + ``interface``) and ``"physics"`` (= ``energetics``)
            still work.
        tcr_only: keep only descriptors the receptor enters. Set this whenever the question being
            asked is about receptors — a peptide- or MHC-only column carries cohort identity.
        invariance: keep one class of :data:`INVARIANCE` — ``"geometric"`` for the docking's
            isometry invariants, ``"topological"`` for the interface surface's homeomorphism
            invariants, ``"compositional"`` for counts over the labelled contact set,
            ``"energetic"`` or ``"categorical"``. Combines with ``family``.

    Returns:
        The matching names, in catalogue order.

    Example:
        >>> descriptors("energetics", tcr_only=True)
        ('Phi_tcr_pep', 'Phi_tcr_mhc', 'Phi_cdr12', 'Phi_cdr3a', 'Phi_cdr3b', 'dPhi_tcr_pep')
        >>> descriptors("physics") == descriptors("energetics")   # retired alias
        True
        >>> descriptors("topology", invariance="topological")
        ('fp_b0_r7', 'fp_b1_r7', 'fp_chi_r7', 'fp_b0_frac_r7', 'fp_b0_r8', 'fp_b1_r8', \
'fp_chi_r8', 'fp_b0_frac_r8')
    """
    if invariance is not None and invariance not in INVARIANCE_CLASSES:
        raise ValueError(
            f"unknown invariance {invariance!r}; expected one of {INVARIANCE_CLASSES}"
        )
    if family is None:
        keep = set(FAMILIES) | {"score"}
    elif family in _FAMILY_ALIASES:
        keep = set(_FAMILY_ALIASES[family])
    elif family in FAMILIES:
        keep = {family}
    else:
        raise ValueError(f"unknown family {family!r}; expected one of {FAMILIES} "
                         f"or an alias {tuple(_FAMILY_ALIASES)}")
    return tuple(n for n, (fam, tcr) in DESCRIPTORS.items()
                 if fam in keep
                 and (tcr or not tcr_only)
                 and (invariance is None or INVARIANCE[n] == invariance))


def _extent(cm) -> float:
    """Distinct TCR residues contacting the pMHC (interface size); default TCR-region selection."""
    import polars as pl
    df = pl.concat([cm.interface("tcr_peptide"), cm.interface("tcr_mhc")])
    nodes = set()
    if df.height:
        for a, i in zip(df["chain.id.from"].to_list(), df["residue.index.from"].to_list()):
            nodes.add((a, i))
    return float(len(nodes))


def _chain_balance(cm) -> float:
    """min(a,b)/(a+b) over TCR:peptide contacts by TCR chain (0.5 = both chains equal, 0 = one only)."""
    tp = cm.interface("tcr_peptide", tcr_regions="all")
    if tp.height == 0:
        return math.nan
    a = b = 0
    for t in tp["chain.type.from"].to_list():
        a += t == "TRA"
        b += t == "TRB"
    return min(a, b) / (a + b) if (a + b) else math.nan


def _interface_symmetry(tp) -> dict[str, float]:
    """CDR3-dominance and TCR chain/loop imbalance from per-loop TCR:peptide contact **counts**.

    ``tp`` is the ``tcr_peptide`` interface table (``tcr_regions="all"``). Unlike ``e_cdr*`` (which are
    interface *energies*), these are pure contact-topology descriptors. Emitted as extra output columns
    (:data:`INTERFACE_SYMMETRY_FEATURES`), not part of :data:`RECOGNITION_FEATURES`.
    """
    import polars as pl
    reg, ch = pl.col("region.type.from"), pl.col("chain.type.from")
    h = lambda f: float(tp.filter(f).height)  # noqa: E731
    n12 = h(reg.is_in(["CDR1", "CDR2"]))                                   # germline CDR1/2 (both chains)
    n3a, n3b = h((reg == "CDR3") & (ch == "TRA")), h((reg == "CDR3") & (ch == "TRB"))
    nA = h(reg.is_in(["CDR1", "CDR2", "CDR3"]) & (ch == "TRA"))            # whole alpha CDRs
    nB = h(reg.is_in(["CDR1", "CDR2", "CDR3"]) & (ch == "TRB"))            # whole beta CDRs
    tot = n12 + n3a + n3b
    return {
        # CDR3 (a+b) share of CDR TCR:peptide contacts -- higher = CDR3-dominated (binder-like; oriented +)
        "cdr3_dominance": (n3a + n3b) / tot if tot else math.nan,
        # |CDR3a - CDR3b| normalised imbalance -- absolute magnitude (direction is tested, not assumed)
        "cdr3_ab_imbalance": abs(n3a - n3b) / (n3a + n3b) if (n3a + n3b) else math.nan,
        # |alpha - beta| whole-CDR contact imbalance, normalised -- absolute magnitude
        "chain_cdr_imbalance": abs(nA - nB) / (nA + nB) if (nA + nB) else math.nan,
    }


def _burial(structure, tcr_ids, pmhc_ids) -> float:
    """Interface ΔSASA = SASA(TCR alone) + SASA(pMHC alone) − SASA(complex) via biopython ShrakeRupley
    (``n_points=100``), reproducing the training-time ``burial``. ΔSASA is an interface quantity, so the
    distal TCR constant domain cancels; computed on a temp PDB of the typed chains."""
    if not tcr_ids or not pmhc_ids:
        return math.nan
    import os
    import tempfile
    from copy import deepcopy

    from Bio.PDB import PDBParser
    from Bio.PDB.Model import Model
    from Bio.PDB.SASA import ShrakeRupley
    from Bio.PDB.Structure import Structure as BioStructure

    from .structure.io import write_pdb

    with tempfile.TemporaryDirectory() as td:
        path = os.path.join(td, "complex.pdb")
        write_pdb(structure, path)
        model = PDBParser(QUIET=True).get_structure("x", path)[0]      # parsed fully into memory
    sr = ShrakeRupley(n_points=100)

    def sasa_of(ids):
        m2 = Model(0)
        for ch in model:
            if ch.id in ids:
                m2.add(deepcopy(ch))
        s2 = BioStructure("t")
        s2.add(m2)
        sr.compute(s2, level="A")
        return sum(a.sasa for ch in m2 for res in ch if res.id[0] == " " for a in res.get_atoms())

    both = set(tcr_ids) | set(pmhc_ids)
    return float((sasa_of(set(tcr_ids)) + sasa_of(set(pmhc_ids))) - sasa_of(both))


def _cdr3_frame_features(structure) -> dict[str, float]:
    """The 18 CDR3-local frame descriptors (:data:`CDR3_FRAME_FEATURES`) for a chain-typed structure.

    Both CDR3 loops are projected onto the pMHC groove frame (see :data:`CDR3_FRAME_FEATURES`). The
    structure must already be chain-typed (``classify_chains``) so its CDR3 regions are populated.
    Undefined terms (no groove frame, missing peptide, or a loop with < 3 Cα) are ``NaN``.
    """
    from .orient.docking import _chain_ca, _groove_frame

    out = {k: math.nan for k in CDR3_FRAME_FEATURES}
    try:
        u, w, n = _groove_frame(structure)
    except Exception:
        return out
    pep = _chain_ca(structure, ("PEPTIDE",))
    if len(pep) < 2:
        return out
    origin = pep.mean(axis=0)
    basis = np.stack([u, w, n])                                        # rows = groove basis
    for loop, ctype in (("cdr3a", "TRA"), ("cdr3b", "TRB")):
        cas = None
        for c in structure.chains:
            if c.chain_type != ctype:
                continue
            for reg in getattr(c, "regions", []) or []:
                if reg.region_type == "CDR3":
                    pts = [r.ca for r in reg.residues if r.ca is not None]
                    if len(pts) >= 3:
                        cas = np.asarray(pts)
                    break
        if cas is None:
            continue
        d = cas.mean(axis=0) - origin
        reach = float(np.linalg.norm(d))
        off = basis @ (d / (reach + 1e-9))
        av = cas[-1] - cas[0]
        ax = basis @ (av / (np.linalg.norm(av) + 1e-9))
        topep = float(np.linalg.norm(cas[:, None, :] - pep[None, :, :], axis=2).min())
        ext = float(np.linalg.norm(cas[-1] - cas[0]))
        for k, v in zip(_CDR3_FRAME_KEYS, (reach, *off, *ax, topep, ext)):
            out[f"{loop}_{k}"] = float(v)
    return out


def recognition_features(source, *, organism: str = "human", potential=None,
                         full: bool = False, annotate: bool = True) -> dict[str, float]:
    """Extract the core recognition vector from a TCR–pMHC structure (path or parsed).

    Returns a dict keyed by :data:`RECOGNITION_FEATURES` (degenerate/undefined terms are ``NaN``):
    docking geometry, per-interface energies (raw ``F`` and poly-alanine ``ΔF``), contact-type
    tallies, interface ΔSASA ``burial``, and the ``mhc_class_bin`` indicator. The structure is
    chain-typed and MHC-annotated in place. :data:`DESCRIPTORS` gives each column's family and
    whether the receptor enters its definition.

    With ``full=True`` the row is extended with the 18 CDR3-frame descriptors
    (:data:`CDR3_FRAME_FEATURES`) — the complete :data:`FULL_FEATURES` vector.
    """
    import polars as pl

    from .annotation import classify_chains
    from .contact_types import contact_type_counts
    from .contactmap import ContactMap
    from .ddg import reference_delta
    from .mhc import annotate_mhc
    from .oracle import _native_peptide
    from .orient.docking import docking_angles
    from .orient.tcrdock_geometry import docking_geometry
    from .pipeline import _interface_energy
    from .potential import mj as _mj
    from .potential import tcren2 as _tcren2
    from .structure import Structure, import_structure

    s = source if isinstance(source, Structure) else import_structure(source)
    if annotate:                                                      # skip if pre-annotated (batch path)
        if all(c.chain_type is None for c in s.chains):
            classify_chains(s, organism=organism, autodetect_species=True)
        annotate_mhc(s)
    tcren_pot = potential or _tcren2()
    mj_pot = _mj()

    cm = ContactMap.from_structure(s)
    native = _native_peptide(s)
    row = {k: math.nan for k in (FULL_FEATURES if full else RECOGNITION_FEATURES)}

    try:                                                              # geometry (docking)
        da = docking_angles(s)
        row["pitch"], row["crossing"] = float(da.incident_angle), float(da.crossing_angle)
        row["crossing_signed"] = float(da.crossing_angle_signed)
    except Exception:
        pass
    try:
        dg = docking_geometry(s)                                     # native TCRdock rigid-body params
        row.update(dock_d=float(dg.d), dock_torsion=float(dg.torsion),
                   dock_tcr_uy=float(dg.tcr_unit_y), dock_tcr_uz=float(dg.tcr_unit_z),
                   dock_mhc_uy=float(dg.mhc_unit_y), dock_mhc_uz=float(dg.mhc_unit_z))
    except (ValueError, KeyError, AttributeError) as exc:
        # Six features (dock_d, dock_torsion, dock_{tcr,mhc}_u{y,z}) stay NaN here. Say why: a bare
        # `except Exception: pass` made an unannotated or unsupported groove indistinguishable from
        # a computed result, and the feature dict has no room for a reason column.
        warnings.warn(f"{s.pdb_id}: docking geometry unavailable ({exc}); "
                      f"the six dock_* features are NaN", RuntimeWarning, stacklevel=2)

    tm = cm.interface("tcr_mhc", tcr_regions="all")                  # interface energetics
    row["Phi_tcr_mhc"] = float(_interface_energy(tm, mj_pot))
    tp = cm.interface("tcr_peptide", tcr_regions="all")
    reg, ch = pl.col("region.type.from"), pl.col("chain.type.from")
    row["Phi_tcr_pep"] = float(_interface_energy(tp, tcren_pot))
    row["Phi_pep_mhc"] = float(_interface_energy(cm.interface("peptide_mhc"), mj_pot))
    row["Phi_cdr12"] = float(_interface_energy(tp.filter(reg.is_in(["CDR1", "CDR2"])), tcren_pot))
    row["Phi_cdr3a"] = float(_interface_energy(tp.filter((reg == "CDR3") & (ch == "TRA")), tcren_pot))
    row["Phi_cdr3b"] = float(_interface_energy(tp.filter((reg == "CDR3") & (ch == "TRB")), tcren_pot))
    if native:
        try:
            row["dPhi_tcr_pep"] = float(reference_delta(cm, native, tcren_pot, interface="tcr_peptide"))
        except Exception:
            pass
        try:
            row["dPhi_pep_mhc"] = float(reference_delta(cm, native, mj_pot, interface="peptide_mhc"))
        except Exception:
            pass

    # The smoothed counterpart of dPhi: the same first difference in sequence, but against the free
    # energy of the residue BACKGROUND rather than of poly-alanine, plus varPhi, the variance of the
    # same local field. varPhi is a second CUMULANT, not a difference of differences -- the only
    # `dd` quantity in the package is ddG, the change in binding free energy upon mutation.
    # Emitted in both directions, because they answer different questions -- the peptide direction
    # varies the peptide with the receptor frozen (the peptide scan), the TCR direction varies the
    # receptor with the peptide frozen (the TCR scan) -- and the TCR direction is split by chain,
    # since a linear model given the two chains apart can form any contrast between them, while one
    # given only their sum cannot. See ddg.smoothed_reference.
    from .ddg import smoothed_reference
    for name, kw in (("pep", {"side": "peptide"}),
                     ("tcr", {"side": "tcr"}),
                     ("tra", {"side": "tcr", "chain": "TRA"}),
                     ("trb", {"side": "tcr", "chain": "TRB"})):
        try:
            sm = smoothed_reference(cm, tcren_pot, **kw)
        except Exception:
            continue
        row[f"dPhi_{name}_soft"] = sm["dPhi"]
        if name in ("pep", "tcr"):
            row[f"varPhi_{name}_soft"] = sm["varPhi"]

    row["extent"] = _extent(cm)                                      # coverage
    row["chain_balance"] = _chain_balance(cm)
    row["n_contacts_tp"] = float(tp.height)
    row["n_pep_contacted"] = float(tp.select("residue.index.to").unique().height if tp.height else 0)
    row["n_contacts_tm"] = float(tm.height)

    # scheme="v1" is pinned, not defaulted: the frozen classifiers below were fitted on these
    # counts, and the current typing (tcren.contact_types "v2") gives different ones.
    ctp = contact_type_counts(cm, "tcr_peptide", scheme="v1")        # contact types
    ctm = contact_type_counts(cm, "tcr_mhc", scheme="v1")
    for t in _CT_TYPES:
        if t != "hydrogen_bond":              # emitted once, as n_hbond (the name Eq. Q uses)
            row[f"ct_tp_{t}"] = float(ctp[f"pairs_{t}"])
        row[f"ct_tm_{t}"] = float(ctm[f"pairs_{t}"])
    row["n_hbond"] = float(ctp["pairs_hydrogen_bond"])

    tcr_ids = [c.chain_id for c in s.chains if c.chain_type in _TCR_TYPES]
    pmhc_ids = [c.chain_id for c in s.chains if c.chain_type is not None and c.chain_type not in _TCR_TYPES]
    row["burial"] = _burial(s, tcr_ids, pmhc_ids)
    row["mhc_class_bin"] = 1.0 if any(getattr(c, "chain_supertype", None) == "MHCII"
                                      for c in s.chains) else 0.0

    if full:                                                          # FramePose CDR3 layer
        row.update(_cdr3_frame_features(s))
    return row


def _stability_clash_columns(s) -> dict[str, float]:
    """Interface steric-clash + TCR:peptide contact-stability descriptors for the recognize table.

    Extra *output* columns, **not** part of :data:`RECOGNITION_FEATURES` or any fitted model: a
    coordinate-only read of forced-pose quality --- steric-clash burden (:func:`tcren.interface_clashes`)
    and contact fragility (:func:`tcren.contact_stability`). NaN where the structure lacks a peptide or
    receptor chain.
    """
    from .clashes import interface_clashes
    from .stability import contact_stability

    out: dict[str, float] = {}
    try:
        cl = interface_clashes(s)
        out["n_clashes"], out["clash_score"] = float(cl.n_clashes), float(cl.clash_score)
    except Exception:  # noqa: BLE001 - no peptide chain etc.
        out["n_clashes"] = out["clash_score"] = math.nan
    try:
        st = contact_stability(s)
        out["exp_lost"] = float(st.exp_lost)
        out["mean_margin"] = float(st.mean_margin)
        out["frac_robust"] = float(st.frac_robust)
    except Exception:  # noqa: BLE001 - no peptide/receptor chain etc.
        out["exp_lost"] = out["mean_margin"] = out["frac_robust"] = math.nan
    return out


def _symmetry_columns(s) -> dict[str, float]:
    """Interface-symmetry extra output columns (:data:`INTERFACE_SYMMETRY_FEATURES`) for the recognize
    table --- CDR3-dominance and α/β contact imbalance from a fresh contact map. NaN on failure."""
    from .contactmap import ContactMap
    try:
        cm = ContactMap.from_structure(s)
        return _interface_symmetry(cm.interface("tcr_peptide", tcr_regions="all"))
    except Exception:  # noqa: BLE001 - no peptide/receptor chain etc.
        return {k: math.nan for k in INTERFACE_SYMMETRY_FEATURES}


def _placement_columns(s) -> dict[str, float]:
    """Receptor-placement extra output columns (:data:`TCR_PLACEMENT_FEATURES`) — the ride height,
    in-plane shift and offset of the CDR centroid over the groove. NaN where the frame is undefined."""
    from .orient.docking import tcr_placement
    try:
        tp = tcr_placement(s)
        return {k: float(getattr(tp, k)) for k in TCR_PLACEMENT_FEATURES}
    except Exception:  # noqa: BLE001 - no groove frame / no receptor chain
        return dict.fromkeys(TCR_PLACEMENT_FEATURES, math.nan)


def _footprint_columns(s, radii=(7.0, 8.0)) -> dict[str, float]:
    """Footprint shape extra output columns — the ``topology`` family (:mod:`tcren.footprint`).

    Rigid-motion invariant, so the structure does not need orienting; it needs only chain typing and
    CDR markup, which the batch annotation has already done by the time this runs."""
    from .footprint import footprint_features
    try:
        return footprint_features(s, radii=radii)
    except Exception:  # noqa: BLE001 - no peptide/receptor chain etc.
        return dict.fromkeys(footprint_topology_features(radii) + FOOTPRINT_SIZE_FEATURES, math.nan)


def _peptide_internal_columns(s) -> dict[str, float]:
    """Intra-peptide extra output columns (:data:`PEPTIDE_INTERNAL_FEATURES`) for the recognize table.

    The peptide's MJ contact energy with itself and the number of those contacts — the term the three
    interface energies leave out. NaN/0 where the structure has no peptide chain."""
    from .contactmap import ContactMap
    from .potential import mj as _mj
    from .scoring import intra_peptide_energy
    try:
        cm = ContactMap.from_structure(s, peptide_internal=True)
        return {"Phi_pep_int": float(intra_peptide_energy(cm, _mj())),
                "n_pep_int": float(cm.peptide_internal.height)}
    except Exception:  # noqa: BLE001 - no peptide chain etc.
        return {k: math.nan for k in PEPTIDE_INTERNAL_FEATURES}


def recognition_table(items, *, organism: str = "human", full: bool = False,
                      threads: int = 1, chunk: int = 64,
                      autodetect_species: bool = True, mechanics: bool = False,
                      include: Sequence[str] | None = None, radii: Sequence[float] = (7.0, 8.0),
                      _mmseqs_threads: int = 0) -> list[dict]:
    """Batched feature (+score) extraction for a whole set of TCR–pMHC structures.

    ``items`` is an iterable of ``(id, structure-or-path)``. The set is annotated with a **single**
    arda call per organism (:func:`tcren.paper.helpers._batch_annotate`) and a **single** mmseqs MHC
    search (:func:`tcren.mhc.annotate_mhc_batch`) — the dataset-scale path that avoids the per-structure
    annotation cost — then :func:`recognition_features` (``full=``) is extracted for each. This emits
    **descriptors only**: the fitted composites and cohort-relative scores that used to ride along
    here were removed in 2.26.0, and scoring is :func:`tcren.reliability.s_free` on the table.
    ``full`` also appends the
    intra-peptide columns :data:`PEPTIDE_INTERNAL_FEATURES` (``Phi_pep_int``, ``n_pep_int``) — the
    peptide's contact energy with itself, which the interface energies omit. Returns one row dict per
    structure (``complex.id`` + features); a structure that fails yields
    ``{"complex.id": id, "error": ...}`` so the batch stays resilient.

    The two stages run **in sequence** and never compete for the machine.

    *Search* is one arda call per organism plus one mmseqs MHC search, each given every core, over
    the whole set. *Featurisation* is where the time actually goes — a 100-pose probe spends 96 s
    there against 2.4 s of arda and 0.9 s of MHC search — and it is pure Python/numpy, so
    ``threads`` > 1 runs it in that many **worker processes**. The flag keeps its name for
    compatibility; it has always meant "how much of this machine may I use".

    It used to mean concurrent *threads* over ``chunk``-sized batches, which was the wrong shape
    twice over: the GIL serialised the 94 % of the work that dominates, and each batch spawned its
    own mmseqs, so N batches asked for N x cores. Sharding the same work across independent
    subprocesses was measured 8x faster, which is what this now does directly.

    ``chunk`` is retained for signature compatibility and is no longer used.

    ``autodetect_species`` searches ``organism`` **and** mouse so a mis-declared cohort is still
    typed correctly. That doubles the annotation cost, so pass ``False`` when the organism is known
    — it halves the mmseqs work and changes nothing else.

    ``mechanics`` appends the :mod:`tcren.mechanics` koff proxies (``n_spring``, ``S_tot``,
    ``K_tens``, ``K_shear``, ``aniso``, ``rupture_force``, ``rupture_work``, ``couple_*``) to the
    same rows. They need the same annotated structure the descriptors do, so computing them here
    costs only their own arithmetic — running ``tcren mechanics`` separately repeats the whole
    parse and both mmseqs searches, and returns a second table keyed differently.
    """
    import os as _os

    from .annotation import classify_chains
    from .annotation.arda_adapter import _import_arda
    from .mhc import annotate_mhc_batch
    from .paper.helpers import _batch_annotate
    from .structure import Structure, import_structure

    items = list(items)
    ids, structs, rows = [], [], []
    for id_, src in items:
        try:
            structs.append(src if isinstance(src, Structure) else import_structure(src))
            ids.append(id_)
        except Exception as exc:  # noqa: BLE001
            rows.append({"complex.id": id_, "error": f"{type(exc).__name__}: {str(exc)[:80]}"})

    if structs:                       # stage 1: one arda call per organism + one MHC search, all cores
        cores = _mmseqs_threads or (_os.cpu_count() or 1)
        orgs = (organism, "mouse") if autodetect_species else (organism,)
        recs = _batch_annotate(structs, _import_arda(), organisms=orgs, threads=cores)
        for i, s in enumerate(structs):
            try:
                classify_chains(s, organism=organism, autodetect_species=autodetect_species,
                                precomputed_records=recs[i])
            except Exception:  # noqa: BLE001 - MHC-only / unannotatable chains stay unset
                pass
        annotate_mhc_batch(structs, threads=cores)

    # stage 2: featurisation, the part that actually costs (94 % of wall time on a 100-pose probe:
    # 96 s against 2.4 s of arda and 0.9 s of MHC search). It is pure Python/numpy, so processes.
    work = [(id_, s, organism, full, mechanics, include, tuple(radii))
            for id_, s in zip(ids, structs)]
    if threads > 1 and len(work) > 1:
        from concurrent.futures import ProcessPoolExecutor
        with ProcessPoolExecutor(max_workers=min(threads, len(work))) as ex:
            rows.extend(ex.map(_featurise_one, work, chunksize=max(1, len(work) // (threads * 4))))
    else:
        rows.extend(_featurise_one(w) for w in work)

    return rows


def _featurise_one(args) -> dict:
    """One structure -> one row. Module-level and self-contained so it pickles to a worker process.

    The structure arrives already annotated: chain typing and the MHC call are batch operations and
    belong to the single search in :func:`recognition_table`, not to a per-structure worker.
    """
    id_, s, organism, full, mechanics, include, radii = args
    if include is not None:
        return _featurise_families(id_, s, organism, include, radii)
    try:
        feats = recognition_features(s, organism=organism, full=full, annotate=False)
        row = {"complex.id": id_, **feats, **_stability_clash_columns(s), **_symmetry_columns(s)}
        if full:                              # the intra-peptide term costs a second contact map
            row.update(_peptide_internal_columns(s))
        if mechanics:
            from .mechanics import interface_mechanics
            row.update(interface_mechanics(s))
        return row
    except Exception as exc:  # noqa: BLE001
        return {"complex.id": id_, "error": f"{type(exc).__name__}: {str(exc)[:80]}"}


def _featurise_families(id_, s, organism: str, include, radii) -> dict:
    """One structure -> one row holding exactly the catalogued descriptors of the requested families.

    Only what is asked for is computed: ``tcren features -i topology`` never builds the energies, and
    ``-i placement`` never runs the spring network. The returned row is filtered against
    :data:`DESCRIPTORS`, so a column exists in the output if and only if the catalogue names it —
    which is what makes the families a partition of the feature table rather than a label on it.
    """
    want = set(include)
    unknown = want - set(FAMILIES)
    if unknown:
        raise ValueError(f"unknown feature families {sorted(unknown)}; expected {FAMILIES}")
    row: dict[str, float] = {}
    try:
        if want & {"placement", "interface", "energetics"}:
            row.update(recognition_features(s, organism=organism, full=True, annotate=False))
            row.update(_symmetry_columns(s), **_peptide_internal_columns(s))
        if want & {"interface", "kinetics"}:                 # clash + contact fragility share a pass
            row.update(_stability_clash_columns(s))
        if "placement" in want:
            row.update(_placement_columns(s))
        if "topology" in want:
            row.update(_footprint_columns(s, radii))
        if "potts" in want:
            from .potts import score_structure
            row.update({k: v for k, v in score_structure(s).items() if k != "pdb.id"})
        if "kinetics" in want:
            from .mechanics import interface_mechanics
            row.update(interface_mechanics(s))
    except Exception as exc:  # noqa: BLE001 - keep the batch alive, one bad structure is one bad row
        return {"complex.id": id_, "error": f"{type(exc).__name__}: {str(exc)[:80]}"}
    keep = {n for n, (fam, _) in DESCRIPTORS.items() if fam in want}
    keep |= {f"fp_{k}_r{r:g}" for r in radii for k in ("b0", "b1", "chi", "b0_frac")} if "topology" in want else set()
    return {"complex.id": id_, **{k: v for k, v in row.items() if k in keep}}


