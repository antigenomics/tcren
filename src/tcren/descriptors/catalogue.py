"""The descriptor catalogue: what every emitted column is, and how to select a subset of them.

**Data and selection only.** Nothing here computes a descriptor, and nothing here imports a module
that does -- the one exception is the two tuples of *names* it splats into :data:`DESCRIPTORS`, which
are strings, not arithmetic. That separation is the point: a caller that only wants to know what a
column means, what it is invariant under or which family it belongs to pays for none of the
structure parsing, contact building or potential loading that computing one would cost.

The layers above this one are :mod:`tcren.descriptors.compute`, which turns a structure into the
values, and :mod:`tcren.descriptors.table`, which runs that over a set. :mod:`tcren.recognition` is
kept as the public name and re-exports all three.
"""
from __future__ import annotations

# The only import in this module, and it is a tuple of *names*: `footprint_topology_features()`
# returns the topology family's column names so `DESCRIPTORS` can splat them, which keeps the
# catalogue and the module that produces those columns from drifting apart. No arithmetic crosses.
from ..topology.footprint import FOOTPRINT_SIZE_FEATURES, footprint_topology_features
from ..topology.literature import LITERATURE_FEATURES

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

#: Where the receptor body sits over the groove (:func:`tcren.docking.tcr_placement`), emitted as extra
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
    # where the receptor body sits over the groove (`tcren.docking.tcr_placement`)
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
    # -- topology: published interface descriptors (`tcren.topology.literature`) -----------------
    # The surface pair reads a channel nothing else here reaches: both faces rasterised as height
    # fields on one shared groove-frame grid, so the gap is their difference cell by cell.
    **{f: ("topology", True) for f in LITERATURE_FEATURES},
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
    # -- the contact map read as a bipartite graph (`tcren.interface_graph`) ----------------------
    # Topological: each is an invariant of the abstract contact graph, unchanged by any deformation
    # that preserves which residue touches which. The component and cycle fractions are literally b0
    # and b1 of the 1-complex, size-normalized; the degree sequence and its correlation are graph
    # invariants in the same sense. None of them reads an Angstrom.
    "g_even_tcr": "topological",
    "g_even_pmhc": "topological",
    "g_comp_frac": "topological",
    "g_alg_conn": "topological",
    "g_cyclo_frac": "topological",
    "g_assort": "topological",
    # Compositional: these read the CDR-loop *labelling* on top of the graph -- which loop reached
    # which partner -- so they are counts over a labelled set, not invariants of the bare complex.
    "g_loop_even": "compositional",
    "g_loop_overlap": "compositional",
    "degree_evenness_tp": "compositional",
    "frac_well_coordinated_tp": "compositional",
    # -- the two surfaces read as height fields, and two graph functionals ------------------------
    # Geometric: every surface quantity is built from Angstrom heights on a metric grid, so a
    # deformation that preserves which residue touches which still moves them. That is the point --
    # the gap is the one channel here that measures space rather than incidence.
    "sc_shape": "geometric",
    "sc_gap_mean": "geometric",
    "sc_gap_sd": "geometric",
    "sc_gap_vol": "geometric",
    "sc_interlock": "geometric",
    "sc_gap_index": "geometric",
    "sc_interlock_frac": "geometric",
    "sc_gap_depth": "geometric",
    "sc_gap_height": "geometric",
    "sc_gap_asym": "geometric",
    "sc_dh": "geometric",
    "sc_cells": "geometric",
    "sc_coverage": "geometric",
    # Compositional: these read the amino-acid labelling painted on the same grid, not its shape.
    "sc_charge": "compositional",
    "sc_phobic": "compositional",
    "sc_charge_prod": "compositional",
    "sc_phobic_prod": "compositional",
    "sc_dcharge": "compositional",
    "sc_dphobic": "compositional",
    # Compositional: contact order reads target SEQUENCE positions and the participation
    # coefficient reads which module an edge lands in -- both are labellings on the contact graph,
    # not invariants of the bare complex.
    "co_pep": "compositional",
    "co_mhc": "compositional",
    "partcoef_tcr": "compositional",
    "partcoef_pmhc": "compositional",
    # -- the Calpha / Cbeta maps read as matrices -------------------------------------------------
    # Geometric, not topological, and for the reason `h0_pers_ent` is: the Gaussian kernel is built
    # from Angstrom distances, so its singular values move under a deformation that leaves the
    # contact set alone. A spectral shape descriptor of a metric object is a metric quantity.
    "m_erank_tp": "geometric",
    "m_gap_tp": "geometric",
    "m_erank_tm": "geometric",
    "m_gap_tm": "geometric",
    "m_face_tp": "geometric",
    "m_face_tm": "geometric",
    # The rank correlation between the two maps is a shape statistic over a sample of residue
    # pairs, not a length -- but it is a correlation between two metric quantities and moves under
    # a deformation that leaves the contact set alone, so it files with them, not with the graph.
    "ca_cb_agreement_tp": "geometric",
    "ca_cb_agreement_tm": "geometric",
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
    # -- the contact map as a bipartite graph -----------------------------------------------------
    "g_even_tcr": ("fraction",
        "Pielou evenness of the contact degrees of the engaged CDR-loop residues, base the engaged "
        "count. 1 when every engaged residue carries the same number of partners."),
    "g_even_pmhc": ("fraction",
        "Pielou evenness of the contact degrees of the engaged pMHC residues, base the engaged count."),
    "g_comp_frac": ("fraction",
        "Connected components of the bipartite contact graph per node; the parameter-free form of "
        "the footprint patch count, needing no Calpha radius."),
    "g_alg_conn": ("ratio",
        "Second-smallest eigenvalue of the normalised Laplacian on the largest contact-graph "
        "component, in [0, 2]. Near 0 when the footprint is about to fall into two patches."),
    "g_cyclo_frac": ("fraction",
        "Contacts beyond a spanning forest over all contacts, (E - V + C) / E, which is the contact "
        "graph's first Betti number made size-free. High when the footprint is interlocked."),
    "g_loop_even": ("fraction",
        "Pielou evenness over the six CDR loops of the number of distinct pMHC residues each "
        "reaches, base 6. Counts partners rather than contacts, so it does not track residue size."),
    "g_loop_overlap": ("fraction",
        "Mean pairwise Jaccard overlap of the engaged CDR loops' pMHC partner sets. High when the "
        "loops crowd onto the same residues instead of partitioning the surface."),
    "g_assort": ("ratio",
        "Degree assortativity of the contact graph: the correlation, over contacts, between the "
        "degrees of the two residues involved."),
    "degree_evenness_tp": ("fraction",
        "Participation ratio of the receptor-side contact degrees across TCR:peptide, in [0, 1]. "
        "Low when a few over-reaching side chains hoard the contact budget."),
    "frac_well_coordinated_tp": ("fraction",
        "Share of contacting receptor residues reaching no more than three peptide residues, the "
        "count a crystal side chain typically makes."),
    # -- the Calpha / Cbeta maps as matrices ------------------------------------------------------
    "m_erank_tp": ("fraction",
        "Effective rank of the CDR-loop x peptide Calpha proximity kernel over its maximum: how "
        "many independent approach modes the interface has. Peptide-length coupled (Spearman "
        "-0.547 on 148 class I crystals) because a longer class I peptide bulges."),
    "m_gap_tp": ("ratio",
        "Second over first singular value of that kernel. Near 0 when the approach is separable "
        "into a loop profile times a peptide profile rather than pairing specific residues."),
    "m_erank_tm": ("fraction",
        "Effective rank fraction of the CDR-loop x MHC-helix kernel. CDR3-length coupled "
        "(Spearman -0.437), so it reads how much loop there is to spread over the helices."),
    "m_gap_tm": ("ratio",
        "Second over first singular value of the CDR-loop x MHC-helix kernel. The least "
        "length-coupled column in the catalogue: -0.012 against CDR3 length, +0.027 against "
        "peptide length, 99.9 per cent of its variance surviving both."),
    "m_face_tp": ("A",
        "Mean Calpha-Calpha minus Cbeta-Cbeta distance over the contacting TCR:peptide residue "
        "pairs. Positive when side chains lean towards each other, negative when the backbones are "
        "the close part and the side chains point away."),
    "m_face_tm": ("A",
        "Mean Calpha-Calpha minus Cbeta-Cbeta distance over the contacting TCR:MHC residue pairs."),
    "ca_cb_agreement_tp": ("ratio",
        "Spearman correlation between the Calpha and Cbeta distance maps over the TCR:peptide "
        "approach shell. High when the side chains track the backbone, as they do in a crystal."),
    "ca_cb_agreement_tm": ("ratio",
        "The same rank correlation across the TCR:MHC approach shell."),
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
 # --- published interface descriptors, `tcren.topology.literature` -----------------------------
 "sc_shape": ("ratio", "Pearson r between the pMHC and TCR height fields over the shared grid; positive is complementary, the receptor riding up where the groove rises. Lawrence & Colman's Sc is the same idea on a dot surface."),
 "sc_charge": ("ratio", "Pearson r between the two charge fields; NEGATIVE is complementary, plus meeting minus."),
 "sc_phobic": ("ratio", "Pearson r between the two Kyte-Doolittle fields; positive is complementary, apolar meeting apolar."),
 "sc_charge_prod": ("ratio", "Mean per-cell product of the two charge fields."),
 "sc_phobic_prod": ("ratio", "Mean per-cell product of the two hydropathy fields."),
 "sc_gap_mean": ("A", "Mean of h(TCR) - h(pMHC) over retained cells. Negative on a real interface: the median cell interdigitates."),
 "sc_gap_sd": ("A", "Spread of the same gap. High when the receptor rests on a few high points rather than meshing."),
 "sc_gap_vol": ("A^3", "Void volume, the gap integrated over the contact plane where it is positive."),
 "sc_interlock": ("A^3", "Interdigitated volume, the gap integrated where it is negative. The larger of the two on a real interface."),
 "sc_gap_index": ("A", "Void volume over retained contact area; the intensive form of the gap-volume channel."),
 "sc_interlock_frac": ("fraction", "Share of retained cells whose gap is negative; the per-structure form of the corpus 71% interdigitation."),
 "sc_gap_depth": ("A", "Mean depth over the interlocked cells alone: how far the receptor reaches in where it does."),
 "sc_gap_height": ("A", "Mean standoff over the void cells alone: how high it stands where it does not mesh."),
 "sc_gap_asym": ("signed fraction", "(void - interlock) / (void + interlock); -1 for a face that only interlocks, +1 for one that only stands off."),
 "sc_dh": ("A", "Mean absolute per-cell height difference between the two faces."),
 "sc_dcharge": ("ratio", "Mean absolute per-cell charge difference between the two faces."),
 "sc_dphobic": ("ratio", "Mean absolute per-cell hydropathy difference between the two faces."),
 "sc_cells": ("count", "Grid cells entering the comparison; bookkeeping, so a low complementarity can be told from a thin one."),
 "sc_coverage": ("fraction", "Retained cells as a share of the occupied pMHC cells in the window; bookkeeping."),
 "co_pep": ("ratio", "Contact order on the peptide: mean sequence separation of the peptide residues one CDR loop reaches, averaged over loops and divided by the peptide's span."),
 "co_mhc": ("ratio", "Contact order on the MHC helices, by the same construction."),
 "partcoef_tcr": ("fraction", "Mean over engaged TCR residues of 1 - sum_s (k_s/k)^2 with the modules peptide and MHC; 0 when every residue reads one target only."),
 "partcoef_pmhc": ("fraction", "The same over engaged pMHC residues with the six CDR loops as modules."),
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
    # -- length coupling, measured on 143 class I Native2026 crystals ---------------------------
    # Two axes, because they catch different columns: peptide length 8-13 and CDR3(alpha+beta)
    # length 19-28. The mechanism is the author's, 2026-09-01. On the peptide axis: class I closes
    # its groove at both ends, so a longer peptide must BULGE, and a bulged peptide sits closer to
    # the receptor and takes more contacts per residue. On the CDR3 axis: a longer loop simply has
    # more residues to spread over the surface it reaches.
    #
    # None of this is a reason to drop a column, and the entries below are not a deprecation. Every
    # one keeps most of its variance after BOTH lengths are regressed out (58 per cent in the worst
    # case), so each carries something length does not; and a length-coupled column is useful
    # BESIDE a length-free one, because a PCA or a linear model can form the contrast that cancels
    # the shared part. What these entries buy is knowing which is which. Reported as
    # (rho against peptide length, rho against CDR3 length, variance surviving both).
    "m_erank_tp": ("suspicious", "peptide-length coupled: -0.547 / -0.105, 58.3 per cent of "
                                 "variance beyond both lengths. It reads the class I bulge."),
    "m_erank_tm": ("suspicious", "CDR3-length coupled: -0.186 / -0.437, 69.4 per cent beyond "
                                 "both. It reads how much loop there is to spread over the helices."),
    "m_gap_tp": ("suspicious", "peptide-length coupled: -0.322 / +0.041, 90.8 per cent beyond both."),
    "ca_cb_agreement_tm": ("suspicious", "coupled to both lengths: -0.349 / -0.334, 70.5 per cent "
                                         "beyond them -- the most length-loaded column here."),
    "degree_evenness_tp": ("suspicious", "peptide-length coupled: -0.450 / -0.006, 88.7 per cent "
                                         "beyond both. It reads the class I bulge."),
    "frac_well_coordinated_tp": ("suspicious", "peptide-length coupled: -0.440 / -0.064, 83.4 per "
                                               "cent beyond both; see degree_evenness_tp."),
    "g_even_tcr": ("suspicious", "peptide-length coupled: -0.368 / -0.061, 91.0 per cent beyond both."),
    "g_loop_even": ("suspicious", "CDR3-length coupled: -0.186 / -0.350, 87.7 per cent beyond both."),
    "g_comp_frac": ("suspicious", "CDR3-length coupled: -0.032 / +0.280, 92.7 per cent beyond both."),
    # Measured 2026-09-02 on the 19,213-structure harmonization corpus, within class I. Of the 23
    # descriptors 2.30.0 added, only these two carry a length: the seven strongest binder
    # separators among them -- the gap channel -- keep 96.9 to 99.6 per cent of their variance
    # beyond both lengths, which is why none of the rest is flagged here.
    "sc_gap_depth": ("suspicious", "peptide-length coupled: +0.366 / +0.124, 71.8 per cent of "
                                   "variance beyond both lengths. A longer class I peptide "
                                   "bulges, and the receptor reaches further in where it does."),
    "co_mhc": ("suspicious", "CDR3-length coupled: +0.011 / +0.345, 88.0 per cent beyond both. "
                             "Contact order divides by the target's span, not the loop's, so a "
                             "longer loop spreads over more helix."),
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
