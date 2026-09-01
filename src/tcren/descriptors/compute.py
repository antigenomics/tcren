"""Structure -> descriptor values: the interface terms this package computes itself.

Everything here reads a :class:`~tcren.contactmap.ContactMap` or a
:class:`~tcren.structure.model.Structure` and returns numbers. The *catalogue* -- what the columns
mean and which family each belongs to -- is :mod:`tcren.descriptors.catalogue`, and the batch
dispatch is :mod:`tcren.descriptors.table`.

The energetics, topology, potts and kinetics families are **not** computed here: they belong to
:mod:`tcren.pipeline` / :mod:`tcren.ddg`, :mod:`tcren.footprint` / :mod:`tcren.interface_graph`,
:mod:`tcren.potts` and :mod:`tcren.mechanics` respectively, and this module calls them. What is left
here is the interface block -- burial, extent, chain balance, the contact-type tallies and the
CDR3-frame placement terms -- which has no other home.

Heavy imports stay function-local so a bare ``import tcren`` remains dependency-light.
"""
from __future__ import annotations

import math
import warnings

import numpy as np

from .catalogue import (
    CDR3_FRAME_FEATURES,
    FULL_FEATURES,
    INTERFACE_SYMMETRY_FEATURES,
    PEPTIDE_INTERNAL_FEATURES,
    RECOGNITION_FEATURES,
    TCR_PLACEMENT_FEATURES,
    _CDR3_FRAME_KEYS,
    _CT_TYPES,
    _EPS,
    _TCR_TYPES,
)

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

    from ..structure.io import write_pdb

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
    from ..docking.angles import _chain_ca, _groove_frame

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

    from ..annotation import classify_chains
    from ..contact_types import contact_type_counts
    from ..contactmap import ContactMap
    from ..energetics.mutation import reference_delta
    from ..mhc import annotate_mhc
    from ..oracle import _native_peptide
    from ..docking.angles import docking_angles
    from ..docking.tcrdock_geometry import docking_geometry
    from ..energetics.scoring import _interface_energy
    from ..potential import mj as _mj
    from ..potential import tcren2 as _tcren2
    from ..structure import Structure, import_structure

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
    from ..energetics.mutation import smoothed_reference
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
    from ..clashes import interface_clashes
    from ..mechanics.stability import contact_stability

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
    from ..contactmap import ContactMap
    try:
        cm = ContactMap.from_structure(s)
        return _interface_symmetry(cm.interface("tcr_peptide", tcr_regions="all"))
    except Exception:  # noqa: BLE001 - no peptide/receptor chain etc.
        return {k: math.nan for k in INTERFACE_SYMMETRY_FEATURES}


def _placement_columns(s) -> dict[str, float]:
    """Receptor-placement extra output columns (:data:`TCR_PLACEMENT_FEATURES`) — the ride height,
    in-plane shift and offset of the CDR centroid over the groove. NaN where the frame is undefined."""
    from ..docking.angles import tcr_placement
    try:
        tp = tcr_placement(s)
        return {k: float(getattr(tp, k)) for k in TCR_PLACEMENT_FEATURES}
    except Exception:  # noqa: BLE001 - no groove frame / no receptor chain
        return dict.fromkeys(TCR_PLACEMENT_FEATURES, math.nan)


def _footprint_columns(s, radii=(7.0, 8.0)) -> dict[str, float]:
    """Footprint shape extra output columns — the ``topology`` family (:mod:`tcren.footprint`).

    Rigid-motion invariant, so the structure does not need orienting; it needs only chain typing and
    CDR markup, which the batch annotation has already done by the time this runs."""
    from ..topology.footprint import footprint_features
    try:
        return footprint_features(s, radii=radii)
    except Exception:  # noqa: BLE001 - no peptide/receptor chain etc.
        return dict.fromkeys(footprint_topology_features(radii) + FOOTPRINT_SIZE_FEATURES, math.nan)


def _peptide_internal_columns(s) -> dict[str, float]:
    """Intra-peptide extra output columns (:data:`PEPTIDE_INTERNAL_FEATURES`) for the recognize table.

    The peptide's MJ contact energy with itself and the number of those contacts — the term the three
    interface energies leave out. NaN/0 where the structure has no peptide chain."""
    from ..contactmap import ContactMap
    from ..potential import mj as _mj
    from ..energetics.scoring import intra_peptide_energy
    try:
        cm = ContactMap.from_structure(s, peptide_internal=True)
        return {"Phi_pep_int": float(intra_peptide_energy(cm, _mj())),
                "n_pep_int": float(cm.peptide_internal.height)}
    except Exception:  # noqa: BLE001 - no peptide chain etc.
        return {k: math.nan for k in PEPTIDE_INTERNAL_FEATURES}
