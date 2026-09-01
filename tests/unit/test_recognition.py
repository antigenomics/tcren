"""Unit tests for the Gaussian BN classifier (pure numpy, synthetic data)."""

from __future__ import annotations

import numpy as np
import pytest

from tcren.recognition import (CDR3_FRAME_FEATURES, DESCRIPTORS, FAMILIES, FULL_FEATURES,
                               INTERFACE_SYMMETRY_FEATURES, RECOGNITION_FEATURES,
                               _interface_symmetry, descriptors)


def test_interface_symmetry_contact_counts():
    import math

    import polars as pl
    # 2 CDR1(TRA) + 1 CDR2(TRB) germline, 3 CDR3a, 5 CDR3b  ->  n12=3, n3a=3, n3b=5; nA=2+3=5, nB=1+5=6
    tp = pl.DataFrame({
        "region.type.from": ["CDR1", "CDR1", "CDR2", "CDR3", "CDR3", "CDR3", "CDR3", "CDR3", "CDR3", "CDR3", "CDR3"],
        "chain.type.from":  ["TRA",  "TRA",  "TRB",  "TRA",  "TRA",  "TRA",  "TRB",  "TRB",  "TRB",  "TRB",  "TRB"],
    })
    s = _interface_symmetry(tp)
    assert set(s) == set(INTERFACE_SYMMETRY_FEATURES)
    assert s["cdr3_dominance"] == pytest.approx(8 / 11)        # (3+5)/(3+3+5), CDR3-dominated -> >0.5
    assert s["cdr3_ab_imbalance"] == pytest.approx(2 / 8)      # |3-5|/8, absolute
    assert s["chain_cdr_imbalance"] == pytest.approx(1 / 11)   # |5-6|/11, absolute
    # empty interface -> all NaN, no divide-by-zero
    empty = pl.DataFrame({"region.type.from": [], "chain.type.from": []}, schema={"region.type.from": pl.Utf8, "chain.type.from": pl.Utf8})
    assert all(math.isnan(v) for v in _interface_symmetry(empty).values())
    # these are extra columns, not part of the frozen model vector
    assert not (set(INTERFACE_SYMMETRY_FEATURES) & set(RECOGNITION_FEATURES))


def _auc(y, score):
    """ROC-AUC via the Mann-Whitney U statistic (scipy only; sklearn isn't in the lean CI env)."""
    from scipy.stats import rankdata
    y = np.asarray(y); npos = int(y.sum()); nneg = len(y) - npos
    r = rankdata(np.asarray(score, float))               # average ranks (tie-safe)
    return (r[y == 1].sum() - npos * (npos + 1) / 2) / (npos * nneg)


def _data(seed=0, n=400, p=5):
    rng = np.random.default_rng(seed)
    base = rng.normal(size=(n, p))
    base[:, 1] += 0.8 * base[:, 0]                       # inject dependence 0 -> 1
    y = rng.integers(0, 2, n)
    m = rng.integers(0, 2, n)
    X = base + y[:, None] * np.array([1.0, 0.8, 0.5, 0.0, -0.6]) + m[:, None] * np.array([0, 0, 0.3, 0.4, 0.0])
    return X, y, m


def test_full_feature_schema():
    assert len(RECOGNITION_FEATURES) == 40
    assert len(CDR3_FRAME_FEATURES) == 18
    assert FULL_FEATURES == RECOGNITION_FEATURES + CDR3_FRAME_FEATURES
    assert len(set(FULL_FEATURES)) == len(FULL_FEATURES) == 58      # no duplicate column names
    # the FramePose strain trio (the forced-pose signal) is present
    assert {"cdr3b_topep", "cdr3b_reach", "cdr3b_ext"} <= set(CDR3_FRAME_FEATURES)
    # every energy is named Phi_*: the potential is a property of the interface, not of the column
    assert not [f for f in FULL_FEATURES if f.startswith(("e_", "tcren_", "mj_", "d_")) and f != "dPhi_tcr_pep"]


# --- the descriptor catalogue --------------------------------------------------------------------------
def test_every_emitted_column_is_catalogued():
    assert set(FULL_FEATURES) <= set(DESCRIPTORS)
    assert set(INTERFACE_SYMMETRY_FEATURES) <= set(DESCRIPTORS)
    for name, (family, tcr) in DESCRIPTORS.items():
        assert family in FAMILIES + ("score",), name
        assert isinstance(tcr, bool), name


def test_peptide_and_mhc_only_descriptors_are_the_ones_we_know():
    """The receptor filter is the whole point: a column the TCR does not enter carries cohort
    identity (epitope, allele), so a model handed one can reach a cohort label without physics.

    ``Phi_pep_int``/``n_pep_int`` belong here too — the peptide's contacts with **itself** are a
    property of the epitope's bound conformation, shared by every TCR that reads it."""
    receptor_free = {"Phi_pep_mhc", "dPhi_pep_mhc", "mhc_class_bin", "Phi_pep_int", "n_pep_int"}
    assert {n for n, (_, tcr) in DESCRIPTORS.items() if not tcr} == receptor_free
    assert not set(descriptors(tcr_only=True)) & receptor_free


def test_descriptors_selector():
    assert descriptors("physics", tcr_only=True) == (
        "Phi_tcr_pep", "Phi_tcr_mhc", "Phi_cdr12", "Phi_cdr3a", "Phi_cdr3b", "dPhi_tcr_pep",
        "dPhi_pep_soft", "varPhi_pep_soft", "dPhi_tcr_soft", "varPhi_tcr_soft",
        "dPhi_tra_soft", "dPhi_trb_soft")
    # the fitted composites are gone: the catalogue holds descriptors only (2.26.0)
    assert not {"q_bind", "s_strain", "p_bind", "p_forced", "p_real", "P_native"} & set(descriptors())
    with pytest.raises(ValueError, match="unknown family"):
        descriptors("score")
    # "energetics" and "physics" are the same family under its current and retired name
    assert descriptors("physics") == descriptors("energetics")
    with pytest.raises(ValueError, match="unknown family"):
        descriptors("thermodynamics")


@pytest.mark.slow
def test_recognition_features_full_end_to_end():
    pytest.importorskip("arda")
    from tcren.paths import reference_structure_path
    from tcren.recognition import recognition_features

    f = recognition_features(reference_structure_path("1ao7"), full=True)
    assert set(f) == set(FULL_FEATURES)
    # the CDR-loop energies partition the TCR:peptide interface energy
    assert f["Phi_tcr_pep"] == pytest.approx(f["Phi_cdr12"] + f["Phi_cdr3a"] + f["Phi_cdr3b"], abs=1e-6)
    assert np.isfinite(f["crossing_signed"]) and abs(f["crossing_signed"]) <= 180.0
    # a real crystal complex has both CDR3 loops engaging the peptide
    assert np.isfinite(f["cdr3b_topep"]) and f["cdr3b_topep"] > 0
    assert np.isfinite(f["cdr3a_reach"]) and f["cdr3a_reach"] > 0


def _example_feats():
    """A plausible real-complex feature row (1ao7-like) keyed by RECOGNITION_FEATURES."""
    return {
        "extent": 26.0, "chain_balance": 0.36, "pitch": 25.0, "crossing": 45.0,
        "crossing_signed": -45.0, "dock_d": 25.0, "dock_torsion": 3.35, "dock_tcr_uy": 0.1,
        "dock_tcr_uz": 0.9, "dock_mhc_uy": 0.2, "dock_mhc_uz": 0.95,
        "Phi_cdr12": 0.2, "Phi_cdr3a": 0.1, "Phi_cdr3b": -0.3,
        "Phi_tcr_pep": -0.5, "Phi_tcr_mhc": -1.5, "Phi_pep_mhc": -2.0,
        "dPhi_tcr_pep": -0.4, "dPhi_pep_mhc": -0.6,
        "dPhi_pep_soft": -1.6, "varPhi_pep_soft": 2.6, "dPhi_tcr_soft": -2.5,
        "varPhi_tcr_soft": 2.1, "dPhi_tra_soft": -1.7, "dPhi_trb_soft": -0.8,
        "n_contacts_tp": 30.0, "n_pep_contacted": 8.0, "n_contacts_tm": 40.0,
        "ct_tp_salt_bridge": 1.0, "ct_tm_salt_bridge": 2.0,
        "ct_tm_hydrogen_bond": 6.0, "ct_tp_aromatic": 1.0, "ct_tm_aromatic": 0.0,
        "ct_tp_hydrophobic": 8.0, "ct_tm_hydrophobic": 10.0, "ct_tp_other": 3.0, "ct_tm_other": 4.0,
        "n_hbond": 5.0, "burial": 1950.0, "mhc_class_bin": 0.0,
    }


def test_recognition_features_names_complete():
    from tcren.recognition import RECOGNITION_FEATURES
    assert len(RECOGNITION_FEATURES) == 40
    assert set(_example_feats()) == set(RECOGNITION_FEATURES)   # the example row covers exactly the model inputs


@pytest.mark.slow
def test_recognition_features_end_to_end():
    pytest.importorskip("arda")                             # annotation only; no _geom C-ext needed
    from pathlib import Path

    from tcren.recognition import RECOGNITION_FEATURES, recognition_features

    pdb = Path(__file__).resolve().parents[1] / "assets" / "pdb" / "1ao7.pdb"
    feats = recognition_features(str(pdb))
    assert set(feats) == set(RECOGNITION_FEATURES)
    assert feats["burial"] > 0 and feats["extent"] > 0          # a real complex has a buried interface


