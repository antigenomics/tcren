"""Unit tests for the Gaussian BN classifier (pure numpy, synthetic data)."""

from __future__ import annotations

import numpy as np
import pytest

from tcren.recognition import (BayesianLogisticRecognizer, CDR3_FRAME_FEATURES, FORCED_POSE_MODEL,
                               DESCRIPTORS, FAMILIES, FULL_FEATURES, GaussianBNClassifier,
                               INTERFACE_SYMMETRY_FEATURES, RECOGNITION_FEATURES, _hill_climb,
                               _interface_symmetry, descriptors, encode_features, forced_pose_score,
                               kit_score)


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


def test_fit_predict_separates_classes():
    X, y, m = _data()
    clf = GaussianBNClassifier([f"x{i}" for i in range(5)], max_parents=2).fit(X, y, m)
    p = clf.predict_proba(X, m)[:, 1]
    assert _auc(y, p) > 0.75


def test_structure_recovers_injected_edge():
    X, y, m = _data()
    clf = GaussianBNClassifier([f"x{i}" for i in range(5)], max_parents=2).fit(X, y, m)
    # some edge between features 0 and 1 (either direction) should be learned
    assert 0 in clf.structure_[1] or 1 in clf.structure_[0]


def test_save_load_roundtrip(tmp_path):
    X, y, m = _data()
    clf = GaussianBNClassifier([f"x{i}" for i in range(5)]).fit(X, y, m)
    f = tmp_path / "bn.json.gz"
    clf.save(f)
    clf2 = GaussianBNClassifier.load(f)
    assert np.allclose(clf.predict_proba(X, m)[:, 1], clf2.predict_proba(X, m)[:, 1])


def test_nan_safe():
    X, y, m = _data()
    clf = GaussianBNClassifier([f"x{i}" for i in range(5)]).fit(X, y, m)
    Xn = X.copy(); Xn[0, 0] = np.nan; Xn[3, 2] = np.nan
    assert np.isfinite(clf.predict_proba(Xn, m)[:, 1]).all()


def test_to_dot_has_class_and_mhc_nodes():
    X, y, m = _data()
    clf = GaussianBNClassifier([f"x{i}" for i in range(5)]).fit(X, y, m)
    dot = clf.to_dot(coef_threshold=0.05)
    assert dot.startswith("digraph BN {") and dot.rstrip().endswith("}")
    assert "y ->" in dot          # class node influences some feature


def test_marginal_over_all_equals_full_llr():
    X, y, m = _data()
    names = [f"x{i}" for i in range(5)]
    clf = GaussianBNClassifier(names).fit(X, y, m)
    full = clf.decision_function(X, m)
    marg = clf.marginal_decision(X, names, m)             # marginalise nothing out
    assert np.corrcoef(full, marg)[0, 1] > 0.999          # same joint log-likelihood ratio


def test_marginal_subset_valid_and_separates():
    X, y, m = _data()
    names = [f"x{i}" for i in range(5)]
    clf = GaussianBNClassifier(names).fit(X, y, m)
    s = clf.marginal_decision(X, ["x0", "x2", "x4"], m)   # keep a subset, marginalise the rest
    assert np.isfinite(s).all() and len(s) == len(y)
    assert _auc(y, s) > 0.6                               # the kept features still carry class signal


def test_hill_climb_empty_on_independent_data():
    rng = np.random.default_rng(1)
    Z = rng.normal(size=(500, 4))                        # independent columns
    struct = _hill_climb(Z, max_parents=2)
    assert sum(len(v) for v in struct.values()) <= 1     # ~no spurious edges


# --- distribution-aware logistic recognizer -----------------------------------------------------------
def test_encode_features_distribution_encodings():
    names = ["extent", "dock_torsion", "chain_balance", "n_hbond", "burial", "dock_tcr_uy"]
    X = np.array([[20.0, 1.5, 0.3, 5.0, 1800.0, 0.1],
                  [10.0, 6.0, 0.5, 2.0, 1200.0, -0.4]])
    Z, enc = encode_features(X, names)
    assert enc == ["extent", "dock_torsion_cos", "dock_torsion_sin", "chain_balance_logit",
                   "burial", "dock_tcr_uy"]                # torsion->cos/sin, balance->logit, n_hbond dropped
    j = enc.index("dock_torsion_cos")
    assert Z[0, j] == pytest.approx(np.cos(1.5))
    k = enc.index("dock_torsion_sin")
    assert Z[0, k] == pytest.approx(np.sin(1.5))
    lg = enc.index("chain_balance_logit")                  # logit(2*0.5) clipped -> large positive
    assert Z[1, lg] > 5.0
    assert Z[:, enc.index("extent")].tolist() == [20.0, 10.0]   # counts stay linear


def test_recognizer_predict_roundtrip_and_nan_safe(tmp_path):
    names = ["extent", "dock_torsion", "chain_balance", "burial"]
    X = np.array([[20.0, 1.5, 0.3, 1800.0], [10.0, 6.0, 0.0, 1200.0], [30.0, 0.3, 0.5, 2200.0]])
    Z, enc = encode_features(X, names)
    mean, sd = Z.mean(0), Z.std(0) + 1e-9                   # realistic train stats -> standardised O(1)
    beta = np.linspace(-1, 1, len(enc))
    rec = BayesianLogisticRecognizer(names, enc, mean, sd, 0.3, beta)
    p = rec.predict_proba(X)[:, 1]
    assert p.shape == (3,) and np.all(np.isfinite(p)) and np.all((p > 0) & (p < 1))
    f = tmp_path / "logit.json.gz"
    rec.save(f)
    assert np.allclose(p, BayesianLogisticRecognizer.load(f).predict_proba(X)[:, 1])
    Xn = X.copy(); Xn[0, 2] = np.nan                       # nan chain_balance -> nan logit -> imputed
    assert np.isfinite(rec.predict_proba(Xn)[:, 1]).all()


# --- full feature vector (core + CDR3-frame) -----------------------------------------------------------
def test_full_feature_schema():
    assert len(RECOGNITION_FEATURES) == 34
    assert len(CDR3_FRAME_FEATURES) == 18
    assert FULL_FEATURES == RECOGNITION_FEATURES + CDR3_FRAME_FEATURES
    assert len(set(FULL_FEATURES)) == len(FULL_FEATURES) == 52      # no duplicate column names
    # the FramePose strain trio (the forced-pose signal) is present
    assert {"cdr3b_topep", "cdr3b_reach", "cdr3b_ext"} <= set(CDR3_FRAME_FEATURES)
    # every energy is named F_*: the potential is a property of the interface, not of the column
    assert not [f for f in FULL_FEATURES if f.startswith(("e_", "tcren_", "mj_", "d_")) and f != "dF_tcr_pep"]


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

    ``F_pep_int``/``n_pep_int`` belong here too — the peptide's contacts with **itself** are a
    property of the epitope's bound conformation, shared by every TCR that reads it."""
    receptor_free = {"F_pep_mhc", "dF_pep_mhc", "mhc_class_bin", "F_pep_int", "n_pep_int"}
    assert {n for n, (_, tcr) in DESCRIPTORS.items() if not tcr} == receptor_free
    assert not set(descriptors(tcr_only=True)) & receptor_free


def test_descriptors_selector():
    assert descriptors("physics", tcr_only=True) == (
        "F_tcr_pep", "F_tcr_mhc", "F_cdr12", "F_cdr3a", "F_cdr3b", "dF_tcr_pep")
    # scores are outputs, excluded unless asked for
    assert "q_bind" not in descriptors()
    assert "q_bind" in descriptors("score", with_scores=True)
    with pytest.raises(ValueError, match="unknown family"):
        descriptors("energetics")


def test_forced_pose_model_shape_and_formula():
    m = FORCED_POSE_MODEL
    assert m["features"] == ("dock_d", "cdr3b_reach", "cdr3b_topep", "cdr3a_ext",
                             "extent_per_ct", "chain_balance")
    assert len(m["coef"]) == 6 and m["cv_auc"] == pytest.approx(0.762)
    feats = {"dock_d": 25.0, "cdr3b_reach": 12.0, "cdr3b_topep": 4.0, "cdr3a_ext": 6.0,
             "extent": 26.0, "n_contacts_tp": 40.0, "chain_balance": 0.4}
    # reproduce the explicit sigmoid (pins the frozen coefficients + extent_per_ct derivation)
    vals = [feats["dock_d"], feats["cdr3b_reach"], feats["cdr3b_topep"], feats["cdr3a_ext"],
            feats["extent"] / feats["n_contacts_tp"], feats["chain_balance"]]
    z = m["intercept"] + sum(w * v for w, v in zip(m["coef"], vals))
    assert forced_pose_score(feats) == pytest.approx(1.0 / (1.0 + np.exp(-z)))
    assert 0.0 <= forced_pose_score(feats) <= 1.0


def test_forced_pose_score_nan_safe():
    assert np.isnan(forced_pose_score({"dock_d": float("nan")}))          # missing features -> NaN
    assert np.isnan(forced_pose_score({"extent": 26.0, "n_contacts_tp": 0.0}))  # zero contacts -> NaN


def test_kit_score_combines_and_orders():
    p_bind = np.array([0.2, 0.5, 0.9, 0.6])
    iptm = np.array([0.6, 0.7, 0.85, 0.4])
    k = kit_score(p_bind, iptm)
    assert k.shape == (4,)
    # equal-weight z-sum: the structure high on both ranks top, low on both ranks bottom
    assert k.argmax() == 2 and k.argmin() == 0
    # translation/scale invariance of each z term (only relative ordering matters)
    assert np.allclose(kit_score(p_bind, iptm), kit_score(p_bind * 10 + 3, iptm - 1))


@pytest.mark.slow
def test_recognition_features_full_end_to_end():
    pytest.importorskip("arda")
    from tcren.paths import reference_structure_path
    from tcren.recognition import recognition_features

    f = recognition_features(reference_structure_path("1ao7"), full=True)
    assert set(f) == set(FULL_FEATURES)
    # the CDR-loop energies partition the TCR:peptide interface energy
    assert f["F_tcr_pep"] == pytest.approx(f["F_cdr12"] + f["F_cdr3a"] + f["F_cdr3b"], abs=1e-6)
    assert np.isfinite(f["crossing_signed"]) and abs(f["crossing_signed"]) <= 180.0
    # a real crystal complex has both CDR3 loops engaging the peptide
    assert np.isfinite(f["cdr3b_topep"]) and f["cdr3b_topep"] > 0
    assert np.isfinite(f["cdr3a_reach"]) and f["cdr3a_reach"] > 0


def test_recognizer_name_mismatch_raises():
    rec = BayesianLogisticRecognizer(["extent", "burial"], ["wrong", "names"],
                                     np.zeros(2), np.ones(2), 0.0, np.zeros(2))
    with pytest.raises(ValueError, match="do not match"):
        rec._design(np.zeros((1, 2)))                      # encoded names don't match the fitted model


# --- structure -> recognition features + P(real) (frozen models are pure numpy; no arda needed) --------
def _example_feats():
    """A plausible real-complex feature row (1ao7-like) keyed by RECOGNITION_FEATURES."""
    return {
        "extent": 26.0, "chain_balance": 0.36, "pitch": 25.0, "crossing": 45.0,
        "crossing_signed": -45.0, "dock_d": 25.0, "dock_torsion": 3.35, "dock_tcr_uy": 0.1,
        "dock_tcr_uz": 0.9, "dock_mhc_uy": 0.2, "dock_mhc_uz": 0.95,
        "F_cdr12": 0.2, "F_cdr3a": 0.1, "F_cdr3b": -0.3,
        "F_tcr_pep": -0.5, "F_tcr_mhc": -1.5, "F_pep_mhc": -2.0, "dF_tcr_pep": -0.4, "dF_pep_mhc": -0.6,
        "n_contacts_tp": 30.0, "n_pep_contacted": 8.0, "n_contacts_tm": 40.0,
        "ct_tp_salt_bridge": 1.0, "ct_tm_salt_bridge": 2.0,
        "ct_tm_hydrogen_bond": 6.0, "ct_tp_aromatic": 1.0, "ct_tm_aromatic": 0.0,
        "ct_tp_hydrophobic": 8.0, "ct_tm_hydrophobic": 10.0, "ct_tp_other": 3.0, "ct_tm_other": 4.0,
        "n_hbond": 5.0, "burial": 1950.0, "mhc_class_bin": 0.0,
    }


def test_recognition_features_names_complete():
    from tcren.recognition import RECOGNITION_FEATURES
    assert len(RECOGNITION_FEATURES) == 34
    assert set(_example_feats()) == set(RECOGNITION_FEATURES)   # the example row covers exactly the model inputs


def test_real_probability_from_frozen_models():
    from tcren.recognition import real_probability
    p = real_probability(_example_feats())                     # loads the shipped logistic + BN
    for k in ("logistic", "bn"):
        assert p[k].shape == (1,) and 0.0 < float(p[k][0]) < 1.0
    p2 = real_probability([_example_feats(), _example_feats()])
    assert p2["logistic"].shape == (2,) and p2["bn"].shape == (2,)


def test_real_probability_nan_safe():
    from tcren.recognition import real_probability
    feats = _example_feats()
    feats["burial"] = float("nan"); feats["dF_tcr_pep"] = float("nan")   # missing terms -> imputed
    p = real_probability(feats)
    assert 0.0 < float(p["logistic"][0]) < 1.0 and 0.0 < float(p["bn"][0]) < 1.0


@pytest.mark.slow
def test_recognition_features_end_to_end():
    pytest.importorskip("arda")                             # annotation only; no _geom C-ext needed
    from pathlib import Path

    from tcren.recognition import RECOGNITION_FEATURES, real_probability, recognition_features

    pdb = Path(__file__).resolve().parents[1] / "assets" / "pdb" / "1ao7.pdb"
    feats = recognition_features(str(pdb))
    assert set(feats) == set(RECOGNITION_FEATURES)
    assert feats["burial"] > 0 and feats["extent"] > 0          # a real complex has a buried interface
    assert 0.0 < float(real_probability(feats)["logistic"][0]) < 1.0


def test_add_cohort_scores_batch_and_error_rows():
    # `recognize --scores` appends the fit-free cohort scores over the whole batch; error rows
    # (structures that failed to parse) must be left untouched.
    from tcren.recognition import _add_cohort_scores
    rows = [{"complex.id": f"c{i}", "burial": 1.0 + i, "n_pep_contacted": 2.0 + i % 3,
             "chain_balance": 0.3 + 0.05 * i, "n_hbond": 5.0 + i, "F_cdr12": 1.0 - 0.2 * i,
             "F_cdr3a": 0.5 + 0.1 * i, "cdr3b_topep": 1.0 + i, "cdr3b_reach": 2.0 - 0.3 * i,
             "extent": 10.0 + i, "n_contacts_tp": 5.0 + i} for i in range(6)]
    rows.append({"complex.id": "bad", "error": "boom"})
    _add_cohort_scores(rows)
    assert all("q_bind" in r and "s_strain" in r for r in rows if "error" not in r)
    assert "q_bind" not in rows[-1]                              # error row untouched
    # cohort scores are within-batch z-like: not all equal, finite
    qs = [r["q_bind"] for r in rows if "error" not in r]
    assert len(set(round(q, 6) for q in qs)) > 1
