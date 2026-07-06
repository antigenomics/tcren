"""Unit tests for the Gaussian BN classifier (pure numpy, synthetic data)."""

from __future__ import annotations

import numpy as np
import pytest

from tcren.recognition import (BayesianLogisticRecognizer, GaussianBNClassifier, _hill_climb,
                               encode_features)


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


def test_recognizer_name_mismatch_raises():
    rec = BayesianLogisticRecognizer(["extent", "burial"], ["wrong", "names"],
                                     np.zeros(2), np.ones(2), 0.0, np.zeros(2))
    with pytest.raises(ValueError, match="do not match"):
        rec._design(np.zeros((1, 2)))                      # encoded names don't match the fitted model


# --- structure -> recognition features + P(real) (frozen models are pure numpy; no arda needed) --------
def _example_feats():
    """A plausible real-complex feature row (1ao7-like) keyed by RECOGNITION_FEATURES."""
    return {
        "extent": 26.0, "e_tcr_mhc": -1.5, "chain_balance": 0.36, "pitch": 25.0, "crossing": 45.0,
        "dock_d": 25.0, "dock_torsion": 3.35, "dock_tcr_uy": 0.1, "dock_tcr_uz": 0.9,
        "dock_mhc_uy": 0.2, "dock_mhc_uz": 0.95, "e_cdr12": 0.2, "e_cdr3a": 0.1, "e_cdr3b": -0.3,
        "F_tcr_pep": -0.5, "F_tcr_mhc": -1.5, "F_pep_mhc": -2.0, "dF_tcr_pep": -0.4, "dF_pep_mhc": -0.6,
        "n_contacts_tp": 30.0, "n_pep_contacted": 8.0, "n_contacts_tm": 40.0,
        "ct_tp_salt_bridge": 1.0, "ct_tm_salt_bridge": 2.0, "ct_tp_hydrogen_bond": 5.0,
        "ct_tm_hydrogen_bond": 6.0, "ct_tp_aromatic": 1.0, "ct_tm_aromatic": 0.0,
        "ct_tp_hydrophobic": 8.0, "ct_tm_hydrophobic": 10.0, "ct_tp_other": 3.0, "ct_tm_other": 4.0,
        "n_hbond": 5.0, "burial": 1950.0, "mhc_class_bin": 0.0,
    }


def test_recognition_features_names_complete():
    from tcren.recognition import RECOGNITION_FEATURES
    assert len(RECOGNITION_FEATURES) == 35
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
