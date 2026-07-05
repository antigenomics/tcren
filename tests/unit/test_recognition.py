"""Unit tests for the Gaussian BN classifier (pure numpy, synthetic data)."""

from __future__ import annotations

import numpy as np
import pytest

from tcren.recognition import GaussianBNClassifier, _hill_climb


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
    from sklearn.metrics import roc_auc_score
    p = clf.predict_proba(X, m)[:, 1]
    assert roc_auc_score(y, p) > 0.75


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


def test_hill_climb_empty_on_independent_data():
    rng = np.random.default_rng(1)
    Z = rng.normal(size=(500, 4))                        # independent columns
    struct = _hill_climb(Z, max_parents=2)
    assert sum(len(v) for v in struct.values()) <= 1     # ~no spurious edges
