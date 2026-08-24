"""P(native): the latent-class Bayes network, and what it replaces.

The claim being pinned is narrow and checkable on synthetic data: with the class node unobserved,
EM recovers the class, and it recovers each channel's SIGN on its own. That second part is the whole
argument for dropping the measured coupling `C*` -- the energy channel's inversion on forced poses
stops being a thing to measure and becomes a coefficient like any other.
"""
from __future__ import annotations

import numpy as np
import polars as pl
import pytest

from tcren.cohort import P_NATIVE_CHANNELS, P_NATIVE_FEATURES, p_native
from tcren.recognition import GaussianBNClassifier

NAMES = ["burial", "n_hbond", "height", "F_tcr_pep", "noise"]
#: burial and n_hbond rise with nativeness, the energy runs BACKWARDS (the forced-pose regime this
#: replaces), height falls, and one column is pure noise that must earn no weight.
TRUE_SIGN = {"burial": +1, "n_hbond": +1, "height": -1, "F_tcr_pep": -1, "noise": 0}


def _cohort(n: int = 400, seed: int = 0):
    rng = np.random.default_rng(seed)
    y = rng.integers(0, 2, n)
    X = np.column_stack([
        2.0 * y + rng.normal(0, 1, n),          # burial
        1.5 * y + rng.normal(0, 1, n),          # n_hbond
        -1.6 * y + rng.normal(0, 1, n),         # height
        -1.8 * y + rng.normal(0, 1, n),         # F_tcr_pep, inverted
        rng.normal(0, 1, n)])                   # noise
    return X, y


def test_em_log_likelihood_is_monotone_with_the_graph_held_fixed():
    m = GaussianBNClassifier(NAMES).fit_em(_cohort()[0])
    ll = m.loglik_
    assert len(ll) > 1
    assert all(b >= a - 1e-9 for a, b in zip(ll, ll[1:])), ll


def test_em_recovers_the_latent_class_without_any_label():
    X, y = _cohort()
    p = GaussianBNClassifier(NAMES).fit_em(X).predict_proba(X)[:, 1]
    assert ((p > 0.5).astype(int) == y).mean() > 0.85


def test_em_learns_each_channel_sign_including_the_inverted_energy():
    """This is what makes the measured coupling C* unnecessary: the sign is fitted, not measured."""
    X, _ = _cohort()
    m = GaussianBNClassifier(NAMES).fit_em(X)
    w = {n: m.nodes_[j]["beta"][-2] for j, n in enumerate(NAMES)}
    for n, sign in TRUE_SIGN.items():
        if sign:
            assert np.sign(w[n]) == sign, (n, w[n])
    assert abs(w["noise"]) < min(abs(w[n]) for n in NAMES if TRUE_SIGN[n])


def test_orientation_is_deterministic_so_the_components_cannot_swap():
    X, _ = _cohort()
    a = GaussianBNClassifier(NAMES).fit_em(X).predict_proba(X)[:, 1]
    b = GaussianBNClassifier(NAMES).fit_em(X).predict_proba(X)[:, 1]
    assert a == pytest.approx(b)
    # ...and it points the way `orient_by` says, not the other way
    m = GaussianBNClassifier(NAMES).fit_em(X, orient_by="burial")
    hi, lo = m.responsibilities_ > 0.5, m.responsibilities_ <= 0.5
    assert X[hi, 0].mean() > X[lo, 0].mean()


def test_anchors_pin_their_rows_and_orient_the_fit():
    X, y = _cohort()
    anchors = {int(i): float(y[i]) for i in range(30)}
    m = GaussianBNClassifier(NAMES).fit_em(X, anchors=anchors)
    for i, v in anchors.items():
        assert m.responsibilities_[i] == pytest.approx(v)
    p = m.predict_proba(X)[:, 1]
    assert ((p > 0.5).astype(int) == y).mean() > 0.85


def test_unsupervised_and_semi_supervised_agree_in_direction():
    X, y = _cohort()
    a = GaussianBNClassifier(NAMES).fit_em(X).predict_proba(X)[:, 1]
    b = GaussianBNClassifier(NAMES).fit_em(
        X, anchors={int(i): float(y[i]) for i in range(30)}).predict_proba(X)[:, 1]
    assert np.corrcoef(a, b)[0, 1] > 0.9


def test_relearning_the_structure_is_off_by_default_and_documented_as_non_monotone():
    m = GaussianBNClassifier(NAMES).fit_em(_cohort()[0], relearn_structure=True, rounds=25)
    assert len(m.loglik_) > 1                          # it runs; monotonicity is not claimed for it
    assert "monotone" in GaussianBNClassifier.fit_em.__doc__


# --- the table-level entry point --------------------------------------------------------------

def _table(n: int = 300):
    X, y = _cohort(n, seed=1)
    t = {n_: X[:, j] for j, n_ in enumerate(NAMES)}
    t.update(chain_balance=X[:, 0] * 0.3, n_pep_contacted=X[:, 1] * 0.4, n_clashes=X[:, 4],
             dock_d=X[:, 2], crossing_signed=X[:, 2] * 0.5, pitch=X[:, 2] * 0.2,
             dock_tcr_uz=X[:, 2] * 0.1, D2_pep24=X[:, 0] * 0.6, fp_b0_frac_r7=-X[:, 0] * 0.2,
             H_cell=X[:, 0] * 0.5, L_canon=X[:, 1] * 0.3, ab_imb=X[:, 4] * 0.1,
             F_tcr_mhc=X[:, 3] * 0.5, dF_tcr_pep=X[:, 3] * 0.7)
    return pl.DataFrame(t), y


def test_p_native_scores_a_polars_table_and_tracks_the_latent_class():
    t, y = _table()
    p = p_native(t)
    assert p.shape == (len(y),)
    assert ((p > 0.5).astype(int) == y).mean() > 0.8


def test_p_native_uses_only_the_requested_channels():
    t, _ = _table()
    _, m = p_native(t, channels=("topology",), return_model=True)
    assert set(m.feature_names) <= set(P_NATIVE_FEATURES["topology"])


def test_p_native_skips_columns_the_caller_did_not_compute():
    t, _ = _table()
    _, m = p_native(t.drop("height", "dock_d"), return_model=True)
    assert "height" not in m.feature_names and "dock_d" not in m.feature_names
    assert "burial" in m.feature_names


def test_p_native_refuses_a_cohort_smaller_than_its_feature_set():
    t, _ = _table()
    with pytest.raises(ValueError, match="more structures than features"):
        p_native(t.head(5))


def test_p_native_refuses_when_no_channel_is_present():
    with pytest.raises(ValueError, match="at least two usable columns"):
        p_native(pl.DataFrame({"burial": np.arange(50.0)}))


def test_default_feature_set_covers_every_channel():
    assert set(P_NATIVE_FEATURES) == set(P_NATIVE_CHANNELS)
    assert all(P_NATIVE_FEATURES[c] for c in P_NATIVE_CHANNELS)
