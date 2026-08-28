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

from tcren.cohort import (P_NATIVE_CHANNELS, P_NATIVE_FEATURES, P_NATIVE_ORIENT, P_NATIVE_POOL,
                          _channel_columns, p_native)
from tcren.recognition import GaussianBNClassifier

NAMES = ["burial", "n_hbond", "height", "neg_energy", "noise"]
#: burial and n_hbond rise with nativeness, the energy runs BACKWARDS (the forced-pose regime this
#: replaces), height falls, and one column is pure noise that must earn no weight.
TRUE_SIGN = {"burial": +1, "n_hbond": +1, "height": -1, "neg_energy": -1, "noise": 0}


def _cohort(n: int = 400, seed: int = 0):
    rng = np.random.default_rng(seed)
    y = rng.integers(0, 2, n)
    X = np.column_stack([
        2.0 * y + rng.normal(0, 1, n),          # burial
        1.5 * y + rng.normal(0, 1, n),          # n_hbond
        -1.6 * y + rng.normal(0, 1, n),         # height
        -1.8 * y + rng.normal(0, 1, n),         # neg_energy, inverted
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
             dock_d=X[:, 2], crossing_signed=X[:, 2] * 0.5, dock_torsion=X[:, 2] * 0.2,
             dock_tcr_uz=X[:, 2] * 0.1, D2_pep24=X[:, 0] * 0.6, fp_b0_frac_r7=-X[:, 0] * 0.2,
             H_cell=X[:, 0] * 0.5, L_canon=X[:, 1] * 0.3, ab_imb=X[:, 4] * 0.1,
             log_z=X[:, 3] * 0.5, log_lik=X[:, 3] * 0.7)
    return pl.DataFrame(t), y


def test_p_native_scores_a_polars_table_and_tracks_the_latent_class():
    t, y = _table()
    p = p_native(t)
    assert p.shape == (len(y),)
    assert ((p > 0.5).astype(int) == y).mean() > 0.8


def test_p_native_uses_only_the_requested_channels():
    t, _ = _table()
    _, m = p_native(t, channels=("topology",), rule="flat", return_model=True)
    assert set(m.feature_names) <= set(P_NATIVE_FEATURES["topology"])


def test_p_native_skips_columns_the_caller_did_not_compute():
    t, _ = _table()
    _, models = p_native(t.drop("height", "dock_d"), return_model=True)
    names = [n for m in models.values() for n in m.feature_names]
    assert "height" not in names and "dock_d" not in names
    assert "burial" in names


def test_p_native_refuses_a_cohort_smaller_than_its_feature_set():
    t, _ = _table()
    with pytest.raises(ValueError, match="more structures than features"):
        p_native(t.head(5))


def test_p_native_refuses_when_no_channel_is_present():
    with pytest.raises(ValueError, match="at least two usable columns"):
        p_native(pl.DataFrame({"burial": np.arange(50.0)}))


def test_default_feature_set_covers_every_channel():
    assert set(P_NATIVE_POOL) == set(P_NATIVE_CHANNELS)
    assert all(_channel_columns(c) for c in P_NATIVE_CHANNELS)


def test_the_pool_partitions_the_feature_families_exactly_once():
    """No descriptor may reach two channels: the log-odds sum would count it twice."""
    fams = [f for c in P_NATIVE_CHANNELS for f in P_NATIVE_POOL[c]]
    assert len(fams) == len(set(fams))
    assert set(fams) == set(P_NATIVE_FEATURES)
    cols = [n for c in P_NATIVE_CHANNELS for n in _channel_columns(c)]
    assert len(cols) == len(set(cols))


def test_geometry_pools_placement_and_interface():
    assert P_NATIVE_POOL["geometry"] == ("placement", "interface")
    assert set(_channel_columns("geometry")) == (set(P_NATIVE_FEATURES["placement"])
                                                 | set(P_NATIVE_FEATURES["interface"]))


def test_the_sum_rule_is_the_sum_of_its_channels():
    """`rule="sum"` must equal the hand-computed log-odds sum, or the equation in the paper is not
    the thing the code computes."""
    t, _ = _table()
    total = p_native(t, rule="sum")
    parts, priors = [], []
    for ch in P_NATIVE_CHANNELS:
        pc, m = p_native(t, channels=(ch,), rule="flat", return_model=True)
        parts.append(np.log(np.clip(pc, 1e-9, 1 - 1e-9) / (1 - np.clip(pc, 1e-9, 1 - 1e-9))))
        priors.append(m.prior_)
    pri = float(np.mean(priors))
    lg = sum(parts) - 2 * np.log(pri / (1 - pri))
    assert np.allclose(total, 1.0 / (1.0 + np.exp(-lg)), atol=1e-9)


def test_T_is_p_native_over_the_topology_channel_alone():
    t, _ = _table()
    assert np.allclose(p_native(t, channels=("topology",)),
                       p_native(t, channels=("topology",), rule="flat"))


def test_energetics_is_oriented_on_the_favourable_direction():
    """Orienting on the wrong end of the energy axis labels the wrong component native.

    Until 2.17.0 the channel read `F_tcr_pep`, a contact-preference sum in which LOWER is more
    favourable, so it was oriented on `-F_tcr_pep`. It now reads the Potts `neg_energy`, where
    HIGHER is more favourable, so the negation is gone — and the test that matters is that the
    orientation feature is a column of the channel and points the favourable way, not which
    spelling it happens to have."""
    orient = P_NATIVE_ORIENT["energetics"]
    assert orient == "neg_energy"
    assert not orient.startswith("-")
    assert orient.lstrip("-") in _channel_columns("energetics")


def test_unknown_rule_is_rejected():
    t, _ = _table()
    with pytest.raises(ValueError, match="rule must be"):
        p_native(t, rule="average")


def test_an_absent_orientation_feature_raises_instead_of_flipping_the_labels():
    """A mixture is identified only up to a label swap, so the orientation feature decides which
    component is called native. When it is missing the old code fell back to the first surviving
    column, whose direction is arbitrary — the 2.17.0 energetics migration turned that into a
    measured sign reversal (Spearman +0.63 became −0.63 against the same reference). It raises now.
    """
    t, _ = _table()
    with pytest.raises(ValueError, match="cannot orient this fit"):
        p_native(t, channels=("geometry",), rule="flat", orient_by="a_column_not_in_the_table")

    # naming a column that IS present still fits, and flipping its sign flips the labelling —
    # which is the reversal the fallback used to introduce without saying so
    up = p_native(t, channels=("geometry",), rule="flat", orient_by="burial")
    down = p_native(t, channels=("geometry",), rule="flat", orient_by="-burial")
    assert np.corrcoef(up, down)[0, 1] < -0.99
