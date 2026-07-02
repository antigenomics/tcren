"""Fast tests for the frozen binder classifier."""

from __future__ import annotations

import pytest

from tcren.binder import BINDER_MODEL, binder_score


def test_model_shape():
    m = BINDER_MODEL
    assert len(m["features"]) == len(m["mu"]) == len(m["sigma"]) == len(m["w"]) == 5
    assert m["denoised_auc"] > 0.9


def test_strong_interface_scores_high_weak_low():
    strong = {"pm_cov_ntcr": 34, "chain_balance": 0.45, "n_hbond": 14, "dSASA": 2400, "pp_combo": 1.6}
    weak = {"pm_cov_ntcr": 18, "chain_balance": 0.10, "n_hbond": 2, "dSASA": 1500, "pp_combo": -1.2}
    ps, pw = binder_score(strong), binder_score(weak)
    assert pw < 0.5 < ps
    assert 0.0 <= pw and ps <= 1.0


def test_score_monotone_in_burial():
    base = {"pm_cov_ntcr": 26, "chain_balance": 0.33, "n_hbond": 7, "dSASA": 1950, "pp_combo": 0.0}
    more = {**base, "dSASA": 2600}
    assert binder_score(more) > binder_score(base)  # more burial -> more binder-like (w_dSASA>0)


@pytest.mark.slow
def test_binder_features_end_to_end():
    pytest.importorskip("arda")
    pytest.importorskip("tcren._geom")
    from tcren.binder import binder_features
    from tcren.paths import reference_structure_path

    feats = binder_features(reference_structure_path("1ao7"))
    assert set(feats) == set(BINDER_MODEL["features"])
    assert feats["pm_cov_ntcr"] > 0 and feats["dSASA"] > 0  # a real TCR-pMHC has a buried interface
    assert 0.0 <= binder_score(feats) <= 1.0
