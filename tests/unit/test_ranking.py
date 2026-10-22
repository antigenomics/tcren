"""The epitope-ranking score: frozen reference, fixed signs, and no query-cohort coupling."""
from __future__ import annotations

import json

import numpy as np
import polars as pl
import pytest

from tcren.ranking import FORMAT, bundled_reference, epitope_rank


def _table(n=4, **overrides):
    base = {"complex.id": [f"m{i}" for i in range(n)],
            "dPhi_tcr_pep": [-6.0, -3.0, 0.0, 2.0][:n],
            "iptm": [0.9, 0.8, 0.6, 0.4][:n]}
    base.update(overrides)
    return pl.DataFrame(base)


def test_bundled_reference_is_the_declared_format_and_carries_its_own_hash():
    ref = bundled_reference()
    assert ref["format"] == FORMAT
    assert len(ref["reference_sha256"]) == 64
    assert [c["name"] for c in ref["components"]] == ["energy", "confidence"]
    # The energy enters negated: lower Phi is the more favourable interface.
    assert {c["name"]: c["direction"] for c in ref["components"]} == {"energy": -1, "confidence": 1}


def test_score_is_the_equal_weight_mean_of_standardized_components():
    ref = bundled_reference()
    t = _table()
    got = epitope_rank(t, ref)
    z = []
    for c in ref["components"]:
        x = np.asarray(t[c["column"]], float)
        z.append(c["direction"] * (x - c["mean"]) / c["sd"])
    assert np.allclose(got["epitope_rank_score"].to_numpy(), np.mean(z, axis=0))
    assert got["n_components"].to_list() == [2] * t.height


def test_a_more_favourable_energy_scores_higher():
    t = _table(n=2, dPhi_tcr_pep=[-8.0, 1.0], iptm=[0.7, 0.7])
    s = epitope_rank(t)["epitope_rank_score"].to_numpy()
    assert s[0] > s[1]


def test_batch_and_single_row_agree_because_nothing_reads_the_cohort():
    t = _table()
    batch = epitope_rank(t)["epitope_rank_score"].to_numpy()
    alone = np.concatenate([epitope_rank(t[i:i + 1])["epitope_rank_score"].to_numpy()
                            for i in range(t.height)])
    assert np.allclose(batch, alone)


def test_unrelated_rows_never_move_an_existing_score():
    t = _table()
    extra = pl.concat([t, _table(n=1).with_columns(pl.lit("other").alias("complex.id"))])
    assert np.allclose(epitope_rank(t)["epitope_rank_score"].to_numpy(),
                       epitope_rank(extra)["epitope_rank_score"].to_numpy()[:t.height])


def test_a_missing_component_is_flagged_and_excluded_rather_than_imputed():
    t = _table(n=2, iptm=[0.8, None])
    got = epitope_rank(t)
    assert got["missing_confidence"].to_list() == [False, True]
    assert got["n_components"].to_list() == [2, 1]
    ref = bundled_reference()
    energy = next(c for c in ref["components"] if c["name"] == "energy")
    expected = energy["direction"] * (t["dPhi_tcr_pep"][1] - energy["mean"]) / energy["sd"]
    assert got["epitope_rank_score"][1] == pytest.approx(expected)


def test_a_row_with_no_observed_component_is_unavailable_not_neutral():
    t = _table(n=1, dPhi_tcr_pep=[None], iptm=[None])
    got = epitope_rank(t)
    assert got["n_components"].to_list() == [0]
    assert got["epitope_rank_score"][0] is None


@pytest.mark.parametrize("table, message", [
    (pl.DataFrame({"dPhi_tcr_pep": [1.0], "iptm": [0.5]}), "complex.id"),
    (_table(n=2).with_columns(pl.Series("complex.id", ["a", "a"])), "unique"),
    (_table(n=1).with_columns(pl.Series("complex.id", [""])), "unique"),
    (_table(n=1).drop("iptm"), "missing component columns"),
    (_table(n=1, dPhi_tcr_pep=[float("inf")]), "infinite"),
])
def test_invalid_input_is_rejected(table, message):
    with pytest.raises(ValueError, match=message):
        epitope_rank(table)


@pytest.mark.parametrize("bad, message", [
    ({"format": "something-else"}, "unsupported reference format"),
    ({"format": FORMAT, "components": []}, "no components"),
    ({"format": FORMAT, "components": [{"name": "e", "column": "dPhi_tcr_pep",
                                        "direction": 0, "mean": 0.0, "sd": 1.0}]}, "direction"),
    ({"format": FORMAT, "components": [{"name": "e", "column": "dPhi_tcr_pep",
                                        "direction": 1, "mean": 0.0, "sd": 0.0}]}, "nonpositive"),
    ({"format": FORMAT, "components": [{"name": "e", "column": "dPhi_tcr_pep"}]}, "missing"),
])
def test_invalid_reference_is_rejected(bad, message):
    with pytest.raises(ValueError, match=message):
        epitope_rank(_table(n=1), bad)


def test_reference_resource_is_valid_json_with_finite_frozen_constants():
    ref = bundled_reference()
    assert json.dumps(ref)
    for c in ref["components"]:
        assert np.isfinite(c["mean"]) and np.isfinite(c["sd"]) and c["sd"] > 0
        assert c["n"] == 729


def test_epitope_routes_use_the_same_raw_component_centering():
    from tcren.screen import CohortKind, bundled_reference as screen_reference, screen_scores

    table = _table(n=3, iptm=[0.8, None, 0.6])
    ranked = epitope_rank(table)
    screened = screen_scores(table, screen_reference(CohortKind.EPITOPE_RANKING), min_cohort=1)
    for name in ("energy", "confidence"):
        np.testing.assert_allclose(ranked[f"z_{name}"], screened[f"z_{name}"], equal_nan=True)


@pytest.mark.parametrize("mean", [float("nan"), float("inf"), -float("inf")])
def test_nonfinite_reference_means_are_rejected(mean):
    reference = bundled_reference()
    reference["components"][0]["mean"] = mean
    with pytest.raises(ValueError, match="nonfinite mean"):
        epitope_rank(_table(), reference)
