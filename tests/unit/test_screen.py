"""The cohort screen: refuses non-cohorts, withholds a score below the floor, never imputes."""
from __future__ import annotations

import numpy as np
import polars as pl
import pytest

from tcren.screen import (
    FLAGS, FORMAT, MIN_COHORT, WITHHOLDING_FLAGS, CohortKind, FlagResponse, apply_flag_response,
    bundled_reference, cohort_flags, cohort_signature, screen_scores, stouffer,
    withhold_reason,
)


def _table(n=MIN_COHORT, **overrides):
    rng = np.random.default_rng(0)
    base = {"complex.id": [f"m{i}" for i in range(n)],
            "iptm": list(rng.uniform(0.6, 0.95, n)),
            "dFT_tcr_pep_cdr": list(rng.normal(0, 2.4, n)),
            "sc_shape": list(rng.uniform(0.3, 0.8, n)),
            "g_cyclo_frac": list(rng.uniform(0.15, 0.35, n)),
            "n_contacts_tp": [20] * n}
    base.update(overrides)
    return pl.DataFrame(base)


# --- cohort detection -------------------------------------------------------------------------

def test_one_receptor_many_peptides_is_epitope_ranking():
    kind, summary = cohort_signature(["AAA", "BBB", "CCC"], [("a", "b")] * 3)
    assert kind is CohortKind.EPITOPE_RANKING
    assert summary == {"n_structures": 3, "n_distinct_peptides": 3, "n_distinct_receptors": 1}


def test_one_peptide_many_receptors_is_tcr_ranking():
    kind, _ = cohort_signature(["AAA"] * 3, [("a", "b"), ("c", "d"), ("e", "f")])
    assert kind is CohortKind.TCR_RANKING


def test_a_mixed_set_is_refused_and_says_what_differs():
    with pytest.raises(ValueError, match="not a cohort: 2 distinct peptides and 2 distinct"):
        cohort_signature(["AAA", "AAA", "BBB", "BBB"],
                         [("a", "b"), ("c", "d"), ("a", "b"), ("c", "d")])


def test_a_set_where_nothing_varies_is_refused():
    with pytest.raises(ValueError, match="nothing varies"):
        cohort_signature(["AAA"] * 4, [("a", "b")] * 4)


@pytest.mark.parametrize("peps, tcrs, message", [
    ([], [], "no structures"),
    (["A"], [], "one peptide and one receptor identity per structure"),
])
def test_malformed_identity_input_is_refused(peps, tcrs, message):
    with pytest.raises(ValueError, match=message):
        cohort_signature(peps, tcrs)


# --- the combination rule ---------------------------------------------------------------------

def test_weighted_stouffer_matches_a_hand_computed_value():
    z = np.array([[1.0, 2.0]])
    w = np.array([0.8, 0.2])
    assert stouffer(z, w)[0] == pytest.approx((0.8 * 1.0 + 0.2 * 2.0) / np.sqrt(0.8**2 + 0.2**2))


def test_weighted_stouffer_reduces_to_plain_stouffer_at_equal_weights():
    z = np.array([[1.0, -0.5, 2.0]])
    assert stouffer(z, np.full(3, 0.3))[0] == pytest.approx(z.sum() / np.sqrt(3))


def test_stouffer_drops_a_missing_component_from_both_sums():
    z = np.array([[1.0, np.nan]])
    assert stouffer(z, np.array([0.8, 0.2]))[0] == pytest.approx(0.8 * 1.0 / 0.8)


def test_stouffer_returns_nan_when_nothing_is_finite():
    assert np.isnan(stouffer(np.array([[np.nan, np.nan]]), np.array([0.5, 0.5]))[0])


@pytest.mark.parametrize("w, message", [
    (np.array([0.5]), "one weight per component"),
    (np.array([-1.0, 2.0]), "non-negative"),
    (np.array([0.0, 0.0]), "not all zero"),
])
def test_invalid_weights_are_rejected(w, message):
    with pytest.raises(ValueError, match=message):
        stouffer(np.array([[1.0, 1.0]]), w)


# --- the cohort floor -------------------------------------------------------------------------

def test_below_the_floor_components_are_emitted_but_no_score():
    got = screen_scores(_table(MIN_COHORT - 1), bundled_reference(CohortKind.TCR_RANKING))
    assert got["screen_score"].null_count() == got.height
    assert got["score_unavailable"][0] == (
        f"cohort of {MIN_COHORT - 1} scorable models is smaller than the minimum of {MIN_COHORT}")
    # the components themselves are still there, which is the point of reporting without scoring
    assert got["z_confidence"].null_count() == 0
    assert got["n_components"].to_list() == [4] * got.height


def test_at_the_floor_a_score_is_emitted():
    got = screen_scores(_table(MIN_COHORT), bundled_reference(CohortKind.TCR_RANKING))
    assert got["screen_score"].null_count() == 0
    assert got["score_unavailable"].null_count() == got.height


# --- scoring semantics ------------------------------------------------------------------------

def test_batch_and_single_row_agree_because_nothing_reads_the_cohort():
    ref = bundled_reference(CohortKind.TCR_RANKING)
    t = _table(MIN_COHORT)
    batch = screen_scores(t, ref)["screen_score"].to_numpy()
    # min_cohort=1 isolates the per-row computation from the floor rule
    alone = np.concatenate([screen_scores(t[i:i + 1], ref, min_cohort=1)["screen_score"].to_numpy()
                            for i in range(t.height)])
    assert np.allclose(batch, alone)


def test_unrelated_rows_never_move_an_existing_score():
    ref = bundled_reference(CohortKind.TCR_RANKING)
    t = _table(MIN_COHORT)
    extra = pl.concat([t, _table(1).with_columns(pl.lit("other").alias("complex.id"))])
    assert np.allclose(screen_scores(t, ref)["screen_score"].to_numpy(),
                       screen_scores(extra, ref)["screen_score"].to_numpy()[:t.height])


def test_default_tcr_ranking_uses_only_confidence_and_topology_at_accepted_ratio():
    ref = bundled_reference(CohortKind.TCR_RANKING)
    t = _table()
    got = screen_scores(t, ref)
    expected = (6 * t["iptm"].to_numpy() + t["g_cyclo_frac"].to_numpy()) / 7
    np.testing.assert_allclose(got["screen_score"].to_numpy(), expected)
    changed = t.with_columns(dFT_tcr_pep_cdr=pl.lit(-100.0), sc_shape=pl.lit(None))
    np.testing.assert_array_equal(screen_scores(changed, ref)["screen_score"].to_numpy(),
                                  got["screen_score"].to_numpy())
    assert got["z_energy"].is_finite().all() and got["z_geometry"].is_finite().all()


def test_reporting_components_cannot_restore_missing_tcr_ranking_channels():
    got = screen_scores(_table(1, iptm=[None], g_cyclo_frac=[None]),
                        bundled_reference(CohortKind.TCR_RANKING), min_cohort=1)
    assert got["screen_score"][0] is None
    assert got["score_unavailable"][0] == "no weighted component observed"


@pytest.mark.parametrize("fraction", [-0.1, 1.1, float("nan")])
def test_invalid_structural_fraction_is_rejected(fraction):
    ref = bundled_reference(CohortKind.TCR_RANKING)
    ref["components"][1]["structural_fraction"] = fraction
    with pytest.raises(ValueError, match="structural_fraction"):
        screen_scores(_table(), ref)


def test_a_missing_component_is_flagged_and_excluded_rather_than_imputed():
    ref = bundled_reference(CohortKind.TCR_RANKING)
    t = _table(MIN_COHORT).with_columns(
        pl.Series("sc_shape", [None] + [0.5] * (MIN_COHORT - 1), dtype=pl.Float64))
    got = screen_scores(t, ref)
    assert got["missing_geometry"].to_list() == [True] + [False] * (MIN_COHORT - 1)
    assert got["n_components"].to_list() == [3] + [4] * (MIN_COHORT - 1)
    assert got["screen_score"][0] is not None


def test_a_row_with_no_observed_component_is_unavailable_not_neutral():
    ref = bundled_reference(CohortKind.TCR_RANKING)
    n = MIN_COHORT
    t = _table(n)
    for c in ("iptm", "dFT_tcr_pep_cdr", "sc_shape", "g_cyclo_frac"):
        t = t.with_columns(pl.Series(c, [None] + t[c].to_list()[1:], dtype=pl.Float64))
    got = screen_scores(t, ref)
    assert got["n_components"][0] == 0
    assert got["screen_score"][0] is None
    assert got["score_unavailable"][0] == "no component observed"


@pytest.mark.parametrize("table, message", [
    (_table(MIN_COHORT).drop("complex.id"), "complex.id"),
    (_table(2).with_columns(pl.Series("complex.id", ["a", "a"])), "unique"),
    (_table(1).with_columns(pl.Series("complex.id", [""])), "unique"),
    (_table(MIN_COHORT).drop("sc_shape"), "missing component columns"),
    (_table(MIN_COHORT).with_columns(
        pl.Series("iptm", [float("inf")] + [0.8] * (MIN_COHORT - 1))), "infinite"),
])
def test_invalid_input_is_rejected(table, message):
    with pytest.raises(ValueError, match=message):
        screen_scores(table, bundled_reference(CohortKind.TCR_RANKING))


@pytest.mark.parametrize("weight", [-0.1, 1.1, float("nan")])
def test_an_out_of_range_weight_is_rejected(weight):
    with pytest.raises(ValueError, match=r"weight must lie in \[0, 1\]"):
        screen_scores(_table(MIN_COHORT), bundled_reference(CohortKind.TCR_RANKING), weight=weight)


@pytest.mark.parametrize("bad, message", [
    ({"format": "something-else"}, "unsupported reference format"),
    ({"format": FORMAT, "components": []}, "no components"),
    ({"format": FORMAT, "components": [{"name": "confidence", "column": "iptm",
                                        "direction": 0, "mean": 0.0, "sd": 1.0}]}, "direction"),
    ({"format": FORMAT, "components": [{"name": "confidence", "column": "iptm",
                                        "direction": 1, "mean": 0.0, "sd": 0.0}]}, "nonpositive"),
    ({"format": FORMAT, "components": [{"name": "confidence", "column": "iptm"}]}, "missing"),
])
def test_invalid_reference_is_rejected(bad, message):
    with pytest.raises(ValueError, match=message):
        screen_scores(_table(MIN_COHORT), bad)


# --- the shipped resources --------------------------------------------------------------------

@pytest.mark.parametrize("task, names", [
    (CohortKind.TCR_RANKING, ["confidence", "energy", "geometry", "topology"]),
    (CohortKind.EPITOPE_RANKING, ["confidence", "energy"]),
])
def test_bundled_references_carry_the_task_component_set(task, names):
    ref = bundled_reference(task)
    assert ref["format"] == FORMAT
    assert len(ref["reference_sha256"]) == 64
    assert [c["name"] for c in ref["components"]] == names
    for c in ref["components"]:
        assert np.isfinite(c["mean"]) and np.isfinite(c["sd"]) and c["sd"] > 0
        assert c["n"] > 0


def test_topology_is_absent_from_the_epitope_reference():
    # its direction is established for receptor ranking and inverted on peptide ranking
    columns = {c["column"] for c in bundled_reference(CohortKind.EPITOPE_RANKING)["components"]}
    assert "g_cyclo_frac" not in columns


def test_the_two_tasks_substitute_different_sides():
    tcr = {c["column"] for c in bundled_reference(CohortKind.TCR_RANKING)["components"]}
    epi = {c["column"] for c in bundled_reference(CohortKind.EPITOPE_RANKING)["components"]}
    assert "dFT_tcr_pep_cdr" in tcr and "dPhi_tcr_pep" not in tcr
    assert "dPhi_tcr_pep" in epi and "dFT_tcr_pep_cdr" not in epi


def test_confidence_carries_the_declared_weight():
    ref = bundled_reference(CohortKind.TCR_RANKING)
    # At weight 1 the bundled TCR score is raw confidence alone.
    got = screen_scores(_table(MIN_COHORT), ref, weight=1.0)
    assert np.allclose(got["screen_score"].to_numpy(), got["iptm"].to_numpy())


@pytest.mark.parametrize("task, expected", [
    # The right weight is task-specific: ipTM is much the strongest input for TCR ranking, while
    # the energy is stronger for epitope ranking, where equal weight reproduces the shipped
    # rank-epitopes result.
    (CohortKind.TCR_RANKING, 6 / 7),
    (CohortKind.EPITOPE_RANKING, 0.5),
])
def test_each_reference_declares_its_own_default_weight(task, expected):
    assert bundled_reference(task)["default_weight"] == expected


def test_the_reference_default_weight_is_used_when_none_is_given():
    ref = bundled_reference(CohortKind.TCR_RANKING)
    t = _table(MIN_COHORT)
    assert np.allclose(screen_scores(t, ref)["screen_score"].to_numpy(),
                       screen_scores(t, ref, weight=ref["default_weight"])["screen_score"].to_numpy())


# --- cohort flags -----------------------------------------------------------------------------

def _confident(n=MIN_COHORT + 10):
    return _table(n).with_columns(pl.Series("iptm", np.linspace(0.84, 0.92, n)))


def _unconfident(n=MIN_COHORT + 10):
    return _table(n).with_columns(pl.Series("iptm", np.linspace(0.60, 0.70, n)))


def test_a_confident_cohort_raises_no_withholding_flag():
    f = cohort_flags(_confident(), bundled_reference(CohortKind.TCR_RANKING))
    assert f["low_confidence"] is False
    assert withhold_reason(f, CohortKind.TCR_RANKING) is None


def test_a_low_confidence_cohort_is_flagged_but_does_not_withhold():
    """The flag reports a regime; it does not change what the user is given.

    It briefly withheld, on the strength of the evaluation panels where it fires on 3 of 29
    cohorts and all 3 rank below chance. On the held-out training pool it fires on 3 of 14 and
    all 3 are good, so acting on it would have suppressed three of the better cohorts.
    """
    ref = bundled_reference(CohortKind.TCR_RANKING)
    t = _unconfident()
    f = cohort_flags(t, ref)
    assert f["low_confidence"] is True
    assert f["cohort_confidence_z"] < ref["cohort_confidence"]["flag_z"]
    assert withhold_reason(f, CohortKind.TCR_RANKING) is None
    got = screen_scores(t, ref)
    assert got["screen_score"].null_count() == 0


def test_no_flag_withholds_on_either_task():
    for task in CohortKind:
        assert WITHHOLDING_FLAGS[task] == ()


def test_an_explicit_suppress_reason_still_withholds():
    # the mechanism stays, so a caller who has their own reason can use it; no flag supplies one
    ref = bundled_reference(CohortKind.TCR_RANKING)
    got = screen_scores(_unconfident(), ref, suppress="caller said so")
    assert got["screen_score"].null_count() == got.height
    assert got["score_unavailable"][0] == "caller said so"


def test_small_cohort_is_flagged_as_well_as_withheld():
    f = cohort_flags(_table(MIN_COHORT - 1), bundled_reference(CohortKind.TCR_RANKING),
                     min_cohort=MIN_COHORT)
    assert f["small_cohort"] is True
    assert f["n_models"] == MIN_COHORT - 1


def test_forced_pose_is_unavailable_unless_the_strain_is_calibrated():
    ref = bundled_reference(CohortKind.TCR_RANKING)
    t = _confident()
    high = np.full(t.height, 5.0)
    # an uncalibrated strain is not comparable to the reference cut, so the flag must abstain
    assert cohort_flags(t, ref, strain=high)["forced_pose"] is None
    assert cohort_flags(t, ref, strain=high, calibrated=True)["forced_pose"] is True
    assert cohort_flags(t, ref, strain=None)["forced_pose"] is None


def test_every_flag_name_is_documented():
    f = cohort_flags(_confident(), bundled_reference(CohortKind.TCR_RANKING))
    assert set(FLAGS) <= set(f)
    for names in WITHHOLDING_FLAGS.values():
        assert set(names) <= set(FLAGS)


def test_cohort_confidence_is_a_z_against_the_reference_distribution():
    """Reported as a distance, not a side of a line: -2.7 and -0.7 are not the same situation."""
    ref = bundled_reference(CohortKind.TCR_RANKING)
    mu, sd = ref["cohort_confidence"]["mean"], ref["cohort_confidence"]["sd"]
    t = _unconfident()
    f = cohort_flags(t, ref)
    assert f["cohort_confidence_z"] == pytest.approx((f["mean_iptm"] - mu) / sd)
    # a confident cohort lands above the flag z, an unconfident one below it
    assert cohort_flags(_confident(), ref)["cohort_confidence_z"] > ref["cohort_confidence"]["flag_z"]


def test_the_flag_z_reproduces_the_quantile_it_replaced():
    # z = -0.6 is the 20th percentile of the reference distribution the earlier hard cut used,
    # so the operating point did not move when the description improved
    c = bundled_reference(CohortKind.TCR_RANKING)["cohort_confidence"]
    assert c["flag_z"] == -0.6
    assert c["mean"] + c["flag_z"] * c["sd"] == pytest.approx(0.7242, abs=5e-4)


def test_disagreement_is_reported_as_a_correlation_not_only_a_flag():
    f = cohort_flags(_confident(), bundled_reference(CohortKind.TCR_RANKING))
    rho = f["confidence_structure_agreement"]
    assert rho != rho or -1.0 <= rho <= 1.0


# --- responses to a flag ----------------------------------------------------------------------

def _scored(table=None):
    ref = bundled_reference(CohortKind.TCR_RANKING)
    t = table if table is not None else _unconfident()
    return t, ref, screen_scores(t, ref)


def test_report_changes_nothing():
    t, ref, sc = _scored()
    out = apply_flag_response(sc, t, ref, response=FlagResponse.REPORT, reason="r")
    assert out.equals(sc)
    assert "flag_response" not in out.columns


def test_reverse_negates_the_score_and_nothing_else():
    t, ref, sc = _scored()
    out = apply_flag_response(sc, t, ref, response=FlagResponse.REVERSE, reason="r")
    assert np.allclose(out["screen_score"].to_numpy(), -sc["screen_score"].to_numpy())
    # the components are untouched, so the evidence behind the response stays visible
    assert np.allclose(out["z_confidence"].to_numpy(), sc["z_confidence"].to_numpy())
    assert out["flag_response"][0] == "reverse"


def test_structure_drops_confidence_entirely():
    t, ref, sc = _scored()
    out = apply_flag_response(sc, t, ref, response=FlagResponse.STRUCTURE, reason="r")
    expected = screen_scores(t, ref, weight=0.0, min_cohort=1)["screen_score"].to_numpy()
    assert np.allclose(out["screen_score"].to_numpy(), expected)


def test_geometry_uses_the_geometry_component_alone():
    t, ref, sc = _scored()
    out = apply_flag_response(sc, t, ref, response=FlagResponse.GEOMETRY, reason="r")
    assert np.allclose(out["screen_score"].to_numpy(), sc["z_geometry"].to_numpy())


def test_geometry_is_refused_when_the_reference_has_no_geometry_component():
    t = _table(MIN_COHORT).select("complex.id", "iptm",
                                  pl.col("dFT_tcr_pep_cdr").alias("dPhi_tcr_pep"))
    ref = bundled_reference(CohortKind.EPITOPE_RANKING)
    with pytest.raises(ValueError, match="no geometry component"):
        apply_flag_response(screen_scores(t, ref), t, ref,
                            response=FlagResponse.GEOMETRY, reason="r")


def test_withhold_nulls_every_score_and_records_the_reason():
    t, ref, sc = _scored()
    out = apply_flag_response(sc, t, ref, response=FlagResponse.WITHHOLD, reason="low_confidence")
    assert out["screen_score"].null_count() == out.height
    assert out["score_unavailable"][0].startswith("withhold applied")


def test_report_is_the_default_response():
    assert FlagResponse("report") is FlagResponse.REPORT
    assert WITHHOLDING_FLAGS[CohortKind.TCR_RANKING] == ()


# --- models the generator never placed on the groove -------------------------------------------

def test_a_model_with_no_interface_is_discarded_not_scored():
    """Zero TCR:peptide contacts means the receptor was docked elsewhere. Every structural
    component for such a model is computed over an empty interface and carries no information."""
    ref = bundled_reference(CohortKind.TCR_RANKING)
    n = MIN_COHORT + 5
    t = _table(n).with_columns(pl.Series("n_contacts_tp", [0] * 3 + [20] * (n - 3)))
    got = screen_scores(t, ref)
    assert got["no_interface"].to_list() == [True] * 3 + [False] * (n - 3)
    assert got["screen_score"][:3].null_count() == 3
    assert got["screen_score"][3:].null_count() == 0
    assert got["score_unavailable"][0].startswith("no TCR:peptide interface")


def test_discarded_models_do_not_count_toward_the_cohort_floor():
    """A cohort that only reaches the floor by including failed models has not reached it."""
    ref = bundled_reference(CohortKind.TCR_RANKING)
    n = MIN_COHORT + 5
    t = _table(n).with_columns(pl.Series("n_contacts_tp", [0] * 8 + [20] * (n - 8)))
    got = screen_scores(t, ref)
    assert got["screen_score"].null_count() == got.height
    reasons = {r for r in got["score_unavailable"].to_list() if r}
    assert any("scorable models is smaller than the minimum" in r for r in reasons)


def test_discarded_models_are_kept_in_the_output():
    # nothing disappears silently; the reader sees which models were dropped and why
    ref = bundled_reference(CohortKind.TCR_RANKING)
    n = MIN_COHORT + 5
    t = _table(n).with_columns(pl.Series("n_contacts_tp", [0] * 3 + [20] * (n - 3)))
    assert screen_scores(t, ref).height == n


def test_the_cohort_reports_how_many_models_have_no_interface():
    ref = bundled_reference(CohortKind.TCR_RANKING)
    n = MIN_COHORT + 5
    t = _table(n).with_columns(pl.Series("n_contacts_tp", [0] * 4 + [20] * (n - 4)))
    f = cohort_flags(t, ref)
    assert f["n_no_interface"] == 4
    assert f["no_interface"] is True
    assert f["frac_no_interface"] == pytest.approx(4 / n)


@pytest.mark.parametrize("response", list(FlagResponse))
def test_flag_responses_preserve_discarded_rows_and_the_cohort_floor(response):
    ref = bundled_reference(CohortKind.TCR_RANKING)
    for n in (MIN_COHORT, MIN_COHORT + 1):
        t = _table(n, n_contacts_tp=[0] + [20] * (n - 1))
        scored = screen_scores(t, ref)
        out = apply_flag_response(scored, t, ref, response=response, reason="low_confidence")
        unavailable = scored["screen_score"].is_null()
        assert out.filter(unavailable)["screen_score"].null_count() == unavailable.sum()
        assert out.filter(unavailable)["score_unavailable"].to_list() == (
            scored.filter(unavailable)["score_unavailable"].to_list())


def test_cohort_flag_uses_the_scorable_floor():
    t = _table(MIN_COHORT, n_contacts_tp=[0] + [20] * (MIN_COHORT - 1))
    assert cohort_flags(t, bundled_reference(CohortKind.TCR_RANKING))["small_cohort"] is True


@pytest.mark.parametrize("mean", [float("nan"), float("inf"), -float("inf")])
def test_nonfinite_reference_means_are_rejected(mean):
    reference = bundled_reference(CohortKind.TCR_RANKING)
    reference["components"][0]["mean"] = mean
    with pytest.raises(ValueError, match="nonfinite mean"):
        screen_scores(_table(), reference)


def test_retired_confidence_coefficient_is_rejected():
    reference = bundled_reference(CohortKind.TCR_RANKING)
    reference["components"][1]["confidence_rho"] = 0.0
    with pytest.raises(ValueError, match="confidence_rho is retired"):
        screen_scores(_table(), reference)


def test_confidence_only_without_iptm_has_explicit_unavailability():
    t = _table(1, iptm=[None])
    result = screen_scores(t, bundled_reference(CohortKind.TCR_RANKING),
                           weight=1.0, min_cohort=1)
    assert result["screen_score"][0] is None
    assert result["score_unavailable"][0] == "no weighted component observed"
    assert np.isfinite(result["z_energy"][0])


def test_template_gate_uses_confidence_only_but_keeps_components():
    t = _table(3, iptm=[0.6, None, 0.9], n_contacts_tp=[20, 20, 0])
    ref = bundled_reference(CohortKind.TCR_RANKING)
    gated = screen_scores(t, ref, weight=0.6, native_template=False, min_cohort=1)
    assert gated["confidence_weight"].to_list() == [1.0] * 3
    assert gated["native_template"].to_list() == [False] * 3
    assert gated["screen_score"][0] == pytest.approx(gated["iptm"][0])
    assert gated["screen_score"][1] is None
    assert gated["screen_score"][2] is None
    assert gated["missing_confidence"].to_list() == [False, True, False]
    assert np.isfinite(gated["z_energy"].to_numpy()).all()
    covered = screen_scores(t, ref, weight=0.6, native_template=True, min_cohort=1)
    plain = screen_scores(t, ref, weight=0.6, min_cohort=1)
    np.testing.assert_allclose(covered["screen_score"].to_numpy(),
                               plain["screen_score"].to_numpy(), equal_nan=True)
    assert covered["screen_score"][1] is not None
    small = screen_scores(t, ref, weight=0.6, native_template=False, min_cohort=3)
    assert small["screen_score"].null_count() == 3


def test_template_metadata_is_explicit_hashed_and_validated(tmp_path):
    from tcren.screen import template_coverage
    import hashlib

    path = tmp_path / "metadata.tsv"
    path.write_text("epitope\tpdb.id\n aaa \tx\nAAA\ty\nBBB\tz\n")
    coverage = template_coverage("aaa", path)
    assert coverage["native_template"] is True
    assert coverage["n_epitopes"] == 2
    assert coverage["metadata_sha256"] == hashlib.sha256(path.read_bytes()).hexdigest()
    assert template_coverage("CCC", path)["native_template"] is False
    path.write_text("epitope\n\n")
    with pytest.raises(ValueError, match="nonempty"):
        template_coverage("AAA", path)


def test_raw_score_bounds_missingness_and_legacy_reference():
    ref = bundled_reference(CohortKind.TCR_RANKING)
    t = _table(4, iptm=[0., 1., None, .7], g_cyclo_frac=[0., 1., .3, None])
    np.testing.assert_allclose(screen_scores(t, ref, min_cohort=1)["screen_score"],
                               [0., 1., .3, .7])
    for column in ("iptm", "g_cyclo_frac"):
        with pytest.raises(ValueError, match=r"must lie in \[0, 1\]"):
            screen_scores(t.with_columns(pl.lit(1.1).alias(column)), ref, min_cohort=1)
    legacy = dict(ref)
    legacy.pop("score_method")
    got = screen_scores(_table(), legacy)
    expected = stouffer(np.column_stack([got["z_confidence"], got["z_topology"]]),
                        np.array([6 / 7, 1 / 7]))
    np.testing.assert_allclose(got["screen_score"], expected)
