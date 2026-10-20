"""Hand-count checks for opportunity-aware, receptor-clustered contrasts."""

import math

import polars as pl
import pytest

from tcren.contacts.class_statistics import class_statistics


def metadata(labels, *, clusters=None, cohorts=None, coverage=None, epitopes=None, dataset="vdjdb"):
    n = len(labels)
    return pl.DataFrame({
        "id": [str(i) for i in range(n)], "label": labels,
        "receptor_cluster": clusters or [str(i) for i in range(n)],
        "dataset": [dataset] * n, "cohort": cohorts or ["case"] * n,
        "epitope": epitopes or ["pep"] * n, "mhc": ["A*01:01"] * n,
        "mhc_class": ["I"] * n, "template_coverage": coverage or ["covered"] * n,
    })


def observations(counts, opportunities, feature="f"):
    return pl.DataFrame({"id": [str(i) for i in range(len(counts))],
                         "feature": [feature] * len(counts), "C": counts, "O": opportunities})


def test_equal_complex_density_counts_and_cluster_variance():
    per, overall = class_statistics(observations([1, 0, 0, 0], [1, 99, 1, 1]), metadata([1, 1, 0, 0]))
    row = per.row(0, named=True)
    assert row["mean_density_binder"] == 0.5
    assert row["sum_C_binder"] == 1
    assert row["sum_O_binder"] == 100
    assert row["sum_density_binder"] == 1
    assert row["n_contact_binder"] == 1
    assert row["prevalence_all_binder"] == 0.5
    assert row["cluster_se"] == pytest.approx(0.5)
    assert row["cluster_df"] == 3
    assert row["residual_df"] == 2
    assert row["cr1_factor"] == 2
    assert row["n_clusters_binder"] == row["n_clusters_nonbinder"] == 2
    assert row["log2_rate_ratio"] == math.inf
    assert row["log2_rate_ratio_reason"] == "nonbinder_zero"
    assert row["set_membership"] == "binder_only"
    assert overall["delta_mean_density"][0] == 0.5


def test_repeated_clusters_and_unequal_class_sizes():
    obs = observations([1, 0, 1, 2], [1, 2, 2, 2])
    per, _ = class_statistics(obs, metadata([1, 0, 0, 0]))
    assert per["mean_density_nonbinder"][0] == 0.5
    assert per["cluster_se"][0] is None
    assert per["n_clusters_binder"][0] == 1
    assert per["inference_status"][0] == "unreplicated_class"
    per, _ = class_statistics(observations([1, 1, 0, 1, 2], [1, 2, 2, 2, 2]),
                              metadata([1, 1, 0, 0, 0]))
    assert per["cluster_se"][0] == pytest.approx(math.sqrt(125 / 864))
    per, _ = class_statistics(observations([1, 0, 0, 0], [1] * 4),
                              metadata([1, 1, 0, 0], clusters=["a", "a", "b", "c"]))
    assert per["n_clusters"][0] == 3
    assert per["cluster_se"][0] is None
    assert per["p_two_sided"][0] is None
    assert per["inference_status"][0] == "unreplicated_class"


def test_absent_opportunities_and_eligible_zero_contacts():
    per, _ = class_statistics(observations([0, 0, 0], [2, 0, 3]), metadata([1, 1, 0]))
    row = per.row(0, named=True)
    assert row["n_models_binder"] == 2
    assert row["n_eligible_binder"] == 1
    assert row["n_no_opportunity_binder"] == 1
    assert row["mean_density_binder"] == 0
    assert row["log2_rate_ratio"] is None
    assert row["log2_rate_ratio_reason"] == "both_zero"
    assert row["set_membership"] == "neither"
    assert row["p_two_sided"] is None
    assert row["ci95_low"] is row["ci95_high"] is None
    assert row["inference_status"] == "unreplicated_class"
    # Omitting an O=0 row has exactly the same estimand and population.
    omitted, _ = class_statistics(observations([0, 0, 0], [2, 0, 3]).filter(pl.col("O") > 0),
                                  metadata([1, 1, 0]))
    assert omitted.equals(per)


@pytest.mark.parametrize("labels,counts,opportunities", [([1, 1], [1, 0], [1, 1]),
                                                         ([1, 0], [1, 0], [1, 0])])
def test_absent_eligible_class(labels, counts, opportunities):
    per, overall = class_statistics(observations(counts, opportunities), metadata(labels))
    assert per["delta_mean_density"][0] is None
    assert per["cluster_se"][0] is None
    assert per["inference_status"][0] == "eligible_class_absent"
    assert overall["n_strata_used"][0] == 0
    assert overall["all_strata_n_models_binder"][0] == labels.count(1)
    assert overall["n_models_binder"][0] == 0


def test_single_cluster_and_negative_infinite_ratio():
    per, _ = class_statistics(observations([0, 1], [1, 1]),
                              metadata([1, 0], clusters=["a", "a"]))
    assert per["log2_rate_ratio"][0] == -math.inf
    assert per["inference_status"][0] == "fewer_than_two_clusters"
    assert per["cluster_se"][0] is None


def test_overall_clusters_sum_before_variance_and_equal_stratum_weights():
    meta = metadata([1, 1, 0, 0] * 2, clusters=["a", "b", "c", "d"] * 2,
                    cohorts=["one"] * 4 + ["two"] * 4)
    per, overall = class_statistics(observations([1, 0, 0, 0, 0, 1, 0, 0], [1] * 8), meta)
    assert per["cluster_se"].to_list() == pytest.approx([0.5] * 2)
    row = overall.row(0, named=True)
    assert row["n_strata_used"] == 2
    assert row["n_clusters"] == 4
    assert row["delta_mean_density"] == 0.5
    assert row["residual_df"] == 4
    assert row["cr1_factor"] == pytest.approx(7 / 3)
    assert row["cluster_se"] == 0
    assert row["inference_status"] == "degenerate_variance"
    assert row["p_two_sided"] is None
    assert row["ci95_low"] is row["ci95_high"] is None
    # A large second case receives the same stratum weight as the small first.
    meta = metadata([1, 0, 1, 1, 0, 0], cohorts=["one"] * 2 + ["two"] * 4)
    _, overall = class_statistics(observations([1, 0, 0, 0, 0, 0], [1] * 6), meta)
    assert overall["mean_density_binder"][0] == 0.5


def test_no_pooling_coverage_tcrvdb_epitopes_or_mhc_alleles_within_strata():
    meta = metadata([1, 0] * 3, coverage=["covered"] * 4 + ["scarce"] * 2,
                    epitopes=["a"] * 2 + ["b"] * 4, dataset="TCRvdb")
    per, overall = class_statistics(observations([1, 0] * 3, [1] * 6), meta)
    assert per.height == overall.height == 3
    changed = meta.with_columns(pl.when(pl.col("id") == "0").then(pl.lit("A*02:01"))
                               .otherwise(pl.col("mhc")).alias("mhc"))
    per, _ = class_statistics(observations([1, 0] * 3, [1] * 6), changed)
    assert per.height == 4


@pytest.mark.parametrize("kind", ["negative", "excess", "fraction", "duplicate", "unknown", "label"])
def test_invalid_inputs_fail(kind):
    obs, meta = observations([1, 0], [1, 1]), metadata([1, 0])
    if kind == "negative":
        obs = obs.with_columns(pl.lit(-1).alias("C"))
    elif kind == "excess":
        obs = obs.with_columns(pl.lit(2).alias("C"))
    elif kind == "fraction":
        obs = obs.with_columns(pl.col("C").cast(pl.Float64))
    elif kind == "duplicate":
        obs = pl.concat([obs, obs.head(1)])
    elif kind == "unknown":
        obs = obs.with_columns(pl.lit("missing").alias("id"))
    else:
        meta = meta.with_columns(pl.lit(2).alias("label"))
    with pytest.raises(ValueError):
        class_statistics(obs, meta)


def test_finite_ratio_and_t_interval():
    from scipy.stats import t

    per, _ = class_statistics(observations([1, 1, 0, 1], [1, 2, 1, 1]),
                              metadata([1, 1, 0, 0]))
    row = per.row(0, named=True)
    assert row["mean_density_binder"] == 0.75
    assert row["mean_density_nonbinder"] == 0.5
    assert row["log2_rate_ratio"] == pytest.approx(math.log2(1.5))
    assert row["log2_rate_ratio_reason"] == "finite"
    assert row["set_membership"] == "shared"
    se = math.sqrt(5 / 16)
    assert row["cluster_se"] == pytest.approx(se)
    assert row["p_two_sided"] == pytest.approx(2 * t.sf(0.25 / se, 3))
    assert row["ci95_low"] == pytest.approx(0.25 - t.ppf(0.975, 3) * se)
    assert row["ci95_high"] == pytest.approx(0.25 + t.ppf(0.975, 3) * se)


@pytest.mark.parametrize("counts", [[1, 0], [0, 0]])
def test_single_model_per_class_never_implies_certainty(counts):
    per, overall = class_statistics(observations(counts, [1, 1]), metadata([1, 0]))
    for frame in (per, overall):
        row = frame.row(0, named=True)
        assert row["delta_mean_density"] == counts[0] - counts[1]
        assert row["n_clusters_binder"] == row["n_clusters_nonbinder"] == 1
        assert row["residual_df"] == 0
        assert not row["class_replication_estimable"]
        assert not row["variance_estimable"]
        for name in ("cluster_se", "t_statistic", "p_two_sided", "ci95_low", "ci95_high"):
            assert row[name] is None


@pytest.mark.parametrize("counts", [[1, 1, 0, 0], [0, 0, 0, 0]])
def test_replicated_constant_density_has_undefined_inference(counts):
    per, overall = class_statistics(observations(counts, [1] * 4), metadata([1, 1, 0, 0]))
    for frame in (per, overall):
        row = frame.row(0, named=True)
        assert row["class_replication_estimable"]
        assert row["variance_estimable"]
        assert row["cluster_se"] == 0
        assert row["inference_status"] == "degenerate_variance"
        for name in ("t_statistic", "p_two_sided", "ci95_low", "ci95_high"):
            assert row[name] is None


def test_overall_unreplicated_stratum_keeps_point_estimate_without_inference():
    meta = metadata([1, 0, 1, 1, 0, 0], cohorts=["one"] * 2 + ["two"] * 4)
    per, overall = class_statistics(observations([1, 0, 1, 0, 0, 1], [1] * 6), meta)
    assert per["class_replication_estimable"].to_list() == [False, True]
    row = overall.row(0, named=True)
    assert row["n_strata_used"] == 2
    assert row["n_strata_unreplicated_class"] == 1
    assert row["delta_mean_density"] == 0.5
    assert row["n_models_binder"] == row["n_models_nonbinder"] == 3
    assert row["n_clusters_binder"] == row["n_clusters_nonbinder"] == 3
    assert row["residual_df"] == 2
    assert not row["class_replication_estimable"]
    assert row["inference_status"] == "unreplicated_class"
    for name in ("cluster_se", "t_statistic", "p_two_sided", "ci95_low", "ci95_high"):
        assert row[name] is None


def test_no_eligible_observations_preserves_population_and_undefined_estimates():
    per, overall = class_statistics(observations([0, 0], [0, 0]), metadata([1, 0]))
    row = per.row(0, named=True)
    assert row["n_no_opportunity_binder"] == row["n_no_opportunity_nonbinder"] == 1
    assert row["n_clusters_binder"] == row["n_clusters_nonbinder"] == 0
    assert row["delta_mean_density"] is None
    assert row["cluster_se"] is row["p_two_sided"] is None
    assert row["inference_status"] == "eligible_class_absent"
    assert overall["n_strata_used"][0] == 0
    assert overall["all_strata_n_models_binder"][0] == 1
