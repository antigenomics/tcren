"""Equal-complex contact-density contrasts with receptor-clustered uncertainty.

Sparse observations contain one row per (id, feature), including every O>0 row
with C=0. O=0,C=0 rows may be omitted. Features are opaque identifiers. The
output-sized stratum/feature grid never expands to model/feature combinations.
"""

from __future__ import annotations

import numpy as np
import polars as pl
from scipy.stats import t as student_t

STRATUM = ["dataset", "cohort", "epitope", "mhc", "mhc_class", "template_coverage"]
OVERALL = ["dataset", "mhc_class", "template_coverage", "overall_epitope"]
COUNTS = ["n_models", "n_eligible", "n_contact", "n_no_opportunity", "sum_C", "sum_O", "sum_density"]


def _inference(frame: pl.DataFrame) -> pl.DataFrame:
    """Vectorized t inference on reduced rows, assuming independent clusters."""
    frame = frame.with_columns(
        (pl.col("mean_density_binder") - pl.col("mean_density_nonbinder"))
        .alias("delta_mean_density"))
    n = (frame["n_eligible_binder"].cast(pl.Int64)
         + frame["n_eligible_nonbinder"].cast(pl.Int64)).to_numpy()
    k = 2 * frame["n_strata_used"].to_numpy() if "n_strata_used" in frame.columns else 2
    residual_df = n - k
    frame = frame.with_columns(pl.Series("residual_df", residual_df))
    delta = frame["delta_mean_density"].cast(pl.Float64).to_numpy()
    g = frame["n_clusters"].to_numpy()
    ss = frame["cluster_if_ss"].to_numpy()
    valid = (np.isfinite(delta) & (g >= 2) & (residual_df > 0)
             & frame["class_replication_estimable"].to_numpy())
    factor = np.full(frame.height, np.nan)
    factor[valid] = g[valid] / (g[valid] - 1) * (n[valid] - 1) / residual_df[valid]
    se = np.full(frame.height, np.nan)
    df = np.full(frame.height, np.nan)
    statistic = np.full(frame.height, np.nan)
    p = np.full(frame.height, np.nan)
    margin = np.full(frame.height, np.nan)
    df[valid] = g[valid] - 1
    se[valid] = np.sqrt(np.maximum(0, factor[valid] * ss[valid]))
    nonzero = valid & (se > 0)
    statistic[nonzero] = delta[nonzero] / se[nonzero]
    p[nonzero] = 2 * student_t.sf(np.abs(statistic[nonzero]), df[nonzero])
    margin[nonzero] = student_t.ppf(0.975, df[nonzero]) * se[nonzero]
    frame = frame.with_columns(pl.Series("variance_estimable", valid))
    frame = frame.with_columns([
        pl.Series(name, values, nan_to_null=True)
        for name, values in (("cluster_se", se), ("cluster_df", df), ("cr1_factor", factor),
                             ("t_statistic", statistic), ("p_two_sided", p),
                             ("ci95_low", delta - margin), ("ci95_high", delta + margin))
    ])
    pos, neg = pl.col("mean_density_binder"), pl.col("mean_density_nonbinder")
    absent = pl.col("delta_mean_density").is_null()
    return frame.with_columns(
        pl.when(absent | ((pos == 0) & (neg == 0))).then(None)
        .otherwise((pos / neg).log(2)).alias("log2_rate_ratio"),
        pl.when(absent).then(pl.lit("eligible_class_absent"))
        .when((pos == 0) & (neg == 0)).then(pl.lit("both_zero"))
        .when(neg == 0).then(pl.lit("nonbinder_zero"))
        .when(pos == 0).then(pl.lit("binder_zero"))
        .otherwise(pl.lit("finite")).alias("log2_rate_ratio_reason"),
        pl.when((pl.col("n_contact_binder") > 0) & (pl.col("n_contact_nonbinder") > 0))
        .then(pl.lit("shared")).when(pl.col("n_contact_binder") > 0)
        .then(pl.lit("binder_only")).when(pl.col("n_contact_nonbinder") > 0)
        .then(pl.lit("nonbinder_only")).otherwise(pl.lit("neither")).alias("set_membership"),
        pl.when(absent).then(pl.lit("eligible_class_absent"))
        .when(pl.col("n_clusters") < 2).then(pl.lit("fewer_than_two_clusters"))
        .when(~pl.col("class_replication_estimable")).then(pl.lit("unreplicated_class"))
        .when(pl.col("residual_df") <= 0).then(pl.lit("no_residual_df"))
        .when(pl.col("cluster_se") == 0).then(pl.lit("degenerate_variance"))
        .otherwise(pl.lit("finite")).alias("inference_status"),
        *[
            pl.when(pl.col(f"n_models_{name}") > 0)
            .then(pl.col(f"n_contact_{name}") / pl.col(f"n_models_{name}"))
            .otherwise(None).alias(f"prevalence_all_{name}")
            for name in ("binder", "nonbinder")
        ], *[
            pl.when(pl.col(f"n_eligible_{name}") > 0)
            .then(pl.col(f"n_contact_{name}") / pl.col(f"n_eligible_{name}"))
            .otherwise(None).alias(f"prevalence_eligible_{name}")
            for name in ("binder", "nonbinder")
        ])


def _cluster_summary(influences: pl.DataFrame, keys: list[str]) -> pl.DataFrame:
    summary = (influences.group_by(keys + ["receptor_cluster"])
               .agg(pl.col("influence").sum())
               .group_by(keys).agg(pl.len().alias("n_clusters"),
                                   pl.col("influence").pow(2).sum().alias("cluster_if_ss")))
    counts = influences.group_by(keys).agg([
        pl.col("receptor_cluster").filter(pl.col("label") == label).n_unique()
        .alias(f"n_clusters_{name}")
        for label, name in ((1, "binder"), (0, "nonbinder"))
    ])
    return summary.join(counts, on=keys)


def class_statistics(
    observations: pl.DataFrame, metadata: pl.DataFrame,
) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Return per-stratum and equal-stratum overall contact-density statistics.

    Metadata has unique id, label (0/1), receptor_cluster and STRATUM columns.
    Observations have unique id/feature, integer C/O with 0<=C<=O. Missing rows
    mean zero opportunities, never an unrecorded eligible zero. Density is C/O;
    its class mean weights eligible complexes equally. log2_rate_ratio compares
    these means without smoothing. Observed set membership is descriptive, not
    evidence that a contact is physically required or forbidden.

    Cluster influence sums use IF=(d-mean)/n in binders and its negative in
    nonbinders. CR1 variance is G/(G-1)*(N-1)/(N-K)*sum(IF_cluster**2),
    where G counts distinct eligible receptor clusters, N counts eligible model
    observations and K is 2 per stratum or twice n_strata_used overall. Analytic
    inference assumes independent clusters and uses an approximate t reference
    with G-1 df, not an exact finite-sample test. Each contributing stratum must
    have at least two distinct eligible receptor clusters in EACH class; no
    variance pooling assumption is added for unreplicated strata. Residual df
    N-K must be positive. Zero variance retains SE=0 but no t, p or CI. There is
    no multiplicity correction.

    Overall means equally weight strata eligible in both classes. Influences
    are weighted before summing across globally shared receptor clusters.
    Overall count fields describe used strata; all_strata_* retain excluded
    populations. Dataset, MHC class and coverage remain separate; TCRvdb
    epitopes also remain separate. Exact alleles remain separate within strata.
    """
    required = ["id", "label", "receptor_cluster", *STRATUM]
    if not set(required).issubset(metadata.columns):
        raise ValueError(f"metadata requires {required}")
    if not {"id", "feature", "C", "O"}.issubset(observations.columns):
        raise ValueError("observations requires id, feature, C, O")
    meta = metadata.select(required)
    obs = observations.select("id", "feature", "C", "O")
    if meta.select(pl.any_horizontal(pl.all().is_null()).any()).item():
        raise ValueError("metadata values must be nonnull")
    if meta["id"].n_unique() != meta.height or not set(meta["label"].unique()).issubset({0, 1}):
        raise ValueError("metadata ids must be unique and labels binary 0/1")
    if obs.select(pl.any_horizontal(pl.all().is_null()).any()).item():
        raise ValueError("observations must be nonnull")
    if any(not obs.schema[c].is_integer() for c in ("C", "O")):
        raise ValueError("C and O must have integer dtype")
    if obs.filter((pl.col("C") < 0) | (pl.col("O") < pl.col("C"))).height:
        raise ValueError("counts require 0 <= C <= O")
    if obs.select("id", "feature").unique().height != obs.height:
        raise ValueError("duplicate id/feature observations")
    if obs.join(meta.select("id"), on="id", how="anti").height:
        raise ValueError("observations contain unknown model ids")
    if obs.schema["feature"] != pl.String:
        raise ValueError("feature must be a string identifier")
    if obs.filter(pl.col("feature").str.len_chars() == 0).height:
        raise ValueError("feature identifiers must be nonempty")
    if obs.is_empty() or meta.is_empty():
        raise ValueError("at least one model and one feature are required")

    keys = [*STRATUM, "feature"]
    grid = meta.select(STRATUM).unique().join(obs.select("feature").unique(), how="cross")
    totals = meta.group_by([*STRATUM, "label"]).agg(pl.len().alias("n_models"))
    eligible = obs.filter(pl.col("O") > 0).join(meta, on="id", how="left", validate="m:1")
    eligible = eligible.with_columns((pl.col("C") / pl.col("O")).alias("density"))
    means = eligible.group_by([*keys, "label"]).agg(
        pl.len().alias("n_eligible"), (pl.col("C") > 0).sum().alias("n_contact"),
        pl.col("C").sum().alias("sum_C"), pl.col("O").sum().alias("sum_O"),
        pl.col("density").sum().alias("sum_density"),
        pl.col("density").mean().alias("mean_density"))
    per = grid
    for label, name in ((1, "binder"), (0, "nonbinder")):
        side = (grid.join(totals.filter(pl.col("label") == label).drop("label"), on=STRATUM, how="left")
                .join(means.filter(pl.col("label") == label).drop("label"), on=keys, how="left")
                .with_columns(pl.col("n_models", "n_eligible", "n_contact", "sum_C", "sum_O", "sum_density").fill_null(0))
                .with_columns((pl.col("n_models") - pl.col("n_eligible")).alias("n_no_opportunity")))
        per = per.join(side.select(keys + [pl.col(c).alias(f"{c}_{name}")
                                          for c in [*COUNTS, "mean_density"]]), on=keys)
    influences = eligible.join(means.select([*keys, "label", "n_eligible", "mean_density"]),
                               on=[*keys, "label"], validate="m:1").with_columns(
        ((pl.col("density") - pl.col("mean_density")) / pl.col("n_eligible")
         * pl.when(pl.col("label") == 1).then(1).otherwise(-1)).alias("influence"))
    per = per.join(_cluster_summary(influences, keys), on=keys, how="left").with_columns(
        pl.col("n_clusters", "n_clusters_binder", "n_clusters_nonbinder").fill_null(0),
        pl.col("cluster_if_ss").fill_null(0.0))
    per = per.with_columns(
        pl.when(pl.col("dataset").str.to_lowercase() == "tcrvdb")
        .then(pl.col("epitope")).otherwise(pl.lit("all")).alias("overall_epitope"),
        ((pl.col("n_eligible_binder") > 0) & (pl.col("n_eligible_nonbinder") > 0)).alias("used"),
        ((pl.col("n_clusters_binder") >= 2) & (pl.col("n_clusters_nonbinder") >= 2))
        .alias("class_replication_estimable"))
    group = [*OVERALL, "feature"]
    count_cols = [f"{c}_{name}" for name in ("binder", "nonbinder") for c in COUNTS]
    overall = per.group_by(group).agg(
        pl.len().alias("n_strata_total"), pl.col("used").sum().alias("n_strata_used"),
        (pl.col("used") & ~pl.col("class_replication_estimable")).sum()
        .alias("n_strata_unreplicated_class"),
        *[pl.col(c).sum().alias(f"all_strata_{c}") for c in count_cols],
        *[pl.col(c).filter(pl.col("used")).sum().alias(c) for c in count_cols],
        *[pl.col(f"mean_density_{name}").filter(pl.col("used")).mean()
          for name in ("binder", "nonbinder")])
    combined_if = (influences.join(per.select(keys + ["overall_epitope", "used"]), on=keys)
                   .filter(pl.col("used"))
                   .join(overall.select(group + ["n_strata_used"]), on=group)
                   .with_columns((pl.col("influence") / pl.col("n_strata_used")).alias("influence")))
    overall = overall.join(_cluster_summary(combined_if, group), on=group, how="left").with_columns(
        pl.col("n_clusters", "n_clusters_binder", "n_clusters_nonbinder").fill_null(0),
        pl.col("cluster_if_ss").fill_null(0.0))
    overall = overall.with_columns(
        ((pl.col("n_strata_used") > 0) & (pl.col("n_strata_unreplicated_class") == 0))
        .alias("class_replication_estimable"))
    return (_inference(per.drop("used").sort(keys)), _inference(overall.sort(group)))
