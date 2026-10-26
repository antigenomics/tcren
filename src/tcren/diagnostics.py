"""Descriptive reductions of recorded structural observations, without fitting."""
from __future__ import annotations

import numpy as np
import polars as pl

from .evaluation import evaluate_scores


def backbone_extrema(rows: pl.DataFrame, keys=("id",)):
    """Retain all tied CA-bond/angle extrema for peptide and receptor regions."""
    keys = list(keys)
    selected = rows.filter(pl.col("metric").is_in(["ca_bond", "ca_angle"]) &
                           pl.col("chain.type").is_in(["TRA", "TRB", "PEPTIDE"]))
    if not selected.height:
        raise ValueError("No peptide/receptor CA bond or angle observations")
    selected = selected.with_columns(pl.when(pl.col("chain.type") == "PEPTIDE")
        .then(pl.lit("peptide")).otherwise(pl.col("region.type") +
            pl.when(pl.col("chain.type") == "TRA").then(pl.lit("alpha")).otherwise(pl.lit("beta"))).alias("region"))
    groups = keys + ["region", "metric", "units"]
    summary = selected.group_by(groups).agg(pl.col("value").min().alias("minimum"),
        pl.col("value").max().alias("maximum"), pl.len().alias("n_windows")).sort(groups)
    extremes = selected.join(summary, on=groups).filter(
        (pl.col("value") == pl.col("minimum")) | (pl.col("value") == pl.col("maximum")))
    features = summary.unpivot(on=["minimum", "maximum"], index=groups + ["n_windows"],
        variable_name="statistic", value_name="value").with_columns(
            pl.concat_str(["region", "metric", "statistic"], separator="__").alias("feature"))
    wide = features.select(*keys, "feature", "value").pivot(on="feature", index=keys, values="value").sort(keys)
    return {"geometry_summary": summary, "extreme_windows": extremes,
            "geometry_features": features, "geometry_wide": wide}


def grouped_summary(rows: pl.DataFrame, columns, groups):
    """Finite counts, distinct counts and linear quartiles, including missing-label groups."""
    groups = list(groups)
    records = []
    for key, part in rows.group_by(groups, maintain_order=True):
        context = dict(zip(groups, key))
        for column in columns:
            values = np.asarray(part[column], float)
            finite = values[np.isfinite(values)]
            records.append({**context, "feature": column, "n_models": part.height,
                "n_finite": len(finite), "n_unique": len(np.unique(finite)),
                "median": float(np.median(finite)) if len(finite) else None,
                "q25": float(np.quantile(finite, .25)) if len(finite) else None,
                "q75": float(np.quantile(finite, .75)) if len(finite) else None})
    return pl.DataFrame(records, infer_schema_length=None)


def context_correlations(rows, features, contexts, groups, keys=("id",), reference="iptm"):
    """Existing pairwise finite geometry/context correlations, stratified by supplied labels.

    Missing grouping labels remain a separate group. Each feature and reference use the
    same finite context/feature/reference rows. This delegates metrics to evaluate_scores.
    """
    groups = list(groups)
    records = []
    for key, part in rows.group_by(groups, maintain_order=True):
        group = dict(zip(groups, key))
        for feature in features:
            for context in contexts:
                result = evaluate_scores(part, {feature: 1}, keys=keys, endpoint=context,
                    kind="continuous", mask="pair", reference=reference)["metrics"]
                pearson = result.filter(pl.col("metric") == "pearson_R").row(0, named=True)
                spearman = result.filter(pl.col("metric") == "spearman_R").row(0, named=True)
                records.append({**group, "feature": feature, "context": context,
                    "n_total": part.height, "n_matched": pearson["n_used"],
                    "pearson_R": pearson["value"], "spearman_rho": spearman["value"],
                    "iptm_R": pearson["baseline"], "iptm_rho": spearman["baseline"]})
    return pl.DataFrame(records, infer_schema_length=None)

