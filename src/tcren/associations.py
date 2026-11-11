"""Descriptive pMHC feature associations with paired score performance.

No fitting, calibration, significance tests, label-selected directions or pooling
of MHC classes/template strata. Recorded scores remain immutable inputs.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import polars as pl
import typer
from scipy.stats import rankdata

STRATA = ["dataset", "task", "mhc_class", "template_coverage"]
GROUPS = STRATA + ["cohort"]


def add_feature_products(features: pl.DataFrame, pairs: list[tuple[str, str]]) -> tuple[pl.DataFrame, list[str]]:
    """Append raw signed per-row products, without labels, centering or fitting.

    Names are ``left__x__right``. Sources must be distinct existing numeric
    columns; nonfinite sources or overflowing products produce missing values.
    These are joint scalar descriptors, not independent interaction estimates
    or spatial co-localization measurements. Units multiply with the sources.
    """
    names, expressions = [], []
    for left, right in pairs:
        name = f"{left}__x__{right}"
        if left == right or name in names or name in features.columns:
            raise ValueError("feature products require distinct sources and unique new names")
        if any(c not in features.columns or not features.schema[c].is_numeric() for c in (left, right)):
            raise ValueError("feature product sources must be existing numeric columns")
        product = pl.col(left).cast(pl.Float64) * pl.col(right).cast(pl.Float64)
        valid = pl.col(left).is_finite() & pl.col(right).is_finite() & product.is_finite()
        expressions.append(pl.when(valid).then(product).otherwise(None).alias(name))
        names.append(name)
    return features.with_columns(expressions), names


def _unique(table, columns, name):
    if any(c not in table.columns for c in columns):
        raise ValueError(f"{name}: required columns {columns}")
    if table.select(pl.any_horizontal(pl.col(columns).is_null()).any()).item():
        raise ValueError(f"{name}: null keys")
    if table.select(columns).is_duplicated().any():
        raise ValueError(f"{name}: duplicate keys")


def _correlation(x, y):
    if len(x) < 3:
        return None, "fewer_than_three_cohorts"
    if np.ptp(x) == 0 or np.ptp(y) == 0:
        return None, "constant"
    return float(np.corrcoef(rankdata(x), rankdata(y))[0, 1]), "ok"


def _metrics(y, x):
    from sklearn.metrics import average_precision_score, precision_recall_curve, roc_auc_score

    names = ("roc_auc", "average_precision", "precision_at_recall_0.1", "precision_at_recall_0.5")
    if len(np.unique(y)) != 2:
        return dict.fromkeys(names)
    p, r, _ = precision_recall_curve(y, x)
    return dict(zip(names, [float(roc_auc_score(y, x)), float(average_precision_score(y, x)),
                           float(p[r >= .1].max()), float(p[r >= .5].max())]))


def associate(results: pl.DataFrame, features: pl.DataFrame, scores: dict[str, int],
              feature_columns: list[str], *, endpoint: str = "y") -> dict[str, pl.DataFrame]:
    """Compare each score with ipTM on its own matched finite rows, then associate.

    ``results`` requires globally unique ``id``, dataset/task/cohort/mhc_class/
    template_coverage, binary or missing endpoint, iptm and explicit score columns. Classes
    are I/II; missing coverage must be explicitly ``unknown``. ``features`` has
    unique id plus selected numeric columns, one pMHC conformation per input row.
    Missing features never remove classification rows. Cohort feature medians are
    computed once using all models, label-blind, with finite counts and quartiles.
    Missing outcomes are excluded only from classification, with membership flags. Correlations use
    one cohort per observation, separately within each of the four strata keys.

    Returns metrics, features, associations and row-level membership tables.
    """
    if not results.height or not scores or not feature_columns:
        raise ValueError("nonempty results, scores and feature columns are required")
    _unique(results, ["id"], "results")
    _unique(features, ["id"], "features")
    required = GROUPS + [endpoint, "iptm", *scores]
    missing = set(required) - set(results.columns)
    if missing or set(feature_columns) - set(features.columns):
        raise ValueError(f"missing required result/feature columns: {sorted(missing)}")
    if any(v not in (-1, 1) for v in scores.values()):
        raise ValueError("each score direction must be +1 or -1")
    if set(scores) & {"pitch", "pitch_angle"}:
        raise ValueError("pitch is excluded as generator-confidence leakage")
    if len(set(feature_columns)) != len(feature_columns) or set(feature_columns) & set(results.columns):
        raise ValueError("feature columns must be distinct and not overlap results")
    if results.select(pl.any_horizontal(pl.col(GROUPS).is_null()).any()).item():
        raise ValueError("null grouping keys; encode unknown coverage explicitly")
    if set(results["mhc_class"]) - {"I", "II"}:
        raise ValueError("mhc_class must be I or II")
    if features.join(results.select("id"), on="id", how="anti").height:
        raise ValueError("features contain unexpected IDs")
    labels = np.asarray(results[endpoint], float)
    if np.isinf(labels).any() or not np.isin(labels[np.isfinite(labels)], [0, 1]).all():
        raise ValueError("endpoint must contain binary labels or missing values")
    confidence = np.asarray(results["iptm"], float)
    if ((confidence[np.isfinite(confidence)] < 0) | (confidence[np.isfinite(confidence)] > 1)).any():
        raise ValueError("finite iptm values must be in [0, 1]")
    joined = results.join(features.select("id", *feature_columns), on="id", how="left", validate="1:1")
    summaries, metrics, membership = [], [], []
    for key, part in joined.sort(GROUPS + ["id"]).group_by(GROUPS, maintain_order=True):
        context = dict(zip(GROUPS, key))
        for feature in feature_columns:
            x = np.asarray(part[feature], float)
            x = x[np.isfinite(x)]
            summaries.append({**context, "feature": feature, "n_rows": part.height,
                "n_feature": len(x), "n_missing": part.height-len(x),
                "median": float(np.median(x)) if len(x) else None,
                "q25": float(np.quantile(x, .25)) if len(x) else None,
                "q75": float(np.quantile(x, .75)) if len(x) else None})
        y, confidence = np.asarray(part[endpoint], float), np.asarray(part["iptm"], float)
        for score, direction in scores.items():
            x = np.asarray(part[score], float) * direction
            ok = np.isfinite(y) & np.isfinite(x) & np.isfinite(confidence)
            membership.append(part.select("id", *GROUPS).with_columns(
                pl.lit(score).alias("score"), pl.Series("used", ok),
                pl.Series("finite_endpoint", np.isfinite(y)),
                pl.Series("finite_score", np.isfinite(x)), pl.Series("finite_iptm", np.isfinite(confidence))))
            current, baseline = _metrics(y[ok], x[ok]), _metrics(y[ok], confidence[ok])
            for metric, value in current.items():
                ref = baseline[metric]
                metrics.append({**context, "endpoint": endpoint, "score": score, "direction": direction,
                    "metric": metric, "n_rows": part.height, "n_used": int(ok.sum()),
                    "n_excluded": int((~ok).sum()), "n_positive": int(y[ok].sum()),
                    "n_negative": int(ok.sum()-y[ok].sum()), "value": value, "iptm_value": ref,
                    "delta": value-ref if value is not None else None,
                    "status": "ok" if value is not None else "fewer_than_two_classes"})
    metrics = pl.DataFrame(metrics, infer_schema_length=None)
    summaries = pl.DataFrame(summaries, infer_schema_length=None)
    association_rows = []
    combined = metrics.join(summaries, on=GROUPS, how="inner")
    keys = STRATA + ["score", "direction", "metric", "feature"]
    for key, part in combined.group_by(keys, maintain_order=True):
        context = dict(zip(keys, key))
        for target in ("value", "iptm_value", "delta"):
            x, y = np.asarray(part["median"], float), np.asarray(part[target], float)
            ok = np.isfinite(x) & np.isfinite(y)
            rho, status = _correlation(x[ok], y[ok])
            association_rows.append({**context, "target": target, "spearman_rho": rho,
                "n_cohorts": int(ok.sum()), "n_cohorts_total": part.height,
                "n_cohorts_missing": int((~ok).sum()), "status": status})
    return {"metrics": metrics, "features": summaries,
            "associations": pl.DataFrame(association_rows, infer_schema_length=None),
            "membership": pl.concat(membership)}


def read_table(path):
    """Read a CSV, TSV or Parquet table without sampling the schema."""
    path = Path(path)
    return pl.read_parquet(path) if path.suffix == ".parquet" else pl.read_csv(
        path, separator="\t" if path.suffix == ".tsv" else ",", infer_schema_length=None)


def association_command(
    results: Path = typer.Option(..., "--results", exists=True, dir_okay=False, help="Recorded scores, labels and grouping keys"),
    features: Path = typer.Option(..., "--features", exists=True, dir_okay=False, help="pMHC features keyed by id"),
    scores: Path = typer.Option(..., "--scores", exists=True, dir_okay=False, help="JSON mapping score columns to fixed +1/-1 directions"),
    out: Path = typer.Option(..., "-o", "--out", help="New output directory; never overwrite results"),
    feature_columns: str = typer.Option(..., "--feature-columns", help="Comma-separated pMHC feature columns"),
    endpoint: str = typer.Option("y", "--endpoint", help="Existing binary endpoint; labels are never redefined"),
    feature_product: list[str] = typer.Option(None, "--feature-product", help="Repeat LEFT:RIGHT for raw signed per-model feature products; no centering or fitting"),
) -> None:
    """Associate pMHC features with recorded score performance versus ipTM.

    RESULTS/FEATURES: CSV, TSV or Parquet. SCORES: JSON {column: +1|-1}.
    OUT must be a new directory. FEATURE_COLUMNS is a comma-separated list.
    Required grouping columns: dataset, task, cohort, mhc_class (I/II),
    template_coverage (use unknown explicitly). No fitting or significance tests.
    """
    from .exports import _sha256
    from . import __version__
    from .provenance import registry_digest

    if out.exists():
        raise typer.BadParameter("output already exists; choose a new directory")
    try:
        directions = json.loads(scores.read_text())
        if not isinstance(directions, dict):
            raise ValueError("scores JSON must be an object mapping columns to directions")
        pairs = []
        for pair in feature_product or []:
            parts = [c.strip() for c in pair.split(":")]
            if len(parts) != 2 or not all(parts):
                raise ValueError("feature products must be LEFT:RIGHT")
            pairs.append(tuple(parts))
        joint, product_columns = add_feature_products(read_table(features), pairs)
        columns = [c.strip() for c in feature_columns.split(",") if c.strip()]
        tables = associate(read_table(results), joint, directions,
                           columns + product_columns, endpoint=endpoint)
    except (ValueError, OSError, pl.exceptions.PolarsError) as exc:
        raise typer.BadParameter(str(exc)) from exc
    out.mkdir(parents=True)
    for name, table in tables.items():
        table.write_csv(out / f"{name}.tsv", separator="\t")
    if pairs:
        joint.write_csv(out / "joint_features.tsv", separator="\t")
    identity = {"inputs": {str(p.resolve()): _sha256(p) for p in (results, features, scores)},
                "feature_columns": ",".join(columns + product_columns), "endpoint": endpoint,
                "feature_products": pairs,
                "directions": directions, "software": {"tcren": __version__, "registry_digest": registry_digest()},
                "feature_summary": "label-blind cohort median; bound-conformation models",
                "association": "descriptive Spearman; one cohort per observation within strata; no tests"}
    (out / "provenance.json").write_text(json.dumps(identity, indent=2, default=str) + "\n")
    typer.echo(f"Wrote {tables['metrics'].height} paired metrics and "
               f"{tables['associations'].height} associations to {out}")
    if not tables["metrics"].filter(pl.col("status") == "ok").height:
        typer.echo("No evaluable matched two-class cohorts", err=True)
        raise typer.Exit(1)
