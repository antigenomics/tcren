"""Established matched ranking/correlation readouts for recorded score tables."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import pearsonr, spearmanr
import typer

from .associations import _metrics, _unique, read_table


def metric_values(y, score, *, kind="binary", pearson_test=False):
    """Reuse inclusive PR-envelope precision and ordinary Pearson/Spearman definitions."""
    if kind == "binary":
        return _metrics(y, score)
    if len(y) < 3 or np.ptp(y) == 0 or np.ptp(score) == 0:
        return dict.fromkeys(["pearson_R", "spearman_R"] + (["pearson_P"] if pearson_test else []))
    correlation = pearsonr(y, score)
    result = {"pearson_R": float(correlation.statistic), "spearman_R": float(spearmanr(y, score).statistic)}
    if pearson_test:
        result["pearson_P"] = float(correlation.pvalue)
    return result


def evaluate_scores(rows, scores, *, groups=(), keys=("id",), endpoint="y", kind="binary",
                    mask="common", reference=None, baseline_scores=None, curves=False,
                    pearson_test=False):
    """Evaluate each declared group, with explicit common or score/reference finite masks.

    Missing outcomes remain in membership but cannot enter metrics. Continuous outcomes
    never depend on binary labels. No cohort pooling, resampling, fitting or new threshold.
    ``baseline_scores`` maps current columns to previous-score columns; their directions
    match the current score. ``reference`` instead compares every arm to one named score.
    """
    groups, keys = list(groups), list(keys)
    baseline_scores = baseline_scores or {}
    reserved = {"endpoint", "score", "direction", "metric", "n_rows", "n_used", "n_excluded",
                "n_positive", "n_negative", "value", "baseline", "delta", "status", "used",
                "finite_outcome", "finite_score", "recall", "precision", "fpr"}
    collisions = reserved.intersection([*groups, *keys])
    if collisions:
        raise ValueError(f"Grouping and identity columns collide with reserved outputs: {sorted(collisions)}")
    if not rows.height or not scores or not keys:
        raise ValueError("Nonempty rows, score directions and identity keys are required")
    if kind not in {"binary", "continuous"} or mask not in {"common", "pair"}:
        raise ValueError("kind must be binary/continuous and mask common/pair")
    if (curves and kind != "binary") or (pearson_test and kind != "continuous"):
        raise ValueError("Curves require binary outcomes; Pearson test requires continuous outcomes")
    if reference and baseline_scores:
        raise ValueError("Choose one reference or per-score baseline columns")
    if set(baseline_scores) - set(scores) or any(type(v) is not int or v not in (-1, 1) for v in scores.values()):
        raise ValueError("Baseline keys must be scored columns; directions must be integer +1/-1")
    required = list(dict.fromkeys([*groups, *keys, endpoint, *scores, *baseline_scores.values(), *([reference] if reference else [])]))
    if set(required) - set(rows.columns):
        raise ValueError(f"Missing columns: {sorted(set(required) - set(rows.columns))}")
    _unique(rows, list(dict.fromkeys([*groups, *keys])), "evaluation rows")
    if groups and rows.select(pl.any_horizontal(pl.col(groups).is_null()).any()).item():
        raise ValueError("Null grouping keys; encode unknown groups explicitly")
    numeric = list(dict.fromkeys([endpoint, *scores, *baseline_scores.values(), *([reference] if reference else [])]))
    rows = rows.with_columns(pl.col(name).cast(pl.Float64) for name in numeric if rows.schema[name] == pl.Null)
    if any(not rows.schema[name].is_numeric() for name in numeric):
        raise ValueError("Outcomes and scores must be numeric columns")
    labels = rows[endpoint].to_numpy()
    if kind == "binary" and not np.isin(labels[np.isfinite(labels)], [0, 1]).all():
        raise ValueError("Finite binary outcomes must be 0 or 1")
    metrics, membership, curve_rows = [], [], []
    ordered = rows.sort(list(dict.fromkeys([*groups, *keys])))
    partitions = ordered.group_by(groups, maintain_order=True) if groups else [((), ordered)]
    for group, part in partitions:
        context = dict(zip(groups, group))
        y = np.asarray(part[endpoint], float)
        finite = {name: np.isfinite(np.asarray(part[name], float)) for name in numeric}
        common = np.logical_and.reduce(list(finite.values()))
        for score, direction in scores.items():
            previous = reference or baseline_scores.get(score)
            used = common if mask == "common" else finite[endpoint] & finite[score]
            if mask == "pair" and previous:
                used &= finite[previous]
            yy = y[used]
            x = direction * np.asarray(part[score], float)[used]
            current = metric_values(yy, x, kind=kind, pearson_test=pearson_test)
            before = metric_values(yy, (scores.get(reference, 1) if reference else direction) *
                                   np.asarray(part[previous], float)[used], kind=kind,
                                   pearson_test=pearson_test) if previous else {}
            n = int(used.sum())
            counts = dict(n_rows=part.height, n_used=n, n_excluded=part.height-n,
                          n_positive=int(yy.sum()) if kind == "binary" else None,
                          n_negative=int((yy == 0).sum()) if kind == "binary" else None)
            for metric, value in current.items():
                baseline = before.get(metric)
                delta = value-baseline if value is not None and baseline is not None else None
                metrics.append({**context, "endpoint": endpoint, "score": score, "direction": direction,
                                "metric": metric, **counts, "value": value, "baseline": baseline,
                                "delta": delta, "status": "ok" if value is not None else "unavailable"})
            membership.append(part.select(list(dict.fromkeys([*groups, *keys]))).with_columns(
                pl.lit(score).alias("score"), pl.Series("used", used),
                pl.Series("finite_outcome", finite[endpoint]), pl.Series("finite_score", finite[score])))
            if curves and len(np.unique(yy)) == 2:
                # One row per distinct descending score: ties enter together, as in the
                # existing publication PR/ROC producer. Include the origin explicitly.
                order = np.argsort(-x, kind="stable")
                xx, yy_ordered = x[order], yy[order]
                ends = np.r_[np.flatnonzero(np.diff(xx) != 0), len(xx)-1]
                tp = np.cumsum(yy_ordered)[ends]
                fp = ends+1-tp
                points = zip(np.r_[0, tp/yy.sum()], np.r_[1, tp/(ends+1)], np.r_[0, fp/(len(yy)-yy.sum())])
                curve_rows.extend({**context, "score": score, **counts, "recall": float(r),
                                   "precision": float(p), "fpr": float(f)} for r, p, f in points)
    output = {"metrics": pl.DataFrame(metrics, infer_schema_length=None), "membership": pl.concat(membership)}
    if curves:
        output["curves"] = pl.DataFrame(curve_rows, infer_schema_length=None) if curve_rows else pl.DataFrame(
            schema={**{g: rows.schema[g] for g in groups}, "score": pl.String, "recall": pl.Float64,
                    "precision": pl.Float64, "fpr": pl.Float64})
    return output


def command(
    results: Path = typer.Option(..., exists=True, dir_okay=False),
    scores: Path = typer.Option(..., exists=True, dir_okay=False, help="JSON {score_column: +1|-1}"),
    output: Path = typer.Option(..., "--output", "-o"),
    groups: str = typer.Option("", help="Comma-separated grouping columns; omitted means one group"),
    keys: str = typer.Option("id", help="Comma-separated row identity columns"),
    endpoint: str = typer.Option("y"),
    kind: str = typer.Option("binary", help="binary or continuous"),
    mask: str = typer.Option("common", help="common or pair finite population"),
    reference: str | None = typer.Option(None),
    baseline_scores: Path | None = typer.Option(None, exists=True, dir_okay=False),
    curves: bool = typer.Option(False, help="Export existing tied-threshold PR/ROC points"),
    pearson_test: bool = typer.Option(False, help="Explicit two-sided Pearson zero-correlation test"),
):
    """Report matched ranking or continuous-outcome metrics without fitting scores."""
    from . import __version__
    from .exports import _sha256

    if output.exists():
        raise ValueError("Output exists; choose a new directory")
    parameters = dict(groups=groups.split(",") if groups else [], keys=keys.split(","), endpoint=endpoint,
                      kind=kind, mask=mask, reference=reference,
                      baseline_scores=json.loads(baseline_scores.read_text()) if baseline_scores else None,
                      curves=curves, pearson_test=pearson_test)
    tables = evaluate_scores(read_table(results), json.loads(scores.read_text()), **parameters)
    output.mkdir(parents=True)
    for name, table in tables.items():
        table.write_csv(output / f"{name}.tsv", separator="\t")
    sources = [results, scores, *([baseline_scores] if baseline_scores else [])]
    (output / "SOURCES.json").write_text(json.dumps({"version": __version__, "parameters": parameters,
        "inputs": {str(p): _sha256(p) for p in sources}, "precision_at_recall": "maximum precision at recall >= target",
        "test_null": "Pearson correlation equals zero" if pearson_test else None}, indent=2)+"\n")
