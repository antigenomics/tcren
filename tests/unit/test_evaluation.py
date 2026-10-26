"""Matched populations, tied thresholds and historical correlation semantics."""
import json

import numpy as np
import polars as pl
import pytest
from scipy.stats import pearsonr, spearmanr
from typer.testing import CliRunner

from tcren.cli import app
from tcren.evaluation import evaluate_scores


def test_ties_missing_masks_and_previous_columns():
    rows = pl.DataFrame({"id": list("abcde"), "cohort": ["a"]*5,
                         "y": [1., 0., 1., 0., None], "a": [1., 1., 0., 0., .5],
                         "b": [1., 0., 1., None, .5], "old:a": [0., 0., 1., 1., .5]})
    result = evaluate_scores(rows, {"a": 1, "b": 1}, groups=["cohort"], mask="pair", curves=True)
    metrics = result["metrics"]
    a = metrics.filter(pl.col("score") == "a")
    assert set(a["n_used"]) == {4}
    assert a.filter(pl.col("metric").is_in(["roc_auc", "average_precision"]))["value"].to_list() == [.5, .5]
    assert result["curves"].filter(pl.col("score") == "a").select("recall", "precision", "fpr").rows() == [(0., 1., 0.), (.5, .5, .5), (1., .5, 1.)]
    common = evaluate_scores(rows, {"a": 1, "b": 1}, baseline_scores={"a": "old:a"})["metrics"]
    assert set(common["n_used"]) == {3}
    assert common.filter((pl.col("score") == "a") & (pl.col("metric") == "roc_auc"))["delta"].item() == -.5
    with pytest.raises(ValueError, match="binary"):
        evaluate_scores(rows.with_columns(pl.lit(2).alias("y")), {"a": 1})
    with pytest.raises(ValueError, match="duplicate"):
        evaluate_scores(pl.concat([rows, rows]), {"a": 1})


def test_continuous_keeps_unlabelled_rows_and_explicit_test():
    rows = pl.DataFrame({"id": list("abcd"), "response": [2., 4., 3., 8.], "a": [1., 2., 3., 4.]})
    result = evaluate_scores(rows, {"a": -1}, endpoint="response", kind="continuous", pearson_test=True)["metrics"]
    values = dict(result.select("metric", "value").iter_rows())
    assert set(result["n_used"]) == {4}
    assert values["pearson_R"] == pytest.approx(pearsonr(rows["response"], -rows["a"]).statistic)
    assert values["spearman_R"] == pytest.approx(spearmanr(rows["response"], -rows["a"]).statistic)
    assert values["pearson_P"] == pytest.approx(pearsonr(rows["response"], -rows["a"]).pvalue)
    constant = evaluate_scores(rows.with_columns(pl.lit(1.).alias("a")), {"a": 1}, endpoint="response", kind="continuous")["metrics"]
    assert constant["value"].null_count() == 2


def test_cli_exports_replayable_metrics_and_refuses_overwrite(tmp_path):
    rows = tmp_path / "rows.tsv"
    pl.DataFrame({"id": ["a", "b", "c", "d"], "y": [1, 0, 1, 0], "s": [1., 1., 0., 0.]}).write_csv(rows, separator="\t")
    directions = tmp_path / "scores.json"
    directions.write_text(json.dumps({"s": 1}))
    out = tmp_path / "evaluated"
    args = ["evaluate", "--results", str(rows), "--scores", str(directions), "--curves", "-o", str(out)]
    result = CliRunner().invoke(app, args)
    assert result.exit_code == 0, result.output
    metrics = pl.read_csv(out / "metrics.tsv", separator="\t")
    np.testing.assert_allclose(metrics["value"], [.5, .5, .5, .5])
    assert json.loads((out / "SOURCES.json").read_text())["test_null"] is None
    assert CliRunner().invoke(app, args).exit_code != 0


def test_reserved_group_and_identity_names_fail_before_overwrite():
    rows = pl.DataFrame({"id": ["a", "b"], "metric": ["m", "m"], "y": [0, 1], "x": [0., 1.]})
    for args in ({"groups": ["metric"]}, {"keys": ["metric"]}):
        with pytest.raises(ValueError, match="reserved outputs"):
            evaluate_scores(rows, {"x": 1}, **args)


def test_all_missing_scores_keep_membership_and_zero_count():
    rows = pl.DataFrame({"id": ["a", "b"], "y": [1, 0], "x": [None, None]})
    result = evaluate_scores(rows, {"x": 1})
    assert result["membership"].height == 2
    assert not result["membership"]["used"].any()
    assert set(result["metrics"]["n_used"]) == {0}
    with pytest.raises(ValueError, match="numeric"):
        evaluate_scores(rows.with_columns(pl.lit("broken").alias("x")), {"x": 1})
