"""Frozen GLM isolation, additive accounting and explicit training exclusions."""
import copy

import numpy as np
import polars as pl
import pytest

from tcren.integral_glm import fit_integral_glm, integral_glm_scores


def sample():
    rng = np.random.default_rng(0)
    n = 80
    y = np.tile([0, 1], n // 2)
    data = pl.DataFrame({"complex.id": [str(i) for i in range(n)],
        "mhc_class_bin": np.zeros(n), "Phi_tcr_pep": rng.normal(size=n),
        "n_contacts_tp": np.full(n, 10), "sc_shape": rng.uniform(-1, 1, n),
        "g_cyclo_frac": rng.uniform(0, 1, n), "iptm": .2 + .6*y})
    labels = data.select("complex.id").with_columns(pl.Series("y", y),
        pl.Series("epitope", ["train_a"]*40+["train_b"]*40),
        pl.Series("receptor_id", [f"r{i}" for i in range(n)]))
    return data, labels


def fit(data, labels, **kwargs):
    return fit_integral_glm(data, labels, excluded_epitopes={"held"},
        excluded_receptors={"held_r"}, provenance={"input_sha256": "a"*64}, **kwargs)


def test_glm_frozen_batch_invariance_missingness_and_additive_logit():
    data, labels = sample()
    model = fit(data, labels.reverse())
    out = integral_glm_scores(data, model)
    assert out.filter(labels["y"] == 1)["integral_glm_score"].min() > .9
    assert out.filter(labels["y"] == 0)["integral_glm_score"].max() < .1
    assert integral_glm_scores(data.head(1), model).equals(out.head(1))
    total = out.select(pl.sum_horizontal(pl.col("^glm_.*_contribution$")).alias("sum"))["sum"]
    np.testing.assert_allclose(total + model["intercept"], out["integral_glm_logit"])
    partial = data.head(1).with_columns(pl.lit(None, dtype=pl.Float64).alias("iptm"))
    assert integral_glm_scores(partial, model)["n_components"][0] == 3
    absent = partial.with_columns(*[pl.lit(None, dtype=pl.Float64).alias(c)
                                   for c in ["Phi_tcr_pep", "sc_shape", "g_cyclo_frac"]])
    assert integral_glm_scores(absent, model)["integral_glm_score"][0] is None
    bad = copy.deepcopy(model)
    bad["scales"][0] = 0
    with pytest.raises(ValueError, match="parameters"):
        integral_glm_scores(data, bad)


@pytest.mark.parametrize("column,value", [("epitope", "held"), ("receptor_id", "held_r"), ("y", 3)])
def test_glm_training_boundary_is_enforced(column, value):
    data, labels = sample()
    with pytest.raises(ValueError):
        fit(data, labels.with_columns(pl.lit(value).alias(column)))


def test_glm_cli_roundtrip_and_frozen_output(tmp_path):
    from typer.testing import CliRunner
    from tcren.cli import app
    data, labels = sample()
    paths = {k: tmp_path/(k+'.tsv') for k in ['features', 'labels', 'epitopes', 'receptors', 'scores']}
    data.write_csv(paths['features'], separator='\t')
    labels.write_csv(paths['labels'], separator='\t')
    pl.DataFrame({'epitope':['held']}).write_csv(paths['epitopes'], separator='\t')
    pl.DataFrame({'receptor_id':['held_r']}).write_csv(paths['receptors'], separator='\t')
    model = tmp_path/'model.json'
    runner = CliRunner()
    args = ['integral-fit', '--features', str(paths['features']), '--labels', str(paths['labels']),
            '--excluded-epitopes', str(paths['epitopes']), '--excluded-receptors', str(paths['receptors']),
            '-o', str(model)]
    result = runner.invoke(app, args)
    assert result.exit_code == 0, result.exception
    before = model.read_bytes()
    assert runner.invoke(app, args).exit_code != 0
    assert model.read_bytes() == before
    result = runner.invoke(app, ['integral-predict', '--features', str(paths['features']),
                                '--model', str(model), '-o', str(paths['scores'])])
    assert result.exit_code == 0, result.exception
    assert pl.read_csv(paths['scores'], separator='\t').height == len(data)
