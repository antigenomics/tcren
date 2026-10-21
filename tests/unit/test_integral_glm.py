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


@pytest.mark.parametrize("pairwise", [False, True])
def test_glm_cli_roundtrip_and_frozen_output(tmp_path, pairwise):
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
    if pairwise:
        args.append('--pairwise')
    result = runner.invoke(app, args)
    assert result.exit_code == 0, result.exception
    before = model.read_bytes()
    assert runner.invoke(app, args).exit_code != 0
    assert model.read_bytes() == before
    result = runner.invoke(app, ['integral-predict', '--features', str(paths['features']),
                                '--model', str(model), '-o', str(paths['scores'])])
    assert result.exit_code == 0, result.exception
    assert pl.read_csv(paths['scores'], separator='\t').height == len(data)


def test_pairwise_glm_keeps_original_support_and_has_no_intercept():
    data, labels = sample()
    model = fit(data, labels, pairwise=True)
    assert model['format'] == 'integral-pairwise-glm-1'
    assert model['n_training'] == 80 and model['n_pairs'] == 800
    assert model['intercept'] == 0
    scores = integral_glm_scores(data, model)
    assert scores.filter(labels['y']==1)['integral_glm_logit'].min() > 0
    assert scores.filter(labels['y']==0)['integral_glm_logit'].max() < 0
    assert integral_glm_scores(data.head(1), model).equals(scores.head(1))


def test_pair_weights_preserve_model_count_and_weight_epitope_support():
    from tcren.integral_glm import _pairwise_design
    design=np.arange(14).reshape(7,2)
    paired,y,w,n=_pairwise_design(design,np.array([0,1,0,0,1,0,1]),np.array(['a']*3+['b']*4))
    assert n==6 and paired.shape==(12,2)
    assert w.sum()==pytest.approx(7)
    assert w[:4].sum()==pytest.approx(7*np.sqrt(3)/(np.sqrt(3)+2))
    assert w[4:].sum()==pytest.approx(14/(np.sqrt(3)+2))
    np.testing.assert_array_equal(paired[:2],-paired[2:4])
    assert y[:4].tolist()==[1,1,0,0]
    with pytest.raises(ValueError,match='one million'):
        _pairwise_design(np.zeros((2001,8)),np.r_[np.zeros(1000),np.ones(1001)],np.full(2001,'a'))
