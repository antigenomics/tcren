"""The score set: the two algebraic identities it rests on, and the reproducibility contract.

Three things are asserted here rather than argued anywhere.

**There is no Jacobian in the back-transform.** A recurring proposal is to correct a posterior
computed in PCA coordinates by ``|det dz/dx|``. For a posterior that factor is constant in ``x``
and appears identically in numerator and denominator, so ``P(c|z) == P(c|x)`` exactly for any
invertible affine map. The Jacobian belongs to the density.

**What loses information under truncation is marginalization, not a determinant**, and for a
Gaussian a marginal is a sub-block of Sigma. Every ``channel_scores`` value is such a sub-block, so
this identity is load-bearing rather than decorative: if it did not hold, the channels would be a
second model wearing the first one's name.

**The frozen coefficients are reproducible from inputs that ship beside them.** This is the test
that distinguishes this read-out from the one the project withdrew.
"""
from __future__ import annotations

import pathlib

import numpy as np
import pytest

from tcren.score import CHANNELS, MODEL_FILE, holdout_model
from tcren.score.model import Joint


def _toy(p: int = 6, n: int = 400, seed: int = 0) -> Joint:
    rng = np.random.default_rng(seed)
    j = Joint(names=[f"x{i}" for i in range(p)])
    X, w = {}, {}
    for c in (0, 1):
        A = rng.normal(size=(p, p))
        X[c] = rng.normal(size=(n, p)) @ A + c * 0.7
        w[c] = np.full(n, 1.0 / n)
    return j.fit(X, w, shrink=False)


def test_a_full_rank_rotation_leaves_the_posterior_unchanged():
    """The Jacobian claim, made testable: a square PCA is a rotation and the log-odds is invariant."""
    j = _toy()
    rng = np.random.default_rng(1)
    X = rng.normal(size=(200, len(j.names)))
    W, jp = j.project(len(j.names))
    assert np.abs(j.log_odds(X) - jp.log_odds(X @ W)).max() < 1e-9


def test_a_sub_block_is_the_marginal_and_not_an_approximation_of_it():
    """`channel_scores` is this identity. Checked against numerical quadrature on a 2D case."""
    j = _toy(p=2, n=600, seed=2)
    rng = np.random.default_rng(3)
    X = rng.normal(size=(40, 2))

    def by_quadrature(x0: float) -> float:
        grid = np.linspace(-25, 25, 4001)
        pts = np.column_stack([np.full_like(grid, x0), grid])
        dens = {}
        for c in (0, 1):
            S, mu = j.cov[c], j.mu[c]
            Z = pts - mu
            q = np.einsum("ij,jk,ik->i", Z, np.linalg.inv(S), Z)
            _, logdet = np.linalg.slogdet(S)
            dens[c] = j.prior[c] * np.exp(-0.5 * (q + logdet))
        return float(np.log(np.trapezoid(dens[1], grid) / np.trapezoid(dens[0], grid)))

    exact = j.log_odds(X, subset=["x0"])
    quad = np.array([by_quadrature(v) for v in X[:, 0]])
    assert np.abs(exact - quad).max() < 1e-8


def test_truncating_a_component_does_change_the_posterior():
    """The other half of the same point: dropping a coordinate is not free, and must not look it."""
    j = _toy()
    rng = np.random.default_rng(4)
    X = rng.normal(size=(200, len(j.names)))
    assert np.abs(j.log_odds(X) - j.log_odds(X, subset=j.names[:3])).max() > 1.0


def test_an_unshrunk_covariance_is_left_alone():
    """`shrink=False` is what the artefact test needs, so it has to really be off."""
    j = _toy()
    assert j.alpha[0] == 0.0 and j.alpha[1] == 0.0


# --------------------------------------------------------------------- the shipped model
#: The whole hold-out descriptor table is 8,292 rows by 147 descriptors, 19 MB, and deliberately
#: does NOT ship: the manifest names the structures and `tcren features` regenerates it. What ships
#: here instead is a **slice** of it, so that every test below runs for anyone who clones the repo
#: rather than skipping on a path only one machine has.
#:
#: 362 structures -- 232 binders and 130 non-binders over 10 epitopes -- carrying all 147 modelling
#: descriptors, every value present and finite, rounded to four significant digits. The row keys are
#: the real `pdb.id` hashes, so the shipped `holdout_manifest` supplies the label, the epitope and
#: the ipTM by an ordinary join and nothing about the fixture is invented. Rebuild it by taking the
#: first N ids per (epitope, class) in `pdb.id` order from a regenerated hold-out table; the counts
#: are chosen so the binder arm clears the 200-row floor `fit_holdout` needs before it estimates the
#: ipTM coordinate at all.
SLICE_DIR = pathlib.Path(__file__).resolve().parents[1] / "assets" / "score"
SLICE_FEATURES = SLICE_DIR / "holdout_slice_features.tsv.gz"
#: What the shipped fitter returned from that slice, frozen: `mu0`, `mu1`, `cov0`, `cov1`, the ipTM
#: covariance row and the transform parameters. Covariances are stored float32 -- 149 x 149 twice is
#: the whole file -- which sets the tolerance the refit is checked at.
SLICE_MODEL = SLICE_DIR / "holdout_slice_model.npz"


def _holdout_features(rows: int | None = 300):
    import gzip

    import polars as pl
    with gzip.open(SLICE_FEATURES, "rb") as fh:
        t = pl.read_csv(fh.read(), separator="\t", infer_schema_length=None)
    return t if rows is None else t.head(rows)


pytestmark_model = pytest.mark.skipif(not MODEL_FILE.exists(),
                                      reason="holdout_model.npz not present in this checkout")


@pytestmark_model
def test_the_shipped_model_loads_and_declares_what_it_was_fitted_on():
    m = holdout_model()
    assert m.n_pos > 1000 and m.n_neg > 200 and m.n_epitopes > 10
    assert len(m.receptor_names) < len(m.joint.names), \
        "the receptor read-out must marginalize out the no-receptor descriptors"
    assert len(m.catalogue_digest) == 64


@pytestmark_model
def test_every_channel_is_populated_and_they_are_not_the_same_number():
    """Five named channels, each a real sub-block. Identical values would mean the index is wrong."""

    m = holdout_model()
    t = _holdout_features()
    ch = m.channel_scores(t)
    assert set(ch) == set(CHANNELS)
    for name, v in ch.items():
        assert np.isfinite(v).sum() > 100, name
    vals = np.array([v[np.isfinite(v)][:50] for v in ch.values()])
    assert np.abs(np.corrcoef(vals) - np.eye(len(ch))).max() > 0.05


@pytestmark_model
def test_scores_do_not_depend_on_what_was_scored_alongside_them():
    """The single-structure claim. Scoring one row must give what scoring 300 gave for that row."""

    from tcren.score import score_table
    t = _holdout_features()
    many = score_table(t)
    ok = np.isfinite(many["binder_score"].to_numpy())
    i = int(np.flatnonzero(ok)[0])
    one = score_table(t[i])
    for c in ("pose_score", "binder_score", "channel_shape", "peptide_score"):
        assert one[c][0] == pytest.approx(many[c][i], rel=1e-12, abs=1e-12), c


def test_refitting_the_committed_slice_reproduces_the_frozen_arrays():
    """The reproducibility contract, asserted on the arrays and run by anyone who clones the repo.

    `P_native` was removed from this project because its coefficients were frozen against a training
    set that no longer existed. The claim that replaces it is that the shipped fitter, handed a
    named input, returns a frozen output -- so the test has to pin the model itself, and a version
    string or a row count is not the model.

    **What is reproduced.** `fit_holdout` is re-run on the committed 362-structure slice (232
    binders, 130 non-binders, 10 epitopes, 147 descriptors) joined to the shipped manifest, and
    every array it emits is compared against `holdout_slice_model.npz`: the two class means, 149
    coordinates each; the two class covariances, 149 x 149 each; the ipTM covariance row, 149; and
    the three per-descriptor transform parameter maps. This is NOT the 8,292-structure fit the wheel
    carries -- that one needs the 19 MB table, and
    `test_the_shipped_model_is_reproducible_from_the_full_hold_out_table` is where it is checked.

    **To what tolerance.** Everything is checked at rtol 1e-6. The binding term is not BLAS:
    `Transformer.fit` takes each Yeo-Johnson lambda from `scipy.stats.yeojohnson_normmax`, whose
    Brent search converges only to its own default tolerance of about 1.5e-8, and a lambda wobble
    that size moves a column mean by about 9e-9. Between macOS arm64 and Linux x86_64 the refit
    means duly differ, by at most 6.0e-9 absolute on 18 of the 149 coordinates and 1.4e-7 relative
    -- so the earlier rtol of 1e-9 on the means asserted a tolerance the fit never promised, passed
    on the machine it was written on and failed on CI. The covariances are stored float32, which
    costs at most 5.9e-8 relative, so 1e-6 is a hundredfold above their storage floor. Every
    tolerance here is orders of magnitude below what a change to the epitope weighting, the
    Yeo-Johnson fit or the Ledoit-Wolf shrinkage would move these arrays by.
    """
    import json
    import tempfile

    from tcren.score import holdout_manifest
    from tcren.score.fit import fit_holdout

    ref = np.load(SLICE_MODEL, allow_pickle=False)
    rmeta = json.loads(str(ref["meta"]))
    with tempfile.TemporaryDirectory() as d:
        out = pathlib.Path(d) / "refit.npz"
        meta = fit_holdout(_holdout_features(rows=None), holdout_manifest(), out=out)
        got = {k: np.asarray(v, float) for k, v in np.load(out, allow_pickle=False).items()
               if k != "meta"}

    # the column set first: a mismatch here means the catalogue moved and the fixture is stale,
    # which would otherwise surface as an unreadable diff over 22,201 covariance entries.
    for k in ("descriptors", "coordinates", "receptor_coordinates"):
        assert meta[k] == rmeta[k], f"{k} moved -- regenerate the slice fixture"

    for k in ("mu0", "mu1", "conf_cov"):
        np.testing.assert_allclose(got[k], np.asarray(ref[k], float),
                                   rtol=1e-6, atol=1e-7, err_msg=k)
    for k in ("cov0", "cov1"):
        np.testing.assert_allclose(got[k], np.asarray(ref[k], float),
                                   rtol=1e-6, atol=1e-9, err_msg=k)
    for k in ("lam", "loc", "scale"):
        assert set(meta[k]) == set(rmeta[k]), k
        np.testing.assert_allclose([meta[k][n] for n in rmeta[k]],
                                   [rmeta[k][n] for n in rmeta[k]],
                                   rtol=1e-6, atol=1e-9, err_msg=k)
    for k in ("prior", "n_pos", "n_neg", "n_epitopes"):
        assert meta[k] == rmeta[k], k
    for k in ("alpha", "conf_mu", "conf_var"):
        assert meta[k] == pytest.approx(rmeta[k], rel=1e-6), k


@pytestmark_model
@pytest.mark.slow
def test_the_shipped_model_is_reproducible_from_the_full_hold_out_table():
    """The same contract at full scale. Opt-in, because its 19 MB input does not ship.

    Regenerate the hold-out descriptor table with `tcren features` over the 8,292 structures
    `holdout_manifest` names, then point `TCREN_HOLDOUT_FEATURES` at it. On the table the shipped
    model was frozen from, the refit is bit-identical: `mu0`, `mu1`, `cov0`, `cov1` and the ipTM
    covariance row all differ by exactly 0, so this is asserted at zero tolerance rather than a
    nominal one.
    """
    import json
    import os
    import tempfile

    import polars as pl

    from tcren.score import holdout_manifest
    from tcren.score.fit import fit_holdout

    src = os.environ.get("TCREN_HOLDOUT_FEATURES", "")
    if not src or not pathlib.Path(src).exists():
        pytest.skip("set TCREN_HOLDOUT_FEATURES to a regenerated hold-out feature table")

    features = pl.read_csv(src, separator="\t", infer_schema_length=None)
    with tempfile.TemporaryDirectory() as d:
        out = pathlib.Path(d) / "refit.npz"
        meta = fit_holdout(features, holdout_manifest(), out=out)
        got = {k: np.asarray(v, float) for k, v in np.load(out, allow_pickle=False).items()
               if k != "meta"}

    shipped = np.load(MODEL_FILE, allow_pickle=False)
    ref = json.loads(str(shipped["meta"]))
    for k in ("mu0", "mu1", "cov0", "cov1", "conf_cov"):
        assert np.array_equal(got[k], np.asarray(shipped[k], float)), k
    for k in ("coordinates", "receptor_coordinates", "prior", "alpha", "n_pos", "n_neg",
              "n_epitopes", "conf_mu", "conf_var", "catalogue_digest"):
        assert meta[k] == ref[k], k


def test_the_artefact_test_finds_a_planted_generator_regularity():
    """A direction the modelled binders are pinned to but real complexes are free on.

    Planted rather than measured, because the point of the test is the criterion, not the corpus:
    if a 100-fold spread ratio on a pinned direction does not flag, the criterion is broken.
    """
    from tcren.reliability import artefact_directions

    rng = np.random.default_rng(0)
    binders = rng.normal(size=(500, 8))
    binders[:, 0] *= 0.01                     # the generator pins this one
    crystals = rng.normal(size=(60, 8))       # physics does not
    r = artefact_directions(binders, crystals)
    assert r["is_artefact"].sum() == 1
    assert r["ratio"][0] > 20
    # a direction nobody pins is not flagged, whatever its spread ratio
    assert not r["is_artefact"][r["sd_binder"] > 0.5].any()


# --------------------------------------------------------------------- the transform layer
def test_every_transform_class_survives_a_legal_input():
    """Nine classes, one per unit family. A transform that quietly emits inf poisons a covariance.

    This was `transform.py`'s `__main__` self-check before the module moved into the package and
    stopped being runnable on its own.
    """
    from tcren.recognition import DESCRIPTORS
    from tcren.score.transform import Transformer, kind

    t = Transformer(names=list(DESCRIPTORS))
    rng = np.random.default_rng(0)
    X = np.empty((200, len(DESCRIPTORS)))
    for j, n in enumerate(DESCRIPTORS):
        X[:, j] = {"logit": rng.uniform(0, 1, 200), "atanh": rng.uniform(-1, 1, 200),
                   "fisher": rng.uniform(-1, 1, 200), "log": rng.uniform(1, 12, 200),
                   "anscombe": rng.integers(0, 40, 200).astype(float),
                   "cossin": rng.uniform(-180, 180, 200)}.get(kind(n), rng.normal(0, 3, 200))
    X[0, 0] = np.nan
    Z = t.fit(X).transform(X)
    assert Z.shape[1] == len(t.out_names())
    assert np.isnan(Z[0, 0]), "a NaN input must stay NaN, never become 0"
    assert np.isfinite(Z[1:, :]).all(), "a legal input produced a non-finite coordinate"
    # a hard 0 and a hard 1 are legal for a fraction and must not go infinite
    frac = next(d for d in DESCRIPTORS if kind(d) == "logit")
    j = list(DESCRIPTORS).index(frac)
    Y = X.copy()
    Y[1, j], Y[2, j] = 0.0, 1.0
    assert np.isfinite(t.transform(Y)[[1, 2], t.out_names().index(frac)]).all()


def test_transforming_a_subset_of_columns_matches_transforming_all_of_them():
    """A feature table missing a whole family is scored by marginalization, and this is why.

    The transform parameters are per descriptor, so restricting the columns is exact rather than
    an approximation. If it were not, every degraded read-out would be silently on a wrong scale.
    """
    from tcren.recognition import DESCRIPTORS
    from tcren.score.transform import Transformer

    names = list(DESCRIPTORS)[:40]
    rng = np.random.default_rng(7)
    X = np.abs(rng.normal(3, 1, size=(150, len(names))))
    t = Transformer(names=names).fit(X)
    full = t.transform(X)
    sub = names[5:25]
    idx = [names.index(n) for n in sub]
    part = t.transform(X[:, idx], names=sub)
    cols = t.out_names()
    for k, n in enumerate(t.out_names(sub)):
        assert np.allclose(part[:, k], full[:, cols.index(n)], equal_nan=True), n


@pytestmark_model
def test_a_table_missing_a_whole_family_is_still_scored():
    """`tcren features -i placement,interface,topology` omits potts and kinetics: 18 coordinates.

    Dropping them must marginalize, not raise -- and must not silently return the same number as
    the full table, because information really was lost.
    """
    from tcren.score import score_table

    t = _holdout_features()
    drop = [c for c in t.columns
            if c in ("log_z", "log_lik", "psi", "n_contacts", "exp_lost")]
    assert drop, "the fixture no longer carries the columns this test drops"
    full = score_table(t)
    thin = score_table(t.drop(drop))
    ok = np.isfinite(full["binder_score"].to_numpy()) & np.isfinite(thin["binder_score"].to_numpy())
    assert ok.sum() > 100
    a, b = full["binder_score"].to_numpy()[ok], thin["binder_score"].to_numpy()[ok]
    assert not np.allclose(a, b), "dropping five coordinates must change the posterior"
    assert np.corrcoef(a, b)[0, 1] > 0.5, "and must not scramble it either"
