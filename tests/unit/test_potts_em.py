"""Monte-Carlo EM over the borderline contact calls — :mod:`tcren.potts.em`.

Two claims are tested here. First, the latent band must reach BOTH ways
around the cutoff, or the fit can only ever remove contacts; that needs `available_pairs(probe=...)`
and the code refuses to run without it. Second, EM must reduce EXACTLY to the pseudolikelihood fit
when nothing is latent, because it is an extension of that fit and not a replacement for it.
"""
from __future__ import annotations

import numpy as np
import polars as pl
import pytest

pytest.importorskip("arda")

from pathlib import Path  # noqa: E402

from tcren.annotation.batch import iter_annotated_set  # noqa: E402
from tcren.mhc import annotate_mhc  # noqa: E402
from tcren.potts import available_pairs, fit_potts  # noqa: E402
from tcren.potts.em import (LATENT_MARGIN, consistency, fit_potts_em, impute,  # noqa: E402
                            latent_mask)

PDB_DIR = Path(__file__).resolve().parents[1] / "assets" / "pdb"
pytestmark = pytest.mark.skipif(not PDB_DIR.exists(), reason="pdb assets missing")


@pytest.fixture(scope="module")
def tables():
    """Both pair tables from ONE annotation pass over ONE set of structure objects.

    Annotating twice and comparing the two results would make the test depend on the annotation
    returning the same structures both times, which is not what is under test here.
    """
    wide, narrow = [], []
    for s in iter_annotated_set(str(PDB_DIR)):
        try:
            annotate_mhc(s)
            w = available_pairs(s, "peptide", probe=5.0 + LATENT_MARGIN)
            n = available_pairs(s, "peptide")
        except Exception:
            continue
        if w.height and n.height:
            wide.append(w)
            narrow.append(n)
    return pl.concat(wide), pl.concat(narrow)


@pytest.fixture(scope="module")
def wide(tables):
    """Pair table with the heavy-atom layer measured past the cutoff — what EM needs."""
    return tables[0]


@pytest.fixture(scope="module")
def narrow(tables):
    """The default table, whose `d_heavy` stops at the cutoff."""
    return tables[1]


def test_probe_widens_d_heavy_without_moving_sigma(wide, narrow):
    """The probe is a measurement radius, not the contact definition."""
    assert wide.height == narrow.height
    assert wide["sigma"].sum() == narrow["sigma"].sum()
    assert wide["d_heavy"].null_count() < narrow["d_heavy"].null_count()


def test_the_default_table_gives_a_one_sided_band_and_is_refused(narrow):
    """Every site in the band is a contact, so EM on it could only ever remove contacts."""
    d = narrow["d_heavy"].fill_null(np.inf).to_numpy()
    band = np.abs(d - 5.0) <= LATENT_MARGIN
    assert band.any()
    assert narrow["sigma"].to_numpy()[band].min() == 1.0
    with pytest.raises(ValueError, match="only reaches inwards"):
        latent_mask(narrow)


def test_the_wide_table_gives_a_two_sided_band(wide):
    mask = latent_mask(wide)
    sigma = wide["sigma"].to_numpy()[mask]
    assert (sigma == 1).any() and (sigma == 0).any()


def test_imputation_touches_only_the_latent_sites(wide):
    """Every confident call is clamped to what the structure did."""
    model = fit_potts(wide)
    stacked, mask, _ = impute(wide, model, draws=3, seed=0)
    from tcren.potts.score import _prepare
    _, _, _, sigma, _, _, _, _ = _prepare(wide, model)
    n = len(sigma)
    for d in range(3):
        got = stacked["sigma"].to_numpy()[d * n:(d + 1) * n]
        assert np.array_equal(got[~mask], sigma[~mask])


def test_draws_are_kept_apart_so_couplings_do_not_span_them(wide):
    """Every coupling edge joins on `pdb.id`; unsuffixed copies would link sites across draws."""
    model = fit_potts(wide)
    stacked, _, _ = impute(wide, model, draws=3, seed=0)
    ids = stacked["pdb.id"].unique().to_list()
    assert all("#" in i for i in ids)
    assert len(ids) == 3 * wide["pdb.id"].n_unique()


def test_no_latent_set_reduces_to_the_pseudolikelihood_fit(wide):
    """`margin=0` makes every site observed, so EM must return the fit it started from."""
    ref = fit_potts(wide)
    got, trace = fit_potts_em(wide, init=ref, margin=0.0, draws=2, steps=2, seed=0)
    assert trace["n_latent"].max() == 0
    assert got.alpha == pytest.approx(ref.alpha, rel=1e-6)
    assert np.allclose(got.h_rec, ref.h_rec, atol=1e-6)
    assert np.allclose(got.h_par, ref.h_par, atol=1e-6)
    assert np.allclose(got.kernel, ref.kernel, atol=1e-6)
    assert np.allclose(got.coupling_array(), ref.coupling_array(), atol=1e-6)


def test_the_same_seed_reproduces_the_fit(wide):
    """Every sampler is seeded, so a run is exactly reproducible however noisy EM's walk is."""
    ref = fit_potts(wide)
    a, _ = fit_potts_em(wide, init=ref, draws=3, steps=2, seed=11)
    b, _ = fit_potts_em(wide, init=ref, draws=3, steps=2, seed=11)
    assert np.array_equal(np.asarray(a.h_rec), np.asarray(b.h_rec))
    assert np.array_equal(a.coupling_array(), b.coupling_array())
    assert a.alpha == b.alpha


def test_a_different_seed_moves_the_fit(wide):
    """If the seed did nothing, the E-step would not be sampling."""
    ref = fit_potts(wide)
    a, _ = fit_potts_em(wide, init=ref, draws=3, steps=2, seed=11)
    b, _ = fit_potts_em(wide, init=ref, draws=3, steps=2, seed=12)
    assert not np.array_equal(np.asarray(a.h_rec), np.asarray(b.h_rec))


def test_counts_report_the_fit_not_the_stacked_copies(wide):
    """The M-step sees `draws` copies of everything; the model must not claim them as data."""
    ref = fit_potts(wide)
    got, _ = fit_potts_em(wide, init=ref, draws=4, steps=1, seed=0)
    assert got.n_structures == ref.n_structures
    assert got.n_sites == ref.n_sites
    assert "EM:" in got.notes


def test_trace_separates_drift_from_jitter(wide):
    """Monte-Carlo EM jitters; `max_shift_vs_init` is what says whether it is also moving."""
    ref = fit_potts(wide)
    _, trace = fit_potts_em(wide, init=ref, draws=3, steps=4, seed=0)
    assert set(trace.columns) >= {"step", "n_latent", "mean_imputed_sigma", "max_shift",
                                  "max_shift_vs_init", "converged"}
    assert trace.height == 4
    assert (trace["n_latent"] > 0).all()


def test_consistency_of_a_model_with_itself_is_zero(wide):
    ref = fit_potts(wide)
    c = consistency(ref, ref)
    assert c["max_abs"].max() == pytest.approx(0.0, abs=1e-12)
    assert c["n_beyond_1"].sum() == 0


def test_consistency_refuses_a_yardstick_it_does_not_have(wide):
    """A block with no spread has no scale; reporting a ratio against ~0 would be noise."""
    ref = fit_potts(wide)
    got, _ = fit_potts_em(wide, init=ref, draws=3, steps=1, seed=0)
    c = consistency(ref, got)
    # one partner class only, so g_class is constant and cannot be scaled
    assert c.filter(pl.col("block") == "g_class")["yardstick"][0] == "degenerate"
    # the coupling kernel is the one block with a real standard error
    assert c.filter(pl.col("block") == "kernel")["yardstick"][0] == "cluster-robust s.e."
