"""Synthetic boundary checks for the October Potts audit; no scientific refits or inputs."""
from dataclasses import replace

import numpy as np
import polars as pl
import pytest
import scipy.sparse as sp

from tcren.potts import PottsModel, ais_log_z, gibbs, mu_star, site_codes
from tcren.potts.em import _eta_vector, consistency, fit_potts_em, latent_mask
from tcren.potts.ensemble import threaded_ensemble
from tcren.potts.fit import irls
from tcren.potts.kernel import edges


def _sites():
    return pl.DataFrame({"pdb.id": ["toy"] * 2, "chain.rec": ["TRB"] * 2,
                         "chain.id.rec": ["T"] * 2, "chain.id.par": ["P"] * 2,
                         "region.rec": ["CDR3"] * 2, "pos.rec": [0, 0],
                         "pos.par": [0, 1], "par.len": [3, 3], "aa.rec": ["A"] * 2,
                         "aa.par": ["A"] * 2, "cls": ["peptide"] * 2,
                         "role.par": ["tcr_facing"] * 2, "d_ca": [5., 5.],
                         "d_heavy": [4.9, 5.1], "sigma": [1., 0.], "probe": [15., 15.]})


def test_physical_mhc_chains_do_not_share_within_region_edges():
    a = _sites().with_columns(pl.lit("mhc").alias("cls"),
                             pl.lit("GROOVE_FLOOR").alias("role.par"))
    b = a.with_columns(pl.lit("Q").alias("chain.id.par"))
    _, _, q = site_codes(pl.concat([a, b]))
    assert q["pchain"].n_unique() == 2
    assert len(edges(q)[0]) == 2
    with pytest.raises(ValueError, match="regenerate legacy"):
        site_codes(a.drop("chain.id.par"))
    # The peptide-only legacy schema remains readable.
    site_codes(_sites().drop("chain.id.par", "chain.id.rec"))


def test_threading_carries_last_kept_chain_matrix(monkeypatch):
    import tcren.potts.ensemble as ensemble
    seen = []

    def fake_gibbs(et, A, cols, rng, *, chains, burn, draws, thin, observer, init):
        seen.append((np.asarray(init).copy(), burn))
        end = np.full((chains, len(et)), len(seen) % 2, dtype=float)
        observer(end)
        return end.mean(0), end.sum(1)

    monkeypatch.setattr(ensemble, "gibbs", fake_gibbs)
    for coupled in (True, False):
        seen.clear()
        threaded_ensemble(_sites(), PottsModel.bundled(),
                          pl.DataFrame({"pdb.id": ["toy", "toy"], "peptide": ["AAA", "GGG"]}),
                          coupled=coupled, chains=2, draws=1, burn=7, candidate_burn=2)
        assert seen[0][1] == 7 and seen[1][1] == 2
        assert seen[1][0].shape == (2, 2)
        assert np.array_equal(seen[1][0], np.ones((2, 2)))


def test_threading_rejects_unrepresented_tail_lengths_and_wrong_arm():
    model = PottsModel.bundled()
    for peptide in ("AA", "AAAA", "AAX"):
        with pytest.raises(ValueError, match="par.len|alphabet"):
            threaded_ensemble(_sites(), model,
                              pl.DataFrame({"pdb.id": ["toy"], "peptide": [peptide]}))
    with pytest.raises(ValueError, match="peptide partner"):
        threaded_ensemble(_sites().with_columns(pl.lit("mhc").alias("cls")), model,
                          pl.DataFrame({"pdb.id": ["toy"], "peptide": ["AAA"]}))


def test_orientation_is_counted_and_compared():
    a = replace(PottsModel.bundled(), g_orient=[0.] * 8)
    b = replace(a, g_orient=[1.] * 8)
    assert a.n_parameters() == PottsModel.bundled().n_parameters() + 8
    assert np.max(np.abs(_eta_vector(a) - _eta_vector(b))) == 1
    row = consistency(a, b).filter(pl.col("block") == "g_orient")
    assert row["max_abs"][0] == 1


def test_latent_mask_uses_probe_provenance_not_contact_counts():
    f = _sites().head(1).with_columns(pl.lit(5.).alias("d_heavy"))
    assert not latent_mask(f, margin=0).any()
    assert latent_mask(f).all()
    assert latent_mask(f.drop("probe")).all()
    with pytest.raises(ValueError, match="requires probe"):
        latent_mask(f.with_columns(pl.lit(5.).alias("probe")))
    with pytest.raises(ValueError, match="steps"):
        fit_potts_em(f, init=PottsModel.bundled(), steps=0)


def test_sampling_domains_matrix_starts_and_bracket_expansion():
    et, A, cols = np.zeros(2), np.zeros((2, 2)), [np.arange(2)]
    rng = np.random.default_rng(0)
    for kwargs in ({"draws": 0}, {"thin": 0}, {"chains": -1}, {"burn": -1},
                   {"init": np.full((2, 2), .5)}):
        with pytest.raises(ValueError):
            gibbs(et, A, cols, rng, chains=2, **kwargs) if "chains" not in kwargs else gibbs(et, A, cols, rng, **kwargs)
    gibbs(et, A, cols, rng, chains=2, burn=0, draws=1, init=np.ones((2, 2)))
    with pytest.raises(ValueError, match="steps"):
        ais_log_z(et, A, cols, rng, steps=0)
    assert mu_star(np.array([0., 1.]), .9, span=.01) == pytest.approx(np.log(9), abs=1e-7)


def test_irls_rejects_invalid_weights_and_unbounded_intercepts():
    X = sp.csr_matrix(np.ones((2, 1)))
    for y, w in ((np.ones(2), np.ones(2)), (np.array([0., 1.]), np.array([-1., 1.])),
                 (np.array([0., 1.]), np.zeros(2))):
        with pytest.raises(ValueError):
            irls(X, y, w, 1.)


def test_available_pairs_retains_chain_specific_groove_lengths():
    from tcren.potts import available_pairs
    from tcren.structure.model import Atom, Chain, RegionMarkup, Residue, Structure

    def chain(key, kind, count, region, y):
        residues = [Residue(i, i + 1, "", "A", "ALA",
                            (Atom("CA", "C", np.array([float(i), y, 0.])),))
                    for i in range(count)]
        return Chain(key, residues, kind, regions=[RegionMarkup(region, 0, count, "A" * count,
                                                                residues)])

    s = Structure("toy", [chain("T", "TRB", 3, "CDR3", 0.),
                          chain("A", "MHCa", 2, "GROOVE_FLOOR", 4.),
                          chain("B", "MHCb", 4, "GROOVE_FLOOR", -4.)])
    pairs = available_pairs(s, partner="mhc")
    assert pairs.filter(pl.col("chain.id.par") == "A")["par.len"].unique().to_list() == [2]
    assert pairs.filter(pl.col("chain.id.par") == "B")["par.len"].unique().to_list() == [4]
    assert pairs["chain.id.rec"].unique().to_list() == ["T"]
    assert pairs["probe"].unique().to_list() == [15.]


def test_sampling_extreme_finite_fields_are_stable():
    et = np.array([-1000., 1000.])
    occ, _ = gibbs(et, np.zeros((2, 2)), [np.arange(2)], np.random.default_rng(0),
                   chains=2, burn=0, draws=1, thin=1)
    assert np.array_equal(occ, [0., 1.])


def test_irls_backtracks_an_oversized_newton_step(monkeypatch):
    solve = np.linalg.solve
    calls = []

    def overshoot_first(A, b):
        calls.append(1)
        return (4 if len(calls) == 1 else 1) * solve(A, b)

    monkeypatch.setattr(np.linalg, "solve", overshoot_first)
    coef, _ = irls(sp.csr_matrix(np.ones((4, 1))), np.array([1., 1., 1., 0.]),
                   np.ones(4), 1.)
    assert coef[0] == pytest.approx(np.log(3), abs=1e-7)


def test_irls_returns_hessian_at_returned_coefficients():
    from scipy.special import expit
    x = np.column_stack([np.ones(10), np.arange(10)])
    y = np.array([0., 0., 0., 1., 0., 1., 1., 1., 1., 1.])
    b, H = irls(sp.csr_matrix(x), y, np.ones(10), 1., tol=.5)
    mu = expit(x @ b)
    want = x.T @ ((mu * (1 - mu))[:, None] * x) + np.diag([0., 1.])
    np.testing.assert_allclose(H, want)


def test_probe_provenance_must_be_finite_and_non_null():
    for value in (None, float("nan"), float("inf")):
        with pytest.raises(ValueError, match="finite and non-null"):
            latent_mask(_sites().with_columns(pl.lit(value, dtype=pl.Float64).alias("probe")))
    with pytest.raises(ValueError, match="draws"):
        gibbs(np.zeros(2), np.zeros((2, 2)), [np.arange(2)], np.random.default_rng(0), draws=True)


def test_partial_gibbs_colours_require_explicit_clamped_start():
    A = np.array([[0., 3.], [3., 0.]])
    et, colours, rng = np.zeros(2), [np.array([0])], np.random.default_rng(0)
    with pytest.raises(ValueError, match="explicitly clamps"):
        gibbs(et, A, colours, rng)
    p, _ = gibbs(et, A, colours, rng, init=np.array([0., 1.]),
                 chains=2, burn=0, draws=1, thin=1)
    assert p[1] == 1.


def test_orientation_uses_compact_ca_indices_and_reports_missing_atoms():
    from tcren.geometry import orientation_cosine
    from tcren.potts.sites import _orientation
    from tcren.structure.model import Atom, Chain, Residue, Structure

    def chain(key, count):
        return Chain(key, [Residue(i, i + 1, "", "A", "ALA", () if i == 1 else
                                  (Atom("CA", "C", np.array([float(i), float(i % 2), 1.])),))
                           for i in range(count)])

    rec, par = chain("T", 5), chain("P", 4)
    frame = pl.DataFrame({"key.tcr.chain": ["T", "T"], "key.par.chain": ["P", "P"],
                          "key.tcr.res": [3, 1], "key.par.res": [3, 2]})
    got = _orientation(Structure("toy", [rec, par]), frame).sort("key.tcr.res")
    expected = orientation_cosine(np.array([r.ca for r in rec.residues if r.ca is not None]),
                                  np.array([r.ca for r in par.residues if r.ca is not None]))[2, 2]
    assert got.filter(pl.col("key.tcr.res") == 3)["cos_cb"][0] == pytest.approx(expected)
    assert np.isnan(got.filter(pl.col("key.tcr.res") == 1)["cos_cb"][0])


def test_orientation_bins_reject_missing_geometry_explicitly():
    model = replace(PottsModel.bundled(), g_orient=[0.] * 8)
    with pytest.raises(ValueError, match="missing Cα geometry"):
        site_codes(_sites().with_columns(pl.lit(float("nan")).alias("cos_cb")), model)
