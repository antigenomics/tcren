"""The ensemble read-outs of :mod:`tcren.potts.ensemble`.

The claim under test is one identity and its two conditions. `Phi(sigma) = sum_a sigma_a phi_a` is
linear, so `<Phi> = sum_a <sigma_a> phi_a` for any joint — but only when `<sigma_a>` is the marginal
of the *coupled* model and is read at the sequence being scored. Every test here pins one half of
that, or pins a quantity the marginals cannot determine at all.
"""
from __future__ import annotations

import numpy as np
import polars as pl
import pytest

from tcren.potential import tcren2
from tcren.potts import PottsModel, available_pairs, contact_probabilities
from tcren.potts.ensemble import (energy_ensemble, map_ensemble, site_energies, threaded_ensemble)
from tcren.potts.sample import exact_log_z
from tcren.potts.score import _one, _prepare, _rng_for

pytest.importorskip("arda")

from pathlib import Path  # noqa: E402

from tcren.annotation import classify_chains  # noqa: E402
from tcren.mhc import annotate_mhc  # noqa: E402
from tcren.structure import parse_structure  # noqa: E402

ASSET = Path(__file__).resolve().parents[1] / "assets" / "pdb" / "1ao7.pdb"
pytestmark = pytest.mark.skipif(not ASSET.exists(), reason="1ao7 asset missing")


@pytest.fixture(scope="module")
def pairs():
    s = parse_structure(str(ASSET), pdb_id="1ao7")
    classify_chains(s)
    annotate_mhc(s)
    return available_pairs(s, "peptide")


@pytest.fixture(scope="module")
def model():
    return PottsModel.bundled()


def test_mean_equals_coupled_marginals_within_monte_carlo_error(pairs, model):
    """`<Phi>` from whole configurations against `sum_a p_a phi_a` from an independent chain.

    `contact_probabilities` starts its chains from the factorised model and `energy_ensemble` from
    the observed map, so these are genuinely different chains sampling the same distribution. Any
    disagreement beyond Monte-Carlo error is a bug in one of them, or a chain that has not mixed.
    """
    ens = energy_ensemble(pairs, model, workers=1, seed=0).row(0, named=True)
    _, _, _, _, _, _, _, q = _prepare(pairs, model)
    phi, _ = site_energies(q, tcren2())
    cp = contact_probabilities(pairs, model, workers=1, seed=0)

    independent = float(cp["p_model"].to_numpy() @ phi)
    se = np.sqrt(ens["phi_var"] / ens["n_draws"])
    assert abs(independent - ens["phi_mean"]) < 10 * se

    # Within one run the identity is exact arithmetic, so it holds to machine precision.
    assert ens["phi_mean"] == pytest.approx(ens["phi_mean_from_marginals"], abs=1e-12)


def test_factorised_marginals_give_a_different_answer(pairs, model):
    """The couplings do work: `p_independent` is not a substitute for `p_model`.

    A test that passed under either marginal would be testing nothing about the coupled model.
    """
    ens = energy_ensemble(pairs, model, workers=1, seed=0).row(0, named=True)
    _, _, _, _, _, _, _, q = _prepare(pairs, model)
    phi, _ = site_energies(q, tcren2())
    cp = contact_probabilities(pairs, model, workers=1, seed=0)
    naive = float(cp["p_independent"].to_numpy() @ phi)
    assert abs(naive - ens["phi_mean"]) > 1e-3


def test_variance_and_free_energy_are_not_the_mean(pairs, model):
    """The two read-outs marginals cannot determine are both non-degenerate."""
    r = energy_ensemble(pairs, model, workers=1, seed=0).row(0, named=True)
    assert r["phi_var"] > 0.0
    # -1/beta log <exp(-beta Phi)> <= <Phi> by Jensen, strictly when the spread is non-zero.
    assert r["phi_free"] < r["phi_mean"]


def test_beta_to_zero_recovers_the_mean(pairs, model):
    """`phi_free -> phi_mean` as `beta -> 0`; the free energy is a soft minimum, not a new scale."""
    r = energy_ensemble(pairs, model, beta=1e-6, workers=1, seed=0).row(0, named=True)
    assert r["phi_free"] == pytest.approx(r["phi_mean"], rel=1e-4)


def test_workers_do_not_change_the_answer(pairs, model):
    """Per-structure results are a function of `(seed, pdb.id)`, never of the table split."""
    a = energy_ensemble(pairs, model, workers=1, seed=0)
    b = energy_ensemble(pairs, model, workers=4, seed=0)
    for col in ("phi_mean", "phi_var", "phi_free", "map_entropy", "hamming_obs"):
        assert a[col].to_list() == b[col].to_list()


def test_anchor_choice_does_not_move_the_stationary_distribution(pairs, model):
    """Starting at the observed map shortens burn-in; it does not sample something else."""
    obs = energy_ensemble(pairs, model, anchor="observed", workers=1, seed=0).row(0, named=True)
    rnd = energy_ensemble(pairs, model, anchor="random", workers=1, seed=0).row(0, named=True)
    se = np.sqrt(obs["phi_var"] / obs["n_draws"])
    assert abs(obs["phi_mean"] - rnd["phi_mean"]) < 10 * se


def test_bad_anchor_is_rejected(pairs, model):
    with pytest.raises(ValueError, match="anchor"):
        energy_ensemble(pairs, model, anchor="nearby", workers=1)


# ----------------------------------------------------------------------------- the configurations

def test_coarse_graining_is_what_makes_a_configuration_recur(pairs, model):
    """Site-level maps never repeat; peptide-position engagement patterns do.

    This is the measurement the default `by="position"` rests on. At 174 sites the state space is
    2**174 and every draw is unique, so a site-level `p_hat` carries no information.
    """
    raw = map_ensemble(pairs, model, by="pair", top_k=5, workers=1, seed=0)
    pos = map_ensemble(pairs, model, by="position", top_k=5, workers=1, seed=0)
    assert raw["n_distinct_total"][0] > 10 * pos["n_distinct_total"][0]
    assert raw["p_hat"].max() < 0.01
    assert pos["p_hat"].max() > 0.10


def test_the_modal_engagement_pattern_is_the_observed_one(pairs, model):
    """On 1ao7 the most probable coarse-grained configuration is the crystal's own."""
    me = map_ensemble(pairs, model, by="position", top_k=5, workers=1, seed=0)
    top = me.sort("rank").row(0, named=True)
    assert top["is_observed"]
    assert top["p_hat"] > 0.2


def test_probabilities_are_a_distribution(pairs, model):
    me = map_ensemble(pairs, model, by="position", top_k=1000, workers=1, seed=0)
    assert me["p_hat"].sum() == pytest.approx(1.0, abs=1e-9)
    assert me["rank"].to_list() == sorted(me["rank"].to_list())


def test_bad_by_is_rejected(pairs, model):
    with pytest.raises(ValueError, match="by must be"):
        map_ensemble(pairs, model, by="cell", workers=1)


# ------------------------------------------------------ the sampler against exact enumeration

def _toy_sites(n_pos: int = 4, n_rec: int = 3) -> pl.DataFrame:
    """A small synthetic available-pair set, few enough sites to enumerate all 2**n states."""
    rng = np.random.default_rng(7)
    aa = list("ACDEFGHIKLMNPQRSTVWY")
    rows = []
    for i in range(n_rec):
        for j in range(n_pos):
            rows.append({"pdb.id": "toy", "aa.rec": aa[rng.integers(20)], "chain.rec": "TRB",
                         "region.rec": "CDR3", "pos.rec": i, "aa.par": aa[rng.integers(20)],
                         "pos.par": j, "par.len": n_pos, "role.par": "tcr_facing",
                         "cls": "peptide", "d_heavy": 4.0, "d_ca": 6.0 + 0.5 * abs(i - j),
                         "sigma": float(rng.integers(2))})
    return pl.DataFrame(rows)


def test_sampled_marginals_match_exact_enumeration(model):
    """12 sites: compare the block-Gibbs marginals against the exact Boltzmann sum."""
    sites = _toy_sites()
    _, starts, e, sigma, E, eoff, kv, _ = _prepare(sites, model)
    et, _, A, cols = _one(0, starts[1], e, sigma, E, eoff, kv, 0)
    n = len(et)
    assert n <= 20

    states = ((np.arange(1 << n)[:, None] >> np.arange(n)) & 1).astype(float)
    logp = states @ et + 0.5 * ((states @ A) * states).sum(1) - exact_log_z(et, A)
    p_exact = np.exp(logp) @ states

    from tcren.potts.sample import gibbs
    p_gibbs, _ = gibbs(et, A, cols, _rng_for(0, "toy"), chains=256, burn=200, draws=400, thin=2)
    assert np.max(np.abs(p_exact - p_gibbs)) < 0.02


def test_exact_configuration_probabilities_on_a_toy_interface(model):
    """`by="pair"` p_hat converges to the exact configuration probabilities when n is small."""
    sites = _toy_sites()
    _, starts, e, sigma, E, eoff, kv, _ = _prepare(sites, model)
    et, _, A, cols = _one(0, starts[1], e, sigma, E, eoff, kv, 0)
    n = len(et)
    states = ((np.arange(1 << n)[:, None] >> np.arange(n)) & 1).astype(float)
    p_exact = np.exp(states @ et + 0.5 * ((states @ A) * states).sum(1) - exact_log_z(et, A))

    me = map_ensemble(sites, model, by="pair", top_k=3, chains=256, burn=200, draws=400, thin=2,
                      workers=1, seed=0).sort("rank")
    # The most probable enumerated state must be the most probable sampled one.
    best = states[int(np.argmax(p_exact))]
    got = np.array([int(c) for c in me["pattern"][0]], dtype=float)
    assert np.array_equal(best, got)
    assert me["p_hat"][0] == pytest.approx(float(p_exact.max()), abs=0.02)


# ------------------------------------------------------------------------------------- threading

def test_threading_re_equilibrates_the_contact_ensemble(pairs, model):
    """The number of contacts moves with the candidate — what a frozen marginal cannot express."""
    cand = pl.DataFrame({"pdb.id": ["1ao7"] * 3,
                         "peptide": ["LLFGYPVYV", "WWWWWWWWW", "GGGGGGGGG"]})
    r = threaded_ensemble(pairs, model, cand, chains=16, draws=25, thin=2, seed=0)
    assert r.height == 3
    assert r["n_contacts_model"].max() - r["n_contacts_model"].min() > 0.5


def test_frozen_marginals_are_a_different_number(pairs, model):
    """`coupled=False` is the shipped approximation; it must not agree with the ensemble mean."""
    cand = pl.DataFrame({"pdb.id": ["1ao7"] * 2, "peptide": ["WWWWWWWWW", "GGGGGGGGG"]})
    hot = threaded_ensemble(pairs, model, cand, coupled=True, chains=16, draws=25, seed=0)
    cold = threaded_ensemble(pairs, model, cand, coupled=False, chains=16, draws=25, seed=0)
    assert np.max(np.abs(hot["phi_mean"].to_numpy() - cold["phi_mean"].to_numpy())) > 1e-6


def test_wrong_length_candidate_is_rejected(pairs, model):
    cand = pl.DataFrame({"pdb.id": ["1ao7"], "peptide": ["LLF"]})
    with pytest.raises(ValueError, match="must equal par.len"):
        threaded_ensemble(pairs, model, cand, chains=8, draws=5)


def test_candidate_outside_the_alphabet_is_rejected(pairs, model):
    cand = pl.DataFrame({"pdb.id": ["1ao7"], "peptide": ["LLFGYPVYX"]})
    with pytest.raises(ValueError, match="alphabet"):
        threaded_ensemble(pairs, model, cand, chains=8, draws=5)


def test_undefined_potential_cells_are_counted_not_swallowed(pairs, model):
    """TCRen2 has no `from = C` row; a site it cannot score contributes zero and is reported."""
    _, _, _, _, _, _, _, q = _prepare(pairs, model)
    phi, undefined = site_energies(q, tcren2())
    n_cys = int((q["aa.rec"] == "C").sum())
    assert int(undefined.sum()) == n_cys
    assert np.all(np.isfinite(phi))
    # and it is reported per structure, not as a whole-set total
    ens = energy_ensemble(pairs, model, workers=1, seed=0)
    assert ens["n_undefined"].to_list() == [n_cys]
