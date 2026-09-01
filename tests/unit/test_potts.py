"""Unit tests for the coupled Potts model over the contact map.

The numerical core is checked against things that admit an exact answer: the partition function of
the uncoupled model is closed-form, and with couplings on a small enough map can be enumerated
outright, so ``ais_log_z`` and ``gibbs`` are compared against truth rather than against themselves.
"""

from __future__ import annotations

import json

import numpy as np
import polars as pl
import pytest

from tcren.potts import (AA, PottsModel, ais_log_z, centred_potential, colour,
                         coupling_matrix, edges, energy, exact_log_z,
                         eta, factorised_log_z, fit_potts, gauge, gibbs, irls,
                         kernel_names, kernel_table, neighbour_counts, peptide_free_energy,
                         score_sites, site_codes)


# --------------------------------------------------------------------------- toy lattice


def _toy_graph(n=13, density=0.25, seed=7):
    """A random coupling graph on ``n`` sites, small enough for exact enumeration."""
    rng = np.random.default_rng(seed)
    eta = rng.normal(-0.5, 1.0, n)
    ea, eb = np.triu_indices(n, 1)
    keep = rng.random(len(ea)) < density
    return eta, ea[keep], eb[keep], rng


def test_factorised_log_z_is_closed_form():
    eta = np.array([-1.0, 0.0, 2.0])
    assert factorised_log_z(eta) == pytest.approx(float(np.log1p(np.exp(eta)).sum()))


def test_exact_log_z_matches_factorised_when_uncoupled():
    eta, *_ = _toy_graph(n=10)
    A = np.zeros((10, 10))
    assert exact_log_z(eta, A) == pytest.approx(factorised_log_z(eta), abs=1e-9)


def test_colouring_leaves_no_edge_inside_a_class():
    _, ea, eb, _ = _toy_graph()
    cols = colour(13, ea, eb)
    A = coupling_matrix(13, ea, eb, np.zeros(len(ea), int), np.array([1.0]))
    for cidx in cols:
        assert not A[np.ix_(cidx, cidx)].any()
    assert sum(len(c) for c in cols) == 13          # every site coloured exactly once


def test_gibbs_reproduces_the_closed_form_when_uncoupled():
    eta, *_, rng = _toy_graph(n=12)
    occ, totals = gibbs(eta, np.zeros((12, 12)), [np.arange(12)], rng,
                        chains=3000, burn=50, draws=150, thin=1)
    assert np.max(np.abs(occ - 1.0 / (1.0 + np.exp(-eta)))) < 0.02
    assert totals.mean() == pytest.approx(float((1 / (1 + np.exp(-eta))).sum()), abs=0.1)


def test_ais_is_exact_when_uncoupled():
    eta, *_, rng = _toy_graph(n=12)
    log_z, ess = ais_log_z(eta, np.zeros((12, 12)), [], rng, particles=32, steps=16)
    assert log_z == pytest.approx(factorised_log_z(eta), abs=1e-12)
    assert ess == 32


@pytest.mark.parametrize("k", [1.5, 0.8, -1.0])
def test_ais_matches_exact_enumeration_with_couplings(k):
    """The real test of the annealing: 2^13 states enumerated, against AIS from beta = 0."""
    eta, ea, eb, rng = _toy_graph(n=13)
    A = coupling_matrix(13, ea, eb, np.zeros(len(ea), int), np.array([k]))
    cols = colour(13, ea, eb)
    exact = exact_log_z(eta, A)
    est = [ais_log_z(eta, A, cols, rng, particles=128, steps=256)[0] for _ in range(3)]
    assert np.max(np.abs(np.array(est) - exact)) < 0.12


def test_gibbs_marginals_match_exact_marginals_with_couplings():
    eta, ea, eb, rng = _toy_graph(n=13)
    A = coupling_matrix(13, ea, eb, np.zeros(len(ea), int), np.array([0.9]))
    cols = colour(13, ea, eb)
    s = ((np.arange(1 << 13)[:, None] >> np.arange(13)) & 1).astype(float)
    lw = s @ eta + 0.5 * ((s @ A) * s).sum(1)
    p = np.exp(lw - lw.max())
    exact = (p / p.sum()) @ s
    occ, _ = gibbs(eta, A, cols, rng, chains=1500, burn=150, draws=250, thin=2)
    assert np.max(np.abs(occ - exact)) < 0.025


def test_energy_sign_convention():
    """Lower energy is more favourable, and a positive coupling lowers it for a co-occurring pair."""
    eta = np.array([0.0, 0.0])
    A = np.array([[0.0, 1.0], [1.0, 0.0]])
    both = energy(np.array([1.0, 1.0]), eta, A)
    one = energy(np.array([1.0, 0.0]), eta, A)
    assert both < one


# --------------------------------------------------------------------------- the fit


def test_irls_solves_its_own_score_equation():
    """At the optimum the penalised gradient is zero — the check that needs no other library."""
    import scipy.sparse as sp
    rng = np.random.default_rng(3)
    X = np.column_stack([np.ones(800), rng.normal(size=(800, 3)),
                         rng.integers(0, 3, (800, 2)).astype(float)])
    truth = np.r_[0.2, 0.5, -0.3, 0.1, 0.4, -0.2]
    y = rng.binomial(1, 1 / (1 + np.exp(-(X @ truth)))).astype(float)
    w = np.ones(800)
    ridge = 0.7
    b, _ = irls(sp.csr_matrix(X), y, w, ridge)
    mu = 1 / (1 + np.exp(-(X @ b)))
    lam = np.full(X.shape[1], ridge)
    lam[0] = 0.0
    assert np.max(np.abs(X.T @ (w * (y - mu)) - lam * b)) < 1e-8


def test_irls_recovers_a_known_logistic_at_zero_ridge():
    import scipy.sparse as sp
    rng = np.random.default_rng(11)
    X = np.column_stack([np.ones(20000), rng.normal(size=(20000, 2))])
    truth = np.array([-0.4, 0.9, -0.6])
    y = rng.binomial(1, 1 / (1 + np.exp(-(X @ truth)))).astype(float)
    b, _ = irls(sp.csr_matrix(X), y, np.ones(20000), 0.0, free=range(3))
    assert np.max(np.abs(b - truth)) < 0.06


def test_gauge_is_zero_sum_and_leaves_eta_unchanged():
    rng = np.random.default_rng(5)
    sizes = [4, 4, 3]
    slices = [slice(1, 5), slice(5, 9), slice(9, 12), slice(12, 28)]
    b = rng.normal(size=28)
    codes = [rng.integers(0, 4, 50), rng.integers(0, 4, 50), rng.integers(0, 3, 50)]
    def eta_of(v):
        out = np.full(50, v[0])
        for ck, s in zip(codes, slices[:-1]):
            out += v[s][ck]
        return out + v[slices[-1]].reshape(4, 4)[codes[0], codes[1]]
    bg, J = gauge(b, slices, free_coupling=True)
    assert np.max(np.abs(eta_of(bg) - eta_of(b))) < 1e-12
    assert abs(J.sum(0)).max() < 1e-12 and abs(J.sum(1)).max() < 1e-12
    for s in slices[:-1]:
        assert abs(bg[s].sum()) < 1e-12


def test_centred_potential_is_double_centred_and_sign_flipped():
    for name in ("tcren2", "mj", "mj1996", "keskin", "tcren"):
        M = centred_potential(name)
        assert M.shape == (20, 20)
        assert np.abs(M.mean(0)).max() < 1e-10
        assert np.abs(M.mean(1)).max() < 1e-10
    with pytest.raises(ValueError):
        centred_potential("not-a-potential")


def test_centred_potential_flips_the_sign_of_the_energy():
    """TCRen is negative-is-favourable; a coupling is a log-odds where positive is more likely."""
    from tcren.potential import tcren2
    from tcren.pose import _double_centred
    raw, index = _double_centred(tcren2())
    M = centred_potential("tcren2")
    i, j = index["L"], index["L"]
    assert M[AA.index("L"), AA.index("L")] == pytest.approx(-float(raw[i, j]))


def test_uncentred_pin_is_the_raw_matrix_and_differs_by_one_body_terms_only():
    """The gauge that lets a pinned model reproduce a referenced contact-map score.

    ``reference_delta`` is a difference of one-body sums, so a double-centred pin re-injects the
    potential's partner-residue column mean scaled by the position's contact count and the identity
    fails. Pinned uncentred the coupling IS the potential, so any linear read-out of the field
    reduces to the potential's own score up to the fitted scale.
    """
    from tcren.potential import tcren2
    M, index = tcren2().as_matrix()
    raw = centred_potential("tcren2", centre=False)
    for a in ("L", "K", "W"):
        assert raw[AA.index(a), AA.index(a)] == pytest.approx(-float(M[index[a], index[a]]))

    # what centring removes is exactly the one-body part, on the cells the potential observes
    keep = [k for k, a in enumerate(AA)
            if a in index and np.isfinite(M[index[a]]).any() and np.isfinite(M[:, index[a]]).any()]
    d = (raw - centred_potential("tcren2", centre=True))[np.ix_(keep, keep)]
    assert np.abs(d - d.mean(1, keepdims=True) - d.mean(0, keepdims=True) + d.mean()).max() < 1e-12
    assert d.mean(0).std() > 0.0            # and it is not a constant, so the gauge matters


def test_pin_centred_reaches_the_fit():
    sites = _synthetic_sites(n_struct=12)
    on = fit_potts(sites, coupling_matrix="tcren2")
    off = fit_potts(sites, coupling_matrix="tcren2", pin_centred=False)
    assert on.n_parameters() == off.n_parameters()
    assert not np.allclose(np.asarray(on.coupling_array()), np.asarray(off.coupling_array()))


# --------------------------------------------------------------------------- sites and kernel


def _synthetic_sites(n_struct=6, n_rec=8, n_par=9, seed=2) -> pl.DataFrame:
    """A grid of sites with a plausible contact pattern — no structures, no annotation."""
    rng = np.random.default_rng(seed)
    rows = []
    for k in range(n_struct):
        for i in range(n_rec):
            for j in range(n_par):
                d = 4.0 + abs(i - 3) + abs(j - 4) + rng.normal(0, 0.5)
                rows.append({
                    "pdb.id": f"s{k}", "aa.rec": AA[rng.integers(0, 20)], "chain.rec": "TRB",
                    "region.rec": "CDR3", "pos.rec": i, "aa.par": AA[rng.integers(0, 20)],
                    "pos.par": j, "role.par": "tcr_facing", "cls": "peptide",
                    "d_heavy": None, "d_ca": float(max(d, 3.5)),
                    "sigma": float(d < 6.5)})
    return pl.DataFrame(rows)


def test_edges_are_symmetric_offsets_within_one_loop():
    sites = _synthetic_sites()
    _, _, q = site_codes(sites)
    ed = edges(q, joint=False)
    assert len(ed) == len(kernel_names(joint=False))
    # (0,+1) links every pair of adjacent partner positions inside one receptor position
    n_adj = ed[0].shape[0]
    assert n_adj == 6 * 8 * 8              # structures x receptor positions x adjacent partners
    for e in ed:
        assert e.ndim == 2 and e.shape[1] == 2
        assert (e[:, 0] != e[:, 1]).all()  # no self-edge


def test_neighbour_counts_count_both_directions():
    sites = _synthetic_sites(n_struct=1, n_rec=3, n_par=3)
    _, _, q = site_codes(sites)
    ed = edges(q, joint=False)
    sigma = np.ones(q.height)
    counts = neighbour_counts(ed, sigma, q.height)
    # a site in the middle of the 3x3 grid has both (0,-1) and (0,+1) neighbours
    mid = q.filter((pl.col("pos.rec") == 1) & (pl.col("pos.par") == 1))["sid"][0]
    assert counts[mid, 0] == 2
    assert counts.shape == (q.height, len(ed))


def test_site_codes_bins_are_global_not_relative():
    """A model fitted on one set must index the same distance bins on another."""
    a = _synthetic_sites(seed=1)
    b = _synthetic_sites(seed=2).with_columns(pl.col("d_ca") + 2.0)
    ca, _, _ = site_codes(a)
    cb, _, _ = site_codes(b)
    assert ca[2].min() != cb[2].min()      # different data -> different bins occupied
    # but the mapping itself is absolute: 7.2 A always lands in bin 14 at DBIN = 0.5
    one = a.head(1).with_columns(pl.lit(7.2).alias("d_ca"))
    assert site_codes(one)[0][2][0] == 14


# --------------------------------------------------------------------------- fit and round-trip


def test_fit_recovers_a_gauge_consistent_model():
    model = fit_potts(_synthetic_sites(n_struct=12), notes="synthetic")
    assert model.n_structures == 12
    assert model.n_contacts > 0
    assert abs(np.sum(model.h_rec)) < 1e-9 and abs(np.sum(model.h_par)) < 1e-9
    J = np.asarray(model.coupling)
    assert np.abs(J.sum(0)).max() < 1e-9 and np.abs(J.sum(1)).max() < 1e-9
    assert len(model.kernel) == len(kernel_names(joint=False))
    assert model.n_parameters() == 1 + 20 + 20 + 31 + 15 + 6 + 2 + 400 + len(model.kernel)


def test_no_couplings_gives_the_factorised_model():
    model = fit_potts(_synthetic_sites(), couplings=False)
    assert model.kernel == []
    # with no couplings the partition function is closed form, so scoring needs no sampling
    from tcren.potts import score_sites
    sc = score_sites(_synthetic_sites(), model, particles=8, steps=4)
    assert np.allclose(sc["log_z"].to_numpy(), sc["log_z0"].to_numpy())


def test_fixed_coupling_matrix_is_one_parameter():
    free = fit_potts(_synthetic_sites(n_struct=12))
    fixed = fit_potts(_synthetic_sites(n_struct=12), coupling_matrix="mj")
    assert fixed.coupling is None and fixed.coupling_matrix_name == "mj"
    assert fixed.beta_matrix is not None
    assert fixed.n_parameters() == free.n_parameters() - 399
    # the free fit has 400 cells to spend, so it cannot fit worse
    assert free.pseudo_loglik >= fixed.pseudo_loglik
    assert np.asarray(fixed.coupling_array()).shape == (20, 20)


def test_model_round_trips_through_json(tmp_path):
    model = fit_potts(_synthetic_sites())
    path = tmp_path / "m.json"
    model.to_json(path)
    back = PottsModel.from_json(path)
    assert back.pseudo_loglik == pytest.approx(model.pseudo_loglik)
    assert np.allclose(back.coupling_array(), model.coupling_array())
    assert back.kernel == model.kernel
    assert json.loads(path.read_text())["alphabet"] == "".join(AA)


def test_kernel_table_is_sorted_by_absolute_z():
    """Strongest first, and a class with no edges in the data (z is NaN) sorts to the bottom."""
    table = kernel_table(fit_potts(_synthetic_sites(n_struct=12)))
    assert table.height == len(kernel_names(joint=False))
    z = np.nan_to_num(table["z"].abs().to_numpy(), nan=-1.0)
    assert np.all(np.diff(z) <= 1e-12)
    assert np.isnan(table["z"].to_numpy()[-1])          # synthetic data has one loop only


# --------------------------------------------------------------------------- scoring


def test_scores_and_probabilities_are_consistent():
    from tcren.potts import contact_probabilities, score_sites
    sites = _synthetic_sites(n_struct=4)
    model = fit_potts(sites)
    sc = score_sites(sites, model, particles=32, steps=64)
    assert sc.height == 4
    assert (sc["log_lik"] <= 0).all()                    # a log probability
    assert (sc["psi"] == sc["log_lik"] / sc["n_sites"]).all()
    assert (sc["log_z"] >= sc["log_z0"] - 1.0).all()     # couplings here are net attractive
    probs = contact_probabilities(sites, model, chains=32, burn=20, draws=40, thin=2)
    assert probs.height == sites.height
    for col in ("p_independent", "p_model", "p_conditional"):
        v = probs[col].to_numpy()
        assert ((v >= 0) & (v <= 1)).all()
    # The intercept is the one unpenalised column, so its score equation is sum(sigma - mu) = 0
    # and the fitted conditional probability averages to the observed contact rate EXACTLY.
    assert probs["p_conditional"].mean() == pytest.approx(sites["sigma"].mean(), abs=1e-9)


def test_sampled_maps_reproduce_the_observed_contact_rate():
    """The generative check, on enough data that the ridge is not the dominant term."""
    from tcren.potts import contact_probabilities
    sites = _synthetic_sites(n_struct=40)
    model = fit_potts(sites)
    probs = contact_probabilities(sites, model, chains=48, burn=60, draws=80, thin=2)
    assert probs["p_model"].mean() == pytest.approx(sites["sigma"].mean(), abs=0.04)


def test_connected_correlations_recover_the_data_they_were_fitted_to():
    """The generative test the pseudolikelihood never optimises: two-point structure."""
    from tcren.potts import connected_correlations
    sites = _synthetic_sites(n_struct=40)
    model = fit_potts(sites)
    cc = connected_correlations(sites, model, chains=48, burn=60, draws=80, thin=2)
    assert cc.height == len(kernel_names(joint=False))
    seen = cc.filter(pl.col("n_edges") > 100)
    assert seen.height >= 6
    r = np.corrcoef(seen["c_data"].to_numpy(), seen["c_model"].to_numpy())[0, 1]
    assert r > 0.8, r


def test_bundled_model_loads_and_is_gauge_consistent():
    model = PottsModel.bundled()
    assert model.n_structures == 362 and model.n_contacts == 7865
    assert model.radius == 15.0 and model.cutoff == 5.0
    assert abs(np.sum(model.h_rec)) < 1e-9
    J = model.coupling_array()
    assert np.abs(J.sum(0)).max() < 1e-9
    # the couplings the crystals show: axial classes positive, off-axis negative
    k = dict(zip(kernel_names(model.joint), model.kernel))
    assert k["K(+1,+0)"] > 0.5 and k["K(+0,+1)"] > 0.5
    assert k["K(+1,+1)"] < -0.5 and k["K(+1,-1)"] < -0.5


# --------------------------------------------------------------------- bound versus unbound


def test_tilt_mean_is_the_plain_mean_at_zero_and_monotone():
    from tcren.potts import tilt_mean
    t = np.array([3.0, 7.0, 7.0, 11.0, 19.0])
    assert abs(tilt_mean(t, 0.0) - t.mean()) < 1e-12
    ms = [tilt_mean(t, m) for m in np.linspace(-4, 4, 40)]
    assert all(b > a for a, b in zip(ms, ms[1:])), "<N>_mu must increase in mu"
    # the tilt saturates at the support
    assert abs(tilt_mean(t, 40.0) - t.max()) < 1e-9
    assert abs(tilt_mean(t, -40.0) - t.min()) < 1e-9


def test_mu_star_round_trips_a_known_tilt():
    """mu_star inverts tilt_mean, which is the whole claim the score rests on."""
    from tcren.potts import mu_star, tilt_mean
    rng = np.random.default_rng(0)
    t = rng.integers(0, 40, 5000).astype(float)
    for mu in (-1.0, -0.25, 0.0, 0.25, 1.0):
        assert abs(mu_star(t, tilt_mean(t, mu)) - mu) < 1e-6
    assert np.isnan(mu_star(t, 1e6)), "outside the support this is extrapolation, not reweighting"


def test_delta_f_empty_against_exact_enumeration():
    """log[P(N>=1)/P(N=0)] = log(Z-1) exactly, because E(empty) = 0."""
    from tcren.potts import delta_f_empty, exact_log_z
    rng = np.random.default_rng(3)
    n = 12
    eta = rng.normal(-1.5, 1.0, n)
    A = np.zeros((n, n))
    iu = np.triu_indices(n, 1)
    A[iu] = rng.normal(0.0, 0.4, len(iu[0]))
    A = A + A.T
    lz = exact_log_z(eta, A)
    # brute force: P(N = 0) is the single empty configuration, whose energy is 0
    s = ((np.arange(1 << n)[:, None] >> np.arange(n)) & 1).astype(float)
    lp = s @ eta + 0.5 * ((s @ A) * s).sum(1) - lz
    p0 = float(np.exp(lp[s.sum(1) == 0]).sum())
    assert abs(delta_f_empty(lz) - np.log((1 - p0) / p0)) < 1e-9


def test_delta_f_threshold_is_the_model_logit_and_z_free():
    from tcren.potts import delta_f_threshold
    t = np.array([1.0, 2.0, 3.0, 4.0])          # P(N>=3) = 1/2
    assert abs(delta_f_threshold(t, 3)) < 1e-12
    assert delta_f_threshold(t, 2) > 0 and delta_f_threshold(t, 4) < 0
    assert delta_f_threshold(t, 0) == float("inf")


def test_count_free_energy_is_a_normalised_histogram():
    from tcren.potts import count_free_energy
    t = np.array([2.0, 3.0, 3.0, 4.0, 4.0, 4.0])
    n, f = count_free_energy(t)
    assert list(n) == [2, 3, 4]
    assert abs(np.exp(-f).sum() - 1.0) < 1e-12
    assert np.argmin(f) == 2, "F(N) is minimal at the modal count"


def test_gibbs_observer_sees_every_kept_draw():
    from tcren.potts import gibbs
    rng = np.random.default_rng(1)
    eta = np.array([0.4, -0.2, 0.9, -1.1])
    A = np.zeros((4, 4))
    seen = []
    _, tot = gibbs(eta, A, [np.arange(4)], rng, chains=16, burn=4, draws=5, thin=2,
                   observer=lambda S: seen.append(S.sum(1).copy()))
    assert len(seen) == 5
    assert np.array_equal(np.concatenate(seen), tot)


def test_bound_unbound_columns_and_identity():
    """df_empty must equal the log(Z-1) of the log_z the same call reports."""
    from tcren.potts import bound_unbound, delta_f_empty, fit_potts
    sites = _synthetic_sites()
    model = fit_potts(sites)
    b = bound_unbound(sites, model, threshold=2, chains=16, burn=20, draws=20, thin=2,
                      particles=16, steps=32, seed=0)
    assert b.height == sites["pdb.id"].n_unique()
    for c in ("df_empty", "df_threshold", "mu_star", "log_z", "n_lo", "n_hi", "n_var"):
        assert c in b.columns
    got = b["df_empty"].to_numpy()
    want = np.array([delta_f_empty(z) for z in b["log_z"].to_numpy()])
    assert np.allclose(got, want, atol=1e-12)


def test_reweighting_matches_direct_simulation_of_the_tilt():
    """<N>_mu by reweighting must equal <N> from a chain actually run at that tilt.

    A LINEAR tilt in N is exactly a constant added to every field, `E - mu N = -(eta + mu).sigma`,
    so the tilted model can be sampled directly and the reweighting identity is checkable rather
    than assumed. (A quadratic tilt is not a field shift; that is what makes it a genuine global
    term.)
    """
    from tcren.potts import colour, coupling_matrix, gibbs, tilt_mean
    rng = np.random.default_rng(11)
    n = 16
    eta = rng.normal(-1.2, 0.8, n)
    _, ea, eb, _ = _toy_graph(n=n, density=0.3, seed=5)
    A = coupling_matrix(n, ea, eb, np.zeros(len(ea), np.int64), np.array([0.6]))
    cols = colour(n, ea, eb)
    _, base = gibbs(eta, A, cols, np.random.default_rng(0), chains=256, burn=200,
                    draws=400, thin=2)
    for mu in (-0.5, 0.0, 0.5):
        _, direct = gibbs(eta + mu, A, cols, np.random.default_rng(1), chains=256, burn=200,
                          draws=400, thin=2)
        assert abs(tilt_mean(base, mu) - direct.mean()) < 0.25, (
            f"reweighting and direct simulation disagree at mu = {mu}")


# --------------------------------------------------------------------------- the sampled outputs


def test_sample_maps_returns_one_row_per_draw_per_structure():
    """`test_sampled_maps_reproduce_the_observed_contact_rate` reads as if it covers this, but it
    calls `contact_probabilities`. `sample_maps` itself had no test."""
    from tcren.potts import sample_maps

    sites = _synthetic_sites(n_struct=3)
    m = fit_potts(sites, couplings=False)
    chains, draws = 8, 20
    d = sample_maps(sites, m, chains=chains, burn=20, draws=draws, thin=2, seed=0)
    assert set(d["pdb.id"].unique()) == {"s0", "s1", "s2"}
    per = d.group_by("pdb.id").len().sort("pdb.id")["len"].to_list()
    assert per == [chains * draws] * 3, "one row per chain per kept draw"
    assert {"pdb.id", "n_contacts_observed", "n_contacts_sampled"} == set(d.columns)
    n = d["n_contacts_sampled"].to_numpy()
    assert (n >= 0).all() and (n <= 8 * 9).all()          # bounded by the available set
    assert n.std() > 0, "the generative check needs the spread, not only the mean"


def test_count_profile_is_a_normalised_landscape_beside_the_observed_counts():
    from tcren.potts import count_profile

    sites = _synthetic_sites(n_struct=4)
    m = fit_potts(sites, couplings=False)
    d = count_profile(sites, m, chains=8, burn=20, draws=30, thin=2, seed=0)
    assert {"n_contacts", "f_model", "f_observed", "count_observed"} == set(d.columns)
    assert np.exp(-d["f_model"].to_numpy()).sum() == pytest.approx(1.0)
    assert d["count_observed"].sum() == 4                 # one observation per structure
    assert d["n_contacts"].is_sorted()


def test_eta_is_the_one_body_log_odds_of_the_site_table():
    """`eta` is only ever reached through `_prepare`; its own contract is a length and a finite."""
    from tcren.potts import eta as eta_fn

    sites = _synthetic_sites(n_struct=2)
    m = fit_potts(sites, couplings=False)
    codes, _, _ = site_codes(sites, m)
    e = eta_fn(codes, m)
    assert e.shape == (sites.height,) and np.isfinite(e).all()


def test_score_sites_emits_the_energy_block_reliability_consumes():
    """`neg_energy` is `s_score`'s Pi block. 2.15.0 emitted only `energy`, the opposite sign, so the
    three-block score was unreachable from the package."""
    sites = _synthetic_sites(n_struct=3)
    m = fit_potts(sites, couplings=False)
    d = score_sites(sites, m, particles=16, steps=32, seed=0)
    assert "neg_energy" in d.columns
    assert np.allclose(d["neg_energy"].to_numpy(), -d["energy"].to_numpy())


def test_log_z_is_a_function_of_the_structure_not_of_its_row_position():
    """The sampled quantities must not depend on how the pair table was assembled.

    Before 2.17.0 every sampler built one `np.random.default_rng(seed)` OUTSIDE its loop over
    structures, so each structure's AIS run consumed the stream the previous structures had
    advanced: the same PDB scored on its own, in a subset, or in a reordered frame came back with a
    different `log Z`. Row order inside a structure mattered too, because it sets the colouring.
    `score.py` now derives the generator from `(seed, pdb.id)` and sorts sites on their own
    identity, so all four of these agree exactly.
    """
    sites = _synthetic_sites(n_struct=4)
    m = fit_potts(sites, couplings=True)
    kw = dict(particles=16, steps=32, seed=0)
    ids = sorted(sites["pdb.id"].unique().to_list())
    ref = score_sites(sites, m, **kw).sort("pdb.id")["log_z"].to_numpy()

    for label, got in (
        ("reversed rows", score_sites(sites.reverse(), m, **kw).sort("pdb.id")),
        ("shuffled rows", score_sites(sites.sample(fraction=1.0, shuffle=True, seed=11), m, **kw)
                          .sort("pdb.id")),
        ("one structure at a time",
         pl.concat([score_sites(sites.filter(pl.col("pdb.id") == i), m, **kw) for i in ids])
           .sort("pdb.id")),
        ("a two-structure subset",
         score_sites(sites.filter(pl.col("pdb.id").is_in(ids[:2])), m, **kw).sort("pdb.id")),
    ):
        n = len(got)
        assert np.allclose(ref[:n] if n < len(ref) else ref, got["log_z"].to_numpy()), label


def test_the_seed_still_selects_the_stream():
    """Per-structure seeding must not silently ignore the seed argument."""
    sites = _synthetic_sites(n_struct=3)
    m = fit_potts(sites, couplings=True)
    a = score_sites(sites, m, particles=16, steps=32, seed=0)["log_z"].to_numpy()
    b = score_sites(sites, m, particles=16, steps=32, seed=1)["log_z"].to_numpy()
    assert not np.allclose(a, b)


def test_splitting_the_work_across_processes_changes_no_number():
    """`workers` is a scheduling knob, never a scientific one.

    Parallelising the per-structure loop is only sound because `_rng_for` derives each structure's
    generator from `(seed, pdb.id)` and `_prepare` sorts sites on their own identity — so a
    structure's numbers do not depend on which chunk it landed in, or on how many chunks there
    were. Measured on 616 TCRvdb structures: score_sites 4.2x and bound_unbound 5.6x at 8
    processes, both with max|diff| exactly 0.
    """
    from tcren.potts import bound_unbound, score_sites

    sites = _synthetic_sites()
    model = fit_potts(sites, couplings=False)
    n = sites["pdb.id"].n_unique()
    workers = max(2, min(4, n // 2))
    for fn, kw in ((score_sites, {"particles": 8, "steps": 4}),
                   (bound_unbound, {"chains": 8, "burn": 4, "draws": 4, "thin": 1,
                                    "particles": 8, "steps": 4})):
        one = fn(sites, model, workers=1, **kw)
        many = fn(sites, model, workers=workers, **kw)
        assert one.equals(many), f"{fn.__name__} depends on the number of workers"


# --------------------------------------------------------------------------- contact_map
# `contact_probabilities` is per residue pair. What an experiment measures is coarser: a CDR loop
# against a peptide position (an MD contact-frequency map), or a peptide position on its own
# (residue importance). `contact_map` closes the pairs onto those grids, and the closure is the
# only arithmetic it adds.


def _map_inputs(n_struct=3):
    sites = _synthetic_sites(n_struct=n_struct)
    return sites, fit_potts(sites, couplings=True)


@pytest.mark.parametrize("by", ["loop", "position"])
def test_contact_map_closes_the_pairs_with_poisson_binomial_p_at_least_one(by):
    """`p_any` must be exactly `1 - prod(1 - p)` over the group, computed here by hand.

    The residues of a loop are distinct pairs with different marginals, so the count of
    simultaneous contacts is Poisson-binomial and has no closed form. Only "at least one" does, and
    it is what a frequency map measures.
    """
    from tcren.potts import contact_map, contact_probabilities

    sites, model = _map_inputs()
    kw = dict(chains=8, burn=4, draws=8, thin=1, seed=0, workers=1)
    pairs = contact_probabilities(sites, model, **kw)
    got = contact_map(sites, model, by=by, **kw)

    keys = ["pdb.id"] + (["region.rec", "pos.par", "aa.par"] if by == "loop"
                         else ["pos.par", "aa.par"])
    assert got.height == pairs.select(keys).n_unique()
    for row in got.iter_rows(named=True):
        grp = pairs
        for k in keys:
            grp = grp.filter(pl.col(k) == row[k])
        p = grp["p_model"].to_numpy()
        assert row["p_any"] == pytest.approx(1.0 - np.prod(1.0 - p), abs=1e-12)
        assert row["p_expected"] == pytest.approx(p.sum(), abs=1e-12)
        assert row["n_pairs"] == len(p)
        assert row["n_observed"] == int(grp["sigma"].sum())
        assert row["observed"] == int(grp["sigma"].max() > 0)


def test_contact_map_by_pair_is_contact_probabilities_unchanged():
    """The passthrough must not quietly reshape anything — the benchmark joins on these columns."""
    from tcren.potts import contact_map, contact_probabilities

    sites, model = _map_inputs()
    kw = dict(chains=8, burn=4, draws=8, thin=1, seed=0, workers=1)
    assert contact_map(sites, model, by="pair", **kw).equals(contact_probabilities(sites, model,
                                                                                   **kw))


def test_a_certain_pair_forces_its_whole_group_to_one_without_nan():
    """`p = 1` sends `log(1 - p)` to -inf, which is the correct answer and must not become `nan`.

    Accumulating in log space is what keeps a twelve-residue loop from underflowing; the price is
    that a saturated pair hits the boundary, so the boundary is pinned here.
    """
    from tcren.potts import score as score_mod

    frame = pl.DataFrame({"pdb.id": ["x", "x", "y", "y"],
                          "region.rec": ["CDR3"] * 4,
                          "pos.par": [1, 1, 2, 2], "aa.par": ["A", "A", "G", "G"],
                          "p_model": [1.0, 0.5, 0.25, 0.5], "sigma": [1.0, 0.0, 0.0, 0.0]})
    # Drive the real aggregation, with the sampler stubbed out: reaching p = 1 through Gibbs takes
    # a contrived interface, whereas the boundary itself is what needs pinning.
    orig = score_mod.contact_probabilities
    try:
        score_mod.contact_probabilities = lambda *a, **k: frame
        got = score_mod.contact_map(frame, None, by="loop").sort("pdb.id")
    finally:
        score_mod.contact_probabilities = orig
    assert got["p_any"][0] == 1.0                       # saturated, not nan
    assert np.isfinite(got["p_any"].to_numpy()).all()
    assert got["p_any"][1] == pytest.approx(1.0 - 0.75 * 0.5)


def test_contact_map_rejects_an_unknown_grouping():
    from tcren.potts import contact_map

    sites, model = _map_inputs(n_struct=1)
    with pytest.raises(ValueError, match="by must be one of"):
        contact_map(sites, model, by="residue", chains=4, burn=2, draws=2, thin=1, workers=1)


def test_contact_map_does_not_depend_on_the_number_of_workers():
    """Same contract as `score_sites`: `workers` is a scheduling knob, never a scientific one."""
    from tcren.potts import contact_map

    sites, model = _map_inputs(n_struct=6)
    kw = dict(by="loop", chains=8, burn=4, draws=8, thin=1, seed=0)
    assert contact_map(sites, model, workers=1, **kw).equals(
        contact_map(sites, model, workers=3, **kw))


# --------------------------------------------------------------------------- peptide_free_energy
#
# `_synthetic_sites` draws `aa.par` per ROW, so one partner position carries several residues --
# which no real peptide does, and which `peptide_free_energy` now rejects. These tests give the
# partner chain one sequence, the way `available_pairs` reads it off a structure.


def _sequence_sites(n_struct=1):
    sites, _ = _map_inputs(n_struct=n_struct)
    seq = {j: AA[(3 * j + 1) % 20] for j in sites["pos.par"].unique()}
    sites = sites.with_columns(pl.col("pos.par").replace_strict(seq).alias("aa.par"))
    return sites, fit_potts(sites, couplings=True)


def test_peptide_free_energy_reproduces_score_sites_log_z0_at_the_observed_sequence():
    """The whole-interface log Z0 is the same number `score_sites` already emits.

    `peptide_free_energy` rebuilds log Z0 as (every site at its own residue) minus (the sites at
    one position) plus (those sites under the candidate). At the residue the structure carries,
    that has to collapse back onto the untouched sum -- otherwise the constant is wrong and every
    reported free energy is offset by it.
    """
    pairs, model = _sequence_sites()
    f = peptide_free_energy(pairs, model)
    ref = float(score_sites(pairs, model)["log_z0"][0])
    obs = f.filter(pl.col("is_observed") == 1)
    assert obs.height == pairs["pos.par"].n_unique()
    assert np.allclose(obs["log_z0"].to_numpy(), ref, atol=1e-9)


def test_peptide_free_energy_is_equimolar_referenced():
    """dF sums to zero over the twenty residues at every position, by construction."""
    pairs, model = _sequence_sites()
    f = peptide_free_energy(pairs, model)
    per_pos = f.group_by("pos.par").agg(pl.col("dF").sum().alias("s"), pl.len().alias("n"))
    assert set(per_pos["n"]) == {20}
    assert np.allclose(per_pos["s"].to_numpy(), 0.0, atol=1e-10)


def test_peptide_free_energy_is_additive_over_positions():
    """log Z0 is a sum over independent sites, so a whole peptide is the sum of its own cells.

    This is what lets one L x 20 table score both a response-matrix cell and a whole library
    peptide; if it ever stops holding, the library arm is silently scoring something else.
    """
    pairs, model = _sequence_sites()
    f = peptide_free_energy(pairs, model)
    wide = {(r["pos.par"], r["aa.par"]): r["dF"] for r in f.iter_rows(named=True)}
    positions = sorted({p for p, _ in wide})
    # take the observed residue everywhere but one position, and check the two routes agree
    obs = {r["pos.par"]: r["aa.par"] for r in f.filter(pl.col("is_observed") == 1).iter_rows(named=True)}
    for a in ("A", "W", "D"):
        swapped = dict(obs) | {positions[0]: a}
        direct = sum(wide[(p, swapped[p])] for p in positions)
        stepwise = sum(wide[(p, obs[p])] for p in positions) - wide[(positions[0], obs[positions[0]])] \
            + wide[(positions[0], a)]
        assert direct == pytest.approx(stepwise)


def test_peptide_free_energy_coupled_is_the_linear_response_of_the_uncoupled_one():
    """d(log Z)/d(eta) = <sigma>, so with the UNCOUPLED marginals the coupled arm must agree with
    a finite difference of `factorised_log_z`. Checked against p_independent rather than p_model,
    because that is the marginal the closed form actually has."""
    pairs, model = _sequence_sites()
    codes, _sizes, _q = site_codes(pairs, model)
    e = eta(codes, model)
    marg = pl.DataFrame({"p_model": 1.0 / (1.0 + np.exp(-e))})
    got = peptide_free_energy(pairs, model, coupled=True, marginals=marg)

    # a small perturbation of one position's residue, scored both ways
    rng = np.random.default_rng(0)
    d = 1e-6 * rng.standard_normal(len(e))
    fd = (factorised_log_z(e + d) - factorised_log_z(e)) / 1e-6
    lr = float(np.sum((1.0 / (1.0 + np.exp(-e))) * d) / 1e-6)
    assert fd == pytest.approx(lr, rel=1e-4)
    assert got.height == 20 * pairs["pos.par"].n_unique()


def test_peptide_free_energy_rejects_marginals_of_the_wrong_length():
    pairs, model = _sequence_sites()
    with pytest.raises(ValueError, match="one row per site"):
        peptide_free_energy(pairs, model, coupled=True,
                            marginals=pl.DataFrame({"p_model": [0.5, 0.5]}))


def test_peptide_free_energy_rejects_a_position_carrying_two_residues():
    """One partner position is one residue. A frame that disagrees has no sequence to thread."""
    pairs, model = _map_inputs(n_struct=1)          # aa.par drawn per row, so positions clash
    with pytest.raises(ValueError, match="more than one residue"):
        peptide_free_energy(pairs, model)
