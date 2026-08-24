"""Cohort-relative scores: shape, orientation, and the crystal-calibration contract."""
import numpy as np
import pytest

from tcren.cohort import (F_TERMS, Q_FEATURES, Q_FEATURES_CORE, Q_FEATURES_GEOM, STRAIN_TERMS, f_score, native_reference, q_score, strain_z, zscore)


def _table(n=120, seed=0, **shift):
    rng = np.random.default_rng(seed)
    cols = ["burial", "n_pep_contacted", "chain_balance", "n_hbond", "F_cdr12", "F_cdr3a",
            "pitch", "F_tcr_mhc", "cdr3b_topep", "cdr3b_reach", "extent"]
    t = {c: rng.normal(size=n) for c in cols}
    t["n_contacts_tp"] = np.abs(rng.normal(size=n)) * 10 + 5
    for k, v in shift.items():
        t[k] = t[k] + v
    return t


def test_zscore_is_standardised():
    z = zscore(np.random.default_rng(1).normal(size=500))
    assert abs(float(np.mean(z))) < 1e-9 and abs(float(np.std(z)) - 1.0) < 1e-9


def test_zscore_rank_is_bounded_and_monotone():
    x = np.array([1.0, 2.0, 3.0, 4.0, np.nan])
    r = zscore(x, method="rank")
    assert np.all(r[:4] >= -1) and np.all(r[:4] <= 1)
    assert np.all(np.diff(r[:4]) > 0)      # order-preserving
    assert np.isnan(r[4])                  # NaN in -> NaN out


def test_native_reference_ships_and_has_q_descriptors():
    ref = native_reference()
    for c in Q_FEATURES_GEOM + ("F_cdr12", "F_cdr3a"):
        assert c in ref and len(ref[c]) > 300      # the ~369 Native2026 crystals
        assert np.isfinite(ref[c]).all()


def test_q_score_defined_for_a_single_structure_via_reference():
    # The deployment path: a cohort-relative Q is undefined for n=1 (sd over one row), but standardizing
    # the four geometry terms against the shipped native reference gives a finite score.
    ref = native_reference()
    one = {c: np.array([float(ref[c][0])]) for c in Q_FEATURES_GEOM}
    q = q_score(one, reference=ref, features=Q_FEATURES_GEOM)
    assert q.shape == (1,) and np.isfinite(q[0])


def test_q_score_default_is_directional_decorrelated_k4():
    import inspect
    sig = inspect.signature(q_score)
    assert sig.parameters["features"].default == Q_FEATURES_GEOM and len(Q_FEATURES_GEOM) == 4
    assert sig.parameters["decorrelate"].default is True
    # default equals the manual z @ Cinv @ 1 against the reference
    rng = np.random.default_rng(5)
    feats = Q_FEATURES_GEOM
    ref = {c: rng.normal(size=400) for c in feats}
    ref["n_pep_contacted"] = ref["burial"] * 0.6 + rng.normal(size=400) * 0.5   # induce C != I
    t = {c: rng.normal(size=50) for c in feats}
    q = q_score(t, reference=ref)
    zc = lambda x, r: (np.asarray(x, float) - np.mean(r)) / np.std(r)
    Z = np.vstack([zc(t[c], ref[c]) for c in feats])
    C = np.cov(np.vstack([zc(ref[c], ref[c]) for c in feats]))
    w = np.linalg.pinv(C) @ np.ones(len(feats))
    assert np.allclose(q, (w[:, None] * Z).sum(0))
    # and it differs from the equal-weight mean because the descriptors are correlated
    assert not np.allclose(q, q_score(t, reference=ref, decorrelate=False))


def test_q_score_reduces_to_equal_weight_mean_when_uncorrelated():
    rng = np.random.default_rng(6)
    ref = {c: rng.normal(size=5000) for c in Q_FEATURES_GEOM}     # independent -> C ~ I
    t = {c: rng.normal(size=300) for c in Q_FEATURES_GEOM}
    qd = q_score(t, reference=ref, decorrelate=True)
    qm = q_score(t, reference=ref, decorrelate=False)
    assert np.corrcoef(qd, qm)[0, 1] > 0.999                      # rank-identical to the mean


def test_native_reference_ranks_like_cohort_relative():
    # Transferability: for a cohort drawn from the native distribution (real inputs match the crystal
    # scales to within ~10%), standardizing against the fixed native reference ranks ~identically to
    # standardizing against the cohort itself. (This holds because the scales match; a cohort on a wildly
    # different scale would not transfer -- which is exactly why the reference is the crystal manifold.)
    ref = native_reference()
    rng = np.random.default_rng(3)
    idx = rng.integers(0, len(ref["burial"]), 200)
    t = {c: ref[c][idx] + rng.normal(0, 0.01 * np.std(ref[c]), 200) for c in Q_FEATURES_GEOM}
    qc = q_score(t, features=Q_FEATURES_GEOM)                               # cohort-relative
    qr = q_score(t, reference=ref, features=Q_FEATURES_GEOM)                # fixed native reference
    rho = np.corrcoef(np.argsort(np.argsort(qc)), np.argsort(np.argsort(qr)))[0, 1]
    assert rho > 0.99


def test_zscore_constant_column_does_not_divide_by_zero():
    assert np.all(zscore(np.full(20, 3.7)) == 0.0)


def test_strain_is_zero_mean_when_self_calibrated():
    # Without a reference the cohort standardises to itself, so the gradient is invisible.
    assert abs(float(np.nanmean(strain_z(_table())))) < 1e-9


def test_strain_reads_positive_against_a_crystal_reference():
    # The contract behind the crystal < AF-real < AF-decoy gradient: a cohort whose CDR3beta
    # reaches further from the peptide must read as MORE strained than the reference.
    crystal = _table(seed=0)
    forced = _table(seed=0, cdr3b_topep=2.0, cdr3b_reach=2.0)
    assert float(np.nanmean(strain_z(forced, reference=crystal))) > 0.5


def test_strain_signs_match_the_documented_physics():
    assert dict(STRAIN_TERMS)["chain_balance"] == -1.0   # less balanced = more forced
    assert dict(STRAIN_TERMS)["cdr3b_reach"] == +1.0     # reaching away = more forced


def test_q_core_drops_n_pep_contacted_and_still_scores():
    t = _table()
    assert "n_pep_contacted" not in Q_FEATURES_CORE and len(Q_FEATURES_CORE) == 4
    # the core score is a valid one-per-row score and differs from the 5-term default
    assert q_score(t, features=Q_FEATURES_CORE).shape == (120,)
    assert not np.allclose(q_score(t, features=Q_FEATURES_CORE), q_score(t))


def test_missing_column_names_itself():
    t = _table()
    del t["n_hbond"]
    with pytest.raises(KeyError, match="n_hbond"):
        q_score(t)
    # and the remedy is given as a library call, not only as a CLI invocation
    with pytest.raises(KeyError, match=r"recognition_table\(items, full=True\)"):
        q_score(t)


def test_q_uses_exactly_five_descriptors():
    assert len(Q_FEATURES) == 5


def test_q_features_geom_is_the_four_geometry_terms():
    # Q_geom = Q_FEATURES without the pp_combo energy contrast
    assert Q_FEATURES_GEOM == ("burial", "n_pep_contacted", "chain_balance", "n_hbond")
    assert "pp_combo" not in Q_FEATURES_GEOM and set(Q_FEATURES_GEOM) < set(Q_FEATURES)


def _table_f(n=120, seed=4):
    """Table with the two F contact-energy terms f_score needs."""
    t = _table(n, seed)
    rng = np.random.default_rng(seed + 1)
    t["F_tcr_pep"] = rng.normal(size=n)
    t["F_tcr_mhc"] = rng.normal(size=n)
    return t


def test_f_score_is_binder_oriented_zscore_of_negated_energy():
    t = _table_f()
    # F = z(-(F_tcr_pep + F_tcr_mhc)): lower raw energy -> higher (more binder-like) score
    expected = zscore(-(np.asarray(t["F_tcr_pep"]) + np.asarray(t["F_tcr_mhc"])))
    assert np.allclose(f_score(t), expected)
    assert F_TERMS == ("F_tcr_pep", "F_tcr_mhc")
    # a structure with lower total energy ranks strictly higher
    assert f_score(t)[np.argmin(np.asarray(t["F_tcr_pep"]) + np.asarray(t["F_tcr_mhc"]))] == f_score(t).max()


def test_q_coupled_r_override_defaults_to_the_measured_coupling():
    """`r=None` must reproduce the published path byte-for-byte; a supplied r replaces it."""
    from scipy.special import erf

    from tcren.cohort import coupling, q_coupled

    rng = np.random.default_rng(7)
    q = rng.normal(size=60)
    e = 0.6 * q + 0.8 * rng.normal(size=60)
    assert np.array_equal(q_coupled(q, e), q_coupled(q, e, r=None))
    # supplying the value the cohort would have measured is a no-op
    assert np.allclose(q_coupled(q, e), q_coupled(q, e, r=coupling(zscore(q), zscore(e))))
    # r = 0 collapses the energy factor to the constant 1/2, leaving geometry alone
    geom = 0.5 * (1.0 + erf(zscore(q) / np.sqrt(2.0)))
    assert np.allclose(q_coupled(q, e, r=0.0), 0.5 * geom)
    # a per-row r is accepted and acts row-wise
    rvec = np.linspace(-0.5, 0.5, 60)
    out = q_coupled(q, e, r=rvec)
    assert out.shape == (60,)
    assert np.isclose(out[10], q_coupled(q, e, r=float(rvec[10]))[10])
