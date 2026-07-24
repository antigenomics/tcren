"""Cohort-relative scores: shape, orientation, and the crystal-calibration contract."""
import numpy as np
import pytest

from tcren.cohort import (F_TERMS, Q_FEATURES, Q_FEATURES_CORE, Q_FEATURES_GEOM, STRAIN_TERMS, f_score,
                          phi_bind, q_f, q_iptm, q_score, strain_z, zscore)


def _table(n=120, seed=0, **shift):
    rng = np.random.default_rng(seed)
    cols = ["burial", "n_pep_contacted", "chain_balance", "n_hbond", "e_cdr12", "e_cdr3a",
            "pitch", "F_tcr_mhc", "cdr3b_topep", "cdr3b_reach", "extent"]
    t = {c: rng.normal(size=n) for c in cols}
    t["n_contacts_tp"] = np.abs(rng.normal(size=n)) * 10 + 5
    for k, v in shift.items():
        t[k] = t[k] + v
    return t


def test_scores_are_one_value_per_row():
    t = _table()
    assert q_score(t).shape == strain_z(t).shape == (120,)
    with pytest.warns(DeprecationWarning):
        assert phi_bind(t).shape == (120,)


def test_zscore_is_standardised():
    z = zscore(np.random.default_rng(1).normal(size=500))
    assert abs(float(np.mean(z))) < 1e-9 and abs(float(np.std(z)) - 1.0) < 1e-9


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


def test_phi_bind_is_deprecated():
    with pytest.warns(DeprecationWarning, match="q_score"):
        phi_bind(_table())


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


def test_q_uses_exactly_five_descriptors():
    assert len(Q_FEATURES) == 5


def test_q_features_geom_is_the_four_geometry_terms():
    # Q_geom = Q_FEATURES without the pp_combo energy contrast
    assert Q_FEATURES_GEOM == ("burial", "n_pep_contacted", "chain_balance", "n_hbond")
    assert "pp_combo" not in Q_FEATURES_GEOM and set(Q_FEATURES_GEOM) < set(Q_FEATURES)


def test_q_iptm_is_zsum_of_iptm_and_q():
    t = _table()
    iptm = np.random.default_rng(3).uniform(0.4, 0.95, size=120)
    combo = q_iptm(t, iptm, features=Q_FEATURES_GEOM)
    assert combo.shape == (120,)
    # exact definition: z(ipTM) + z(Q_geom)
    expected = zscore(iptm) + zscore(q_score(t, features=Q_FEATURES_GEOM))
    assert np.allclose(combo, expected)
    # a missing ipTM falls back to z(Q) for that row, so the whole vector stays finite and rankable
    iptm2 = iptm.copy(); iptm2[0] = np.nan
    out = q_iptm(t, iptm2, features=Q_FEATURES_GEOM)
    assert np.isfinite(out).all()
    assert np.isclose(out[0], zscore(q_score(t, features=Q_FEATURES_GEOM))[0])
    # all ipTM missing -> rank by the model geometry alone (== z(Q))
    allnan = np.full(120, np.nan)
    assert np.allclose(q_iptm(t, allnan, features=Q_FEATURES_GEOM),
                       zscore(q_score(t, features=Q_FEATURES_GEOM)))


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


def test_q_f_is_zq_plus_signed_zf():
    t = _table_f()
    zq = zscore(q_score(t, features=Q_FEATURES_GEOM))
    assert np.allclose(q_f(t), zq + f_score(t))                 # sign=+1 default -> z(Q)+z(F)
    assert np.allclose(q_f(t, sign=-1.0), zq - f_score(t))      # forced-pose form z(Q)-z(F)
    # the two signs differ exactly by 2 z(F)
    assert np.allclose(q_f(t) - q_f(t, sign=-1.0), 2 * f_score(t))
