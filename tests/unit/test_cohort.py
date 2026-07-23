"""Cohort-relative scores: shape, orientation, and the crystal-calibration contract."""
import numpy as np
import pytest

from tcren.cohort import Q_FEATURES, STRAIN_TERMS, phi_bind, q_score, strain_z, zscore


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
    for fn in (q_score, phi_bind, strain_z):
        assert fn(t).shape == (120,)


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


def test_phi_bind_adds_the_docking_axis_to_q():
    t = _table()
    assert not np.allclose(phi_bind(t), q_score(t))


def test_missing_column_names_itself():
    t = _table()
    del t["n_hbond"]
    with pytest.raises(KeyError, match="n_hbond"):
        q_score(t)


def test_q_uses_exactly_five_descriptors():
    assert len(Q_FEATURES) == 5
