"""S_free and the frozen calibration: the n=1 guarantee and the reload path. 2026-08-28"""
import json
from importlib import resources

import numpy as np
import pytest

from tcren import reliability as rel
from tcren.cohort import Q_FEATURES_GEOM, q_score


@pytest.fixture(scope="module")
def ref():
    return rel.reliability_reference()


def test_reference_is_complete_in_every_block_descriptor(ref):
    """np.cov propagates a single NaN into every cell, so the covariance -- and with it Q and T --
    is only defined if the shipped reference is complete-case in all nine descriptors."""
    for c in (*Q_FEATURES_GEOM, *rel.T_FEATURES_TOPO):
        assert c in ref, c
        assert np.isfinite(ref[c]).all(), f"{c} carries a non-finite value"


def test_s_free_is_defined_for_a_single_structure(ref):
    """The property `cohort.p_native` cannot have: it refits per call and raises when rows <=
    features, so it is undefined for one row."""
    one = {k: v[:1] for k, v in ref.items()}
    s = rel.s_free(one, energy=one["neg_energy"])
    assert s.shape == (1,) and np.isfinite(s[0])
    assert np.isfinite(rel.s_free(one)[0])          # and without the energy term


def test_block_native_mean_is_zero_and_sd_is_not_one(ref):
    """Var(Q) = s' C^-1 s, which is not 1 -- so the outer transform is a DIVIDE, and the centring
    it is usually written with is a no-op. Both halves are asserted because the manuscript says so."""
    m = rel.moments()["blocks"]
    for name, feats, sgn in (("Q", Q_FEATURES_GEOM, None),
                             ("T", rel.T_FEATURES_TOPO, rel.T_SIGNS)):
        v = q_score(ref, reference=ref, features=feats, signs=sgn)
        assert abs(np.nanmean(v)) < 1e-9, f"{name} native mean is not 0"
        assert abs(m[name]["sd"] - np.nanstd(v)) < 1e-6
        assert abs(m[name]["sd"] - 1.0) > 0.1, f"{name} sd is ~1; the divide would be a no-op"


def test_signs_change_the_direction():
    """q_score's `signs` is what lets the topology block carry a term that runs the other way."""
    ref = rel.reliability_reference()
    a = q_score(ref, reference=ref, features=rel.T_FEATURES_TOPO, signs=rel.T_SIGNS)
    b = q_score(ref, reference=ref, features=rel.T_FEATURES_TOPO)
    assert not np.allclose(a, b)
    with pytest.raises(ValueError, match="signs has"):
        q_score(ref, reference=ref, features=rel.T_FEATURES_TOPO, signs=(1.0, 1.0))


def test_frozen_calibration_reloads_exactly():
    """The shipped JSON must reproduce the fitted probabilities on reload, or the number a user
    reads is not the number that was validated."""
    raw = json.loads((resources.files("tcren.data") / "reliability_moments.json").read_text())
    for link, c in raw["calibration"].items():
        s = np.linspace(-4, 4, 17)
        want = 1.0 / (1.0 + np.exp(-(c["slope"] * s + c["intercept"])))
        assert np.allclose(rel.p_binder(s, link=link), want, atol=1e-12)
    with pytest.raises(KeyError, match="no frozen link"):
        rel.p_binder(np.zeros(3), link="nope")


def test_p_binder_is_a_probability(ref):
    p = rel.p_binder(rel.s_free(ref, energy=ref["neg_energy"]))
    ok = np.isfinite(p)
    assert ok.sum() > 300 and ((p[ok] >= 0) & (p[ok] <= 1)).all()


def test_af_band_clamps_rather_than_extrapolating():
    b = rel.af_band([-99.0, 0.85, 99.0], reference="binder_bm|ipTM")
    assert b[0]["band"] == 0 and b[-1]["band"] == 9
    assert 0.0 <= b[1]["p_nonbinder"] <= 1.0
    assert b[-1]["ci_lo"] <= b[-1]["p_nonbinder"] <= b[-1]["ci_hi"]
    with pytest.raises(KeyError, match="no band table"):
        rel.af_band([0.8], reference="nope")


def test_catalogues_are_non_empty():
    assert rel.available_links() and rel.available_bands()


def test_t_score_wires_its_own_default_reference(ref):
    """`t_score` is reached transitively through `s_free`, so its own features/signs wiring was
    never exercised: pointing it at the wrong descriptors would only show up as a value drift."""
    direct = rel.t_score(ref)
    hand = q_score(ref, reference=ref, features=rel.T_FEATURES_TOPO, signs=rel.T_SIGNS)
    assert np.allclose(direct, hand, equal_nan=True)
    assert rel.t_score(ref, reference=ref).shape == (len(direct),)


def test_pi_frozen_matches_the_shipped_moments():
    """The constant and the JSON key are two copies of one fact and could drift apart silently."""
    m = rel.moments()
    assert m["pi_frozen"] == rel.PI_FROZEN
    assert rel.PI_FROZEN in m["blocks"]
    assert {"blocks", "calibration", "af_bands", "pi_frozen"} <= set(m)


def test_pi_frozen_is_the_column_potts_emits():
    """The blocker of 2.15.0: `s_free`'s energy block named a column no tcren function wrote.

    Asserted against the columns the scorers actually return, not against their source text: the
    literal moved out of `score_sites` when the per-structure loop was parallelised in 2.17.0, and
    a source grep called that a regression when nothing had changed about the output.
    """
    from test_potts import _synthetic_sites

    from tcren.potts import bound_unbound, fit_potts, score_sites

    sites = _synthetic_sites()
    model = fit_potts(sites, couplings=False)
    for fn, kw in ((score_sites, {"particles": 8, "steps": 4}),
                   (bound_unbound, {"chains": 8, "burn": 4, "draws": 4, "thin": 1,
                                    "particles": 8, "steps": 4})):
        assert rel.PI_FROZEN in fn(sites, model, workers=1, **kw).columns, fn.__name__


def test_inversion_flag_is_defined_for_one_structure_and_needs_the_energy(ref):
    """The forced-pose detector reads the energy block against the two shape blocks, so with no
    energy term there is nothing to invert and it must say so rather than return a number."""
    e = np.asarray(ref[rel.PI_FROZEN], float)
    f = rel.inversion_flag(ref, energy=e)
    assert f.shape == (len(e),) and np.isfinite(f).sum() > 300
    # genuine crystals are the reference manifold, so the flag sits near 0 on them by construction
    assert abs(np.nanmean(f)) < 0.5, np.nanmean(f)
    one = {k: [v[0]] for k, v in ref.items()} if isinstance(ref, dict) else ref[:1]
    assert np.isfinite(rel.inversion_flag(one, reference=ref, energy=e[:1])).all()
    assert np.isnan(rel.inversion_flag(ref)).all(), "no energy must not silently score"


def test_screening_yield_is_the_cut_and_never_a_nan_measurement():
    s = np.arange(100.0)
    y = rel.screening_yield(s, 0.1, prevalence=0.2)
    assert y["n_tested"] == 10 and y["threshold"] == 90.0
    assert y["expected_hits"] == pytest.approx(2.0)
    # enrichment needs labels this function does not have; it must be absent, not NaN
    assert "enrichment" not in y
    assert "expected_hits" not in rel.screening_yield(s, 0.1)
    assert rel.screening_yield(s, 1.0)["n_tested"] == 100
    with pytest.raises(ValueError, match="budget"):
        rel.screening_yield(s, 1.5)
