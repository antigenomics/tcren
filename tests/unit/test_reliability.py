"""S_free and the frozen calibration: the n=1 guarantee and the reload path. 2026-08-28"""
import json
from importlib import resources

import numpy as np
import polars as pl
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
    """The property the discarded cohort posterior could not have: it refit per call and raised when rows <=
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


def _corr_inputs():
    """The native crystals, their energy and contact count -- everything a correction reads."""
    ref = rel.reliability_reference()
    return ref, np.asarray(ref[rel.PI_FROZEN], float), np.asarray(ref["n_contacts"], float)


def test_every_frozen_correction_resolves_and_an_unknown_one_raises():
    assert rel.available_corrections() == ["binder_bm|ipTM", "binder_bm|pLDDT",
                                           "tcrvdb|ipTM", "tcrvdb|pLDDT"]
    ref, e, n = _corr_inputs()
    for k in rel.available_corrections():
        out = rel.correct_confidence(ref, np.full(len(e), 0.85), reference=k, energy=e, contacts=n)
        assert np.isfinite(out["p_corrected"]).any()
    with pytest.raises(KeyError, match="no frozen correction"):
        rel.correct_confidence(ref, np.full(len(e), 0.85), reference="nope|ipTM")


def test_the_decomposition_reconstructs_the_corrected_probability_exactly():
    """`delta_logit` is what the structure added, so the two probabilities differ by exactly it.

    This is the property the CLI reports on, and the reason the correction is readable at all: a
    caller can see whether a number moved because of the generator or because of the coordinates.
    """
    ref, e, n = _corr_inputs()
    r = rel.correct_confidence(ref, np.full(len(e), 0.85), energy=e, contacts=n)
    lo_c = np.log(r["p_confidence"] / (1 - r["p_confidence"]))
    lo_k = np.log(r["p_corrected"] / (1 - r["p_corrected"]))
    ok = np.isfinite(lo_c) & np.isfinite(lo_k)
    assert np.allclose((lo_k - lo_c)[ok], r["delta_logit"][ok])


def test_a_better_structure_is_corrected_upwards_at_the_same_confidence():
    """Hold the generator's confidence fixed and vary only the coordinates."""
    ref, e, n = _corr_inputs()
    conf = np.full(len(e), 0.85)
    r = rel.correct_confidence(ref, conf, energy=e, contacts=n)
    s = r["s_free"]
    ok = np.isfinite(s) & np.isfinite(r["p_corrected"])
    good, bad = s[ok] >= np.nanpercentile(s[ok], 75), s[ok] <= np.nanpercentile(s[ok], 25)
    assert r["p_corrected"][ok][good].mean() > r["p_corrected"][ok][bad].mean()
    # and the confidence-only reading cannot tell them apart, because it never saw the structure
    assert np.allclose(r["p_confidence"][ok][good].mean(), r["p_confidence"][ok][bad].mean())


def test_the_correction_is_defined_for_a_single_structure():
    """The property the whole score family rests on: no cohort, no refit, n = 1 is enough."""
    ref, e, n = _corr_inputs()
    one = {k: np.asarray(v)[:1] for k, v in ref.items()}
    r = rel.correct_confidence(one, np.array([0.9]), energy=e[:1], contacts=n[:1])
    assert r["p_corrected"].shape == (1,)
    assert np.isfinite(r["p_corrected"]).all()
    # and it agrees with the same row scored inside the whole set
    whole = rel.correct_confidence(ref, np.full(len(e), 0.9), energy=e, contacts=n)
    assert np.allclose(r["p_corrected"][0], whole["p_corrected"][0])


def test_the_contact_term_drops_out_rather_than_being_imputed():
    """No contact count is a smaller model, not a wrong one -- and it says so with NaN."""
    ref, e, n = _corr_inputs()
    conf = np.full(len(e), 0.85)
    with_n = rel.correct_confidence(ref, conf, energy=e, contacts=n)
    without = rel.correct_confidence(ref, conf, energy=e, contacts=None)
    assert np.isnan(without["n_contacts"]).all()
    assert not np.allclose(with_n["delta_logit"], without["delta_logit"])
    # dropping a term must not introduce NaN of its own: what is undefined without the contact
    # count is exactly what was already undefined with it (rows missing a descriptor block)
    assert np.array_equal(np.isfinite(without["p_corrected"]), np.isfinite(with_n["p_corrected"]))


def test_a_table_whose_contact_count_is_not_the_potts_one_is_refused():
    """Regression, 2026-08-29. `tcren features -i placement,interface,topology` used to emit the
    footprint's CDR-loop tally under `n_contacts` -- 66 on 1ao7 where the Potts engaged-pair count
    is 29 -- and the frozen correction standardizes that column against the Potts population. The
    contact term read z = +7.2 for a native interface instead of +1.2 and the corrected probability
    moved, with no error, no warning and no NaN to notice. The footprint tally is `n_loop_contacts`
    now, and a table that still calls it `n_contacts` has to be refused rather than corrected.
    """
    ref, e, n = _corr_inputs()
    conf = np.full(len(e), 0.85)
    topo_only = {k: v for k, v in ref.items() if k not in rel._POTTS_MARKERS}
    assert "n_contacts" in topo_only, "the column the old topology pass wrote, under the old name"
    with pytest.raises(ValueError, match="potts"):
        rel.correct_confidence(topo_only, conf, contacts=topo_only["n_contacts"])
    # the same refusal through a polars frame, which is what `tcren diagnose` hands it
    frame = pl.DataFrame({k: np.asarray(v, float) for k, v in topo_only.items()})
    with pytest.raises(ValueError, match="n_loop_contacts"):
        rel.correct_confidence(frame, conf, contacts=frame["n_contacts"].to_numpy())
    # it is the provenance that is refused and not the missing energy: drop the column and the
    # same table corrects fine, with the contact term dropped rather than mis-standardized
    dropped = {k: v for k, v in topo_only.items() if k != "n_contacts"}
    out = rel.correct_confidence(dropped, conf)
    assert np.isfinite(out["p_corrected"]).any() and np.isnan(out["n_contacts"]).all()
    # and a table that does carry the potts columns is untouched by the guard
    assert np.isfinite(rel.correct_confidence(ref, conf, energy=e, contacts=n)["p_corrected"]).any()
