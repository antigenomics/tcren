"""The self-consistent loop — :mod:`tcren.potts.selfcon`.

One claim carries this module and it is testable exactly rather than approximately: at
``margin = 0`` the crystal's distance read pins every site, so the contact posterior IS the observed
map and the probability-weighted re-derivation must reproduce the plain derivation **to the last
bit**. That is the check that the two halves of the loop are wired to the same contacts. Everything
else here is downstream of it.

The second claim is that a positive margin actually moves mass off the crystal. A loop whose
posterior is the crystal map for every margin would converge beautifully and mean nothing.
"""
from __future__ import annotations

import numpy as np
import polars as pl
import pytest

pytest.importorskip("arda")

from pathlib import Path  # noqa: E402

from tcren.annotation.batch import iter_annotated_set  # noqa: E402
from tcren.mhc import annotate_mhc  # noqa: E402
from tcren.potential import derive_tcren, tcren2  # noqa: E402
from tcren.potts import available_pairs, fit_potts  # noqa: E402
from tcren.potts.model import PottsModel, centred_potential  # noqa: E402
from tcren.potts.selfcon import (LOGIT_CLIP, contact_posterior,  # noqa: E402
                                 fit_selfconsistent, observation_field, weighted_contacts)

PDB_DIR = Path(__file__).resolve().parents[1] / "assets" / "pdb"
pytestmark = pytest.mark.skipif(not PDB_DIR.exists(), reason="pdb assets missing")

CUTOFF = 5.0
MARGIN = 0.5


@pytest.fixture(scope="module")
def sites():
    """One pair table with the heavy-atom layer measured past the cutoff, from one pass."""
    frames = []
    for s in iter_annotated_set(str(PDB_DIR)):
        try:
            annotate_mhc(s)
            f = available_pairs(s, "peptide", probe=CUTOFF + 3 * MARGIN)
        except Exception:  # noqa: BLE001 - an unusable asset is skipped, as everywhere here
            continue
        if f.height:
            frames.append(f)
    if not frames:
        pytest.skip("no usable pair tables in the assets")
    return pl.concat(frames)


@pytest.fixture(scope="module")
def model():
    return PottsModel.bundled()


# -- the observation field ---------------------------------------------------------------------


def test_observation_field_is_monotone_and_saturates():
    f = pl.DataFrame({"d_heavy": [0.0, 3.0, 4.5, 5.0, 5.5, 8.0, None]})
    v = observation_field(f, cutoff=CUTOFF, margin=MARGIN)
    assert np.all(np.diff(v[:-1]) < 0), "the field must fall as the pair gets further apart"
    assert v[3] == pytest.approx(0.0), "a pair exactly at the cutoff is a coin toss"
    assert v[0] == pytest.approx((CUTOFF - 0.0) / MARGIN), "the field is the scaled distance"
    assert v[-1] == pytest.approx(-LOGIT_CLIP), "an unmeasured pair is a confident non-contact"
    # It saturates only where the scaled distance would run past the clip -- a sharp read of a
    # pair far outside the cutoff, which is where an unclipped logit would overflow the sigmoid.
    sharp = observation_field(pl.DataFrame({"d_heavy": [0.0, 50.0]}), cutoff=CUTOFF, margin=0.01)
    assert sharp[0] == pytest.approx(LOGIT_CLIP)
    assert sharp[1] == pytest.approx(-LOGIT_CLIP)


def test_zero_margin_is_the_hard_indicator(sites):
    """The limit, not an approximation: the field becomes +-clip on the observed indicator."""
    v = observation_field(sites, cutoff=CUTOFF, margin=0.0)
    assert set(np.unique(v)) <= {-LOGIT_CLIP, LOGIT_CLIP}
    assert np.array_equal(v > 0, sites["sigma"].to_numpy().astype(bool))


def test_a_wider_margin_releases_more_of_the_interface(sites):
    """A field that never leaves the crystal makes the whole loop a no-op, so assert it moves."""
    tight = observation_field(sites, cutoff=CUTOFF, margin=0.1)
    loose = observation_field(sites, cutoff=CUTOFF, margin=2.0)
    n_tight = int(np.sum(np.abs(tight) < 3.0))
    n_loose = int(np.sum(np.abs(loose) < 3.0))
    assert n_loose > n_tight > 0


# -- the posterior and the bridge to the derivation ---------------------------------------------


def test_zero_margin_posterior_is_the_observed_map(sites, model):
    post = contact_posterior(sites, model, margin=0.0, workers=1)
    assert np.array_equal(post["p_post"].to_numpy(), post["sigma"].to_numpy().astype(float))


def test_posterior_moves_off_the_crystal_at_a_positive_margin(sites, model):
    post = contact_posterior(sites, model, margin=MARGIN, chains=16, draws=20, seed=0, workers=1)
    p = post["p_post"].to_numpy()
    assert np.all((p >= 0.0) & (p <= 1.0))
    assert ((p > 0.05) & (p < 0.95)).sum() > 0, "no site is uncertain; the margin does nothing"


def test_weighted_contacts_emits_the_schema_the_derivation_reads(sites, model):
    post = contact_posterior(sites, model, margin=0.0, workers=1)
    cw = weighted_contacts(post)
    assert cw.columns == ["pdb.id", "residue.aa.from", "residue.aa.to", "p_post"]
    derive_tcren(cw, weight_col="p_post")   # raises if the schema is wrong


def test_the_reduction_gate_probability_weighting_reproduces_the_plain_derivation(sites, model):
    """At margin 0 the soft tally IS the hard tally, cell for cell."""
    post = contact_posterior(sites, model, margin=0.0, workers=1)
    soft = derive_tcren(weighted_contacts(post), weight_col="p_post")
    hard = derive_tcren(post.filter(pl.col("sigma") == 1).select(
        "pdb.id",
        pl.col("aa.rec").alias("residue.aa.from"),
        pl.col("aa.par").alias("residue.aa.to")))
    j = (soft.matrix.rename({"value": "soft"})
         .join(hard.matrix.rename({"value": "hard"}), on=["residue.aa.from", "residue.aa.to"]))
    assert j.height == soft.matrix.height == hard.matrix.height
    assert np.array_equal(j["soft"].to_numpy(), j["hard"].to_numpy())
    assert post["p_post"].sum() == pytest.approx(float(post["sigma"].sum()), abs=0)


# -- pinning the coupling to a matrix that was never shipped ------------------------------------


def test_centred_potential_takes_a_potential_object():
    assert np.allclose(centred_potential(tcren2()), centred_potential("tcren2"))
    assert np.allclose(centred_potential(tcren2(), centre=False),
                       centred_potential("tcren2", centre=False))


def test_a_pinned_model_round_trips_and_counts_one_parameter(sites, tmp_path):
    m = fit_potts(sites, coupling_matrix=tcren2(), pin_centred=False)
    assert m.coupling_pinned and m.coupling_matrix_name is None
    assert m.beta_matrix is not None
    free = fit_potts(sites)
    assert m.n_parameters() == free.n_parameters() - 399, "a pin is one scale, not 400 cells"
    p = tmp_path / "pinned.json"
    m.to_json(p)
    back = PottsModel.from_json(p)
    assert back.coupling_pinned
    assert np.allclose(back.coupling_array(), m.coupling_array())


def test_a_named_pin_still_counts_one_and_keeps_its_name(sites):
    m = fit_potts(sites, coupling_matrix="tcren2", pin_centred=False)
    assert m.coupling_matrix_name == "tcren2" and not m.coupling_pinned
    assert m.coupling is None
    assert m.n_parameters() == fit_potts(sites).n_parameters() - 399


# -- the loop ------------------------------------------------------------------------------------


def test_undamped_at_zero_margin_reaches_an_exact_fixed_point(sites):
    """The posterior cannot move, so the second iterate must equal the first, exactly."""
    pot, model, trace = fit_selfconsistent(sites, margin=0.0, damping=1.0, steps=4, workers=1)
    assert bool(trace["converged"][-1])
    assert trace.height == 2, "undamped, the loop is at its fixed point after one step"
    assert trace["max_abs_de"][-1] == 0.0
    hard = derive_tcren(sites.filter(pl.col("sigma") == 1).select(
        "pdb.id",
        pl.col("aa.rec").alias("residue.aa.from"),
        pl.col("aa.par").alias("residue.aa.to")))
    j = (pot.matrix.rename({"value": "loop"})
         .join(hard.matrix.rename({"value": "plain"}), on=["residue.aa.from", "residue.aa.to"]))
    assert np.allclose(j["loop"].to_numpy(), j["plain"].to_numpy())


def test_damping_makes_the_cell_movement_decay(sites):
    _, _, trace = fit_selfconsistent(sites, margin=0.0, damping=0.5, steps=5, tol=0.0, workers=1)
    de = trace["max_abs_de"].to_numpy()
    assert np.all(np.diff(de) < 0), "a damped fixed-point iteration must contract"
    assert de[1] / de[0] == pytest.approx(0.5, rel=1e-6), "the rate IS the damping"


def test_rms_movement_is_bounded_by_the_max_and_contracts_with_it(sites):
    """`rms_de` is the same movement averaged rather than maximised, and it is why it is there.

    `max_abs_de` is a maximum over 380 cells, so its Monte-Carlo floor is far higher than the
    r.m.s. of the same 380 numbers -- measured on Native2026, raising draws 100 -> 400 moves the
    max's plateau 0.00391 -> 0.00217 (ratio 0.553 against the 0.500 that quartering the variance
    predicts), which is a floor no draw budget drives to zero. The r.m.s. is what to read when the
    question is whether the matrix is still moving. Both must contract under damping.
    """
    _, _, trace = fit_selfconsistent(sites, margin=0.0, damping=0.5, steps=4, tol=0.0, workers=1)
    mx, rms = trace["max_abs_de"].to_numpy(), trace["rms_de"].to_numpy()
    assert np.all(rms > 0) and np.all(rms <= mx + 1e-12), "an r.m.s. cannot exceed the max"
    assert np.all(np.diff(rms) < 0), "the r.m.s. must contract with the iteration"
    assert rms[1] / rms[0] == pytest.approx(0.5, rel=1e-6), "damping sets the rate for both"


def test_a_positive_margin_moves_the_matrix_and_records_the_uncertain_mass(sites):
    pot, model, trace = fit_selfconsistent(sites, margin=MARGIN, steps=2, chains=16, draws=20,
                                           workers=1, notes="assets")
    assert trace["mass_from_uncertain"][0] > 0.0
    assert trace["contact_mass"][0] > 0.0
    assert model.coupling_pinned and np.isfinite(model.beta_matrix)
    assert "self-consistent" in model.notes
    assert pot.name == "TCRen2-EM" and pot.matrix.height == tcren2().matrix.height


def test_steps_must_be_at_least_one(sites):
    with pytest.raises(ValueError, match="steps must be at least 1"):
        fit_selfconsistent(sites, steps=0)
