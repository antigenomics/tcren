"""The length-generalising Ca-distance prior — :mod:`tcren.potts.geom`.

The prior exists so a receptor:peptide pair can be placed without a template of the right shape.
Two things therefore have to hold: it must emit a table every other `tcren.potts` function accepts
unchanged, and reading a surface back at a length it was never fitted on must work rather than
silently return the training lengths.
"""
from __future__ import annotations

import numpy as np
import polars as pl
import pytest

pytest.importorskip("arda")

from pathlib import Path  # noqa: E402

from tcren.annotation.batch import iter_annotated_set  # noqa: E402
from tcren.geometry import reach_max  # noqa: E402
from tcren.mhc import annotate_mhc  # noqa: E402
from tcren.potts import PottsModel, available_pairs, eta, site_codes  # noqa: E402
from tcren.potts.geom import (REGIONS_MODELLED, GeomPrior, feasible,  # noqa: E402
                              fit_geom_prior, loops_from_structure, predicted_distances,
                              predicted_sites, region_peptide_dmap)
from tcren.potts.sites import SITE_COLUMNS  # noqa: E402

PDB_DIR = Path(__file__).resolve().parents[1] / "assets" / "pdb"
pytestmark = pytest.mark.skipif(not PDB_DIR.exists(), reason="pdb assets missing")


@pytest.fixture(scope="module")
def structures():
    out = []
    for s in iter_annotated_set(str(PDB_DIR)):
        try:
            annotate_mhc(s)
        except Exception:
            pass
        out.append(s)
    return out


@pytest.fixture(scope="module")
def prior(structures):
    return fit_geom_prior(structures, grid=12, notes="assets")


def test_all_three_loops_are_modelled(prior):
    """CDR3 alone is 69.3% of TCR:peptide contacts; a CDR3-only prior caps its own recall."""
    got = {k.split("|")[0].split(":")[1] for k in prior.keys()}
    assert got == set(REGIONS_MODELLED)


def test_surfaces_are_read_back_at_the_lengths_asked_for(prior):
    """The graphon coordinate is the whole length-generalisation mechanism."""
    for L, M in ((8, 9), (17, 15), (5, 11)):
        m, sd = predicted_distances(prior, "TRB", "CDR3", loop_len=L, pep_len=M)
        assert m.shape == (L, M) and sd.shape == (L, M)
        assert np.all(np.isfinite(m)) and np.all(sd >= 0)


def test_a_length_never_seen_still_returns_a_surface(prior, structures):
    """The point of the prior: a CDR3 of 25 has no template here and must still be placeable."""
    seen = set()
    for s in structures:
        for k, v in loops_from_structure(s).items():
            if k == "TRB:CDR3":
                seen.add(len(v))
    assert 25 not in seen
    m, _ = predicted_distances(prior, "TRB", "CDR3", loop_len=25, pep_len=9)
    assert m.shape == (25, 9) and np.all(np.isfinite(m))


def test_predicted_sites_are_accepted_by_the_potts_model(prior, structures):
    """The design constraint: the emitted schema is the one `available_pairs` emits."""
    s = next(x for x in structures if x.pdb_id == "1ao7")
    loops = loops_from_structure(s)
    pep = next(c.sequence() for c in s.chains if c.chain_type == "PEPTIDE")
    pred = predicted_sites(prior, loops, pep)
    assert pred.height > 0
    assert set(SITE_COLUMNS).issubset(pred.columns)
    m = PottsModel.bundled()
    codes, _, _ = site_codes(pred.drop("d_ca_sd", "feasible"), m)
    e = eta(codes, m)
    assert len(e) == pred.height and np.all(np.isfinite(e))


def test_sigma_is_a_placeholder_and_says_so(prior, structures):
    """Nothing here observed a contact; a caller reading `sigma` as data is reading a zero."""
    s = next(x for x in structures if x.pdb_id == "1ao7")
    pep = next(c.sequence() for c in s.chains if c.chain_type == "PEPTIDE")
    pred = predicted_sites(prior, loops_from_structure(s), pep)
    assert pred["sigma"].sum() == 0
    assert pred["d_heavy"].null_count() == pred.height


def test_the_radius_bounds_what_is_emitted(prior, structures):
    s = next(x for x in structures if x.pdb_id == "1ao7")
    pep = next(c.sequence() for c in s.chains if c.chain_type == "PEPTIDE")
    wide = predicted_sites(prior, loops_from_structure(s), pep, radius=15.0)
    tight = predicted_sites(prior, loops_from_structure(s), pep, radius=8.0)
    assert tight.height < wide.height
    assert tight["d_ca"].max() <= 8.0


def test_predictions_track_the_real_distances(prior, structures):
    """In-sample, but it must at least be correlated or nothing downstream can work."""
    s = next(x for x in structures if x.pdb_id == "1ao7")
    pep = next(c.sequence() for c in s.chains if c.chain_type == "PEPTIDE")
    pred = predicted_sites(prior, loops_from_structure(s), pep, radius=15.0)
    obs = available_pairs(s, "peptide").filter(pl.col("region.rec").is_in(list(REGIONS_MODELLED)))
    key = ["chain.rec", "region.rec", "pos.rec", "pos.par"]
    j = obs.join(pred.select(*key, pl.col("d_ca").alias("d_pred")), on=key, how="inner")
    assert j.height > 50
    r = float(np.corrcoef(j["d_ca"].to_numpy(), j["d_pred"].to_numpy())[0, 1])
    assert r > 0.6


def test_feasibility_uses_the_closed_form_reach(prior):
    """A distance no loop of that length can span is flagged, and the bound is `reach_max`."""
    neck = prior.neck_mean["TRB:CDR3"]
    assert feasible(prior, "TRB", "CDR3", 14, 5.0)
    far = reach_max(3, neck, prior.bond) + neck / 2.0 + 50.0
    assert not feasible(prior, "TRB", "CDR3", 3, far)


def test_region_dmap_reaches_every_modelled_region(structures):
    s = next(x for x in structures if x.pdb_id == "1ao7")
    for region in REGIONS_MODELLED:
        d = region_peptide_dmap(s, "TRB", region)
        assert d is not None and d.ndim == 2 and min(d.shape) >= 2


def test_round_trips_through_json(prior, tmp_path):
    p = tmp_path / "geom.json"
    prior.to_json(p)
    back = GeomPrior.from_json(p)
    assert back.keys() == prior.keys()
    assert back.grid == prior.grid
    np.testing.assert_allclose(np.asarray(back.mean[prior.keys()[0]]),
                               np.asarray(prior.mean[prior.keys()[0]]))


def test_an_unseen_mhc_class_falls_back_rather_than_failing(prior):
    """Class II lives only in the crystals; a class the fit never saw must not be a hard error."""
    m, _ = prior.surface("TRB", "CDR3", "MHC-something-else")
    assert m.shape == (prior.grid, prior.grid)


def test_fitting_on_nothing_is_an_error_not_an_empty_prior():
    with pytest.raises(ValueError, match="no structure yielded"):
        fit_geom_prior([])


# --- the shipped prior -----------------------------------------------------------------------
#
# `fit_geom_prior` above is fitted on the 17 bundled crystals, which is enough to test the machinery
# and far too few to ship. These check the artifact in the wheel: the object that lets a contact map
# be scored with no structure of the complex at all.

@pytest.fixture(scope="module")
def shipped():
    return GeomPrior.bundled()


def test_the_prior_ships_and_loads_by_name(shipped):
    from tcren.potts import BUNDLED_PRIOR
    assert BUNDLED_PRIOR == "geom_tcr_peptide"
    assert shipped.grid == 12
    # The same 359 alpha-beta Native2026 crystals TCRen2 is derived from: 12 not alpha-beta and
    # 3 withdrawn, of 374. A different count means the prior and the potential have drifted apart.
    assert shipped.n_structures == 359


def test_both_MHC_classes_carry_a_surface(shipped):
    """Without the groove annotation every crystal reports class I and the split collapses."""
    assert len(shipped.keys()) == 12
    for chain in ("TRA", "TRB"):
        for region in REGIONS_MODELLED:
            for cls in ("MHCI", "MHCII"):
                assert f"{chain}:{region}|{cls}" in shipped.mean
    assert shipped.n["TRB:CDR3|MHCI"] == 265
    assert shipped.n["TRB:CDR3|MHCII"] == 94


def test_a_prediction_scores_under_the_shipped_potts_model(shipped, structures):
    """Sequence in, contact probability out — no coordinate of the complex is read."""
    from tcren.potts import contact_probabilities
    s = next(x for x in structures if x.pdb_id == "1ao7")
    pep = next(c.sequence() for c in s.chains if c.chain_type == "PEPTIDE")
    sites = predicted_sites(shipped, loops_from_structure(s), pep)
    assert sites.height > 50
    probs = contact_probabilities(sites.drop("d_ca_sd", "feasible"),
                                  PottsModel.bundled("potts_tcr_peptide_sc"), workers=1)
    assert probs["p_model"].is_between(0.0, 1.0).all()


def test_a_non_standard_residue_is_dropped_not_crashed_on(shipped):
    """`X` has no code in the model's alphabet, and two Native2026 crystals carry one."""
    sites = predicted_sites(shipped, {"TRB:CDR3": "ASSXGELF"}, "GILXFVFTL")
    assert sites.height
    assert "X" not in set(sites["aa.rec"]) and "X" not in set(sites["aa.par"])


@pytest.fixture(scope="module")
def ens():
    """The shipped ENSEMBLE prior: one surface per contributing crystal, not their mean."""
    from tcren.potts import BUNDLED_ENSEMBLE, GeomPrior
    return GeomPrior.bundled(BUNDLED_ENSEMBLE)


def test_the_ensemble_prior_ships_and_keeps_one_surface_per_crystal(ens):
    assert ens.surfaces is not None
    assert ens.grid == 32
    for chain in ("TRA", "TRB"):
        for region in REGIONS_MODELLED:
            for cls in ("MHCI", "MHCII"):
                k = f"{chain}:{region}|{cls}"
                assert ens.surfaces[k].shape == (ens.n[k], 32, 32)
    assert ens.n["TRB:CDR3|MHCI"] == 265
    assert ens.n["TRB:CDR3|MHCII"] == 94


def test_the_ensemble_readback_is_one_geometry_per_member(ens):
    """Each slice is one crystal carried to the target lengths, not a per-cell average."""
    from tcren.potts import predicted_ensemble
    d = predicted_ensemble(ens, "TRB", "CDR3", loop_len=13, pep_len=9, mhc_class="MHCI")
    assert d.shape == (ens.n["TRB:CDR3|MHCI"], 13, 9)
    # the mean of the members is the mean surface; no member is that mean
    assert not np.allclose(d.min(axis=0), d.mean(axis=0))


def test_averaging_over_geometries_beats_averaging_the_surfaces(ens, shipped, structures):
    """The whole point: the mean surface loses the contact count, the ensemble keeps it.

    Measured over 364 held-out crystals the ensemble expects 24.62 contacts against 21.02
    observed and the mean surface 9.80; this asserts the direction and the rough scale on one
    structure, which is what a unit test can carry.
    """
    from tcren.potts import contact_probabilities, predict_contact_map
    m = PottsModel.bundled("potts_tcr_peptide_sc")
    s = next(x for x in structures if x.pdb_id == "1ao7")
    pep = next(c.sequence() for c in s.chains if c.chain_type == "PEPTIDE")
    loops = loops_from_structure(s)

    flat = predicted_sites(shipped, loops, pep)
    flat = contact_probabilities(flat.drop("d_ca_sd", "feasible"), m, workers=1)
    ens_map = predict_contact_map(ens, loops, pep, m, members=24, workers=1)

    assert ens_map.height
    assert ens_map["p_model"].is_between(0.0, 1.0).all()
    assert ens_map["p_avail"].is_between(0.0, 1.0).all()
    assert ens_map["p_model"].sum() > 2.0 * flat["p_model"].sum()


def test_the_npz_round_trip_keeps_the_stack(tmp_path, ens):
    from tcren.potts import GeomPrior
    p = tmp_path / "g.npz"
    ens.to_npz(p)
    back = GeomPrior.from_npz(p)
    assert sorted(back.surfaces) == sorted(ens.surfaces)
    # stored as hundredths of an angstrom, so equality holds to the model's own resolution
    for k, v in ens.surfaces.items():
        assert np.allclose(back.surfaces[k], v, atol=0.01)
