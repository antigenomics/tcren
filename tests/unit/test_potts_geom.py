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
