"""The all-vs-all Cα table: every TCR read against every pMHC in the canonical frame.

The one test that matters is the diagonal --- a pseudo-complex whose TCR and pMHC come from the
*same* crystal must reproduce that crystal's own Cα distances. It validates the whole vectorised
engine at once, because nothing in the broadcast treats the diagonal specially.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import polars as pl
import pytest

from tcren.potts.cross import (CROSS_COLUMNS, MIN_CA, PAR_COLUMNS, REC_COLUMNS, _BUCKET,
                               cross_pairs, partner_residues, receptor_residues, residue_tables,
                               write_cross_pairs)
from tcren.potts.geom import region_peptide_dmap

PDB_DIR = Path(__file__).resolve().parents[1] / "assets" / "pdb"
pytestmark = pytest.mark.skipif(not PDB_DIR.exists(), reason="pdb assets missing")


@pytest.fixture(scope="module")
def structures():
    pytest.importorskip("arda")
    from tcren.annotation.batch import iter_annotated_set
    from tcren.mhc.regions import annotate_mhc_batch

    loaded = list(iter_annotated_set(PDB_DIR))
    annotate_mhc_batch(loaded)          # without this every structure reports MHCI
    return loaded


@pytest.fixture(scope="module")
def tables(structures):
    return residue_tables(structures)


@pytest.fixture(scope="module")
def cross(tables):
    return cross_pairs(*tables)


def test_the_native_diagonal_reproduces_the_crystals_own_distances(structures, cross):
    """The gate. If this passes, the broadcast is right for every off-diagonal pairing too."""
    by_id = {s.pdb_id: s for s in structures}
    worst, checked = 0.0, 0
    native = cross.filter(pl.col("is_native"))
    for (pdb_id, chain, region), group in native.group_by(
            ["pdb.id.rec", "chain.rec", "region.rec"]):
        observed = region_peptide_dmap(by_id[pdb_id], chain, region)
        if observed is None:
            continue
        got = group.sort("pos.rec", "pos.par")["d_ca"].to_numpy().reshape(observed.shape)
        worst = max(worst, float(np.abs(got - observed).max()))
        checked += got.size
    assert checked > 3000, f"only {checked} native cells checked"
    # float32 on a ~30 Å distance is ~3.6e-6 Å; anything larger is a real disagreement.
    assert worst < 1e-5, f"native pairing differs from the crystal by {worst:.2e} Å"


def test_every_pairing_is_present(tables):
    """Row count is the product of the two tables, bucket by bucket, at an unbounded radius."""
    rec, par = tables
    expected = 0
    for (cls, _chain, _region, _loop_len), r in rec.group_by(
            ["mhc_class", "chain.rec", "region.rec", "len.rec"]):
        expected += len(r) * len(par.filter(pl.col("mhc_class") == cls))
    assert len(cross_pairs(rec, par, min_dist=0.0)) == expected


def test_the_clash_filter_drops_no_row_a_real_complex_produced(tables):
    """3.35 Å sits under the 3.366 Å closest approach any solved complex makes."""
    unfiltered = cross_pairs(*tables, min_dist=0.0)
    filtered = cross_pairs(*tables)
    assert filtered["d_ca"].min() >= MIN_CA
    native_before = unfiltered.filter(pl.col("is_native"))
    native_after = filtered.filter(pl.col("is_native"))
    assert len(native_after) == len(native_before), (
        f"the clash filter lost {len(native_before) - len(native_after)} native rows")
    # It is a per-row filter: no pseudo-complex is removed wholesale.
    def pairings(f):
        return set(f.select("pdb.id.rec", "pdb.id.par").unique().iter_rows())
    assert pairings(filtered) == pairings(unfiltered)


def test_an_mhc_class_is_never_crossed(cross, tables):
    """A class-I TCR on a class-II groove is not a pose, and the bucketing must forbid it."""
    rec, par = tables
    rec_class = dict(rec.select("pdb.id", "mhc_class").unique().iter_rows())
    par_class = dict(par.select("pdb.id", "mhc_class").unique().iter_rows())
    mixed = cross.filter(
        pl.col("pdb.id.rec").replace_strict(rec_class)
        != pl.col("pdb.id.par").replace_strict(par_class))
    assert mixed.is_empty(), f"{len(mixed)} rows cross an MHC class"


def test_the_declared_sort_order_holds(cross):
    """The fit reshapes a bucket straight out of the rows, so the order is part of the contract."""
    assert cross.columns == list(CROSS_COLUMNS)
    keys = [*_BUCKET, "pdb.id.rec", "pdb.id.par", "pos.rec", "pos.par"]
    assert cross.select(keys).equals(cross.select(keys).sort(keys))


def test_a_radius_selects_a_strict_subset(tables):
    bounded = cross_pairs(*tables, radius=15.0)
    assert bounded["d_ca"].max() <= 15.0
    assert 0 < len(bounded) < len(cross_pairs(*tables))


def test_dropping_the_diagonal_drops_exactly_the_native_rows(tables, cross):
    without = cross_pairs(*tables, include_self=False)
    assert not without["is_native"].any()
    assert len(without) == len(cross) - int(cross["is_native"].sum())


def test_the_residue_tables_carry_their_declared_schema(tables):
    rec, par = tables
    assert rec.columns == list(REC_COLUMNS)
    assert par.columns == list(PAR_COLUMNS)
    # len.rec is the resolved region length: every row of a loop agrees on it, and it equals the
    # number of rows that loop contributed.
    counts = rec.group_by("pdb.id", "chain.rec", "region.rec").agg(
        pl.len().alias("rows"), pl.col("len.rec").n_unique().alias("distinct"),
        pl.col("len.rec").first().alias("declared"))
    assert (counts["distinct"] == 1).all()
    assert (counts["rows"] == counts["declared"]).all()
    assert set(par["role.par"].unique()) <= {"anchor", "tcr_facing"}


def test_glycine_is_flagged_rather_than_dropped(tables):
    """Missing a Cβ is a missing covariate, not a reason to lose a backbone position."""
    rec, _ = tables
    glycine = rec.filter(pl.col("aa.rec") == "G")
    assert len(glycine) > 0, "no glycine in the fixture loops"
    assert not glycine["has_cb"].any()
    # The placeholder is the Cα itself, so the row still carries a position.
    assert np.allclose(glycine.select("cb_x", "cb_y", "cb_z").to_numpy(),
                       glycine.select("ca_x", "ca_y", "ca_z").to_numpy())


def test_cos_cb_points_at_the_peptide_or_away_from_it(cross):
    """+1 is a side chain aimed at the partner, -1 opposed, 0 where there is no side chain."""
    with_cb = cross.filter(pl.col("has_cb"))
    assert with_cb["cos_cb"].min() >= -1.0 - 1e-5
    assert with_cb["cos_cb"].max() <= 1.0 + 1e-5
    assert (cross.filter(~pl.col("has_cb"))["cos_cb"].abs() < 1e-9).all()
    # A residue in contact is far more often turned toward the peptide than away from it.
    near = with_cb.filter(pl.col("d_ca") < 8.0)["cos_cb"].mean()
    far = with_cb.filter(pl.col("d_ca") > 20.0)["cos_cb"].mean()
    assert near > far, f"cos_cb {near:.3f} near vs {far:.3f} far — the sign convention is inverted"


def test_a_structure_missing_a_peptide_is_dropped_from_both_tables(structures):
    from tcren.structure.model import PEPTIDE_TYPE

    stripped = structures[0].__class__(
        structures[0].pdb_id,
        [c for c in structures[0].chains if c.chain_type != PEPTIDE_TYPE],
        structures[0].complex_species, structures[0].cell_type)
    assert partner_residues(stripped).is_empty()
    assert not receptor_residues(stripped).is_empty()      # the loops are still there
    with pytest.raises(ValueError, match="CDR loop and a peptide"):
        residue_tables([stripped])


def test_the_parquet_round_trips(tables, tmp_path):
    rows, parts = write_cross_pairs(*tables, tmp_path / "cross", radius=15.0, part_rows=5_000)
    assert parts > 1, "part_rows did not split the output"
    back = pl.read_parquet(tmp_path / "cross" / "*.parquet")
    assert len(back) == rows
    assert back.columns == list(CROSS_COLUMNS)
    assert back.equals(cross_pairs(*tables, radius=15.0))


def test_a_rewrite_replaces_the_old_parts(tables, tmp_path):
    """A smaller second run must not leave the first run's parts behind to be read as data."""
    out = tmp_path / "cross"
    write_cross_pairs(*tables, out, part_rows=5_000)
    _, parts = write_cross_pairs(*tables, out, radius=8.0, part_rows=5_000)
    assert len(list(out.glob("part-*.parquet"))) == parts


def test_an_impossible_radius_is_refused_rather_than_written_empty(tables):
    with pytest.raises(ValueError, match="no receptor:peptide pair"):
        cross_pairs(*tables, radius=1.0, min_dist=0.0)
