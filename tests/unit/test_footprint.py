"""Footprint shape — coverage entropy and topology.

The invariants worth pinning are the ones that make the numbers comparable across structures:
the partitions must have the size they claim, the diversity measures must be maximal on a uniform
footprint and fall when it concentrates, the topology must count patches and holes correctly on
shapes whose answer is known by hand, and every feature must be invariant under rigid motion --
that last one is what lets the module skip canonical orientation.
"""
from __future__ import annotations

import numpy as np
import polars as pl
import pytest

from tcren.footprint import (
    CELL_LOOPS,
    FOOTPRINT_FEATURES,
    _diversity,
    _flag_betti,
    _gf2_rank,
    _h0_persistence_entropy,
    _selfcheck,
    cell_counts,
    footprint_batch,
    footprint_features,
    footprint_score,
)

from tcren.structure.model import PEPTIDE_TYPE, Atom, Chain, RegionMarkup, Residue, Structure

_THREE = {"L": "LEU", "D": "ASP", "K": "LYS", "W": "TRP"}


def _chain(cid, ctype, n, origin, aas, loops, rng):
    """A chain of `n` residues scattered about `origin`, optionally marked up with CDR1/2/3."""
    res = []
    for i in range(n):
        x, y, z = np.asarray(origin, float) + rng.normal(scale=1.8, size=3)
        aa = aas[i % len(aas)]
        res.append(Residue(i, i + 1, "", aa, _THREE[aa],
                           (Atom("CA", "C", np.array([x, y, z])),
                            Atom("CB", "C", np.array([x + 0.7, y, z])))))
    c = Chain(cid, res, chain_type=ctype)
    if loops:
        span = n // 3
        c.regions = [RegionMarkup(r, i * span, (i + 1) * span - 1,
                                  "".join(x.aa for x in res[i * span:(i + 1) * span]),
                                  res[i * span:(i + 1) * span])
                     for i, r in enumerate(("CDR1", "CDR2", "CDR3"))]
    return c


def _full_complex():
    """Peptide at z=0, MHC below it, both receptor chains above with all three CDRs marked up."""
    rng = np.random.default_rng(0)
    return Structure("synth_full", [
        _chain("C", PEPTIDE_TYPE, 9, [0.0, 0.0, 0.0], "DKW", False, rng),
        _chain("D", "TRA", 18, [-2.0, 0.0, 4.0], "LDK", True, rng),
        _chain("B", "TRB", 18, [2.0, 0.0, 4.0], "LDW", True, rng),
        _chain("A", "MHCa", 24, [0.0, 0.0, -4.0], "LDW", False, rng)])


def test_selfcheck_passes():
    _selfcheck()


# --- the partition -------------------------------------------------------------------------

def test_cell_counts_partitions_into_loops_targets_and_peptide_bands():
    t = cell_counts(_full_complex())
    assert set(t["loop"].to_list()) <= set(CELL_LOOPS)
    assert set(t["target"].to_list()) == {"pep", "mhc"}
    # the peptide side must actually reach all three bands, or the 24-cell partition is a 12-cell
    # one wearing a hat -- which is exactly what a null `pos.to` used to produce
    assert set(t.filter(pl.col("target") == "pep")["band"].to_list()) == {"pN", "pM", "pC"}
    assert t.filter(pl.col("target") == "mhc")["band"].unique().to_list() == ["mhc"]
    assert t["n"].sum() > 0


def test_peptide_banding_makes_the_finer_partition_finer():
    r = footprint_features(_full_complex())
    # D2 is the effective number of engaged cells, so splitting cells can only raise it
    assert r["D2_pep24"] > r["D2_cell"] > 0
    assert r["D2_cell"] >= r["D2_loop"]


def test_every_declared_feature_is_present_and_finite():
    r = footprint_features(_full_complex())
    missing = [k for k in FOOTPRINT_FEATURES if k not in r]
    assert not missing, missing
    assert all(np.isfinite(r[k]) for k in FOOTPRINT_FEATURES), \
        {k: r[k] for k in FOOTPRINT_FEATURES if not np.isfinite(r[k])}


def test_contact_totals_add_up():
    r = footprint_features(_full_complex())
    assert r["n_contacts"] == r["n_pep_contacts"] + r["n_mhc_contacts"]


def test_unannotated_mhc_is_flagged_not_silently_empty():
    """The MHC pass must run after chain typing, or six of the twelve cells are unreachable.

    `classify_chains` types an MHC chain generically as "MHC"; `interface("tcr_mhc")` matches the
    supertype `annotate_mhc` assigns. Skipping it emptied the MHC half of the partition with no
    error -- p_germ_mhc fell from ~0.78 to ~0.06 and H_cell was normalised by ln 12 over a
    partition half of which could never be occupied.
    """
    s = _full_complex()
    for c in s.chains:
        if c.chain_type == "MHCa":
            c.chain_type = "MHC"                       # what classify_chains leaves behind
    with pytest.warns(RuntimeWarning, match="MHC chains are not annotated"):
        t = cell_counts(s)
    assert "mhc" not in set(t["target"].to_list())     # and this is exactly why it warns


def test_annotated_mhc_reaches_the_mhc_cells():
    t = cell_counts(_full_complex())
    assert t.filter(pl.col("target") == "mhc")["n"].sum() > 0


# --- the diversity measures ----------------------------------------------------------------

def test_uniform_composition_is_exactly_one_and_maximal():
    d = _diversity(np.full(12, 5.0), 12)
    assert d["H"] == pytest.approx(1.0)
    assert d["D1"] == pytest.approx(12.0)
    assert d["D2"] == pytest.approx(12.0)
    assert d["J"] == pytest.approx(1.0)


def test_concentrating_the_footprint_lowers_every_diversity():
    even = _diversity(np.full(12, 5.0), 12)
    skew = _diversity(np.array([49.0] + [1.0] * 11), 12)
    for k in ("H", "D1", "D2", "J"):
        assert skew[k] < even[k], (k, skew[k], even[k])


def test_hill_order_two_discounts_rare_cells_more_than_order_one():
    # one dominant cell plus a long tail of singletons: D2 must fall further below D1
    n = np.array([60.0] + [1.0] * 11)
    d = _diversity(n, 12)
    assert d["D2"] < d["D1"]


def test_diversity_of_an_empty_footprint_is_null_not_zero():
    d = _diversity(np.zeros(12), 12)
    assert all(np.isnan(v) for v in d.values()), d


# --- topology ------------------------------------------------------------------------------

def test_flag_complex_counts_a_hole_only_while_the_diagonal_is_unspanned():
    square = np.array([[0, 0, 0], [2, 0, 0], [2, 2, 0], [0, 2, 0]], float)
    assert _flag_betti(square, 2.5) == (1.0, 1.0)   # sides joined, diagonal not: one hole
    assert _flag_betti(square, 3.0) == (1.0, 0.0)   # diagonal joined: triangles fill the hole
    assert _flag_betti(square, 1.0) == (4.0, 0.0)   # nothing joined: four patches


def test_a_fragmented_footprint_has_more_patches():
    one = np.array([[0, 0, 0], [1, 0, 0], [2, 0, 0], [3, 0, 0]], float)
    two = np.array([[0, 0, 0], [1, 0, 0], [50, 0, 0], [51, 0, 0]], float)
    assert _flag_betti(one, 1.5)[0] == 1.0
    assert _flag_betti(two, 1.5)[0] == 2.0


def test_gf2_rank_sees_a_dependency_that_only_holds_mod_two():
    # the three rows sum to zero over GF(2) but are independent over the rationals
    assert _gf2_rank(np.array([[1, 1, 0], [0, 1, 1], [1, 0, 1]], np.uint8)) == 2
    assert _gf2_rank(np.eye(4, dtype=np.uint8)) == 4


def test_h0_persistence_entropy_is_maximal_for_evenly_spaced_points():
    even = np.array([[i, 0, 0] for i in range(6)], float)
    clumped = np.array([[0, 0, 0], [0.1, 0, 0], [0.2, 0, 0], [9, 0, 0], [9.1, 0, 0], [20, 0, 0]],
                       float)
    assert _h0_persistence_entropy(even) == pytest.approx(1.0)
    assert _h0_persistence_entropy(clumped) < _h0_persistence_entropy(even)


# --- the property that lets us skip canonical orientation -----------------------------------

def test_every_feature_is_invariant_under_rigid_motion():
    """No canonical frame is needed: rotate and translate the whole complex, get the same row."""
    s = _full_complex()
    theta = 0.7
    R = np.array([[np.cos(theta), -np.sin(theta), 0.0],
                  [np.sin(theta), np.cos(theta), 0.0], [0.0, 0.0, 1.0]])
    shift = np.array([13.0, -5.0, 2.5])
    moved = Structure(s.pdb_id, [
        Chain(c.chain_id,
              [Residue(r.seq_index, r.pdb_index, r.insertion_code, r.aa, r.resname,
                       tuple(Atom(name=a.name, element=a.element,
                                  coord=(R @ np.asarray(a.coord, float)) + shift)
                             for a in r.atoms))
               for r in c.residues],
              chain_type=c.chain_type) for c in s.chains])
    for c_new, c_old in zip(moved.chains, s.chains):
        if c_old.regions:
            by = {r.seq_index: r for r in c_new.residues}
            c_new.regions = [type(reg)(region_type=reg.region_type,
                                       start_seq_index=reg.start_seq_index,
                                       end_seq_index=reg.end_seq_index, sequence=reg.sequence,
                                       residues=[by[r.seq_index] for r in reg.residues])
                             for reg in c_old.regions]

    a, b = footprint_features(s), footprint_features(moved)
    assert set(a) == set(b)
    for k in a:
        assert a[k] == pytest.approx(b[k], rel=1e-9, abs=1e-9), (k, a[k], b[k])


# --- batch + cohort score --------------------------------------------------------------------

def test_batch_over_structures_yields_one_row_each():
    s = _full_complex()
    t = footprint_batch([s, s])
    assert t.height == 2
    assert "pdb.id" in t.columns
    assert all(k in t.columns for k in FOOTPRINT_FEATURES)


def test_batch_of_nothing_is_an_empty_frame_not_an_error():
    assert footprint_batch([]).height == 0


def test_footprint_score_is_the_sum_of_two_standardised_channels():
    t = pl.DataFrame({"pdb.id": list("abcd"),
                      "D2_pep24": [4.0, 8.0, 12.0, 16.0],
                      "fp_b0_frac_r7": [0.4, 0.3, 0.2, 0.1]})
    out = footprint_score(t)
    assert out["fp_score"].to_list() == pytest.approx(
        (out["z_coverage"] + out["z_patch"]).to_list())
    # both channels point the same way: more effective cells and fewer patches score higher
    assert out["fp_score"].to_list() == sorted(out["fp_score"].to_list())
    assert out["z_coverage"].mean() == pytest.approx(0.0, abs=1e-12)


def test_one_missing_structure_does_not_flatten_a_whole_channel():
    """polars propagates NaN through mean/std, so a single contact-free row used to null the column,
    `fill_nan(0.0)` flattened it to zero, and the score silently became its other channel alone."""
    t = pl.DataFrame({"pdb.id": list("abcde"),
                      "D2_pep24": [4.0, 8.0, 12.0, 16.0, float("nan")],
                      "fp_b0_frac_r7": [0.4, 0.3, 0.2, 0.1, 0.5]})
    out = footprint_score(t)
    good = out["z_coverage"].to_numpy()[:4]
    assert not np.allclose(good, 0.0), good          # the four measured rows still carry the channel
    assert out["z_coverage"].to_numpy()[4] == 0.0    # ...and the missing one contributes nothing
    assert np.all(np.diff(good) > 0)


def test_footprint_score_standardises_within_group():
    t = pl.DataFrame({"pdb.id": list("abcd"), "epitope": ["X", "X", "Y", "Y"],
                      "D2_pep24": [4.0, 8.0, 104.0, 108.0],
                      "fp_b0_frac_r7": [0.4, 0.3, 0.4, 0.3]})
    out = footprint_score(t, group="epitope")
    # the two cohorts sit at wildly different D2 scales; within-cohort z removes the offset
    assert out["z_coverage"].to_list()[:2] == pytest.approx(out["z_coverage"].to_list()[2:])
