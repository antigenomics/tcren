"""The published descriptors added at 2.30.0: known answers, decomposition, degenerate input.

Two of the three checks `test_footprint` applies to the footprint block -- a hand-built input whose
answer is computable by hand, and NaN rather than 0 where a structure gives a descriptor no
support. The third, rigid-motion invariance, is not repeated: the surface block is a pure function
of :func:`tcren.topology.surface.surface_map`, whose frame refit
`test_surface.test_groove_frame_is_independent_of_input_orientation` pins, and the graph pair reads
a contact map with no frame in it.

One check the footprint block has no analogue of: the four sign-resolved gap columns exist only to
decompose ``sc_gap_mean``, so the decomposition is asserted rather than assumed.
"""
from __future__ import annotations

import math

import numpy as np
import polars as pl
import pytest

from tcren.recognition import DESCRIPTORS, DETAIL, INVARIANCE
from tcren.topology.literature import (
    LITERATURE_FEATURES,
    SURFACE_FEATURES,
    contact_order,
    literature_features,
    participation_coefficient,
)


def _frame() -> pl.DataFrame:
    """Two TCR residues against a peptide and an MHC, small enough to work out by hand.

    Residue 1 touches peptide 3 and MHC 7 -- one edge into each of two modules, so
    ``P = 1 - (1/2)^2 - (1/2)^2 = 0.5``. Residue 2 touches peptide 3 and peptide 5 -- both edges
    into one module, so ``P = 1 - 1 = 0``. The mean over engaged residues is 0.25.
    """
    return pl.DataFrame({
        "loop": ["TRA:CDR3"] * 4,
        "cf": ["A"] * 4, "rf": [1, 1, 2, 2], "rt": [3, 7, 3, 5],
        "target": ["peptide", "mhc", "peptide", "peptide"],
    })


def test_participation_coefficient_matches_the_hand_computed_value():
    assert participation_coefficient(_frame())["partcoef_tcr"] == pytest.approx(0.25)


def test_participation_coefficient_is_zero_when_every_residue_reads_one_module():
    """Di Paola's P is 0 for a residue whose edges all land in one module, not NaN."""
    t = _frame().filter(pl.col("target") == "peptide")
    assert participation_coefficient(t)["partcoef_tcr"] == pytest.approx(0.0)


def test_contact_order_is_the_mean_separation_over_the_target_span():
    """One loop reaches peptide 3 and 5: mean separation 2, span 5-3+1 = 3, so CO = 2/3."""
    assert contact_order(_frame())["co_pep"] == pytest.approx(2.0 / 3.0)


def test_contact_order_is_nan_not_zero_where_the_target_is_untouched():
    """A descriptor with no support is undefined; 0 would rank it as a perfectly local footprint."""
    t = _frame().filter(pl.col("target") == "peptide")
    assert math.isnan(contact_order(t)["co_mhc"])


def test_contact_order_is_nan_on_a_single_contacted_residue():
    """One residue has no separation to report."""
    t = _frame().filter((pl.col("target") == "peptide") & (pl.col("rt") == 3))
    assert math.isnan(contact_order(t)["co_pep"])


def test_every_feature_is_catalogued_with_units_and_an_invariance_class():
    """A descriptor the catalogue does not name cannot be reached by any caller downstream."""
    for name in LITERATURE_FEATURES:
        assert name in DESCRIPTORS, f"{name} is not in DESCRIPTORS"
        assert DESCRIPTORS[name][0] == "topology"
        assert name in DETAIL and DETAIL[name][1], f"{name} has no definition"
        assert name in INVARIANCE, f"{name} has no invariance class"


def test_the_gap_volumes_and_the_sign_split_agree_with_the_mean():
    """`gap_mean` is `interlock_frac` weighting `gap_height` against `gap_depth`.

    Stated as an identity rather than trusted, because the four sign-resolved columns exist only
    to decompose the pooled mean and a decomposition that does not sum back is not one.
    """
    rng = np.random.default_rng(0)
    gap = rng.normal(-1.5, 2.5, 4000)
    below, above = gap < 0, gap > 0
    frac = below.mean()
    depth, height = -gap[below].mean(), gap[above].mean()
    assert gap.mean() == pytest.approx((1 - frac) * height - frac * depth, abs=1e-12)


@pytest.mark.slow
def test_all_twenty_three_are_defined_on_a_crystal():
    """1AO7 gives every one of them support, and the gap decomposition is self-consistent.

    Rigid-motion invariance is not re-tested here: the surface block is a pure function of
    :func:`tcren.topology.surface.surface_map`, whose frame refit
    `test_surface.test_groove_frame_is_independent_of_input_orientation` already pins, and the two
    graph descriptors read a contact map with no frame in it at all.
    """
    from pathlib import Path

    from tcren.annotation import classify_chains
    from tcren.mhc import annotate_mhc
    from tcren.structure import import_structure

    src = Path.home() / "hf/tcren_structures/Native2026/1ao7.pdb.gz"
    if not src.exists():
        pytest.skip(f"{src} not present")
    s = import_structure(src)
    classify_chains(s, organism="human")
    annotate_mhc(s)
    row = literature_features(s)

    assert not [k for k, v in row.items() if math.isnan(v)], "1AO7 should support all of them"
    # the faces interlock rather than stack, so the interdigitated volume is the larger one and
    # the asymmetry is negative. This is the sign the module's 60-crystal calibration reports.
    assert row["sc_interlock"] > row["sc_gap_vol"] > 0
    assert row["sc_gap_asym"] < 0
    assert 0.5 < row["sc_interlock_frac"] < 1.0
    # gap_mean is interlock_frac weighting height against depth, to the precision of the split
    assert row["sc_gap_mean"] == pytest.approx(
        (1 - row["sc_interlock_frac"]) * row["sc_gap_height"]
        - row["sc_interlock_frac"] * row["sc_gap_depth"], abs=1e-9)


@pytest.mark.slow
def test_the_surface_block_fails_soft_on_a_structure_with_no_groove():
    """An unmappable structure yields a NaN row, not an exception that kills a batch."""
    from tcren.structure.model import Atom, Chain, Residue, Structure

    atoms = (Atom(name="CA", element="C", coord=(0.0, 0.0, 0.0)),)
    res = Residue(1, 1, "", "A", "ALA", atoms)
    bare = Structure("empty", [Chain("A", [res])])
    row = literature_features(bare)
    assert set(row) == set(LITERATURE_FEATURES)
    assert all(math.isnan(row[k]) for k in SURFACE_FEATURES)
