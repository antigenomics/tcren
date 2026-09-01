"""What each descriptor is invariant under, and why geometry and topology are different questions.

Geometry is the study of properties preserved by distance-preserving transformations; topology is
the study of properties preserved by continuous deformation. The catalogue records which is which,
so a block built from a family cannot silently be a block of counts.
"""

from __future__ import annotations

import pytest

from tcren.cohort import Q_FEATURES_GEOM
from tcren.recognition import DESCRIPTORS, INVARIANCE, INVARIANCE_CLASSES, descriptors
from tcren.reliability import T_FEATURES_TOPO


def test_every_descriptor_is_classified_exactly_once():
    assert set(INVARIANCE) == set(DESCRIPTORS)
    assert set(INVARIANCE.values()) <= set(INVARIANCE_CLASSES)


def test_placement_is_metric_throughout():
    """Distances, angles and direction cosines in the groove frame: isometry invariants, all."""
    assert {INVARIANCE[d] for d in descriptors("placement")} == {"geometric"}


def test_the_betti_block_is_the_topological_one():
    """Only the Betti numbers, the Euler characteristic and their normalized forms qualify."""
    assert descriptors(invariance="topological") == (
        "fp_b0_r7", "fp_b1_r7", "fp_chi_r7", "fp_b0_frac_r7",
        "fp_b0_r8", "fp_b1_r8", "fp_chi_r8", "fp_b0_frac_r8",
    )


def test_persistence_entropy_is_metric_not_topological():
    """The H0 barcode's bar lengths are the MST's edge lengths in angstroms.

    Persistent homology is a metric construction, so the entropy of a length distribution is not a
    homeomorphism invariant however topological the barcode's ranks are.
    """
    assert INVARIANCE["h0_pers_ent"] == "geometric"


def test_the_topology_family_is_mostly_compositional():
    """The diversity measures read the labelling of the 12-/24-cell partition, not the shape."""
    topo = descriptors("topology")
    counts = [d for d in topo if INVARIANCE[d] == "compositional"]
    invariants = [d for d in topo if INVARIANCE[d] == "topological"]
    assert len(counts) == 20
    assert len(invariants) == 8
    assert len(counts) > len(invariants)


def test_the_shipped_blocks_are_dominated_by_counts():
    """The state this classification exists to make visible, pinned so a change is deliberate.

    ``Q`` is named for interface geometry and carries one continuous quantity of four -- an area --
    with no angle, distance or height in it at all. ``T`` is named for shape and carries one
    topological invariant of five. Seven of the nine terms across both blocks are counts over the
    labelled contact set, which is why the two correlate more than their names suggest they should.
    """
    q = [INVARIANCE[f] for f in Q_FEATURES_GEOM]
    t = [INVARIANCE[f] for f in T_FEATURES_TOPO]
    assert q.count("geometric") == 1 and q.count("compositional") == 3
    assert t.count("topological") == 1 and t.count("compositional") == 4
    assert (q + t).count("compositional") == 7


def test_invariance_composes_with_family_and_receptor_filters():
    metric_placement = descriptors("placement", invariance="geometric", tcr_only=True)
    assert set(metric_placement) == set(descriptors("placement"))
    assert descriptors("interface", invariance="geometric") == ("burial", "clash_score")


def test_an_unknown_invariance_class_raises():
    with pytest.raises(ValueError, match="unknown invariance"):
        descriptors(invariance="isometric")
