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
    """The Betti numbers of the two complexes, and nothing that reads an angstrom.

    Two complexes now qualify, not one. The `fp_*` block is the flag complex built on the contacted
    pMHC Calpha atoms at a chosen radius; the `g_*` block is the bipartite *contact* graph, whose
    component and cycle fractions are that 1-complex's own b0 and b1, size-normalized and needing no
    radius at all. The degree evenness and assortativity terms sit here for the same reason: the
    degree sequence of a graph is a graph invariant, unchanged by any deformation that preserves
    which residue touches which.
    """
    assert descriptors(invariance="topological") == (
        "g_even_tcr", "g_even_pmhc", "g_comp_frac", "g_alg_conn", "g_cyclo_frac", "g_assort",
        "fp_b0_r7", "fp_b1_r7", "fp_chi_r7", "fp_b0_frac_r7",
        "fp_b0_r8", "fp_b1_r8", "fp_chi_r8", "fp_b0_frac_r8",
    )


def test_the_map_block_is_metric_not_topological():
    """The Calpha-against-Cbeta comparison reads angstroms, so it is geometric.

    Same standard as `h0_pers_ent`: `m_face_*` is a mean of a difference of two distances, and
    `ca_cb_agreement_*` correlates two metric maps -- both move under a deformation that leaves the
    contact set alone, so neither is a homeomorphism invariant however shape-like it reads.
    """
    for d in ("m_erank_tp", "m_gap_tp", "m_erank_tm", "m_gap_tm",
              "m_face_tp", "m_face_tm", "ca_cb_agreement_tp", "ca_cb_agreement_tm"):
        assert INVARIANCE[d] == "geometric", d


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
    assert len(counts) == 24
    assert len(invariants) == 14
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
