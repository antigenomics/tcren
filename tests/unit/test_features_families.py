"""The descriptor catalogue as a partition, and the family-selected feature pass.

Two things are worth pinning here. The families must *partition* the catalogue -- every emitted
column belongs to exactly one, and asking for a family returns exactly its members -- because that
is what lets a cross-channel independence claim mean anything. And selecting families must actually
select: `-i topology` should not pay for the energies, and it must return the same numbers the
dedicated topology command returns.
"""
from __future__ import annotations

import pytest

from tcren.footprint import FOOTPRINT_SIZE_FEATURES, footprint_topology_features
from tcren.recognition import (
    _FAMILY_ALIASES,
    DESCRIPTORS,
    FAMILIES,
    TCR_PLACEMENT_FEATURES,
    descriptors,
)


def test_every_descriptor_has_a_known_family():
    bad = {n: fam for n, (fam, _) in DESCRIPTORS.items() if fam not in set(FAMILIES)}
    assert not bad, bad


def test_families_partition_the_catalogue():
    total = sum(len(descriptors(f)) for f in FAMILIES)
    assert total == len(DESCRIPTORS)
    # and no name is claimed twice
    seen = [n for f in FAMILIES for n in descriptors(f)]
    assert len(seen) == len(set(seen))


def test_retired_family_names_still_resolve():
    assert descriptors("geometry") == tuple(
        n for n in DESCRIPTORS if DESCRIPTORS[n][0] in _FAMILY_ALIASES["geometry"])
    assert set(descriptors("geometry")) == set(descriptors("placement")) | set(descriptors("interface"))
    assert descriptors("physics") == descriptors("energetics")


def test_unknown_family_is_an_error_naming_the_options():
    with pytest.raises(ValueError, match="unknown family"):
        descriptors("topolgy")


def test_placement_holds_the_translational_terms_no_angle_can_see():
    p = set(descriptors("placement"))
    assert set(TCR_PLACEMENT_FEATURES) <= p
    # the pose angles belong here too, and the contact chemistry does not
    assert {"crossing_signed", "dock_d", "dock_torsion"} <= p
    assert not p & {"burial", "n_hbond", "Phi_tcr_pep"}


def test_contact_counts_are_interface_size_not_topology():
    """The shape channel must not carry the interface's size, or the two correlate by construction."""
    topo = set(descriptors("topology"))
    assert not topo & set(FOOTPRINT_SIZE_FEATURES)
    assert set(FOOTPRINT_SIZE_FEATURES) <= set(descriptors("interface"))
    assert set(footprint_topology_features()) == topo


def test_the_engaged_pair_count_belongs_to_potts_and_to_no_other_family():
    """`n_contacts` is the Potts count of available pairs that engaged, and nothing else.

    Regression, 2026-08-29: the footprint wrote its CDR-loop tally under the same name, so
    `-i placement,interface,topology` emitted one quantity and `-i ...,potts` the other, silently,
    under one column. The Potts moments are what any energy term must be standardized against.
    """
    assert DESCRIPTORS["n_contacts"] == ("potts", True)
    assert "n_contacts" in descriptors("potts")
    for f in set(FAMILIES) - {"potts"}:
        assert "n_contacts" not in descriptors(f), f
    assert "n_loop_contacts" in FOOTPRINT_SIZE_FEATURES and "n_contacts" not in FOOTPRINT_SIZE_FEATURES


def test_energetics_is_the_energy_channel_and_nothing_else():
    e = set(descriptors("energetics"))
    assert e == {"Phi_tcr_pep", "Phi_tcr_mhc", "Phi_cdr12", "Phi_cdr3a", "Phi_cdr3b",
                 "dPhi_pep_soft", "varPhi_pep_soft", "dPhi_tcr_soft", "varPhi_tcr_soft",
                 "dPhi_tra_soft", "dPhi_trb_soft",
                 "dPhi_tcr_pep", "Phi_pep_mhc", "dPhi_pep_mhc", "Phi_pep_int"}


def test_tcr_only_drops_the_cohort_identity_columns():
    """A pMHC-only column is shared by every receptor of an epitope, so it leaks cohort identity."""
    assert "mhc_class_bin" in descriptors("interface")
    assert "mhc_class_bin" not in descriptors("interface", tcr_only=True)
    assert "Phi_pep_mhc" not in descriptors("energetics", tcr_only=True)


def test_scores_are_excluded_from_every_family_by_default():
    for f in FAMILIES:
        assert not set(descriptors(f)) & {"S", "Q", "T"}
