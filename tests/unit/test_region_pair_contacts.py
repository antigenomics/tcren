"""Unit tests for region-pair contact summaries (synthetic structure, no mmseqs)."""

from __future__ import annotations

import numpy as np

from tcren.project2d import region_pair_contacts, region_pair_summary
from tcren.structure.model import Atom, Chain, RegionMarkup, Residue, Structure


def _atom(name: str, xyz) -> Atom:
    return Atom(name=name, element=name[0], coord=np.asarray(xyz, dtype=float))


def _synthetic() -> Structure:
    """An Asp (peptide) and a Lys (TRB CDR3) placed for one salt-bridge contact."""
    asp = Residue(0, 1, "", "D", "ASP",
                  (_atom("CA", (0, 0, 0)), _atom("CB", (1, 0, 0)), _atom("OD1", (2, 0, 0))))
    lys = Residue(0, 100, "", "K", "LYS",
                  (_atom("CA", (2, 0, 6)), _atom("CB", (2, 0, 5)), _atom("NZ", (2, 0, 3))))
    pep = Chain("P", [asp], chain_type="PEPTIDE")
    pep.regions = [RegionMarkup("PEPTIDE", 0, 0, "D", [asp])]
    trb = Chain("E", [lys], chain_type="TRB")
    trb.regions = [RegionMarkup("CDR3", 0, 0, "K", [lys])]
    return Structure(pdb_id="syn", chains=[pep, trb])


def test_region_pair_contacts_closest_classifies_and_orders():
    df = region_pair_contacts(_synthetic(), kind="closest")
    assert df.height == 1
    row = df.row(0, named=True)
    # Canonical ordering puts 'peptide' before 'trb'.
    assert (row["complex_chain_1"], row["region_1"]) == ("peptide", "peptide")
    assert (row["complex_chain_2"], row["region_2"]) == ("trb", "cdr3")
    assert row["contact_type"] == "salt_bridge"
    assert row["min_dist"] == 3.0


def test_region_pair_summary_bond_breakdown_and_thresholds():
    s = _synthetic()
    closest = region_pair_summary(s, kind="closest")
    assert closest.height == 1
    r = closest.row(0, named=True)
    assert r["n_contacts"] == 1 and r["n_salt_bridge"] == 1 and r["n_hydrogen_bond"] == 0
    # Cβ (5.1 Å ≤ 8) and Cα (6.3 Å ≤ 12) representative-atom contacts each give one pair.
    assert region_pair_summary(s, kind="cb")["n_contacts"].sum() == 1
    assert region_pair_summary(s, kind="ca")["n_contacts"].sum() == 1
