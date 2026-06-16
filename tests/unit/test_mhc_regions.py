"""Unit tests for MHC groove partitioning (no network / mmseqs needed)."""

from __future__ import annotations

from pathlib import Path

from tcren.mhc.regions import partition_chain
from tcren.structure import parse_structure

REPO = Path(__file__).resolve().parents[2]
PDB_DIR = REPO / "data" / "PDB_structures"


def test_class_i_groove_regions():
    chain = parse_structure(PDB_DIR / "1ao7.pdb").chain("A")
    regions = {r.region_type: r for r in partition_chain(chain, "MHCI", "MHCa")}
    assert set(regions) == {"HELIX_A1", "HELIX_A2", "GROOVE_FLOOR"}

    # HELIX_A1 should map to the α1 helix (~mature residues 57-84).
    a1_pdb = [r.pdb_index for r in regions["HELIX_A1"].residues]
    assert min(a1_pdb) >= 50 and max(a1_pdb) <= 90
    # HELIX_A2 to the α2 helix (~138-180).
    a2_pdb = [r.pdb_index for r in regions["HELIX_A2"].residues]
    assert min(a2_pdb) >= 130 and max(a2_pdb) <= 185
    # Floor sits below both helices (spans α1 and α2 sheet strands).
    assert len(regions["GROOVE_FLOOR"].residues) > 50


def test_class_ii_two_helices():
    s = parse_structure(PDB_DIR / "4ozg.pdb")
    a = {r.region_type for r in partition_chain(s.chain("A"), "MHCII", "MHCa")}
    b = {r.region_type for r in partition_chain(s.chain("B"), "MHCII", "MHCb")}
    assert "HELIX_A1" in a and "GROOVE_FLOOR" in a
    assert "HELIX_B1" in b and "GROOVE_FLOOR" in b


def test_b2m_has_no_groove():
    # B2M is not a groove chain; partitioning returns nothing.
    chain = parse_structure(PDB_DIR / "1ao7.pdb").chain("B")
    assert partition_chain(chain, "MHCI", "B2M") == []
