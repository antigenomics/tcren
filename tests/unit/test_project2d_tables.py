"""Unit tests for the complementarity-map tables (fast — no arda)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from tcren.contacts.geometry import all_atom_contacts
from tcren.project2d.frame import ProjectionResult
from tcren.project2d.tables import classify_contact, contacts_table, residue_markup_table
from tcren.structure import parse_structure
from tcren.structure.model import RegionMarkup

REPO = Path(__file__).resolve().parents[2]
PDB_DIR = REPO / "data" / "PDB_structures"


def test_classify_contact_cases():
    # salt bridge: Asp OD1 ↔ Lys NZ within 4 Å
    assert classify_contact("D", "K", "OD1", "NZ", 3.2) == "salt_bridge"
    assert classify_contact("K", "E", "NZ", "OE2", 3.8) == "salt_bridge"
    # hydrogen bond: N/O atoms ≤ 3.5 Å (not a charged pair)
    assert classify_contact("S", "T", "OG", "OG1", 2.9) == "hydrogen_bond"
    # aromatic: two aromatic residues, ring carbons
    assert classify_contact("F", "Y", "CZ", "CE1", 4.2) == "aromatic"
    # hydrophobic: carbon–carbon, non-aromatic
    assert classify_contact("L", "V", "CD1", "CG1", 4.0) == "hydrophobic"
    # polar: an N/O involved but not an H-bond distance
    assert classify_contact("S", "L", "OG", "CB", 4.5) == "polar"


def test_contacts_table_threshold_validation():
    s = parse_structure(PDB_DIR / "5m01.pdb")
    with pytest.raises(ValueError):
        contacts_table(s, threshold=2.9)
    with pytest.raises(ValueError):
        contacts_table(s, threshold=6.1)


def test_contacts_table_matches_geometry():
    s = parse_structure(PDB_DIR / "5m01.pdb")
    raw = all_atom_contacts(s, cutoff=5.0)
    ct = contacts_table(s, threshold=5.0)
    assert ct.height == raw.height  # same contact set, just renamed + classified
    assert set(ct.columns) == {
        "structure_id", "structure_chain_1", "structure_chain_2",
        "residue_index_1", "residue_index_2", "aa_index_1", "aa_index_2",
        "min_dist", "contact_type", "backbone_1", "backbone_2",
    }
    assert ct["min_dist"].max() <= 5.0
    assert ct["contact_type"].is_in(
        ["salt_bridge", "hydrogen_bond", "aromatic", "hydrophobic", "polar", "other"]
    ).all()


def test_residue_markup_schema_and_keys():
    # Manually type a parsed structure (no arda) to exercise the markup table.
    s = parse_structure(PDB_DIR / "5m01.pdb")
    s.chain("G").chain_type = "TRA"
    s.chain("G").regions = [
        RegionMarkup("CDR3", s.chain("G").residues[90].seq_index,
                     s.chain("G").residues[95].seq_index, "XXXXXX",
                     s.chain("G").residues[90:96])
    ]
    s.chain("P").chain_type = "PEPTIDE"
    mk = residue_markup_table(s)
    assert mk.columns[:8] == [
        "structure_id", "structure_chain", "complex_chain", "complex_region",
        "residue_index", "aa_index", "aa_len", "aa",
    ]
    g = mk.filter(mk["structure_chain"] == "G")
    # aa_index == seq_index (0..N-1); residue_index == pdb numbering; aa_len == chain length
    assert g["aa_index"].to_list() == list(range(g.height))
    assert g["aa_len"].unique().to_list() == [len(s.chain("G").residues)]
    assert g.filter(g["complex_region"] == "cdr3").height == 6
    assert mk.filter(mk["complex_chain"] == "peptide").height == len(s.chain("P").residues)


def test_residue_markup_joins_projection_coords():
    s = parse_structure(PDB_DIR / "5m01.pdb")
    s.chain("P").chain_type = "PEPTIDE"
    keys = [("P", r.seq_index) for r in s.chain("P").residues]
    coords = np.arange(len(keys) * 3, dtype=float).reshape(-1, 3)
    proj = ProjectionResult(keys=keys, coords3d=coords, frame="pca")
    mk = residue_markup_table(s, proj)
    pep = mk.filter(mk["complex_chain"] == "peptide").sort("aa_index")
    assert pep["u"].null_count() == 0
    assert pep["x"][0] == pytest.approx(0.0) and pep["v"][1] == pytest.approx(4.0)
