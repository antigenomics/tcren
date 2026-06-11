"""Unit tests for PDB parsing and the structure data model."""

from __future__ import annotations

from pathlib import Path

import pytest

from tcren.structure import parse_structure

REPO = Path(__file__).resolve().parents[2]
PDB_DIR = REPO / "data" / "PDB_structures"


@pytest.fixture(scope="module")
def s5m01():
    return parse_structure(PDB_DIR / "5m01.pdb")


def test_seq_index_is_sequential_per_chain(s5m01):
    for chain in s5m01.chains:
        indices = [r.seq_index for r in chain.residues]
        assert indices == list(range(len(chain.residues)))


def test_pdb_numbering_preserved(s5m01):
    # Author numbering is kept verbatim and is independent of seq_index.
    pep = s5m01.chain("P")
    assert pep.residues[0].seq_index == 0
    assert all(isinstance(r.pdb_index, int) for r in pep.residues)


def test_sequence_roundtrip(s5m01):
    pep = s5m01.chain("P")
    assert pep.sequence() == "".join(r.aa for r in pep.residues)
    assert len(pep.sequence()) == len(pep.residues)


def test_non_standard_atom_residue_kept_as_x():
    # 5jhd peptide (chain C) begins with the AMN cap, an ATOM record mir keeps as 'X'.
    s = parse_structure(PDB_DIR / "5jhd.pdb")
    pep = s.chain("C")
    assert pep.residues[0].aa == "X"
    assert pep.residues[0].resname == "AMN"


def test_hetatm_modified_residue_skipped():
    # 6v0y peptide (chain C) contains CIR (citrulline) HETATM residues that mir drops.
    s = parse_structure(PDB_DIR / "6v0y.pdb")
    pep = s.chain("C")
    assert "CIR" not in {r.resname for r in pep.residues}
    # The kept residues stay contiguously indexed despite the dropped HETATMs.
    assert [r.seq_index for r in pep.residues] == list(range(len(pep.residues)))


def test_waters_excluded(s5m01):
    for chain in s5m01.chains:
        assert "HOH" not in {r.resname for r in chain.residues}


def test_altloc_atoms_all_retained():
    # 9nmx has a disordered Ser whose alternate conformer is needed for one contact;
    # keeping all altlocs means at least one atom name appears more than once.
    s = parse_structure(PDB_DIR / "9nmx.pdb")
    has_duplicate_atom_name = False
    for chain in s.chains:
        for res in chain.residues:
            names = [a.name for a in res.atoms]
            if len(names) != len(set(names)):
                has_duplicate_atom_name = True
                break
    assert has_duplicate_atom_name
