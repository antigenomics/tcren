"""Unit tests for PDB parsing and the structure data model."""

from __future__ import annotations

from pathlib import Path

import pytest

from tcren.structure import parse_structure

REPO = Path(__file__).resolve().parents[2]
PDB_DIR = REPO / "tests" / "assets" / "pdb"


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


def test_repr_summarises_instead_of_dumping_every_atom(s5m01):
    # The dataclass repr expands every atom and its coordinate (~0.5 MB), which floods a notebook
    # cell and makes any error message that interpolates a structure unreadable.
    r = repr(s5m01)
    assert len(r) < 300
    assert "5m01" in r and f"P:?({len(s5m01.chain('P').residues)})" in r  # unclassified: type '?'
    assert len(repr(s5m01.chain("P"))) < 100


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


# --- B-factors ---------------------------------------------------------------------------------
# parse_structure drops them on purpose -- they are not part of the geometry the package reasons
# about -- so mean_bfactor is the one way to read a generated model's own pLDDT back off disk
# without writing a second PDB parser. The CPL benchmark used to carry exactly that parser.

def test_mean_bfactor_reads_the_column_and_respects_the_chain():
    import math

    from tcren import mean_bfactor

    p = PDB_DIR / "1ao7.pdb"
    whole, chain_a = mean_bfactor(p), mean_bfactor(p, "A")
    assert 1.0 < whole < 200.0                       # a plausible crystallographic B
    assert chain_a != whole                          # one chain is not the whole file
    assert math.isnan(mean_bfactor(p, "Z"))          # a chain that is not in the file


# --- the vectorised PDB fast path ------------------------------------------------------------
@pytest.mark.parametrize("asset", sorted(p.name for p in PDB_DIR.glob("*")))
def test_fast_pdb_path_matches_biopython_exactly(asset, monkeypatch):
    """The fast path is exact or it does not run: every asset must parse identically to Biopython.

    Biopython's PDBParser is 86% of the wall clock of a dataset-scale pass through tcren, so ATOM
    records are sliced as one uint8 array instead. That is only admissible while the two agree
    atom for atom, which is what this asserts across every shipped structure.
    """
    import numpy as np

    from tcren.structure import io

    path = PDB_DIR / asset
    fast = io.parse_structure(path)
    monkeypatch.setattr(io, "_parse_pdb_fast", lambda *a, **k: None)
    ref = io.parse_structure(path)

    assert [c.chain_id for c in fast.chains] == [c.chain_id for c in ref.chains]
    for cf, cr in zip(fast.chains, ref.chains):
        assert len(cf.residues) == len(cr.residues)
        for rf, rr in zip(cf.residues, cr.residues):
            assert (rf.seq_index, rf.pdb_index, rf.insertion_code, rf.aa, rf.resname) == \
                   (rr.seq_index, rr.pdb_index, rr.insertion_code, rr.aa, rr.resname)
            assert len(rf.atoms) == len(rr.atoms)
            for af, ar in zip(rf.atoms, rr.atoms):
                assert (af.name, af.element) == (ar.name, ar.element)
                assert np.array_equal(af.coord, ar.coord)


def test_fast_pdb_path_declines_a_file_with_no_element_column():
    """A blank element column is Biopython's inference problem -- bail out, never guess."""
    from tcren.structure.io import _parse_pdb_fast

    line = "ATOM      1  N   MET A   1      10.000  10.000  10.000  1.00 20.00"
    assert _parse_pdb_fast((line.ljust(80) + "\n").encode(), "x", True) is None
