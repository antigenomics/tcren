"""Unit tests for covalent (single-chain) peptide detection and splitting.

No fused-peptide constructs exist in the bundled datasets, so a single-chain pMHC is
assembled synthetically (peptide + Gly/Ser linker + real MHC heavy chain) to exercise
the alignment-based detection and split.
"""

from __future__ import annotations

from pathlib import Path

from tcren.mhc.linker import (
    check_against_mhc,
    detect_linked_peptide,
    split_linked_peptides,
)
from tcren.structure import parse_structure
from tcren.structure.model import Atom, Chain, Residue, Structure

REPO = Path(__file__).resolve().parents[2]
MHC_CHAIN = parse_structure(REPO / "tests" / "assets" / "cgene" / "1ao7_full.pdb").chain("A")

PEPTIDE = "GILGFVFTL"
LINKER = "GGGGSGGGGS"


def _res(i: int, aa: str) -> Residue:
    return Residue(
        seq_index=i, pdb_index=i + 1, insertion_code="", aa=aa, resname="ALA",
        atoms=(Atom("CA", "C", __import__("numpy").zeros(3)),),
    )


def _single_chain_construct() -> Chain:
    """peptide + GS-linker + real MHC heavy chain, as one fused chain."""
    fused_aa = PEPTIDE + LINKER
    residues = [_res(i, aa) for i, aa in enumerate(fused_aa)]
    offset = len(residues)
    for j, r in enumerate(MHC_CHAIN.residues):
        residues.append(
            Residue(offset + j, r.pdb_index, r.insertion_code, r.aa, r.resname, r.atoms)
        )
    return Chain(chain_id="M", residues=residues, chain_type="MHCa")


def test_alignment_check_flags_extension():
    chain = _single_chain_construct()
    check = check_against_mhc(chain)
    assert check.is_mhc
    assert check.best_ref == "MHCI|MHCa"
    assert check.n_term_extra == len(PEPTIDE) + len(LINKER)


def test_detect_strips_linker_and_returns_peptide():
    chain = _single_chain_construct()
    seg = detect_linked_peptide(chain)
    assert seg is not None
    assert "".join(r.aa for r in seg) == PEPTIDE


def test_split_creates_peptide_chain():
    chain = _single_chain_construct()
    structure = Structure(pdb_id="synthetic", chains=[chain])
    split = split_linked_peptides(structure)
    assert split == ["M"]
    peptide_chains = [c for c in structure.chains if c.chain_type == "PEPTIDE"]
    assert len(peptide_chains) == 1
    assert peptide_chains[0].sequence() == PEPTIDE
    # The MHC chain no longer carries the fused peptide/linker.
    assert structure.chain("M").sequence().startswith(MHC_CHAIN.sequence()[:20])


def test_conventional_mhc_chain_not_split():
    # A plain MHC chain (no fused peptide) must be left untouched.
    chain = Chain(chain_id="A", residues=list(MHC_CHAIN.residues), chain_type="MHCa")
    assert detect_linked_peptide(chain) is None
    structure = Structure(pdb_id="conv", chains=[chain])
    assert split_linked_peptides(structure) == []
