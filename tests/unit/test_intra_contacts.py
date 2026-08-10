"""Unit tests for the intra-chain contact scope and the peptide-internal wrapper.

``all_atom_contacts`` has always masked to inter-chain residue pairs, so no score in the
package can see an interaction a chain makes with itself. ``scope`` opens that mask. The
first test is the guard that matters: the default must stay byte-identical to the legacy
output, because every interface energy in the package is built on it.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from tcren.contacts.geometry import all_atom_contacts, peptide_internal_contacts
from tcren.structure.io import parse_structure
from tcren.structure.model import Atom, Chain, Residue, Structure

ASSET = Path(__file__).resolve().parents[1] / "assets" / "pdb" / "1ao7.pdb"


def _res(seq, coords):
    atoms = tuple(
        Atom(name=f"X{i}", element="C", coord=np.asarray(c, dtype=float))
        for i, c in enumerate(coords)
    )
    return Residue(
        seq_index=seq, pdb_index=seq, insertion_code="", aa="A", resname="ALA", atoms=atoms
    )


def _two_chain_structure():
    # Chain A is a 4-residue run laid out so that residues 0 and 3 touch (a turn), while
    # 0-1 and 1-2 touch only because they are sequence neighbours. Chain B sits 2 A off
    # residue 0, giving one inter-chain pair.
    a = Chain(
        chain_id="A",
        chain_type="PEPTIDE",
        residues=[
            _res(0, [(0.0, 0.0, 0.0)]),
            _res(1, [(3.0, 0.0, 0.0)]),
            _res(2, [(6.0, 0.0, 0.0)]),
            _res(3, [(3.0, 3.0, 0.0)]),
        ],
    )
    b = Chain(chain_id="B", chain_type="TRA", residues=[_res(0, [(0.0, 2.0, 0.0)])])
    return Structure(pdb_id="t", chains=[a, b])


def test_default_scope_is_byte_identical_to_legacy():
    """The no-regression guard: every interface score depends on this staying put."""
    structure = parse_structure(ASSET)
    legacy = all_atom_contacts(structure, cutoff=5.0)
    scoped = all_atom_contacts(structure, cutoff=5.0, scope="inter")
    assert legacy.equals(scoped)
    assert (legacy["chain.id.from"] != legacy["chain.id.to"]).all()


def test_intra_scope_returns_only_same_chain_pairs():
    structure = _two_chain_structure()
    intra = all_atom_contacts(structure, cutoff=5.0, scope="intra")
    assert intra.height > 0
    assert (intra["chain.id.from"] == intra["chain.id.to"]).all()
    # A residue never contacts itself.
    assert (intra["residue.index.from"] != intra["residue.index.to"]).all()


def test_all_scope_is_the_union():
    structure = _two_chain_structure()
    inter = all_atom_contacts(structure, cutoff=5.0, scope="inter")
    intra = all_atom_contacts(structure, cutoff=5.0, scope="intra")
    both = all_atom_contacts(structure, cutoff=5.0, scope="all")
    assert both.height == inter.height + intra.height


def test_bad_scope_raises():
    with pytest.raises(ValueError, match="scope must be"):
        all_atom_contacts(_two_chain_structure(), scope="peptide")


def test_peptide_internal_drops_sequence_neighbours():
    structure = _two_chain_structure()
    contacts = peptide_internal_contacts(structure, cutoff=5.0, min_seq_sep=3)
    sep = (contacts["residue.index.to"] - contacts["residue.index.from"]).abs()
    assert (sep >= 3).all()
    # Only the 0-3 turn survives; the neighbour pairs 0-1, 1-2, 2-3 are covalent geometry.
    assert contacts.height == 1
    assert contacts["residue.index.from"][0] == 0
    assert contacts["residue.index.to"][0] == 3


def test_intra_scope_on_a_deposited_structure_stays_within_chains():
    structure = parse_structure(ASSET)
    intra = all_atom_contacts(structure, cutoff=5.0, scope="intra")
    assert intra.height > 0
    assert (intra["chain.id.from"] == intra["chain.id.to"]).all()
    assert (intra["residue.index.from"] != intra["residue.index.to"]).all()


def test_peptide_internal_is_restricted_to_the_peptide_chain():
    structure = _two_chain_structure()
    contacts = peptide_internal_contacts(structure, cutoff=5.0)
    peptide = next(c for c in structure.chains if c.chain_type == "PEPTIDE")
    assert contacts.height > 0
    assert (contacts["chain.id.from"] == peptide.chain_id).all()
    assert (contacts["chain.id.to"] == peptide.chain_id).all()
    assert "n_atom_contacts" in contacts.columns


def test_peptide_internal_without_a_peptide_chain_is_empty():
    a = Chain(chain_id="A", chain_type="TRA", residues=[_res(0, [(0.0, 0.0, 0.0)])])
    structure = Structure(pdb_id="t", chains=[a])
    contacts = peptide_internal_contacts(structure)
    assert contacts.height == 0
    assert "chain.id.from" in contacts.columns
