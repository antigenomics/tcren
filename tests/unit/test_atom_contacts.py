"""Unit tests for the opt-in atom-atom contact count (``n_atom_contacts``).

The legacy ``all_atom_contacts`` keeps one row per residue pair (the closest heavy-atom
pair). With ``count_atoms=True`` an extra ``n_atom_contacts`` column carries the number
of heavy-atom pairs within the cutoff for that residue pair. These tests pin three
properties: the count is always ``>= 1``, it is ``>=`` the residue-pair row count
(trivially 1 per pair), and ``count_atoms=False`` is byte-identical to the legacy output.
"""

from __future__ import annotations

import numpy as np

from tcren.contacts.geometry import all_atom_contacts
from tcren.structure.model import Atom, Chain, Residue, Structure


def _res(seq, aa, coords):
    atoms = tuple(
        Atom(name=f"X{i}", element="C", coord=np.asarray(c, dtype=float))
        for i, c in enumerate(coords)
    )
    return Residue(
        seq_index=seq, pdb_index=seq, insertion_code="", aa=aa, resname="ALA", atoms=atoms
    )


def _multi_atom_structure():
    # Residue A:0 has three atoms; B:0 has two atoms. Several inter-residue atom pairs
    # fall within 5 A, so n_atom_contacts for the (A:0, B:0) pair must exceed 1.
    a = Chain("A", [_res(0, "G", [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (2.0, 0.0, 0.0)])])
    b = Chain("B", [_res(0, "K", [(3.0, 0.0, 0.0), (4.0, 0.0, 0.0)])])
    return Structure("multi", [a, b])


def test_count_atoms_off_is_byte_exact():
    s = _multi_atom_structure()
    legacy = all_atom_contacts(s, cutoff=5.0)
    assert "n_atom_contacts" not in legacy.columns
    # Adding the column then dropping it must reproduce the legacy frame exactly.
    counted = all_atom_contacts(s, cutoff=5.0, count_atoms=True)
    assert counted.drop("n_atom_contacts").equals(legacy)


def test_n_atom_contacts_at_least_one():
    s = _multi_atom_structure()
    df = all_atom_contacts(s, cutoff=5.0, count_atoms=True)
    assert df.height >= 1
    assert (df["n_atom_contacts"] >= 1).all()


def test_n_atom_contacts_at_least_residue_pair_count():
    # The residue-pair row count is always 1 per kept pair, so every count is >= 1; and a
    # multi-atom residue pair must record strictly more than that single closest pair.
    s = _multi_atom_structure()
    df = all_atom_contacts(s, cutoff=5.0, count_atoms=True)
    # All five atom-atom pairs (3 x 2) are within 5 A: 0-3,0-4,1-3,1-4,2-3,2-4 distances
    # are 3,4,2,3,1,2 -> all <= 5, so the single residue pair has n_atom_contacts == 6.
    assert df.height == 1  # one residue pair
    assert df["n_atom_contacts"][0] == 6
    assert df["n_atom_contacts"][0] >= df.height


def test_count_matches_brute_force():
    s = _multi_atom_structure()
    df = all_atom_contacts(s, cutoff=5.0, count_atoms=True)
    # Brute-force count of heavy-atom pairs within cutoff for the single residue pair.
    a0 = s.chains[0].residues[0].atoms
    b0 = s.chains[1].residues[0].atoms
    brute = sum(
        1
        for x in a0
        for y in b0
        if np.linalg.norm(x.coord - y.coord) <= 5.0
    )
    assert df["n_atom_contacts"][0] == brute
