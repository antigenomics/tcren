"""Unit tests for contact geometry against a brute-force reference."""

from __future__ import annotations

import numpy as np

from tcren.contacts.geometry import all_atom_contacts
from tcren.structure.model import Atom, Chain, Residue, Structure


def _res(chain_seq, aa, coords):
    atoms = tuple(
        Atom(name=f"X{i}", element="C", coord=np.asarray(c, dtype=float))
        for i, c in enumerate(coords)
    )
    return Residue(seq_index=chain_seq, pdb_index=chain_seq, insertion_code="", aa=aa, resname="ALA", atoms=atoms)


def _toy_structure():
    # Two chains; residue coordinates chosen to straddle the 5 A cutoff.
    chain_a = Chain(
        "A",
        [
            _res(0, "G", [(0.0, 0.0, 0.0)]),
            _res(1, "A", [(10.0, 0.0, 0.0)]),
        ],
    )
    chain_b = Chain(
        "B",
        [
            _res(0, "K", [(4.0, 0.0, 0.0)]),  # 4 A from A:0  -> contact
            _res(1, "L", [(0.0, 6.0, 0.0)]),  # 6 A from A:0  -> no contact
            _res(2, "M", [(5.0, 0.0, 0.0)]),  # exactly 5 A from A:0 -> contact (<=)
        ],
    )
    return Structure("toy", [chain_a, chain_b])


def _brute_force(structure, cutoff):
    pairs = set()
    chains = structure.chains
    for ci, c1 in enumerate(chains):
        for c2 in chains[ci + 1 :]:
            for r1 in c1.residues:
                for r2 in c2.residues:
                    dmin = min(
                        np.linalg.norm(a1.coord - a2.coord)
                        for a1 in r1.atoms
                        for a2 in r2.atoms
                    )
                    if dmin <= cutoff:
                        pairs.add((c1.chain_id, r1.seq_index, c2.chain_id, r2.seq_index))
    return pairs


def test_contacts_match_brute_force():
    s = _toy_structure()
    df = all_atom_contacts(s, cutoff=5.0)
    got = {
        (r["chain.id.from"], r["residue.index.from"], r["chain.id.to"], r["residue.index.to"])
        for r in df.iter_rows(named=True)
    }
    assert got == _brute_force(s, 5.0)


def test_cutoff_is_inclusive():
    s = _toy_structure()
    df = all_atom_contacts(s, cutoff=5.0)
    # A:0 <-> B:2 is at exactly 5.0 A and must be present.
    hit = df.filter(
        (df["chain.id.from"] == "A")
        & (df["residue.index.from"] == 0)
        & (df["chain.id.to"] == "B")
        & (df["residue.index.to"] == 2)
    )
    assert hit.height == 1
    assert abs(hit["dist"][0] - 5.0) < 1e-9


def test_no_intra_chain_contacts():
    s = _toy_structure()
    df = all_atom_contacts(s, cutoff=100.0)
    assert (df["chain.id.from"] != df["chain.id.to"]).all()


def test_closest_atom_per_residue_pair():
    # Two atoms on one residue; the reported distance must be the minimum.
    a = Chain("A", [_res(0, "G", [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)])])
    b = Chain("B", [_res(0, "K", [(3.0, 0.0, 0.0)])])
    df = all_atom_contacts(Structure("t", [a, b]), cutoff=5.0)
    assert df.height == 1
    assert abs(df["dist"][0] - 2.0) < 1e-9  # 3.0 - 1.0
