"""Fast synthetic tests for interface clash detection (tcren.clashes)."""

from __future__ import annotations

import numpy as np
import pytest

from tcren.clashes import ClashReport, has_clashes, interface_clashes
from tcren.structure.model import PEPTIDE_TYPE, Atom, Chain, Residue, Structure


def _res(i, resname, aa, atoms):
    return Residue(i, i + 1, "", aa, resname, tuple(atoms))


def _atom(name, el, xyz):
    return Atom(name, el, np.asarray(xyz, float))


def _complex(pep_atoms, partners):
    """partners: list of (chain_id, chain_type, [atoms])."""
    pep = Chain("C", [_res(0, "PHE", "F", pep_atoms)], chain_type=PEPTIDE_TYPE)
    chains = [pep]
    for cid, ctype, atoms in partners:
        chains.append(Chain(cid, [_res(0, "TYR", "Y", atoms)], chain_type=ctype))
    return Structure("synth", chains)


def test_overlapping_atoms_clash():
    # peptide C at origin, MHC O at 2.0 Å → overlap = 1.70 + 1.52 − 2.0 = 1.22 (severe)
    s = _complex([_atom("CZ", "C", [0, 0, 0])], [("D", "MHCa", [_atom("OH", "O", [2.0, 0, 0])])])
    rep = interface_clashes(s)
    assert isinstance(rep, ClashReport)
    assert rep.n_clashes == 1
    assert rep.n_severe == 1
    assert rep.max_overlap == pytest.approx(1.70 + 1.52 - 2.0)
    assert rep.by_partner == {"MHCa": 1}
    assert rep.clashing
    assert rep.worst[0].partner_chain_type == "MHCa"


def test_separated_atoms_clash_free():
    s = _complex([_atom("CZ", "C", [0, 0, 0])], [("D", "MHCa", [_atom("OH", "O", [50.0, 0, 0])])])
    rep = interface_clashes(s)
    assert rep.n_clashes == 0
    assert rep.max_overlap == 0.0
    assert not rep.clashing
    assert not has_clashes(s)


def test_tolerance_boundary():
    # C–C at 3.0 Å: overlap = 3.40 − 3.0 = 0.40; not a clash at tol 0.4, is one at tol 0.3
    s = _complex([_atom("CA", "C", [0, 0, 0])], [("D", "MHCa", [_atom("CA", "C", [3.0, 0, 0])])])
    assert interface_clashes(s, tolerance=0.4).n_clashes == 0
    assert interface_clashes(s, tolerance=0.3).n_clashes == 1


def test_per_partner_breakdown():
    # one clash to TRA, one to MHCa
    s = _complex(
        [_atom("CZ", "C", [0, 0, 0])],
        [("A", "TRA", [_atom("CB", "C", [2.5, 0, 0])]),
         ("D", "MHCa", [_atom("OH", "O", [2.0, 0, 0])])],
    )
    rep = interface_clashes(s)
    assert rep.by_partner == {"TRA": 1, "MHCa": 1}
    assert rep.n_clashes == 2
    assert rep.worst[0].overlap >= rep.worst[1].overlap  # sorted worst-first


def test_hydrogens_ignored():
    # an H right on top of the peptide atom must not count (heavy-atom clashes only)
    s = _complex([_atom("CZ", "C", [0, 0, 0])], [("D", "MHCa", [_atom("H", "H", [0.1, 0, 0])])])
    assert interface_clashes(s).n_clashes == 0


def test_no_peptide_chain_raises():
    s = Structure("x", [Chain("D", [_res(0, "TYR", "Y", [_atom("OH", "O", [0, 0, 0])])], chain_type="MHCa")])
    with pytest.raises(ValueError, match="no peptide chain"):
        interface_clashes(s)


def test_intra_peptide_pairs_not_counted():
    # two peptide atoms overlapping each other, no partner → clash-free (only cross pairs count)
    pep = Chain("C", [_res(0, "PHE", "F", [_atom("CZ", "C", [0, 0, 0]), _atom("CE", "C", [1.0, 0, 0])])],
                chain_type=PEPTIDE_TYPE)
    mhc = Chain("D", [_res(0, "TYR", "Y", [_atom("OH", "O", [50, 0, 0])])], chain_type="MHCa")
    assert interface_clashes(Structure("x", [pep, mhc])).n_clashes == 0
