"""Substituting a residue must move atoms, not just letters.

The virtual ΔΔG path re-indexes a mutant sequence over the wild type's contact map, so a contact
that exists only because a long side chain reaches across survives the mutation on paper. For an
alanine target that is checkable against exact geometry: alanine's heavy atoms are backbone + Cβ
and nothing else, so no rotamer, no relaxation and no choice enter.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from tcren.annotation import classify_chains
from tcren.contactmap import ContactMap
from tcren.ddg import alanine_scan, ddg, reference_delta
from tcren.potential import tcren2
from tcren.refine.substitute import substitute_peptide, virtual_cb
from tcren.structure import parse_structure

PDB = Path(__file__).resolve().parents[1] / "assets" / "pdb" / "1ao7.pdb"
ALANINE_ATOMS = {"N", "CA", "C", "O", "CB"}


@pytest.fixture(scope="module")
def complex_():
    s = parse_structure(str(PDB))
    classify_chains(s)
    return s


@pytest.fixture(scope="module")
def peptide(complex_):
    return next(c.sequence() for c in complex_.chains if c.chain_type == "PEPTIDE")


def _pep_chain(structure):
    return next(c for c in structure.chains if c.chain_type == "PEPTIDE")


# --- geometry -------------------------------------------------------------------------


def test_virtual_cb_reproduces_the_crystallographic_one(complex_):
    """Built from N/Cα/C alone, against every real Cβ in the complex."""
    d = []
    for chain in complex_.chains:
        for r in chain.residues:
            at = {a.name: a.coord for a in r.atoms}
            if r.aa == "G" or not {"N", "CA", "C", "CB"} <= at.keys():
                continue
            d.append(np.linalg.norm(virtual_cb(at["N"], at["CA"], at["C"]) - at["CB"]))
    assert len(d) > 200
    assert np.median(d) < 0.15
    assert np.percentile(d, 99) < 0.5


def test_virtual_cb_sits_at_a_real_bond_length(complex_):
    r = next(r for c in complex_.chains for r in c.residues
             if {"N", "CA", "C"} <= {a.name for a in r.atoms})
    at = {a.name: a.coord for a in r.atoms}
    assert np.linalg.norm(virtual_cb(at["N"], at["CA"], at["C"]) - at["CA"]) == pytest.approx(
        1.53, abs=0.05)


# --- what the substitution produces ------------------------------------------------------


def test_poly_alanine_is_exactly_alanine(complex_, peptide):
    """Every residue must end up with alanine's five heavy atoms — no more, no fewer."""
    out = substitute_peptide(complex_, "A" * len(peptide))
    for r in _pep_chain(out).residues:
        assert r.aa == "A"
        assert {a.name for a in r.atoms} == ALANINE_ATOMS


def test_a_glycine_position_gets_a_cbeta_built(complex_, peptide):
    """1ao7's LLFGYPVYV has a glycine at position 4, which has no Cβ to keep."""
    i = peptide.index("G")
    assert not any(a.name == "CB" for a in _pep_chain(complex_).residues[i].atoms)
    out = substitute_peptide(complex_, "A" * len(peptide))
    assert any(a.name == "CB" for a in _pep_chain(out).residues[i].atoms)


def test_mutating_to_glycine_removes_the_cbeta(complex_, peptide):
    out = substitute_peptide(complex_, "G" * len(peptide))
    assert all({a.name for a in r.atoms} == ALANINE_ATOMS - {"CB"}
               for r in _pep_chain(out).residues)


def test_the_backbone_does_not_move(complex_, peptide):
    before = {a.name: a.coord for a in _pep_chain(complex_).residues[0].atoms}
    after = {a.name: a.coord for a in
             _pep_chain(substitute_peptide(complex_, "A" * len(peptide))).residues[0].atoms}
    for name in ("N", "CA", "C", "O"):
        assert np.allclose(before[name], after[name])


# --- what it does to the contact map and the score ----------------------------------------


def test_poly_alanine_loses_the_contacts_its_side_chains_were_making(complex_, peptide):
    """Truncating at Cβ can only remove contacts, and here it removes many."""
    native = ContactMap.from_structure(complex_).interface("tcr_peptide").height
    ala = ContactMap.from_structure(
        substitute_peptide(complex_, "A" * len(peptide))).interface("tcr_peptide").height
    assert ala < native
    assert ala / native < 0.75


def test_the_structural_reference_differs_from_the_virtual_one(complex_, peptide):
    cm, pot = ContactMap.from_structure(complex_), tcren2()
    virtual = reference_delta(cm, peptide, pot)
    structural = reference_delta(cm, peptide, pot, structure=complex_)
    assert abs(structural - virtual) > 0.5


def test_a_self_substitution_is_not_free_structurally(complex_, peptide):
    """Rebuilding the native sequence truncates it to Cβ, so it is NOT the identity.

    Pinned because it is the trap: the structural path is a comparison against backbone + Cβ, and
    reading it as "the same peptide, moved" would be wrong.
    """
    cm, pot = ContactMap.from_structure(complex_), tcren2()
    assert ddg(cm, peptide, peptide, pot) == pytest.approx(0.0, abs=1e-12)
    assert ddg(cm, peptide, peptide, pot, structure=complex_) != pytest.approx(0.0, abs=1e-6)


def test_a_non_peptide_interface_is_still_exactly_zero(complex_, peptide):
    cm, pot = ContactMap.from_structure(complex_), tcren2()
    assert ddg(cm, peptide, "A" * len(peptide), pot,
               interface="tcr_mhc", structure=complex_) == 0.0


def test_alanine_scan_zeroes_a_position_that_is_already_alanine(complex_):
    """A no-op substitution must cost nothing, structurally as well as virtually."""
    pep = _pep_chain(complex_).sequence()
    seq = "A" + pep[1:]
    s2 = substitute_peptide(complex_, seq)
    scan = alanine_scan(ContactMap.from_structure(s2), seq, tcren2(), structure=s2)
    assert scan.filter(scan["pos"] == 0)["ddG"].item() == pytest.approx(0.0, abs=1e-12)


def test_alanine_scan_reports_one_row_per_position(complex_, peptide):
    scan = alanine_scan(ContactMap.from_structure(complex_), peptide, tcren2(),
                        structure=complex_)
    assert scan.height == len(peptide)
    assert scan["wt_aa"].to_list() == list(peptide)
