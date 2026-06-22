"""Backbone-preserving peptide substitution."""

from __future__ import annotations

import numpy as np
import pytest

from tcren.refine import substitute_peptide
from tcren.structure.model import Atom, Chain, Residue, Structure


def _peptide_residue(i, aa, resname):
    # backbone (+ Cβ for non-Gly) + one extra side-chain atom (CG) that substitution must drop
    names = ["N", "CA", "C", "O"] + ([] if aa == "G" else ["CB", "CG"])
    atoms = tuple(Atom(n, n[0], np.array([float(i), j, 0.0])) for j, n in enumerate(names))
    return Residue(seq_index=i, pdb_index=i + 1, insertion_code="", aa=aa, resname=resname, atoms=atoms)


def _structure(seq="AGW"):
    three = {"A": "ALA", "G": "GLY", "W": "TRP"}
    pep = Chain("C", [_peptide_residue(i, a, three[a]) for i, a in enumerate(seq)],
                chain_type="PEPTIDE")
    other = Chain("A", [_peptide_residue(0, "L", "LEU")], chain_type="TRA")
    return Structure("test", [pep, other])


def test_roundtrip_preserves_backbone_drops_sidechain():
    s = _structure("AGW")
    out = substitute_peptide(s, "AGW")  # identity swap
    pep = out.chain("C")
    for orig, new in zip(s.chain("C").residues, pep.residues):
        orig_bb = {a.name: a.coord for a in orig.atoms if a.name in {"N", "CA", "C", "O", "CB"}}
        new_bb = {a.name: a.coord for a in new.atoms}
        assert set(new_bb) == set(orig_bb)                      # CG dropped, backbone+CB kept
        for n in orig_bb:
            assert np.allclose(new_bb[n], orig_bb[n])           # coords preserved


def test_substitution_changes_identity_and_glycine_drops_cb():
    s = _structure("AAA")
    out = substitute_peptide(s, "GAA")
    res0 = out.chain("C").residues[0]
    assert res0.aa == "G" and res0.resname == "GLY"
    assert {a.name for a in res0.atoms} == {"N", "CA", "C", "O"}  # no Cβ for glycine


def test_other_chains_untouched():
    s = _structure("AGW")
    out = substitute_peptide(s, "WWW")
    assert out.chain("A") is s.chain("A")                        # TCR chain shared, unchanged


def test_length_mismatch_and_bad_aa_raise():
    s = _structure("AGW")
    with pytest.raises(ValueError, match="length mismatch"):
        substitute_peptide(s, "AG")
    with pytest.raises(ValueError, match="non-standard"):
        substitute_peptide(s, "AGX")
    with pytest.raises(ValueError, match="no PEPTIDE chain"):
        substitute_peptide(Structure("x", [Chain("A", [], chain_type="TRA")]), "AAA")
