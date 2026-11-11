"""Scene contacts preserve witnesses, residue identity and missing CA diagnostics."""
import numpy as np
import pytest

from tcren.complementarity.frame import CanonicalFrame
from tcren.complementarity.scene import build_scene
from tcren.structure.model import Atom, Chain, Residue, Structure


def residue(index, x, *, ca=True):
    atoms = [Atom("N", "N", np.array([x, 0., 0.]))]
    if ca:
        atoms.append(Atom("CA", "C", np.array([x, 1., 0.])))
    return Residue(index, 42, "A" if index else "", "G", "GLY", tuple(atoms))


FRAME = CanonicalFrame(np.eye(3), np.zeros(3), "ref", "apex")


def test_all_atom_witnesses_and_insertion_identity():
    structure = Structure("test", [Chain("P", [residue(0, 0)], "PEPTIDE"),
                                   Chain("Q", [residue(1, 3)], "PEPTIDE")])
    scene = build_scene(structure, FRAME, cutoff=4)
    assert len(scene["atom_contacts"]) == 4
    assert scene["residue_contacts"] == [dict(source=0, target=1, contact_count=4)]
    assert scene["residues"][0]["contact_count"] == 4
    assert scene["residues"][0]["contact_degree"] == 1
    assert scene["residues"][0]["atom_count"] == 2
    assert scene["atoms"][0]["contact_count"] == 2
    assert scene["atoms"][0]["partner_residue_count"] == 1
    assert "xyz" not in scene["residues"][0]
    assert scene["residues"][1]["insertion"] == "A"
    assert scene["residues"][1]["number"] == 42


def test_zero_contacts_and_missing_ca_retained():
    structure = Structure("test", [Chain("P", [residue(0, 0, ca=False)], "PEPTIDE"),
                                   Chain("Q", [residue(0, 30)], "PEPTIDE")])
    scene = build_scene(structure, FRAME)
    assert scene["residues"][0]["ca_atom"] is None
    assert scene["residues"][0]["contact_count"] == 0
    assert scene["atom_contacts"] == scene["residue_contacts"] == []
    assert scene["warnings"]


def test_trace_does_not_bridge_coordinate_gap():
    structure = Structure("test", [Chain("P", [residue(0, 0), residue(1, 30)], "PEPTIDE")])
    assert build_scene(structure, FRAME)["trace"] == []


@pytest.mark.parametrize("cutoff", [0, -1, float("nan"), float("inf")])
def test_invalid_cutoff(cutoff):
    with pytest.raises(ValueError):
        build_scene(Structure("empty", []), FRAME, cutoff=cutoff)


def test_empty_scene():
    scene = build_scene(Structure("empty", []), FRAME)
    assert scene["atoms"] == scene["residues"] == []
