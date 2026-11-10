import numpy as np

from tcren.complementarity.frame import CanonicalFrame
from tcren.complementarity.scene import build_scene
from tcren.structure.model import Atom, Chain, Residue, Structure

FRAME = CanonicalFrame(np.eye(3), np.zeros(3), "ref", "apex")


def test_disulfide_is_covalent_and_not_a_noncovalent_contact():
    def residue(x):
        return Residue(0, 1, "", "C", "CYS", (Atom("SG", "S", np.array([x, 0., 0.])),))
    structure = Structure("disulfide", [Chain("A", [residue(0)], "TRA"),
                                        Chain("B", [residue(2.05)], "TRB")])
    scene = build_scene(structure, FRAME)
    assert len(scene["bonds"]) == 1
    assert scene["bonds"][0]["kind"] == "disulfide"
    assert not scene["atom_contacts"]
    assert [a["covalent_degree"] for a in scene["atoms"]] == [1, 1]
    assert [a["contact_count"] for a in scene["atoms"]] == [0, 0]


def test_named_connectivity_rejects_close_unbonded_sidechain_atoms():
    atoms = tuple(Atom(name, "C", np.array(xyz, dtype=float)) for name, xyz in
                  [("CG1", [0, 0, 0]), ("CG2", [1.5, 0, 0]), ("CB", [.75, 1.3, 0])])
    structure = Structure("valine", [Chain("A", [Residue(0, 1, "", "V", "VAL", atoms)], "TRA")])
    scene = build_scene(structure, FRAME)
    pairs = {(edge["source"], edge["target"]) for edge in scene["bonds"]}
    assert pairs == {(0, 2), (1, 2)}
    assert (0, 1) not in pairs


def test_framework_and_unannotated_atoms_are_retained():
    residue = Residue(0, 1, "", "G", "GLY", (Atom("CA", "C", np.zeros(3)),))
    scene = build_scene(Structure("framework", [Chain("A", [residue], "TRA")]), FRAME)
    assert scene["residues"][0]["region"] == "UNANNOTATED"
    assert scene["residues"][0]["ca_atom"] == 0
