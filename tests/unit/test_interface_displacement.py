import numpy as np
import pytest

from tcren.structure.model import Atom, Chain, RegionMarkup, Residue, Structure


def test_pull_loses_original_cycle_and_leaves_coordinates_unchanged():
    from tcren.mechanics.displacement import interface_displacement

    def residue(i, x, y):
        return Residue(i, i, "", "A", "ALA", (Atom("CA", "C", np.array([x, y, 0.])),))

    left = [residue(0, 0, -1), residue(1, 0, 1)]
    right = [residue(0, 4.5, -1), residue(1, 4.5, 1)]
    pose = Structure("pose", [
        Chain("T", left, "TRA", "TRAB", regions=[RegionMarkup("CDR3", 0, 1, "AA", left)]),
        Chain("P", right, "PEPTIDE", "PEPTIDE"),
    ])
    rows = {row.direction: row for row in interface_displacement(pose)
            if row.interface == "tcr_peptide"}
    assert rows["rest"].edges == 4
    assert rows["rest"].beta1 == 1
    assert rows["pull"].edges == rows["pull"].beta1 == 0
    assert rows["pull"].retained == 0
    assert rows["pull"].elastic_work > 0
    assert rows["shear+1"].retained == rows["shear-1"].retained
    assert rows["rest"].elastic_work == pytest.approx(0)
    np.testing.assert_array_equal(pose.chain("T").residues[0].ca, [0, -1, 0])
