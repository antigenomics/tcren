import numpy as np
import pytest

from tcren.contacts.opportunities import native_opportunities
from tcren.structure.model import Atom, Chain, Residue, Structure


def _residue(i, *atoms):
    return Residue(i, i, "", "A", "ALA", tuple(
        Atom(name, "C", np.array([x, 0., 0.])) for name, x in atoms
    ))


def test_native_opportunities_is_a_directed_union_of_shell_and_contacts():
    structure = Structure("pose", [
        Chain("Z", [_residue(0, ("CA", 0), ("CG", 10))], "TRA"),
        Chain("A", [
            _residue(0, ("CA", 3), ("CB", 4)),
            _residue(1, ("CA", -8)),
            _residue(2, ("CA", 20), ("CG", 11)),
            _residue(3, ("CG", 11)),
        ], "PEPTIDE"),
        Chain("M", [_residue(0, ("CA", 6))], "B2M"),
        Chain("X", [_residue(0, ("CA", 6))], "IGH"),
    ])
    t = native_opportunities(structure)
    assert set(t["interface"]) == {"tcr_peptide", "tcr_mhc", "peptide_mhc"}
    assert "X" not in t["chain.id.from"]
    rows = t.filter(t["interface"] == "tcr_peptide").sort("residue.index.to")
    assert rows["chain.id.from"].to_list() == ["Z"] * 3
    assert rows["residue.index.to"].to_list() == [0, 1, 2]
    assert rows["within_shell"].to_list() == [True, True, False]
    assert rows["contacted"].to_list() == [True, False, True]
    assert rows["cb_distance"].to_list() == pytest.approx([4, 8, 20])
    assert rows["min_dist"].to_list() == [3, None, 1]
    assert rows["n_atom_pairs"].to_list() == [2, 0, 1]
    assert rows["n_heavy_from"].to_list() == [2] * 3


def test_native_opportunities_empty_input_preserves_columns():
    rows = native_opportunities(Structure("empty", []))
    assert not rows.height
    assert {"contacted", "within_shell", "atom.from", "interface"} <= set(rows.columns)
