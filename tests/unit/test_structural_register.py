"""Backbone correspondence is independent of motif, direction and author numbering."""
import copy

import numpy as np
import pytest

from tcren.refine.structural_register import select_peptide, structural_register
from tcren.structure.model import Atom, Chain, Residue, Structure


def test_register_and_fused_selection(monkeypatch):
    xyz = np.array([[i * 3., (i % 2) * .8, (i % 3) * .2] for i in range(11)])
    residues = [Residue(i, i - 24, "", "A", "ALA", (Atom("CA", "C", x),)) for i, x in enumerate(xyz)]
    mhc = [Chain("A", [], "MHCa", "MHCII"), Chain("B", [], "MHCb", "MHCII")]
    ref = Structure("ref", mhc + [Chain("P", residues[1:10], "PEPTIDE", "PEPTIDE")])
    query = Structure("query", copy.deepcopy(mhc) + [Chain("P", list(residues), "PEPTIDE", "PEPTIDE")])
    anchors = np.eye(3)
    monkeypatch.setattr("tcren.refine.structural_register._matched_anchors", lambda *a: (anchors, anchors))
    result = structural_register(query, ref, 1)
    assert result["core_start"][0] == 2 and result["direction"][0] == 1
    assert result["core_rmsd_A"][0] < 1e-12 and result["runner_up_margin_A"][0] > 0
    query.chains[-1].residues.reverse()
    assert structural_register(query, ref, 1)["direction"][0] == -1
    query.chains[-1].residues = residues[:6]
    with pytest.raises(ValueError, match="no complete"):
        structural_register(query, ref, 1)
    # Auxiliary peptide must not mask a fused antigen. Keep physical source indices.
    query.chains[1].residues = residues
    select_peptide(query, "B", -24, -16)
    assert query.by_type("PEPTIDE")[0].residues == residues[:9]
    assert len(query.chains[1].residues) == 2
    assert query.chains[2].chain_type is None


def test_modified_geometry_does_not_become_arginine(tmp_path):
    from tcren.structure import parse_structure
    pdb = tmp_path / "modified.pdb"
    pdb.write_text("HETATM    1  CA  CIR P  74       1.000   2.000   3.000  1.00 20.00           C  \nEND\n")
    s = parse_structure(pdb, keep_modified=True)
    assert s.chains[0].sequence() == "X"
    assert s.chains[0].residues[0].resname == "CIR"
