"""Tests for peptide register / forced-pose detection and correction (tcren.refine.register)."""

from __future__ import annotations

import math

import numpy as np
import pytest

from tcren.refine.register import RegisterReport, check_register, fix_register
from tcren.structure.model import PEPTIDE_TYPE, Atom, Chain, Residue, Structure


def _res(i, resname, aa, atoms):
    return Residue(i, i + 1, "", aa, resname, tuple(atoms))


def _atom(name, el, xyz):
    return Atom(name, el, np.asarray(xyz, float))


def _pep(seq, base=(0.0, 0.0, 0.0)):
    """A peptide chain, one Cα-only residue per letter (spread out so it is clash-free internally)."""
    residues = [
        _res(i, "ALA", c, [_atom("CA", "C", (base[0] + 4.0 * i, base[1], base[2]))])
        for i, c in enumerate(seq)
    ]
    return Chain("C", residues, chain_type=PEPTIDE_TYPE)


def _mhc(atoms):
    return Chain("D", [_res(0, "GLU", "E", atoms)], chain_type="MHCa")


def test_no_reference_register_undetermined():
    # clean pose, no reference → register cannot be decided but clashes are reported
    s = Structure("m", [_pep("SIINFEKL"), _mhc([_atom("OE1", "O", [500, 0, 0])])])
    rep = check_register(s)
    assert isinstance(rep, RegisterReport)
    assert rep.wrong_register is None
    assert math.isnan(rep.anchor_rmsd)
    assert rep.clashes.n_clashes == 0
    assert not rep.suspect
    assert "no reference" in rep.reason


def test_severe_clash_makes_pose_suspect_without_reference():
    # an MHC atom jammed into the peptide → severe clash → suspect even absent a reference
    s = Structure("m", [_pep("SIINFEKL"), _mhc([_atom("OE1", "O", [0.0, 0, 0])])])
    rep = check_register(s)
    assert rep.wrong_register is None  # still can't call register without a reference
    assert rep.clashes.n_severe >= 1
    assert rep.suspect  # ... but the geometry is flagged


def test_fix_register_length_mismatch_raises():
    model = Structure("m", [_pep("SIINFEKL"), _mhc([_atom("OE1", "O", [50, 0, 0])])])       # 8-mer
    template = Structure("t", [_pep("GILGFVFTL"), _mhc([_atom("OE1", "O", [50, 0, 0])])])    # 9-mer
    with pytest.raises(ValueError, match="length mismatch"):
        fix_register(model, template)


def _annotated_1ao7():
    pytest.importorskip("arda")
    from tcren.annotation import classify_chains
    from tcren.mhc import annotate_mhc
    from tcren.paths import reference_structure_path
    from tcren.structure import parse_structure

    s = parse_structure(reference_structure_path("1ao7"), pdb_id="1ao7")
    classify_chains(s, organism="human")
    annotate_mhc(s)
    return s


@pytest.mark.slow
def test_check_register_reference_real_structure():
    """A structure checked against itself is in-register: anchor-Cα RMSD ≈ 0, not wrong."""
    s = _annotated_1ao7()
    rep = check_register(s, s)
    assert rep.wrong_register is False
    assert rep.anchor_rmsd == pytest.approx(0.0, abs=1e-6)
    assert rep.anchors  # class-I anchors were resolved (P2 + PΩ)


@pytest.mark.slow
def test_fix_register_real_structure():
    """fix_register re-threads via model_peptide (needs the compiled refine kernel)."""
    pytest.importorskip("tcren._refine")  # engine kernel; built in CI, may be absent in a bare dev venv
    s = _annotated_1ao7()
    res = fix_register(s, s, engine="dope", n_steps=50)
    assert res.structure is not None
