"""Fast synthetic tests for TCR:peptide contact stability / fragility (tcren.stability)."""

from __future__ import annotations

import numpy as np
import pytest

from tcren.stability import StabilityReport, contact_stability
from tcren.structure.model import PEPTIDE_TYPE, Atom, Chain, Residue, Structure


def _res(i, resname, aa, atoms):
    return Residue(i, i + 1, "", aa, resname, tuple(atoms))


def _atom(name, el, xyz):
    return Atom(name, el, np.asarray(xyz, float))


def _complex(pep_atoms_per_res, tcr_atoms_per_res, tcr_type="TRB"):
    """pep_atoms_per_res / tcr_atoms_per_res: list (per residue) of [atoms]."""
    pep = Chain("C", [_res(i, "GLY", "G", ats) for i, ats in enumerate(pep_atoms_per_res)],
                chain_type=PEPTIDE_TYPE)
    tcr = Chain("B", [_res(i, "TYR", "Y", ats) for i, ats in enumerate(tcr_atoms_per_res)],
                chain_type=tcr_type)
    return Structure("synth", [pep, tcr])


def test_margins_and_fragility():
    # one TCR residue at origin; peptide residues at 3.0 (robust), 4.5 (fragile), 6.0 (beyond cutoff)
    s = _complex([[_atom("CA", "C", [3.0, 0, 0])], [_atom("CA", "C", [4.5, 0, 0])],
                  [_atom("CA", "C", [6.0, 0, 0])]], [[_atom("CA", "C", [0, 0, 0])]])
    rep = contact_stability(s)
    assert isinstance(rep, StabilityReport)
    assert rep.n_contacts == 2  # the 6 Å residue is not a contact
    assert rep.mean_margin == pytest.approx(1.25)
    assert rep.frac_robust == pytest.approx(0.5)
    assert rep.frac_marg_lt1 == pytest.approx(0.5)
    assert rep.exp_lost == pytest.approx(0.25)


def test_closest_atom_defines_the_contact():
    # a residue pair's distance is its *closest* heavy-atom pair, not any atom pair
    s = _complex([[_atom("N", "N", [4.9, 0, 0]), _atom("CA", "C", [3.0, 0, 0])]],
                 [[_atom("CA", "C", [0, 0, 0])]])
    rep = contact_stability(s)
    assert rep.n_contacts == 1
    assert rep.mean_margin == pytest.approx(2.0)  # from the 3.0 Å atom, not the 4.9 Å one


def test_cutoff_and_delta_are_tunable():
    s = _complex([[_atom("CA", "C", [4.5, 0, 0])]], [[_atom("CA", "C", [0, 0, 0])]])
    assert contact_stability(s, cutoff=4.0).n_contacts == 0  # 4.5 Å now beyond cutoff
    assert contact_stability(s, delta=2.0).frac_robust == 0.0  # margin 0.5 < delta 2.0


def test_native_matches_numpy_reference(monkeypatch):
    import tcren.stability as stability

    pep = [[_atom("CA", "C", [d, 0, 0]), _atom("CB", "C", [d, 1.0, 0])] for d in (3.0, 3.8, 4.6)]
    tcr = [[_atom("CA", "C", [0, 0, 0])], [_atom("CA", "C", [0, 2.0, 0])]]
    s = _complex(pep, tcr)
    native = contact_stability(s)
    monkeypatch.setattr(stability, "_geom", None)  # force the numpy reference path
    reference = contact_stability(s)
    assert native.n_contacts == reference.n_contacts > 0
    assert native.mean_margin == pytest.approx(reference.mean_margin)
    assert native.exp_lost == pytest.approx(reference.exp_lost)
    assert native.frac_robust == pytest.approx(reference.frac_robust)


def test_missing_chains_raise():
    pep_only = Structure("x", [Chain("C", [_res(0, "GLY", "G", [_atom("CA", "C", [0, 0, 0])])],
                                       chain_type=PEPTIDE_TYPE)])
    with pytest.raises(ValueError, match="no receptor"):
        contact_stability(pep_only)
    tcr_only = Structure("x", [Chain("B", [_res(0, "TYR", "Y", [_atom("CA", "C", [0, 0, 0])])],
                                     chain_type="TRB")])
    with pytest.raises(ValueError, match="no peptide"):
        contact_stability(tcr_only)
