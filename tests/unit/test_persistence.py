"""Persistence agrees with the existing MST entropy and resolves two loops."""
import builtins

import numpy as np
import pytest

from tcren.topology.persistence import _persistence_summary
from tcren.topology.footprint import _h0_persistence_entropy


def test_h0_mst_and_two_circles():
    pytest.importorskip("gudhi")
    theta = np.arange(16)*2*np.pi/16
    circle = np.column_stack([np.cos(theta), np.sin(theta), np.zeros(16)])
    points = np.concatenate([circle, circle+[6, 0, 0]])
    row = _persistence_summary(points, 2.)
    assert row['pers1_bars'] == 2
    assert row['pers0_infinite'] == 2
    row = _persistence_summary(points, 12.)
    assert row['pers0_entropy'] == pytest.approx(_h0_persistence_entropy(points), abs=1e-12)
    assert row['pers0_bars'] == 31
    assert row['pers0_infinite'] == 1


def test_missing_gudhi_fails_explicitly(monkeypatch):
    original = builtins.__import__
    def guarded(name, *args, **kwargs):
        if name == 'gudhi':
            raise ImportError('missing gudhi')
        return original(name, *args, **kwargs)
    monkeypatch.setattr(builtins, '__import__', guarded)
    with pytest.raises(ImportError, match='topology'):
        _persistence_summary(np.zeros((3, 3)), 12.)


def test_empty_and_invalid_inputs():
    pytest.importorskip("gudhi")
    assert _persistence_summary(np.empty((0, 3)), 12.)['pers1_bars'] == 0
    with pytest.raises(ValueError):
        _persistence_summary(np.zeros((3, 3)), float('nan'))


def test_ca_cb_missing_atoms_and_glycine_only_fallback(monkeypatch):
    pytest.importorskip("gudhi")
    import tcren.topology.persistence as module
    from tcren.structure.model import Atom, Chain, Residue, Structure
    residues = [Residue(0, 1, '', 'G', 'GLY', (Atom('CA', 'C', np.array([0., 0., 0.])),)),
                Residue(1, 2, '', 'A', 'ALA', (Atom('CA', 'C', np.array([1., 0., 0.])),)),
                Residue(2, 3, '', 'A', 'ALA', (Atom('CA', 'C', np.array([2., 0., 0.])),
                                               Atom('CB', 'C', np.array([2., 1., 0.]))))]
    structure = Structure('atoms', [Chain('A', residues)])
    parts = {'tcr_peptide': {(('A', 0), ('A', 1)), (('A', 1), ('A', 2))}, 'tcr_mhc': set(), 'peptide_mhc': set()}
    monkeypatch.setattr(module, '_recognition_contact_edges', lambda cm: parts)
    row = module.persistence_features(structure, contact_map=object())
    assert row['points'] == 3
    assert row['points_cb'] == 2
    assert row['missing_atoms_cb'] == 1
    assert row['gly_ca_fallback_cb'] == 1
    cb = module.persistence_features(structure, atom='CB', contact_map=object())
    assert cb['points'] == 2
    assert cb['points_ca'] == 3
