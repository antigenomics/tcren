"""Fast unit tests for the native interface-geometry kernel (tcren._geom)."""

from __future__ import annotations

import numpy as np
import pytest

_geom = pytest.importorskip("tcren._geom")

R, PROBE = 1.7, 1.4
_ISO = 4 * np.pi * (R + PROBE) ** 2  # SASA of a lone atom


def test_sasa_isolated_atom():
    s = _geom.shrake_rupley(np.array([[0.0, 0.0, 0.0]]), np.array([R]), PROBE, 960)
    assert abs(s[0] - _ISO) / _ISO < 0.01  # lone atom → full sphere


def test_sasa_far_atoms_additive():
    s = _geom.shrake_rupley(np.array([[0.0, 0.0, 0.0], [50.0, 0.0, 0.0]]),
                            np.array([R, R]), PROBE, 960)
    assert abs(s.sum() - 2 * _ISO) / (2 * _ISO) < 0.01


def test_sasa_close_atoms_buried():
    s = _geom.shrake_rupley(np.array([[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]]),
                            np.array([R, R]), PROBE, 960)
    assert s.sum() < 2 * _ISO  # mutual burial reduces accessible area


def test_sasa_converges_in_n_points():
    lo = _geom.shrake_rupley(np.array([[0.0, 0.0, 0.0]]), np.array([R]), PROBE, 100)[0]
    hi = _geom.shrake_rupley(np.array([[0.0, 0.0, 0.0]]), np.array([R]), PROBE, 3000)[0]
    assert abs(lo - _ISO) / _ISO < 0.05 and abs(hi - _ISO) / _ISO < 0.005


def test_interface_hbonds_counts_pairs_within_cutoff():
    donors = np.array([[0.0, 0.0, 0.0], [10.0, 0.0, 0.0]])
    acceptors = np.array([[3.0, 0.0, 0.0], [3.4, 0.0, 0.0]])
    assert _geom.interface_hbonds(donors, acceptors, 3.5) == 2  # both within 3.5 of donor 0
    assert _geom.interface_hbonds(donors, acceptors, 2.0) == 0


def test_contact_descriptors_size_and_balance():
    tra = np.array([[0.0, 0.0, 0.0]])
    trb = np.array([[8.0, 0.0, 0.0]])
    pep = np.array([[3.0, 0.0, 0.0], [5.0, 0.0, 0.0]])  # atom0 near TRA, atom1 near TRB
    mhc = np.zeros((0, 3))
    cd = _geom.contact_descriptors(tra, np.array([0], np.int32), trb, np.array([0], np.int32),
                                   pep, mhc, 5.0, 4.5)
    assert cd["pm_cov_ntcr"] == 2  # both TCR residues engage the peptide
    assert cd["n_pep_near_tra"] == 1 and cd["n_pep_near_trb"] == 1
    assert cd["chain_balance"] == pytest.approx(0.5)  # min(1,1)/(1+1)
