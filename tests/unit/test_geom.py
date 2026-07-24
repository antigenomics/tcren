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


def test_interface_clashes_overlap_and_tolerance():
    a = np.array([[0.0, 0.0, 0.0]])           # radius 1.7
    b = np.array([[2.0, 0.0, 0.0],            # d=2.0 → overlap 1.4 (clash)
                  [3.3, 0.0, 0.0],            # d=3.3 → overlap 0.1 (below tol 0.4)
                  [50.0, 0.0, 0.0]])          # far, no overlap
    ra = np.array([1.7]); rb = np.array([1.7, 1.7, 1.7])
    d = _geom.interface_clashes(a, ra, b, rb, 0.4)
    assert list(d["i"]) == [0] and list(d["j"]) == [0]        # only the 2.0 Å pair clashes
    assert d["overlap"][0] == pytest.approx(1.4)
    # loosen the tolerance below 0.1 → the 3.3 Å pair now also clashes
    d2 = _geom.interface_clashes(a, ra, b, rb, 0.05)
    assert sorted(d2["j"]) == [0, 1]


def test_contact_stability_margins_and_fragility():
    tcr = np.array([[0.0, 0.0, 0.0]])
    pep = np.array([[3.0, 0.0, 0.0],   # margin 2.0 (robust)
                    [4.5, 0.0, 0.0],   # margin 0.5 (fragile)
                    [6.0, 0.0, 0.0]])  # beyond 5 Å cutoff → not a contact
    st = _geom.contact_stability(pep, np.array([0, 1, 2], np.int32),
                                 tcr, np.array([0], np.int32), 5.0, 1.0)
    assert st["n5"] == 2
    assert st["mean_margin"] == pytest.approx(1.25)          # (2.0 + 0.5) / 2
    assert st["frac_robust"] == pytest.approx(0.5)
    assert st["frac_marg_lt1"] == pytest.approx(0.5)
    assert st["exp_lost"] == pytest.approx(0.25)             # 0 + clip((1-0.5)/2)
