"""Fast synthetic tests for the cross-map pose-consistency descriptors (tcren.pose)."""

from __future__ import annotations

import numpy as np
import pytest

from tcren.pose import (
    POSE_FEATURES,
    _double_centred,
    _pair_j,
    _selfcheck,
    _spearman,
    c_score,
    pose_consistency,
    pose_native_reference,
)
from tcren.structure.model import PEPTIDE_TYPE, Atom, Chain, Residue, Structure


def _atom(name, el, xyz):
    return Atom(name, el, np.asarray(xyz, float))


def _res(i, resname, aa, xyz, cb_offset=1.0):
    """A residue with a Cα at ``xyz`` and a Cβ ``cb_offset`` further along +x."""
    x, y, z = xyz
    return Residue(i, i + 1, "", aa, resname, (
        _atom("CA", "C", [x, y, z]), _atom("CB", "C", [x + cb_offset, y, z]),
    ))


def _complex(pep, tcr):
    """``pep``/``tcr``: list of ``(aa, xyz, cb_offset)``."""
    three = {"A": "ALA", "L": "LEU", "K": "LYS", "D": "ASP", "W": "TRP", "Y": "TYR"}
    pc = Chain("C", [_res(i, three[a], a, xyz, off) for i, (a, xyz, off) in enumerate(pep)],
               chain_type=PEPTIDE_TYPE)
    tc = Chain("B", [_res(i, three[a], a, xyz, off) for i, (a, xyz, off) in enumerate(tcr)],
               chain_type="TRB")
    return Structure("synth", [pc, tc])


class _Pot:
    """A directional 3-letter potential: (L,D) favourable, (L,K) neutral, (L,W) repulsive."""

    def as_matrix(self):
        idx = {"L": 0, "D": 1, "K": 2, "W": 3}
        m = np.array([
            [0.0, -3.0, 0.0, 3.0],
            [-3.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0],
            [3.0, 0.0, 0.0, 0.0],
        ])
        return m, idx


def _line(pep_specs):
    """One TCR Leu at the origin; peptide residues strung out along +x at the given distances."""
    tcr = [("L", (0.0, 0.0, 0.0), -1.0)]  # Cβ points back toward -x, i.e. away from the peptide
    pep = [(aa, (d, 0.0, 0.0), off) for aa, d, off in pep_specs]
    return _complex(pep, tcr)


def test_double_centring_zeroes_rows_and_columns_of_an_asymmetric_matrix():
    class _P:
        def as_matrix(self):
            return np.array([[1.0, 2.0, 9.0], [3.0, 0.0, 1.0], [5.0, 4.0, 2.0]]), {"A": 0, "B": 1, "C": 2}

    j, idx = _double_centred(_P())
    assert np.allclose(j.sum(axis=0), 0.0)
    assert np.allclose(j.sum(axis=1), 0.0)
    assert idx == {"A": 0, "B": 1, "C": 2}


def test_double_centring_tolerates_an_unobserved_cell():
    class _P:
        def as_matrix(self):
            return np.array([[1.0, np.nan], [3.0, 2.0]]), {"A": 0, "B": 1}

    j, _ = _double_centred(_P())
    assert np.isnan(j[0, 1]) and np.isfinite(j[1, 0])


def test_pair_j_is_nan_outside_the_alphabet():
    j, idx = _double_centred(_Pot())
    assert np.isnan(_pair_j(["L"], ["Z"], j, idx)[0])
    assert np.isfinite(_pair_j(["L"], ["D"], j, idx)[0])


def test_spearman_guards():
    assert np.isnan(_spearman(np.array([1.0, 2.0]), np.array([1.0, 2.0])))   # n < 3
    assert np.isnan(_spearman(np.ones(5), np.arange(5.0)))                   # constant axis
    assert _spearman(np.arange(5.0), np.arange(5.0)) == pytest.approx(1.0)


def test_c_local_is_positive_when_the_tight_contact_is_the_complementary_one():
    # favourable (L,D) tight at 3.0 A; neutral (L,K) at 4.0; repulsive (L,W) loose at 4.9
    s = _line([("D", 3.0, 1.0), ("K", 4.0, 1.0), ("W", 4.9, 1.0)])
    d = pose_consistency(s, potential=_Pot())
    assert d["n_contacts"] == 3
    assert d["c_local"] == pytest.approx(1.0)
    assert d["e_tight_minus_loose"] > 0
    assert d["margin_energy_slope"] > 0


def test_c_local_inverts_on_a_forced_pose():
    # the same chemistry with the distance order reversed: the repulsive pair is now the tight one
    s = _line([("D", 4.9, 1.0), ("K", 4.0, 1.0), ("W", 3.0, 1.0)])
    d = pose_consistency(s, potential=_Pot())
    assert d["c_local"] == pytest.approx(-1.0)
    assert d["e_tight_minus_loose"] < 0


def test_frac_close_favourable_counts_only_the_below_median_contacts():
    s = _line([("D", 3.0, 1.0), ("K", 4.0, 1.0), ("W", 4.9, 1.0)])
    d = pose_consistency(s, potential=_Pot())
    # median distance is 4.0; the two pairs at or below it are (L,D) J<0 and (L,K) J>0 after
    # centring, so exactly half of the close contacts are favourable
    assert 0.0 <= d["frac_close_favourable"] <= 1.0


def test_sidechain_toward_is_positive_when_cbeta_points_at_the_partner():
    toward = pose_consistency(_line([("D", 4.0, -1.0)]), potential=_Pot())["sidechain_toward"]
    away = pose_consistency(_line([("D", 4.0, +1.0)]), potential=_Pot())["sidechain_toward"]
    # the peptide Cβ at -1.0 leans back toward the TCR, so Cβ-Cβ closes relative to Cα-Cα
    assert toward > away


def test_descriptors_are_nan_not_zero_when_there_is_no_interface():
    s = _complex([("D", (60.0, 0.0, 0.0), 1.0)], [("L", (0.0, 0.0, 0.0), 1.0)])
    d = pose_consistency(s, potential=_Pot())
    assert d["n_contacts"] == 0
    assert all(np.isnan(d[k]) for k in ("c_local", "e_tight_minus_loose", "margin_energy_slope"))


def test_too_few_contacts_gives_nan_rather_than_a_made_up_correlation():
    s = _line([("D", 3.0, 1.0), ("W", 4.0, 1.0)])  # 2 contacts, below the 3 a rho needs
    d = pose_consistency(s, potential=_Pot())
    assert d["n_contacts"] == 2 and np.isnan(d["c_local"])


def test_bundled_reference_is_loadable_and_covers_every_feature():
    ref = pose_native_reference()
    for f in POSE_FEATURES:
        assert f in ref and len(ref[f]) > 100
        assert np.isfinite(ref[f]).all()


def test_c_score_is_defined_for_a_single_row_and_orders_correctly():
    ref = pose_native_reference()
    median = {k: [float(np.median(ref[k]))] for k in POSE_FEATURES}
    poor = {k: [float(np.median(ref[k]) - 2 * np.std(ref[k]))] for k in POSE_FEATURES}
    assert np.isfinite(c_score(median)[0])                 # one structure, no cohort
    assert c_score(median)[0] > c_score(poor)[0]           # higher = more crystal-like


def test_selfcheck_runs():
    _selfcheck()
