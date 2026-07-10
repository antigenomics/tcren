"""Omega-loop geometry: the three-condition criterion, Frenet invariants, block layouts."""
import numpy as np
import pytest

from tcren.loops import (
    NECK_RANGE, block_layouts, find_junctions, frenet, is_omega_loop, kabsch_rmsd,
    omega_stats, structural_block_position,
)


def _helix(n, radius=2.3, rise=1.5, turn=100.0):
    """An ideal alpha-helical C-alpha trace: constant curvature and torsion."""
    t = np.arange(n) * np.deg2rad(turn)
    return np.stack([radius * np.cos(t), radius * np.sin(t), rise * np.arange(n)], axis=1)


def _closed_loop(n, radius=4.0):
    """A planar arc whose termini come back close together: an omega-like loop."""
    t = np.linspace(0, 1.75 * np.pi, n)
    return np.stack([radius * np.cos(t), radius * np.sin(t), np.zeros(n)], axis=1)


# ---------------------------------------------------------------- junction location

def test_find_junctions_locates_c_to_fgxg():
    seq = "GVTQTPKFQVLKTGQSMTLQCAQDMNHEYMSWYRQDPGMGLRLIHYSVGAGITDQGEVPNGYNVSRSTTEDFPLRLLSAAPSQTSVYFCASRPGLAGGRPEQYFGPGTRLTVT"
    ca = np.zeros((len(seq), 3))
    js = find_junctions(seq, ca)
    assert len(js) == 1
    j = js[0]
    assert j.seq == "CASRPGLAGGRPEQYF"       # AIRR junction: anchors included
    assert j.cdr3 == "ASRPGLAGGRPEQY"        # IMGT CDR3: two residues shorter
    assert seq[j.cys] == "C" and seq[j.fw] == "F"


def test_find_junctions_rejects_implausible_lengths():
    assert find_junctions("CFGAG", np.zeros((5, 3))) == []          # too short
    assert find_junctions("AAAAFGAG", np.zeros((8, 3))) == []       # no Cys
    with pytest.raises(ValueError):
        find_junctions("CASSFGQG", np.zeros((3, 3)))                # seq/coord mismatch


# ---------------------------------------------------------------- the three conditions

def test_omega_stats_reports_all_three_conditions():
    loop = _closed_loop(12)
    s = omega_stats(loop)
    assert s["n_residues"] == 12
    assert s["length_ok"] and s["termini_close"] and s["compact"]
    assert 0 < s["ratio"] < 2 / 3
    assert is_omega_loop(loop)


def test_an_extended_segment_is_not_an_omega_loop():
    """The compactness rule (d_end < 2/3 * d_max) is what rejects a straight chain: its
    termini ARE the most distant pair, so ratio == 1. This is exactly why CDR1 fails."""
    straight = np.stack([np.arange(8) * 3.8, np.zeros(8), np.zeros(8)], axis=1)
    s = omega_stats(straight)
    assert s["ratio"] == pytest.approx(1.0)
    assert s["length_ok"]              # 8 residues: passes (a)
    assert not s["termini_close"]      # 26.6 A: fails (b)
    assert not s["compact"]            # fails (c)
    assert not is_omega_loop(straight)


def test_compactness_alone_can_reject_a_short_close_segment():
    """A segment can be the right length AND have termini under 10 A apart, yet never
    reverse direction. Two conditions are not enough; the third is load-bearing. This is
    the shape CDR1 has -- its termini are the most distant pair (ratio == 1.00)."""
    # a shallow 6-residue arc: termini 8 A apart, and they are the farthest pair
    ca = np.array([[0, 0, 0], [1.6, 0.5, 0], [3.2, 0.8, 0], [4.8, 0.8, 0],
                   [6.4, 0.5, 0], [8.0, 0, 0]], dtype=float)
    s = omega_stats(ca)
    assert s["length_ok"]                        # passes (a): 6 residues
    assert s["termini_close"]                    # passes (b): 8 A < 10 A
    assert s["ratio"] == pytest.approx(1.0)      # ...but the termini ARE the longest span
    assert not s["compact"]                      # fails (c)
    assert not is_omega_loop(ca)


def test_relax_length_admits_long_cdr_h3():
    long_loop = _closed_loop(22)
    assert not is_omega_loop(long_loop)                      # 22 > 16
    assert is_omega_loop(long_loop, relax_length=True)


def test_neck_range_is_the_measured_band_not_folklore():
    assert NECK_RANGE == (5.0, 7.5)   # 4.5 A would select ~0% of real junctions


# ---------------------------------------------------------------- Frenet

def test_frenet_is_constant_along_an_ideal_helix():
    kappa, tau = frenet(_helix(12))
    assert len(kappa) == 10 and len(tau) == 9
    assert np.allclose(kappa, kappa[0], atol=1e-6)
    assert np.allclose(tau, tau[0], atol=1e-6)


def test_frenet_is_invariant_under_rigid_motion():
    """The whole point of the descriptor: no superposition needed."""
    ca = _helix(10)
    theta = 0.7
    rot = np.array([[np.cos(theta), -np.sin(theta), 0],
                    [np.sin(theta), np.cos(theta), 0],
                    [0, 0, 1]])
    moved = ca @ rot.T + np.array([13.0, -4.0, 2.5])
    k1, t1 = frenet(ca)
    k2, t2 = frenet(moved)
    assert np.allclose(k1, k2, atol=1e-6)
    assert np.allclose(t1, t2, atol=1e-6)


def test_frenet_torsion_flips_sign_under_reflection():
    ca = _helix(10)
    mirrored = ca * np.array([1.0, 1.0, -1.0])
    _, t1 = frenet(ca)
    _, t2 = frenet(mirrored)
    assert np.allclose(t1, -t2, atol=1e-6)   # chirality is captured


def test_frenet_marks_inflection_points_nan():
    """Torsion is undefined where consecutive tangents are collinear."""
    ca = np.array([[0, 0, 0], [1, 0, 0], [2, 0, 0], [3, 0, 0], [3, 1, 0]], dtype=float)
    _, tau = frenet(ca)
    assert np.isnan(tau).any()


def test_frenet_needs_four_points():
    with pytest.raises(ValueError):
        frenet(np.zeros((3, 3)))


# ---------------------------------------------------------------- block layouts

def test_block_layouts_match_seqtree_gapblock_indexing():
    # equal length: one layout, the identity correspondence
    assert block_layouts(4, 4)[0] == [(0, 0), (1, 1), (2, 2), (3, 3)]
    assert len(block_layouts(4, 4)) == 5   # positions 0..4, all giving the same pairing

    # query longer by one: block position i skips query residue i
    lay = block_layouts(4, 3)
    assert len(lay) == 4                          # i in [0, 3] inclusive at both ends
    assert lay[0] == [(1, 0), (2, 1), (3, 2)]     # leading block
    assert lay[3] == [(0, 0), (1, 1), (2, 2)]     # trailing block

    # ref longer by one: the offset moves to the ref index
    lay = block_layouts(3, 4)
    assert lay[0] == [(0, 1), (1, 2), (2, 3)]
    assert lay[3] == [(0, 0), (1, 1), (2, 2)]


def test_kabsch_rmsd_is_zero_after_rigid_motion():
    a = _closed_loop(9)
    theta = 1.1
    rot = np.array([[np.cos(theta), 0, -np.sin(theta)], [0, 1, 0], [np.sin(theta), 0, np.cos(theta)]])
    assert kabsch_rmsd(a, a @ rot.T + 7.0) == pytest.approx(0.0, abs=1e-8)
    with pytest.raises(ValueError):
        kabsch_rmsd(a, a[:-1])


def test_structural_block_position_recovers_a_planted_insertion():
    """Build a loop, insert one residue at a known position, and check the geometry finds it."""
    base = _closed_loop(11)
    insert_at = 5
    extra = (base[insert_at - 1] + base[insert_at]) / 2 + np.array([0.0, 0.0, 3.0])
    grown = np.insert(base, insert_at, extra, axis=0)
    best, rmsd, all_rmsd = structural_block_position(grown, base)
    assert len(all_rmsd) == len(base) + 1
    assert abs(best - insert_at) <= 1
    assert rmsd < 0.5
