"""Omega-loop geometry: the three-condition criterion, Frenet invariants, block layouts."""
import numpy as np
import pytest

from tcren.loops import (
    NECK_RANGE, block_layouts, cb_orientation, find_junctions, frenet, frenet_frame, gap_runs,
    is_omega_loop, is_single_block, kabsch_rmsd, omega_stats, ramachandran, structural_align,
    structural_block_position, virtual_cb,
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


# ---------------------------------------------------------------- model-independent alignment

def test_gap_runs_and_single_block_predicate():
    assert gap_runs("MMMM") == []
    assert gap_runs("MMDDMM") == [("D", 2, 2)]
    assert gap_runs("MMDMMIM") == [("D", 2, 1), ("I", 5, 1)]
    assert is_single_block("MMMM")            # gapless
    assert is_single_block("MMDDMM")          # one block
    assert not is_single_block("MMDMMIM")     # two blocks -- what a single-block model cannot express
    assert not is_single_block("MDMMMDM")     # same sequence, two blocks


def test_structural_align_is_the_identity_on_a_rigidly_moved_copy():
    a = _closed_loop(12)
    theta = 0.9
    rot = np.array([[np.cos(theta), -np.sin(theta), 0],
                    [np.sin(theta), np.cos(theta), 0], [0, 0, 1]])
    pairs, rmsd, ops = structural_align(a, a @ rot.T + np.array([5.0, -3.0, 2.0]))
    assert pairs == [(i, i) for i in range(12)]
    assert ops == "M" * 12
    assert rmsd == pytest.approx(0.0, abs=1e-6)
    assert is_single_block(ops)


def test_structural_align_recovers_a_planted_insertion_as_one_block():
    base = _closed_loop(12)
    insert_at = 5
    extra = (base[insert_at - 1] + base[insert_at]) / 2 + np.array([0.0, 0.0, 3.0])
    grown = np.insert(base, insert_at, extra, axis=0)
    pairs, rmsd, ops = structural_align(grown, base)
    assert is_single_block(ops), ops
    (op, start, length), = gap_runs(ops)
    assert op == "D" and length == 1          # the extra residue of `grown` is unmatched
    assert abs(start - insert_at) <= 1
    assert rmsd < 0.5


def test_structural_align_can_return_two_blocks():
    """The oracle must be able to disagree with the single-block model, or it proves nothing."""
    base = _closed_loop(14)
    # remove one residue near each end: the true correspondence needs TWO gap blocks
    trimmed = np.delete(base, [4, 9], axis=0)
    _, _, ops = structural_align(base, trimmed)
    assert not is_single_block(ops), f"expected two blocks, got {ops}"
    assert len(gap_runs(ops)) == 2


def test_structural_align_needs_six_residues():
    with pytest.raises(ValueError):
        structural_align(np.zeros((5, 3)), np.zeros((8, 3)))


# ---------------------------------------------------------------- Cbeta on the Frenet frame

def _rot(theta, axis=2):
    c, s = np.cos(theta), np.sin(theta)
    m = np.eye(3)
    i, j = [(1, 2), (0, 2), (0, 1)][axis]
    m[i, i] = m[j, j] = c
    m[i, j], m[j, i] = -s, s
    return m


def test_frenet_frame_is_orthonormal_and_right_handed():
    t, n, b = frenet_frame(_helix(10))
    assert t.shape == n.shape == b.shape == (8, 3)
    for v in (t, n, b):
        assert np.allclose(np.linalg.norm(v, axis=1), 1.0)
    assert np.allclose(np.einsum("ij,ij->i", t, n), 0, atol=1e-9)
    assert np.allclose(np.einsum("ij,ij->i", t, b), 0, atol=1e-9)
    assert np.allclose(np.einsum("ij,ij->i", n, b), 0, atol=1e-9)
    assert np.allclose(np.cross(n, b), t, atol=1e-9)   # right-handed (t, n, b)


def test_cb_orientation_is_rigid_motion_invariant():
    """The point of expressing the side chain in the loop's own frame: no superposition needed."""
    ca = _helix(10)
    cb = ca + np.array([0.6, 0.9, 1.2])          # a fixed offset stands in for a side chain
    p1, a1 = cb_orientation(ca, cb)
    rot = _rot(0.8)
    p2, a2 = cb_orientation(ca @ rot.T + 11.0, cb @ rot.T + 11.0)
    assert np.allclose(p1, p2, atol=1e-8)
    assert np.allclose(a1, a2, atol=1e-8)


def test_cb_azimuth_flips_under_reflection():
    """Handedness lives in the azimuth, just as it lives in the Frenet torsion."""
    ca = _helix(10)
    cb = ca + np.array([0.6, 0.9, 1.2])
    mirror = np.array([1.0, 1.0, -1.0])
    _, a1 = cb_orientation(ca, cb)
    _, a2 = cb_orientation(ca * mirror, cb * mirror)
    assert np.allclose(a1, -a2, atol=1e-8)


def test_cb_orientation_rejects_mismatched_shapes():
    with pytest.raises(ValueError):
        cb_orientation(_helix(8), _helix(7))


def test_virtual_cb_reproduces_ideal_tetrahedral_geometry():
    """Bond length ~1.53 A and both N-CA-CB and C-CA-CB near the tetrahedral 110 degrees."""
    n = np.array([[-1.458, 0.0, 0.0]])
    ca = np.array([[0.0, 0.0, 0.0]])
    c = np.array([[0.55, 1.42, 0.0]])
    cb = virtual_cb(n, ca, c)
    assert np.linalg.norm(cb - ca) == pytest.approx(1.53, abs=0.05)

    def ang(u, v):
        u, v = u / np.linalg.norm(u), v / np.linalg.norm(v)
        return np.degrees(np.arccos(np.clip(u @ v, -1, 1)))

    assert ang((n - ca)[0], (cb - ca)[0]) == pytest.approx(110.5, abs=2.0)
    assert ang((c - ca)[0], (cb - ca)[0]) == pytest.approx(110.5, abs=2.0)
    # CB is well off the N-CA-C plane, on the side the construction fixes. Whether that side is
    # the one real L-amino acids use is checked against observed CB in scripts/cb_contacts.py.
    normal = np.cross((n - ca)[0], (c - ca)[0])
    assert abs(normal @ (cb - ca)[0]) > 1.0


def test_ramachandran_matches_a_hand_computed_dihedral():
    """phi_1 = dihedral(C0, N1, CA1, C1). Chosen so b0, b1, b2 are the three unit axes: -90 deg."""
    nn = np.array([[2.0, 0.0, 0.0], [0.0, 0.0, 0.0], [1.0, 1.0, 1.0]])
    ca = np.array([[2.0, 1.0, 0.0], [0.0, 1.0, 0.0], [1.0, 2.0, 1.0]])
    c = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 1.0], [1.0, 2.0, 2.0]])
    phi, psi = ramachandran(nn, ca, c)
    assert len(phi) == len(psi) == 1
    assert phi[0] == pytest.approx(-90.0, abs=1e-6)


def test_ramachandran_is_rigid_invariant_and_flips_under_reflection():
    rng = np.random.default_rng(0)
    ca = np.cumsum(rng.normal(size=(9, 3)), axis=0)
    nn = ca + rng.normal(scale=0.5, size=(9, 3))
    c = ca + rng.normal(scale=0.5, size=(9, 3))
    p1, s1 = ramachandran(nn, ca, c)

    rot, shift = _rot(0.6, axis=1), np.array([2.0, -1.0, 4.0])
    p2, s2 = ramachandran(nn @ rot.T + shift, ca @ rot.T + shift, c @ rot.T + shift)
    assert np.allclose(p1, p2, atol=1e-8) and np.allclose(s1, s2, atol=1e-8)

    m = np.array([1.0, 1.0, -1.0])            # dihedrals are chiral: reflection negates them
    p3, s3 = ramachandran(nn * m, ca * m, c * m)
    assert np.allclose(p1, -p3, atol=1e-8) and np.allclose(s1, -s3, atol=1e-8)

    with pytest.raises(ValueError):
        ramachandran(nn[:3], ca, c)
