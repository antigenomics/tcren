"""Unit tests for the TCRdock docking-geometry primitives (pure geometry, no structure/mmseqs)."""

from __future__ import annotations

import numpy as np
import pytest

from tcren.orient.tcrdock_geometry import (
    DockingConsensus,
    DockingGeometry,
    _dihedral,
    _encode,
    _kabsch,
    _rotation_axis,
    _stub_from_three_points,
    _symmetry_stub,
    consensus_distribution,
    docking_pvalue,
    docking_zscore,
)


def test_kabsch_recovers_known_rotation():
    rng = np.random.default_rng(0)
    mov = rng.normal(size=(8, 3))
    theta = 0.7
    R = np.array([[np.cos(theta), -np.sin(theta), 0], [np.sin(theta), np.cos(theta), 0], [0, 0, 1]])
    fix = (R @ mov.T).T + np.array([3.0, -1.0, 2.0])
    out = _kabsch(fix, mov)
    assert np.allclose(out, fix, atol=1e-8)


def test_rotation_axis_180_degrees():
    # 180° rotation about x maps frame I to diag(1,-1,-1); axis should be ±x
    axes1 = np.eye(3)
    axes2 = np.diag([1.0, -1.0, -1.0])
    n = _rotation_axis(axes1, axes2)
    assert abs(abs(n[0]) - 1.0) < 1e-6 and abs(n[1]) < 1e-6 and abs(n[2]) < 1e-6


def test_symmetry_stub_recovers_axis_and_is_orthonormal():
    # build two halves related by a 180° rotation about the x-axis
    rng = np.random.default_rng(1)
    a = rng.normal(size=(4, 3)) + np.array([0.0, 2.0, 1.0])
    b = (np.diag([1.0, -1.0, -1.0]) @ a.T).T          # Rx(180) of half A
    coords = np.vstack([a, b])
    stub = _symmetry_stub(coords, point_towards=np.array([5.0, 0.0, 0.0]))
    axes = stub["axes"]
    # x-axis parallel to the true symmetry axis, pointing toward point_towards (+x)
    assert axes[0][0] > 0.99
    # orthonormal frame
    assert np.allclose(axes @ axes.T, np.eye(3), atol=1e-6)
    # origin at the midpoint of the two half-COMs
    assert np.allclose(stub["origin"], 0.5 * (a.mean(0) + b.mean(0)), atol=1e-6)


def test_stub_from_three_points_orthonormal():
    axes, origin = _stub_from_three_points(np.array([1.0, 0, 0]), np.array([0.0, 0, 0]), np.array([0.0, 1, 0]))
    assert np.allclose(axes @ axes.T, np.eye(3), atol=1e-9)


@pytest.mark.parametrize("sign", [+1.0, -1.0])
def test_dihedral_right_angle(sign):
    # classic 90° dihedral setup
    p1 = np.array([1.0, 0.0, 0.0])
    p2 = np.array([0.0, 0.0, 0.0])
    p3 = np.array([0.0, 0.0, 1.0])
    p4 = np.array([0.0, sign * 1.0, 1.0])
    ang = np.degrees(_dihedral(p1, p2, p3, p4))
    assert abs(abs(ang) - 90.0) < 1e-6


def test_docking_geometry_dataclass_roundtrip():
    dg = DockingGeometry(29.1, 3.76, -0.01, 0.11, 0.24, 0.03)
    d = dg.to_dict()
    assert set(d) == {"d", "torsion", "tcr_unit_y", "tcr_unit_z", "mhc_unit_y", "mhc_unit_z"}
    assert d["d"] == pytest.approx(29.1)


# ---- docking-canonicality consensus (Mahalanobis Z-score) --------------------------------------------

def _synthetic_geoms(n, rng):
    """A tight cluster of plausible class-I docking geometries around a canonical centre."""
    out = []
    for _ in range(n):
        out.append(DockingGeometry(
            d=30.0 + rng.normal(0, 1.0),
            torsion=3.3 + rng.normal(0, 0.1),
            tcr_unit_y=0.05 + rng.normal(0, 0.05), tcr_unit_z=0.10 + rng.normal(0, 0.05),
            mhc_unit_y=0.20 + rng.normal(0, 0.05), mhc_unit_z=0.02 + rng.normal(0, 0.05)))
    return out


def test_encode_dimension_and_torsion_wrap():
    # torsion is sin-cos encoded, so 0 and 2*pi encode identically (no artificial cliff)
    a = _encode(DockingGeometry(30.0, 0.0, 0.1, 0.1, 0.1, 0.1))
    b = _encode(DockingGeometry(30.0, 2 * np.pi, 0.1, 0.1, 0.1, 0.1))
    assert a.shape == (7,)
    assert np.allclose(a, b, atol=1e-9)


def test_consensus_center_is_closer_than_outlier():
    rng = np.random.default_rng(3)
    cons = consensus_distribution(_synthetic_geoms(200, rng), mhc_class="I")
    center = DockingGeometry(30.0, 3.3, 0.05, 0.10, 0.20, 0.02)      # ~ the consensus mean
    outlier = DockingGeometry(45.0, 3.3, 0.05, 0.10, 0.20, 0.02)     # far in d (non-canonical)
    z_center = docking_zscore(center, cons)
    z_outlier = docking_zscore(outlier, cons)
    assert z_center >= 0.0
    assert z_outlier > z_center + 3.0                                # much less canonical
    # p-value is the mirror: canonical high, outlier ~0
    assert docking_pvalue(center, cons) > docking_pvalue(outlier, cons)
    assert 0.0 <= docking_pvalue(outlier, cons) < 0.01


def test_consensus_dict_roundtrip():
    rng = np.random.default_rng(4)
    cons = consensus_distribution(_synthetic_geoms(50, rng))
    back = DockingConsensus.from_dict(cons.to_dict())
    g = DockingGeometry(31.0, 3.2, 0.0, 0.1, 0.2, 0.0)
    assert docking_zscore(g, cons) == pytest.approx(docking_zscore(g, back), rel=1e-9)
    assert back.n == cons.n and back.mhc_class == cons.mhc_class


def test_consensus_needs_enough_points():
    rng = np.random.default_rng(5)
    with pytest.raises(ValueError):
        consensus_distribution(_synthetic_geoms(5, rng))            # <= 7 params -> underdetermined
