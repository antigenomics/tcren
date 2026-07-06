"""Unit tests for the TCRdock docking-geometry primitives (pure geometry, no structure/mmseqs)."""

from __future__ import annotations

import numpy as np
import pytest

from tcren.orient.tcrdock_geometry import (
    DockingGeometry,
    _dihedral,
    _kabsch,
    _rotation_axis,
    _stub_from_three_points,
    _symmetry_stub,
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
