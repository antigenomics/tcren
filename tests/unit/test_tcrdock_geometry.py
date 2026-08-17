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


# --- end-to-end: the class-I and class-II grooves -------------------------------------------------------
# The MHC core is chain-typed by hand here so these run without mmseqs; only the TCR framework markup
# (arda) is imported. The class-II halves live on separate chains (α1 on MHCa, β1 on MHCb) but resolve to
# the same six strand offsets as class I -- see tcrdock_geometry.CORE_OFFSETS_0X.
_MHC_TYPING = {
    "1ao7": {"A": ("MHCa", "MHCI"), "B": ("B2M", "MHCI"), "C": ("PEPTIDE", None),
             "D": ("TRA", None), "E": ("TRB", None)},
    "4ozg": {"A": ("MHCa", "MHCII"), "B": ("MHCb", "MHCII"), "J": ("PEPTIDE", None),
             "G": ("TRA", None), "H": ("TRB", None)},
    "6v0y": {"A": ("MHCa", "MHCII"), "B": ("MHCb", "MHCII"), "C": ("PEPTIDE", None),
             "D": ("TRA", None), "E": ("TRB", None)},
}


def _typed(pdb_id):
    from pathlib import Path

    from tcren.annotation import annotate_tcr_chains
    from tcren.structure import parse_structure

    pdb_dir = Path(__file__).resolve().parents[1] / "assets" / "pdb"
    s = parse_structure(pdb_dir / f"{pdb_id}.pdb")
    for c in s.chains:
        if c.chain_id in _MHC_TYPING[pdb_id]:
            c.chain_type, c.chain_supertype = _MHC_TYPING[pdb_id][c.chain_id]
    annotate_tcr_chains(s)
    return s


def test_class_i_core_keeps_all_six_pairs():
    from tcren.orient.tcrdock_geometry import _mhc_core_ca

    pytest.importorskip("arda")
    core = _mhc_core_ca(_typed("1ao7"))
    assert core.shape == (12, 3)          # 6 pairs, unchanged from the pre-class-II behaviour


@pytest.mark.parametrize("pdb_id", ["4ozg", "6v0y"])
def test_class_ii_core_pairs_both_chains(pdb_id):
    from tcren.orient.tcrdock_geometry import _MIN_CORE_PAIRS, _mhc_core_ca

    pytest.importorskip("arda")
    core = _mhc_core_ca(_typed(pdb_id))
    assert core is not None, "class-II β-sheet core not located"
    assert core.shape[0] >= 2 * _MIN_CORE_PAIRS and core.shape[0] % 2 == 0


@pytest.mark.slow
@pytest.mark.parametrize("pdb_id", ["1ao7", "4ozg", "6v0y"])
def test_docking_geometry_runs_on_both_mhc_classes(pdb_id):
    from tcren.orient.tcrdock_geometry import docking_geometry

    pytest.importorskip("arda")
    g = docking_geometry(_typed(pdb_id))
    # Class II must land in the same physical range as class I, not merely "not raise".
    assert 20.0 < g.d < 45.0
    assert 0.0 <= g.torsion < 2 * np.pi
    for u in (g.tcr_unit_y, g.tcr_unit_z, g.mhc_unit_y, g.mhc_unit_z):
        assert -1.0 <= u <= 1.0


@pytest.mark.slow
def test_class_i_docking_geometry_is_unchanged():
    from tcren.orient.tcrdock_geometry import docking_geometry

    pytest.importorskip("arda")
    g = docking_geometry(_typed("1ao7"))
    # Pinned against the pre-class-II implementation (verified bit-identical, not merely close).
    assert g.d == pytest.approx(31.301818408808607, abs=1e-9)
    assert np.degrees(g.torsion) == pytest.approx(192.05654000750124, abs=1e-9)
