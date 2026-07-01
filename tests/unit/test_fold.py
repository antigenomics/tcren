"""Fast tests for the open-source fold layer: anchors, the CCD kernel, engines, and end-to-end."""

from __future__ import annotations

import numpy as np
import pytest

from tcren.refine import anchors as A

_fold = pytest.importorskip("tcren._fold")


# --- anchor prediction (pure stdlib, no deps) ---------------------------------------------------

def test_class1_anchors_p2_pomega():
    d = A.decompose("LLFGYPVYV")  # 9-mer -> class I
    assert d.mhc_class == "MHCI"
    assert d.anchors == (1, 8)  # 0-based P2 and PΩ
    assert d.presentation == "XLXXXXXXV"
    assert d.tcr_facing == "LXFGYPVYX"


def test_class1_anchor_indices_scale_with_length():
    assert A.anchor_indices("SIINFEKL", "MHCI") == (1, 7)  # 8-mer
    assert A.anchor_indices("GLCTLVAML", "MHCI") == (1, 8)  # 9-mer


def test_class2_register_inference():
    # Long peptide -> class II; a P1-hydrophobic register is preferred.
    d = A.decompose("PKYVKQNTLKLAT")
    assert d.mhc_class == "MHCII"
    assert len(d.anchors) == 4  # P1/P4/P6/P9 of the chosen core
    s = d.anchors[0]
    assert d.anchors == (s, s + 3, s + 5, s + 8)


def test_class2_short_peptide_no_core():
    assert A.anchor_indices("SHORT", "MHCII") == ()  # < 9 residues


# --- CCD kernel (C++) ---------------------------------------------------------------------------

def _bent_chain(n=10, step=3.8, seed=0):
    rng = np.random.default_rng(seed)
    v = np.array([1.0, 0.3, 0.0]); v /= np.linalg.norm(v)
    pts = [np.zeros(3)]
    for i in range(n - 1):
        ang = 0.7 * (-1) ** i + 0.2 * rng.standard_normal()
        c, s = np.cos(ang), np.sin(ang)
        v = np.array([c * v[0] - s * v[1], s * v[0] + c * v[1], 0.1 * rng.standard_normal()])
        v /= np.linalg.norm(v)
        pts.append(pts[-1] + step * v)
    return np.asarray(pts)


def test_ccd_reaches_reachable_target_and_preserves_bond_lengths():
    coords = _bent_chain()
    n = len(coords)
    bonds = np.array([[i, i + 1] for i in range(n - 1)], dtype=np.int32)
    moving = np.array([n - 1], dtype=np.int32)
    weights = np.array([1.0])

    # A target the tip can reach by rotating bond 0 (so it is exactly reachable).
    o = coords[0]; ax = coords[1] - o; ax /= np.linalg.norm(ax)
    th = 0.6; c, s = np.cos(th), np.sin(th)
    v = coords[n - 1] - o
    tip = (o + v * c + np.cross(ax, v) * s + ax * np.dot(ax, v) * (1 - c)).reshape(1, 3)

    closed, rmsd, iters = _fold.ccd_close(
        np.ascontiguousarray(coords), bonds, moving, np.ascontiguousarray(tip), weights, 2000, 1e-3)
    closed = np.asarray(closed)
    assert np.linalg.norm(closed[n - 1] - tip[0]) < 1e-2  # reached the target
    seg = np.linalg.norm(np.diff(closed, axis=0), axis=1)
    assert np.allclose(seg, 3.8, atol=1e-6)  # rigid linkage: bond lengths preserved


def test_ccd_rejects_out_of_range_indices():
    coords = _bent_chain(n=6)
    bonds = np.array([[0, 1]], dtype=np.int32)
    weights = np.array([1.0])
    target = np.array([[0.0, 0.0, 0.0]])
    # moving index past the end of coords -> ValueError, not a segfault/UB.
    with pytest.raises(ValueError):
        _fold.ccd_close(np.ascontiguousarray(coords), bonds,
                        np.array([99], dtype=np.int32), target, weights, 10, 1e-3)
    # targets shorter than moving -> ValueError.
    with pytest.raises(ValueError):
        _fold.ccd_close(np.ascontiguousarray(coords), bonds,
                        np.array([2, 3], dtype=np.int32), target, np.ones(2), 10, 1e-3)


def test_ccd_monotone_decrease_multi_anchor():
    coords = _bent_chain(n=12)
    n = len(coords)
    bonds = np.array([[i, i + 1] for i in range(n - 1)], dtype=np.int32)
    moving = np.array([5, n - 1], dtype=np.int32)
    targets = coords[moving] + np.array([[0.5, 0.5, 0.0], [0.5, -0.5, 0.3]])
    weights = np.ones(2)

    def rmsd0():
        d = coords[moving] - targets
        return float(np.sqrt((d ** 2).sum(1).mean()))

    _closed, rmsd, _it = _fold.ccd_close(
        np.ascontiguousarray(coords), bonds, moving, np.ascontiguousarray(targets), weights, 1000, 1e-3)
    assert rmsd <= rmsd0() + 1e-9  # CCD never increases the anchor RMSD


# --- engine registry ----------------------------------------------------------------------------

def test_engine_registry():
    from tcren.refine.engines import ENGINES, available_engines, get_engine

    assert set(ENGINES) == {"dope", "ccd", "openmm", "promod3"}
    assert "ccd" in available_engines()  # _fold built (importorskip above)
    with pytest.raises(KeyError):
        get_engine("nope")


def test_optional_engines_unavailable_raise_not_import_error():
    from tcren.refine.engines import EngineUnavailable, OpenMMEngine, ProMod3Engine

    for eng in (OpenMMEngine(), ProMod3Engine()):
        if not eng.available():
            with pytest.raises(EngineUnavailable):
                eng.run(None, A.decompose("LLFGYPVYV"))


# --- end-to-end (needs arda/mmseqs to annotate; slow) -------------------------------------------

@pytest.mark.slow
def test_model_peptide_dope_recovers_native_pose():
    pytest.importorskip("arda")
    from tcren.annotation import classify_chains
    from tcren.mhc import annotate_mhc
    from tcren.paths import reference_structure_path
    from tcren.refine import model_peptide, peptide_rmsd
    from tcren.structure import parse_structure

    s = parse_structure(reference_structure_path("1ao7"), pdb_id="1ao7")
    classify_chains(s, organism="human")
    annotate_mhc(s)

    res = model_peptide(s, engine="dope", seed=1, n_steps=1500)
    rm = peptide_rmsd(res.structure, s, anchors=res.anchors)
    assert rm.backbone_rmsd < 1.0  # native pose is a DOPE optimum
