"""Unit tests for the PyMOL figure layer (:mod:`tcren.viz.pymol`).

The point of this module is that none of the *decisions* need PyMOL: the triad geometry, the
screen-space label placement and the corner compositing are ordinary Python, so they are tested
here without a renderer. Only the two tests at the bottom shell out to PyMOL, and they skip when
it is absent.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest

from tcren.viz.pymol import (
    CANONICAL_AXES,
    CORNERS,
    PALETTES,
    composite,
    gizmo_cgo,
    gizmo_scene,
    groove_scene,
    interface_scene,
    label_points,
    overlay_scene,
    render,
    VIEW_SIDE,
    VIEW_TOP,
)

REPO = Path(__file__).resolve().parents[2]
CANON = REPO / "data" / "Canonical2026"
_CYLINDER_LEN, _CONE_LEN = 14, 17
_STRIDE = _CYLINDER_LEN + _CONE_LEN


# --- the axis frame ---------------------------------------------------------------------------

def test_canonical_axes_name_the_frame_in_orient():
    """The labels must stay tied to what tcren.orient.frame actually does."""
    x, y, z = CANONICAL_AXES
    assert (x.letter, y.letter, z.letter) == ("x", "y", "z")
    assert "PC3" in x.definition and "PC2" in y.definition and "PC1" in z.definition
    assert "C-terminus" in y.definition, "y is signed toward the peptide C-terminus"
    assert "TCR" in z.definition and "−z" in z.definition, "the MHC sits at −z"
    # Short forms have to fit a small tile.
    assert all(len(a.short) <= 5 for a in CANONICAL_AXES)


# --- triad geometry ---------------------------------------------------------------------------

def test_each_arm_runs_along_its_own_axis_from_the_origin():
    cgo, tips = gizmo_cgo()
    assert len(cgo) == 3 * _STRIDE
    for i in range(3):
        start = cgo[i * _STRIDE + 1: i * _STRIDE + 4]
        end = cgo[i * _STRIDE + 4: i * _STRIDE + 7]
        assert start == [0.0, 0.0, 0.0]
        assert end[i] > 0
        assert all(abs(end[k]) < 1e-12 for k in range(3) if k != i)
        assert tips[i][i] > end[i], "the arrowhead extends past the shaft"


def test_head_is_wider_than_the_shaft():
    """Otherwise it renders as a spike, not an arrow."""
    cgo, _ = gizmo_cgo()
    assert cgo[_CYLINDER_LEN + 7] > cgo[7] * 2


def test_thin_by_default():
    """'Thin' is the whole brief: the shaft must be a small fraction of the arm."""
    cgo, tips = gizmo_cgo()
    assert cgo[7] / tips[0][0] < 0.05


@pytest.mark.parametrize("kwargs", [
    {"palette": "neon"}, {"arm": 0}, {"radius": -1}, {"head_length": 1.0}, {"head_radius": 0},
])
def test_gizmo_cgo_rejects_nonsense(kwargs):
    with pytest.raises(ValueError):
        gizmo_cgo(**kwargs)


def test_palettes_are_rgb_triples_in_unit_range():
    for name, colours in PALETTES.items():
        assert len(colours) == 3, name
        assert all(len(c) == 3 and all(0.0 <= v <= 1.0 for v in c) for c in colours), name


# --- label placement --------------------------------------------------------------------------

def test_labels_push_outward_past_their_tips():
    tips = gizmo_cgo()[1]
    pts = label_points([1, 0, 0, 0, 1, 0, 0, 0, 1], tips)
    assert pts[0][0] > tips[0][0] and pts[1][1] > tips[1][1]


def test_head_on_axis_label_falls_back_instead_of_landing_on_the_origin():
    """Under VIEW_SIDE the y axis points at the viewer and has no projected direction.

    Anchoring it at its tip would drop the label onto the origin, on top of the other two.
    """
    tips = gizmo_cgo()[1]
    pts = label_points(VIEW_SIDE, tips)
    assert pts[1][0] < 0 and pts[1][2] < 0, pts[1]
    # ... and it must not coincide with either of the others.
    for other in (pts[0], pts[2]):
        assert sum((pts[1][i] - other[i]) ** 2 for i in range(3)) > 1.0


def test_offset_controls_how_far_out_labels_sit():
    tips = gizmo_cgo()[1]
    near = label_points(VIEW_TOP, tips, offset=1.2)
    far = label_points(VIEW_TOP, tips, offset=3.0)

    def gap(p, t):
        return sum((p[i] - t[i]) ** 2 for i in range(3)) ** 0.5

    assert gap(far[0], tips[0]) > gap(near[0], tips[0])


# --- scene generation -------------------------------------------------------------------------

def test_gizmo_scene_names_the_axes_and_fixes_the_camera():
    scene = gizmo_scene(VIEW_SIDE)
    assert "load_cgo" in scene and "set_view" in scene and "zoom" in scene
    for axis in CANONICAL_AXES:
        assert f"label={axis.short!r}" in scene
    assert f"label={CANONICAL_AXES[2].label!r}" in gizmo_scene(VIEW_SIDE, short_labels=False)


def test_gizmo_scene_needs_a_full_rotation():
    with pytest.raises(ValueError):
        gizmo_scene([1, 0, 0])


@pytest.mark.parametrize("factory,args", [
    (overlay_scene, (["1ao7", "1bd2"], "/tmp")),
    (groove_scene, ("1ao7", "/tmp")),
    (interface_scene, ("1ao7", "/tmp", {"TRA": [26, 27], "TRB": [95]})),
])
def test_presets_load_view_and_zoom(factory, args):
    body = factory(*args)
    assert "cmd.load" in body and "cmd.set_view" in body and "cmd.zoom" in body


def test_overlay_caps_how_many_structures_it_draws():
    """Past a handful an overlay stops being readable, so the preset truncates."""
    body = overlay_scene([f"p{i}" for i in range(40)], "/tmp")
    assert body.count("cmd.load") == 8
    assert overlay_scene([f"p{i}" for i in range(40)], "/tmp", limit=3).count("cmd.load") == 3


def test_interface_scene_selects_the_given_cdr_residues():
    body = interface_scene("1ao7", "/tmp", {"TRA": [26, 27], "TRB": [95]})
    assert "chain A and resi 26+27" in body
    assert "chain B and resi 95" in body


def test_interface_scene_survives_a_chain_with_no_cdr_residues():
    """An unannotatable chain must not produce `resi ` and abort the whole render."""
    body = interface_scene("1ao7", "/tmp", {"TRA": [], "TRB": []})
    assert "resi 0" in body and "resi \n" not in body


# --- corner compositing -----------------------------------------------------------------------

def _ink_bbox(path):
    """Bounding box of the dark tile against the white base.

    Not alpha: the base is opaque, as a real render is, so every pixel would qualify.
    """
    from PIL import Image
    im = Image.open(path).convert("L")
    px = im.load()
    pts = [(x, y) for y in range(im.height) for x in range(im.width) if px[x, y] < 128]
    assert pts, f"no dark pixels in {path}"
    xs, ys = [p[0] for p in pts], [p[1] for p in pts]
    return min(xs), min(ys), max(xs), max(ys)


@pytest.fixture
def tiles(tmp_path):
    from PIL import Image
    base = tmp_path / "base.png"
    tile = tmp_path / "tile.png"
    Image.new("RGBA", (400, 300), (255, 255, 255, 255)).save(base)
    Image.new("RGBA", (60, 60), (0, 0, 0, 255)).save(tile)
    return base, tile


def test_every_corner_lands_inside_the_frame(tiles, tmp_path):
    base, tile = tiles
    from PIL import Image
    for name in CORNERS:
        out = tmp_path / f"{name}.png"
        composite(base, tile, out, corner=name, scale=0.2, margin=0.02)
        im = Image.open(out)
        x0, y0, x1, y1 = _ink_bbox(out)
        assert 0 <= x0 and x1 < im.width, name
        assert 0 <= y0 and y1 < im.height, name


def test_corners_sit_on_the_side_they_are_named_for(tiles, tmp_path):
    base, tile = tiles
    got = {}
    for name in CORNERS:
        out = tmp_path / f"{name}.png"
        composite(base, tile, out, corner=name, scale=0.2, margin=0.02)
        got[name] = _ink_bbox(out)
    assert got["bottom-left"][0] < got["bottom-right"][0]
    assert got["top-left"][0] < got["top-right"][0]
    # Image y runs downward, so "bottom" has the larger row index.
    assert got["bottom-left"][1] > got["top-left"][1]
    assert got["bottom-right"][1] > got["top-right"][1]


def test_scale_sets_the_tile_width(tiles, tmp_path):
    base, tile = tiles
    for scale in (0.1, 0.3):
        out = tmp_path / f"s{scale}.png"
        composite(base, tile, out, corner="bottom-left", scale=scale, margin=0.0)
        x0, _, x1, _ = _ink_bbox(out)
        assert abs((x1 - x0 + 1) - 400 * scale) <= 2, scale


@pytest.mark.parametrize("kwargs", [{"corner": "middle"}, {"scale": 0.0}, {"margin": 0.7}])
def test_composite_rejects_nonsense(tiles, tmp_path, kwargs):
    base, tile = tiles
    with pytest.raises(ValueError):
        composite(base, tile, tmp_path / "x.png", **kwargs)


# --- the renderer itself (needs PyMOL) ----------------------------------------------------------

pymol_only = pytest.mark.skipif(shutil.which("pymol") is None, reason="pymol not on PATH")


@pytest.mark.slow
@pymol_only
def test_render_writes_a_png_with_the_gizmo_in_the_corner(tmp_path):
    pid = "1ao7"
    if not (CANON / f"{pid}.pdb.gz").exists():
        pytest.skip("data/Canonical2026 not fetched; run `tcren fetch-data`")
    from PIL import Image

    plain = render(groove_scene(pid, CANON), tmp_path / "plain.png", size=(320, 320), gizmo=False)
    withax = render(groove_scene(pid, CANON), tmp_path / "ax.png", size=(320, 320),
                    corner="bottom-left")
    assert Image.open(plain).size == (320, 320) == Image.open(withax).size
    # The gizmo tile is a temporary and must not be left behind.
    assert not (tmp_path / "ax_gizmo.png").exists()

    # The corner the gizmo went into must differ from the same render without it; the opposite
    # corner must not. Compare in the bottom-left and top-right 25% squares.
    # Two independent ray traces of the same scene differ by a few grey levels of antialiasing
    # jitter, so the opposite corner is checked for "no ink was added" rather than for bit
    # equality -- a gizmo is worth tens of levels over an 80 px square, jitter is worth single
    # digits.
    import numpy as np
    a = np.asarray(Image.open(plain).convert("L"), dtype=float)
    b = np.asarray(Image.open(withax).convert("L"), dtype=float)
    n = 80
    here = np.abs(a[-n:, :n] - b[-n:, :n]).max()
    there = np.abs(a[:n, -n:] - b[:n, -n:]).max()
    assert here > 20, f"bottom-left should have gained the gizmo (got {here})"
    assert there < 12, f"top-right should carry no gizmo ink (got {there})"
    assert here > there * 3, (here, there)


@pytest.mark.slow
@pymol_only
def test_probe_rotation_reports_the_scene_rotation(tmp_path):
    from tcren.viz.pymol import probe_rotation
    rot = probe_rotation(f"cmd.set_view({list(VIEW_TOP) + [0, 0, -100, 0, 0, 0, 50, 150, -20]!r})")
    assert len(rot) == 9
    assert rot == pytest.approx(list(VIEW_TOP), abs=1e-5)
