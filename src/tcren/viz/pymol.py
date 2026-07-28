"""Publication-ready PyMOL renders of canonically-oriented TCR–pMHC complexes.

A render of an oriented structure is only interpretable if the reader can tell which way the
canonical frame points, so the centrepiece here is the **axis gizmo** — a thin, arrow-headed triad
in a corner of the image, turning with the camera, naming what each direction means.

**The frame it draws.** :mod:`tcren.orient.frame` puts every structure into one frame by PCA of the
reference complex, and :data:`CANONICAL_AXES` is that frame written out in words:

=====  =========================  ================  ============================
axis   definition in code         figure label      equivalent in the literature
=====  =========================  ================  ============================
``x``  PC3, the thin axis         ``groove width``  ``PC2_MHC``, groove width
``y``  PC2, +y to peptide C-term  ``peptide N→C``   ``PC1_MHC``, groove long axis
``z``  PC1, +z toward the TCR     ``pMHC→TCR``      ``PC3_MHC``, groove normal
=====  =========================  ================  ============================

The principal-component *numbers* differ from the docking-geometry literature (SwiftTCR, TCR3d)
because those fit the MHC groove alone while :mod:`tcren.orient.frame` fits the whole complex, in
which the MHC→TCR direction carries the most variance. The three directions are the same three;
only their ranking differs. Naming them for what they are is the point — ``pMHC→TCR`` is readable
in a figure, ``z`` is not.

**How the gizmo is placed.** In its own render pass, not by projecting a world position into the
corner. PyMOL's orthoscopic viewport does not span the world height that ``field_of_view`` and the
camera distance imply — measured on a real scene it is out by about a quarter — so a gizmo placed
by that arithmetic lands off-frame. Rendering the triad alone on a transparent background and
compositing it at pixel coordinates makes its size and position exact by construction, and leaves
nothing to reverse-engineer. It also guarantees the molecule can never occlude it.

**Why the geometry lives here.** PyMOL runs under its own interpreter and cannot import ``tcren``,
so the scripts this module emits are lists of literals and every decision behind them is ordinary
testable Python.

Example:
    >>> from tcren.viz.pymol import render, overlay_scene       # doctest: +SKIP
    >>> render(overlay_scene(["1ao7"], "data/Canonical2026"), "overlay.png")   # doctest: +SKIP
"""

from __future__ import annotations

import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

__all__ = [
    "Axis", "CANONICAL_AXES", "CORNERS", "PALETTES", "CHAIN_COLOURS",
    "gizmo_cgo", "label_points", "gizmo_scene", "probe_rotation", "render", "composite",
    "overlay_scene", "groove_scene", "interface_scene",
]


@dataclass(frozen=True, slots=True)
class Axis:
    """One canonical axis: its letter, what it means, and how to label it in a figure."""

    letter: str
    label: str
    short: str
    definition: str


#: The canonical frame of :mod:`tcren.orient.frame`, in figure-ready words. ``short`` is for the
#: corner gizmo, where anything longer than a few characters collides with the arrows.
CANONICAL_AXES: tuple[Axis, Axis, Axis] = (
    Axis("x", "groove width", "width",
         "PC3, the thin axis; across the cleft, α1↔α2 helix separation"),
    Axis("y", "peptide N→C", "N→C",
         "PC2, the groove/peptide axis, signed toward the peptide C-terminus"),
    Axis("z", "pMHC→TCR", "TCR",
         "PC1, the MHC→TCR long axis, signed toward the TCR; the MHC sits at −z"),
)

#: Corner anchors, as the sign of (x, y) in image space with y up.
CORNERS: dict[str, tuple[int, int]] = {
    "bottom-left": (-1, -1), "bottom-right": (1, -1),
    "top-left": (-1, 1), "top-right": (1, 1),
}

#: ``mono`` — one restrained grey for all three arrows, told apart by their labels. The default,
#: because the structure already spends colour on chains and an orientation gizmo that competes
#: with the molecule makes a worse figure. ``okabe-ito`` is the colourblind-safe triad for when the
#: axes themselves are the subject.
PALETTES: dict[str, tuple[tuple[float, float, float], ...]] = {
    "mono": ((0.20, 0.20, 0.20),) * 3,
    "okabe-ito": ((0.835, 0.369, 0.000), (0.000, 0.620, 0.451), (0.000, 0.447, 0.698)),
}

#: Chain roles after ``tcren orient``, and the colours used throughout these scenes.
CHAIN_COLOURS: dict[str, tuple[str, str]] = {
    "A": ("Vα", "marine"), "B": ("Vβ", "orange"), "C": ("peptide", "yellow"),
    "D": ("MHC α", "grey70"), "E": ("MHC β / β2m", "grey60"),
}

_CYLINDER, _CONE = 9.0, 27.0        # PyMOL CGO opcodes; see pymol.cgo


def gizmo_cgo(
    *,
    arm: float = 10.0,
    radius: float = 0.30,
    head_length: float = 0.30,
    head_radius: float = 2.8,
    palette: str = "mono",
) -> tuple[list[float], list[tuple[float, float, float]]]:
    """The axis triad as CGO, at the world origin, in world units.

    Absolute size does not matter: the triad is rendered on its own and scaled to the gizmo tile,
    so only the *proportions* here are visible. What they control is how the figure reads — a thin
    shaft with a distinct arrowhead, rather than the fat default axes that dominate a panel.

    Args:
        arm: Arrow length.
        radius: Shaft radius. Thin relative to ``arm`` is the intent.
        head_length: Arrowhead length as a fraction of ``arm``.
        head_radius: Arrowhead base radius as a multiple of ``radius``.
        palette: A key of :data:`PALETTES`.

    Returns:
        ``(cgo, tips)`` — the flat CGO float list, and the three arrow tips, where labels go.

    Raises:
        ValueError: If ``palette`` is unknown, or a size is non-positive.

    Example:
        >>> cgo, tips = gizmo_cgo()
        >>> tips[1]
        (0.0, 10.0, 0.0)
    """
    if palette not in PALETTES:
        raise ValueError(f"unknown palette {palette!r}; choose from {sorted(PALETTES)}")
    if min(arm, radius, head_length, head_radius) <= 0:
        raise ValueError("arm, radius, head_length and head_radius must all be positive")
    if head_length >= 1.0:
        raise ValueError("head_length is a fraction of arm and must be < 1")

    cgo: list[float] = []
    tips: list[tuple[float, float, float]] = []
    for i, (r, g, b) in enumerate(PALETTES[palette]):
        d = [0.0, 0.0, 0.0]
        d[i] = 1.0
        joint = [d[k] * arm * (1.0 - head_length) for k in range(3)]
        tip = [d[k] * arm for k in range(3)]
        cgo += [_CYLINDER, 0.0, 0.0, 0.0, *joint, radius, r, g, b, r, g, b]
        cgo += [_CONE, *joint, *tip, radius * head_radius, 0.0, r, g, b, r, g, b, 1.0, 1.0]
        tips.append(tuple(tip))
    return cgo, tips


def label_points(rotation, tips, *, offset: float = 1.42, head_on: float = 0.25):
    """Where each axis label goes, pushed clear of its arrow **in screen space**.

    Anchoring a label to its tip in world coordinates fails exactly where it matters: an axis
    pointing at the viewer foreshortens to a dot, and its label lands on top of the origin and the
    other two labels. Pushing outward along the axis's *projected* direction keeps every label
    clear whatever the camera does, and an axis too head-on to have a projected direction is
    pushed down-left instead, where it reads as belonging to the dot at the origin.

    Args:
        rotation: The nine floats of the camera rotation, column-major as PyMOL stores them.
        tips: The three arrow tips from :func:`gizmo_cgo`.
        offset: How far past the tip to sit, as a multiple of the arm length.
        head_on: Projected length below which an axis counts as pointing at the viewer.

    Returns:
        Three world-space points, rounded for embedding in a script.

    Example:
        >>> pts = label_points([1,0,0, 0,1,0, 0,0,1], [(10,0,0), (0,10,0), (0,0,10)])
        >>> pts[0][0] > 10          # pushed out past the tip
        True
    """
    import math as _m

    rot = [float(v) for v in rotation]
    out = []
    for k, tip in enumerate(tips):
        arm = _m.sqrt(sum(c * c for c in tip)) or 1.0
        sx, sy = rot[k * 3], rot[k * 3 + 1]          # this axis, projected onto the screen
        p = _m.hypot(sx, sy)
        if p < head_on:                              # head-on: no direction to push along
            sx, sy = -0.7071, -0.7071                # down-left, clear of the other two
        else:
            sx, sy = sx / p, sy / p
        # A screen-space push (sx, sy, 0) becomes this world vector under the inverse rotation;
        # `rot` is column-major, so the transpose that inverts it is rot[i * 3 + j].
        push = [(rot[i * 3] * sx + rot[i * 3 + 1] * sy) * arm * (offset - 1.0) for i in range(3)]
        out.append([round(tip[i] + push[i], 4) for i in range(3)])
    return out


def gizmo_scene(
    rotation: list[float],
    *,
    axes: tuple[Axis, Axis, Axis] = CANONICAL_AXES,
    short_labels: bool = True,
    label_size: float = 26.0,
    label_colour: str = "gray20",
    label_offset: float = 1.42,
    **cgo_kwargs,
) -> str:
    """A PyMOL scene holding only the triad, seen under ``rotation``.

    Args:
        rotation: The first nine floats of the main scene's ``cmd.get_view()`` — the same camera
            rotation, so the triad reports the orientation the molecule is actually drawn in.
        axes: The axes to draw, in x, y, z order.
        short_labels: Use :attr:`Axis.short` (fits a small tile) rather than :attr:`Axis.label`.
        label_size: PyMOL ``label_size``. Large, because the tile is later scaled down.
        label_colour: Any PyMOL colour name.
        label_offset: Where the label sits along the arm, as a multiple of its length.
        **cgo_kwargs: Forwarded to :func:`gizmo_cgo`.

    Returns:
        A PyMOL script body. Labels are pseudoatoms because CGO has no text primitive; they
        ray-trace like any other label.

    Raises:
        ValueError: If ``rotation`` is not nine floats.
    """
    if rotation is None or len(rotation) != 9:
        raise ValueError("rotation must be the first nine floats of cmd.get_view()")
    cgo, tips = gizmo_cgo(**cgo_kwargs)
    lines = [
        "from pymol.cgo import CYLINDER, CONE",
        f"cmd.load_cgo({cgo!r}, 'gizmo')",
        f"cmd.set('label_size', {label_size})",
        f"cmd.set('label_color', {label_colour!r})",
        "cmd.set('label_font_id', 7)",
    ]
    for axis, pos in zip(axes, label_points(rotation, tips, offset=label_offset)):
        text = axis.short if short_labels else axis.label
        lines.append(f"cmd.pseudoatom('lab_{axis.letter}', pos={pos!r}, label={text!r})")
    # The same rotation as the molecule, framed on the triad alone. `zoom` on everything keeps the
    # labels inside the tile; without it a long label is clipped at the edge.
    lines += [
        f"cmd.set_view({list(rotation) + [0.0, 0.0, -60.0, 0.0, 0.0, 0.0, 20.0, 100.0, -20.0]!r})",
        "cmd.zoom('all', buffer=1.5, complete=1)",
    ]
    return "\n".join(lines)


def composite(base_png, tile_png, out_png, *, corner: str = "bottom-left",
              scale: float = 0.19, margin: float = 0.025):
    """Paste the gizmo tile into a corner of the render, preserving alpha.

    Args:
        base_png: The molecule render.
        tile_png: The gizmo render, transparent outside the arrows.
        out_png: Where to write. May be ``base_png``.
        corner: A key of :data:`CORNERS`.
        scale: Tile width as a fraction of the base image width.
        margin: Inset from the edges, as a fraction of the base image width.

    Returns:
        The path written.

    Raises:
        ValueError: If ``corner`` is unknown or ``scale``/``margin`` leave no room.
    """
    from PIL import Image

    if corner not in CORNERS:
        raise ValueError(f"unknown corner {corner!r}; choose from {sorted(CORNERS)}")
    if not 0 < scale <= 1 or not 0 <= margin < 0.5:
        raise ValueError("scale must be in (0, 1] and margin in [0, 0.5)")

    base = Image.open(base_png).convert("RGBA")
    tile = Image.open(tile_png).convert("RGBA")
    w = max(1, int(round(base.width * scale)))
    tile = tile.resize((w, max(1, round(w * tile.height / tile.width))), Image.LANCZOS)
    pad = int(round(base.width * margin))
    sx, sy = CORNERS[corner]
    x = pad if sx < 0 else base.width - tile.width - pad
    y = base.height - tile.height - pad if sy < 0 else pad     # image y runs downward
    base.alpha_composite(tile, (x, y))
    base.save(out_png)
    return Path(out_png)


_HEADER = "\n".join([
    "from pymol import cmd",
    # Off first, and for the session. PyMOL re-frames the camera on every new object by default,
    # so loading a CGO or a label pseudoatom would silently zoom the scene onto it.
    "cmd.set('auto_zoom', 0)",
    "cmd.bg_color('white')",
    "cmd.set('ray_opaque_background', 0)",
    "cmd.set('orthoscopic', 1)",
    "cmd.set('ray_shadows', 0)",
    "cmd.set('antialias', 2)",
    "cmd.set('cartoon_sampling', 12)",
    "cmd.set('specular', 0.15)",
    "cmd.set('ambient', 0.18)",
    "cmd.set('direct', 0.55)",
    "cmd.set('reflect', 0.35)",
    # Flat, even lighting and no fog, so panels of the same figure are comparable rather than
    # each being shaded by how deep its own bounding box happens to be.
    "cmd.set('depth_cue', 0)",
    "cmd.set('ray_trace_fog', 0)",
])


def _run(script: str, *, pymol_bin: str | None = None) -> str:
    """Execute a script body under a headless PyMOL; return its stdout."""
    exe = pymol_bin or shutil.which("pymol")
    if exe is None:
        raise RuntimeError("pymol not found on PATH; install it or pass pymol_bin=")
    with tempfile.NamedTemporaryFile("w", suffix=".py", delete=False) as fh:
        fh.write(_HEADER + "\n" + script + "\n")
        path = fh.name
    try:
        done = subprocess.run([exe, "-cq", path], check=True, capture_output=True, text=True)
    finally:
        Path(path).unlink(missing_ok=True)
    return done.stdout


def probe_rotation(scene: str, *, pymol_bin: str | None = None) -> list[float]:
    """Run ``scene`` and report the camera rotation it ends up with.

    The gizmo must be drawn under the same rotation as the molecule, and that is not known until
    the scene has loaded, turned and zoomed — ``cmd.zoom`` leaves the rotation alone but
    ``cmd.turn`` and ``cmd.orient`` do not, so reading it back beats assuming it.

    Args:
        scene: A PyMOL script body that loads structures and sets the view.
        pymol_bin: Override the ``pymol`` executable.

    Returns:
        The first nine floats of ``cmd.get_view()``.
    """
    out = _run(scene + "\nprint('ROT', *['%.6f' % v for v in cmd.get_view()[:9]])",
               pymol_bin=pymol_bin)
    for line in reversed(out.splitlines()):
        if line.startswith("ROT "):
            return [float(v) for v in line.split()[1:]]
    raise RuntimeError(f"no view reported by pymol; output was:\n{out[-2000:]}")


def render(
    scene: str,
    png: str | Path,
    *,
    size: tuple[int, int] = (1200, 1200),
    dpi: int = 300,
    gizmo: bool = True,
    corner: str = "bottom-left",
    gizmo_scale: float = 0.19,
    gizmo_margin: float = 0.025,
    pymol_bin: str | None = None,
    **gizmo_kwargs,
) -> Path:
    """Ray-trace ``scene`` to ``png``, with the canonical-frame gizmo in a corner.

    Args:
        scene: A PyMOL script body: load structures, style them, set and zoom the view.
        png: Output path.
        size: ``(width, height)`` in pixels. 1200 px at 300 dpi is a 4-inch figure panel.
        dpi: Written into the PNG so a document places it at a known physical size.
        gizmo: Draw the axis triad. Turn it off where the frame is not the point.
        corner: Which corner the gizmo goes in; a key of :data:`CORNERS`.
        gizmo_scale: Gizmo tile width as a fraction of the image width.
        gizmo_margin: Gizmo inset from the edges, as a fraction of the image width.
        pymol_bin: Override the ``pymol`` executable.
        **gizmo_kwargs: Forwarded to :func:`gizmo_scene` / :func:`gizmo_cgo`.

    Returns:
        The path written.

    Example:
        >>> render('cmd.load("1ao7.pdb")\\ncmd.show("cartoon")', "fig.png")   # doctest: +SKIP
    """
    png = Path(png)
    _run(f"{scene}\ncmd.ray({size[0]}, {size[1]})\ncmd.png(r'{png}', dpi={dpi})",
         pymol_bin=pymol_bin)
    if not png.exists():
        raise RuntimeError(f"pymol did not write {png}")
    if not gizmo:
        return png
    rotation = probe_rotation(scene, pymol_bin=pymol_bin)
    tile_px = max(64, int(size[0] * gizmo_scale * 2))       # oversampled, then downscaled sharp
    tile = png.with_name(png.stem + "_gizmo.png")
    _run(f"{gizmo_scene(rotation, **gizmo_kwargs)}\ncmd.ray({tile_px}, {tile_px})"
         f"\ncmd.png(r'{tile}', dpi={dpi})", pymol_bin=pymol_bin)
    try:
        composite(png, tile, png, corner=corner, scale=gizmo_scale, margin=gizmo_margin)
    finally:
        tile.unlink(missing_ok=True)
    return png


# --- scene presets ----------------------------------------------------------------------------
# Chain roles after `tcren orient`: A=Vα, B=Vβ, C=peptide, D=MHCα, E=MHCβ/β2m.

_STYLE = "\n".join([
    'cmd.hide("everything")',
    'cmd.show("cartoon")',
] + [f'cmd.color("{c}", "chain {ch}")' for ch, (_, c) in CHAIN_COLOURS.items()])

#: Looking down +x — the thin axis — so the image plane is the groove axis by the docking normal.
#: This is the view in which the crossing and incident angles are visible as angles.
VIEW_SIDE = [0., 1., 0., 0., 0., 1., 1., 0., 0.]
#: Looking down +z, the docking normal: the groove from where the TCR sits.
VIEW_TOP = [1., 0., 0., 0., 1., 0., 0., 0., 1.]


def _scene(loads: str, body: str, rotation: list[float], zoom: str) -> str:
    view = list(rotation) + [0., 0., -300., 0., 0., 0., 100., 520., -20.]
    return f"{loads}\n{body}\ncmd.set_view({view!r})\n{zoom}"


def overlay_scene(ids, canon_dir, *, limit: int = 8, transparency: float = 0.55) -> str:
    """Superpose a set of oriented structures, seen side-on.

    Args:
        ids: PDB ids present in ``canon_dir`` as ``<id>.pdb.gz``.
        canon_dir: The canonical (oriented) structure directory.
        limit: Draw at most this many; past roughly eight the overlay stops being readable.
        transparency: Cartoon transparency, so the spread of the ensemble shows through.

    Returns:
        A PyMOL scene body for :func:`render`.
    """
    ids = list(ids)[:limit]
    loads = "\n".join(f'cmd.load(r"{Path(canon_dir) / f"{p}.pdb.gz"}", "{p}")' for p in ids)
    body = f'{_STYLE}\ncmd.set("cartoon_transparency", {transparency})'
    return _scene(loads, body, VIEW_SIDE, 'cmd.zoom("all", buffer=6, complete=1)')


def groove_scene(pid, canon_dir, *, surface: bool = False) -> str:
    """One complex from above the groove: peptide as sticks in the MHC cleft.

    The layout histo.fyi uses for its structure pages — a pale MHC with the peptide threaded
    along the cleft, which is the most legible way to show what is presented.

    Args:
        pid: PDB id.
        canon_dir: The canonical (oriented) structure directory.
        surface: Add a translucent molecular surface over the MHC ribbon.

    Returns:
        A PyMOL scene body for :func:`render`.
    """
    loads = f'cmd.load(r"{Path(canon_dir) / f"{pid}.pdb.gz"}", "m")'
    body = "\n".join([
        'cmd.hide("everything")',
        'cmd.show("cartoon", "chain D+E")',
        'cmd.color("grey80", "chain D+E")',
        'cmd.color("salmon", "chain D+E and ss H")',      # the groove helices
        'cmd.color("palecyan", "chain D+E and ss S")',    # the β-sheet floor
        'cmd.show("sticks", "chain C")',
        'cmd.color("yellow", "chain C")',
        'cmd.util.cnc("chain C")',
        'cmd.set("stick_radius", 0.25)',
    ] + ([
        'cmd.show("surface", "chain D+E")',
        'cmd.set("transparency", 0.45)',
        'cmd.set("surface_quality", 1)',
    ] if surface else []))
    return _scene(loads, body, VIEW_TOP, 'cmd.zoom("chain C", buffer=12, complete=1)')


def interface_scene(pid, canon_dir, cdr_resi) -> str:
    """The recognition interface: peptide plus the CDR loops that touch it.

    Args:
        pid: PDB id.
        canon_dir: The canonical (oriented) structure directory.
        cdr_resi: ``{"TRA": [...], "TRB": [...]}`` PDB residue numbers, as from annotating the
            pre-orientation structure — orientation preserves numbering.

    Returns:
        A PyMOL scene body for :func:`render`.
    """
    ra = "+".join(str(i) for i in cdr_resi.get("TRA", [])) or "0"
    rb = "+".join(str(i) for i in cdr_resi.get("TRB", [])) or "0"
    sel = f'(chain A and resi {ra}) or (chain B and resi {rb})'
    loads = f'cmd.load(r"{Path(canon_dir) / f"{pid}.pdb.gz"}", "m")'
    body = "\n".join([
        'cmd.hide("everything")',
        'cmd.show("cartoon", "chain D+E")',
        'cmd.color("grey80", "chain D+E")',
        'cmd.set("cartoon_transparency", 0.6, "chain D+E")',
        'cmd.show("sticks", "chain C")',
        'cmd.color("yellow", "chain C")',
        'cmd.util.cnc("chain C")',
        f'cmd.show("cartoon", "chain A and resi {ra}")',
        f'cmd.color("marine", "chain A and resi {ra}")',
        f'cmd.show("cartoon", "chain B and resi {rb}")',
        f'cmd.color("orange", "chain B and resi {rb}")',
        f'cmd.show("sticks", "{sel}")',
        'cmd.set("stick_radius", 0.16)',
    ])
    return _scene(loads, body, VIEW_SIDE, f'cmd.zoom("chain C or ({sel})", buffer=8, complete=1)')


def _selfcheck() -> None:
    """Triad geometry, corner arithmetic and the guards — no PyMOL, no structures needed."""
    cgo, tips = gizmo_cgo()
    assert len(tips) == 3
    # Two primitives per axis. CYLINDER is 14 floats (opcode, 2 points, radius, 2 colours) and
    # CONE is 17 (opcode, 2 points, 2 radii, 2 colours, 2 cap flags).
    assert len(cgo) == 3 * (14 + 17), len(cgo)
    assert cgo[0] == _CYLINDER and cgo[14] == _CONE

    # Each arm starts at the origin and runs along its own axis, and no other.
    for i in range(3):
        start, end = cgo[i * 31 + 1: i * 31 + 4], cgo[i * 31 + 4: i * 31 + 7]
        assert start == [0.0, 0.0, 0.0], start
        assert end[i] > 0, (i, end)
        assert all(abs(end[k]) < 1e-12 for k in range(3) if k != i), (i, end)
        assert tips[i][i] > end[i], "the head must extend past the shaft"

    # The head is wider than the shaft, or it does not read as an arrow at all.
    shaft_r, head_r = cgo[7], cgo[14 + 7]
    assert head_r > shaft_r * 2, (shaft_r, head_r)
    assert gizmo_cgo(arm=20.0)[1][0][0] == 20.0

    for kw in ({"palette": "neon"}, {"arm": 0}, {"radius": -1}, {"head_length": 1.0}):
        try:
            gizmo_cgo(**kw)
        except ValueError:
            pass
        else:
            raise AssertionError(f"expected ValueError for {kw}")

    # Labels push outward along each axis's projected direction. Under VIEW_SIDE the y axis points
    # at the viewer, so it has no projected direction and must fall back to the down-left corner
    # rather than piling onto the origin with the other two.
    tips = gizmo_cgo()[1]
    face_on = label_points([1, 0, 0, 0, 1, 0, 0, 0, 1], tips)
    assert face_on[0][0] > tips[0][0] and face_on[1][1] > tips[1][1], face_on
    side = label_points(VIEW_SIDE, tips)
    assert side[1][0] < 0 and side[1][2] < 0, f"head-on label must fall back down-left: {side[1]}"
    def _gap(pt, tip):
        return sum((pt[i] - tip[i]) ** 2 for i in range(3)) ** 0.5

    far = label_points(VIEW_SIDE, tips, offset=2.5)
    assert _gap(far[0], tips[0]) > _gap(side[0], tips[0]), "offset must push the label further out"

    scene = gizmo_scene(VIEW_SIDE)
    assert "load_cgo" in scene and "pseudoatom" in scene and "set_view" in scene
    for axis in CANONICAL_AXES:
        assert axis.short in scene, axis
    assert CANONICAL_AXES[2].label in gizmo_scene(VIEW_SIDE, short_labels=False)
    try:
        gizmo_scene([1, 0, 0])
    except ValueError:
        pass
    else:
        raise AssertionError("expected ValueError for a short rotation")

    # Corner arithmetic: each corner must land the tile on the right side, inside the frame.
    from PIL import Image
    base, tile = Image.new("RGBA", (400, 300)), Image.new("RGBA", (80, 80), (0, 0, 0, 255))
    import tempfile as _tf
    with _tf.TemporaryDirectory() as d:
        bp, tp = Path(d) / "b.png", Path(d) / "t.png"
        base.save(bp); tile.save(tp)
        seen = {}
        for name in CORNERS:
            out = Path(d) / f"{name}.png"
            composite(bp, tp, out, corner=name, scale=0.2, margin=0.02)
            px = Image.open(out).convert("RGBA")
            ink = [(x, y) for y in range(px.height) for x in range(px.width)
                   if px.getpixel((x, y))[3] > 0]
            xs, ys = [p[0] for p in ink], [p[1] for p in ink]
            seen[name] = (min(xs), min(ys), max(xs), max(ys))
            assert min(xs) >= 0 and max(xs) < px.width, name
            assert min(ys) >= 0 and max(ys) < px.height, name
        assert seen["bottom-left"][0] < seen["bottom-right"][0], seen
        assert seen["bottom-left"][1] > seen["top-left"][1], seen
        assert seen["top-right"][0] > seen["top-left"][0], seen
        for bad in ({"corner": "middle"}, {"scale": 0}, {"margin": 0.9}):
            try:
                composite(bp, tp, Path(d) / "x.png", **bad)
            except ValueError:
                pass
            else:
                raise AssertionError(f"expected ValueError for {bad}")

    for fn, args in ((overlay_scene, (["1ao7"], "/tmp")), (groove_scene, ("1ao7", "/tmp")),
                     (interface_scene, ("1ao7", "/tmp", {"TRA": [1, 2], "TRB": [3]}))):
        body = fn(*args)
        assert "cmd.load" in body and "cmd.set_view" in body and "cmd.zoom" in body, fn.__name__
    assert "resi 1+2" in interface_scene("1ao7", "/tmp", {"TRA": [1, 2], "TRB": [3]})
    assert len(overlay_scene([str(i) for i in range(40)], "/tmp").split("cmd.load")) - 1 == 8

    print("pymol viz self-check OK: arms run along their own axes with wider heads, the gizmo "
          "scene names the canonical frame, all four corners land inside the frame, and the "
          "three scene presets load, view and zoom")


if __name__ == "__main__":
    _selfcheck()
