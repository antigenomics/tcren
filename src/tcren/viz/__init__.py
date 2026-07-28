"""Visualization: 2D complementarity-map SVG, 3D pocket/CDR overlay, PyMOL figure renders."""

from .palette import color_for
from .pymol import CANONICAL_AXES, groove_scene, interface_scene, overlay_scene, render
from .svg2d import render_complementarity_map

__all__ = [
    "render_complementarity_map", "color_for", "view_pocket_cdr",
    "render", "overlay_scene", "groove_scene", "interface_scene", "CANONICAL_AXES",
]


def view_pocket_cdr(*args, **kwargs):
    """Lazy proxy to :func:`tcren.viz.pocket3d.view_pocket_cdr` (needs py3Dmol)."""
    from .pocket3d import view_pocket_cdr as _impl

    return _impl(*args, **kwargs)
