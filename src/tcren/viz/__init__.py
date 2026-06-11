"""Visualization: 2D complementarity-map SVG + 3D pocket/CDR overlay."""

from .palette import color_for
from .svg2d import render_complementarity_map

__all__ = ["render_complementarity_map", "color_for", "view_pocket_cdr"]


def view_pocket_cdr(*args, **kwargs):
    """Lazy proxy to :func:`tcren.viz.pocket3d.view_pocket_cdr` (needs py3Dmol)."""
    from .pocket3d import view_pocket_cdr as _impl

    return _impl(*args, **kwargs)
