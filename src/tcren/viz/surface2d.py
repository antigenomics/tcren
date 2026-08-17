"""Direct SVG builder for pMHC surface-topology maps (:mod:`tcren.surface`).

One ``<rect>`` per grid cell, coloured by the channel being shown, with the cell's height,
chemistry and owning region carried as ``data-*`` attributes and a ``<title>`` tooltip — same
contract as :mod:`tcren.viz.svg2d`: the SVG is a figure *and* a queryable artifact. Pure string
building, no dependencies.

Two ramps, chosen by what the channel means rather than by taste. ``h`` is a magnitude above the
groove floor, so it gets a sequential viridis-like ramp. ``phobic`` and ``charge`` have a
meaningful zero (neutral / uncharged), so they get a diverging ramp **centred on zero**, not on
the data range — a range-fitted diverging ramp paints the least-hydrophobic cell of an all-greasy
surface as if it were hydrophilic.
"""

from __future__ import annotations

import numpy as np

# Sequential ramp for heights: viridis, subsampled. Perceptually uniform and colourblind-safe.
_SEQUENTIAL = ("#440154", "#472d7b", "#3b528b", "#2c728e", "#21918c",
               "#28ae80", "#5ec962", "#addc30", "#fde725")
# Diverging ramp for signed channels: ColorBrewer RdBu, reversed (blue = negative, red = positive).
_DIVERGING = ("#2166ac", "#4393c3", "#92c5de", "#d1e5f0", "#f7f7f7",
              "#fddbc7", "#f4a582", "#d6604d", "#b2182b")

#: Channels whose zero is meaningful, so the ramp is centred rather than range-fitted.
SIGNED_CHANNELS = ("phobic", "charge")

_OUTLINE = {"peptide": "#D55E00", "mhc_helix_a1": "#41ab5d", "mhc_helix_a2": "#006d2c",
            "mhc_helix_b1": "#fec44f", "mhc_floor": "#d9d9d9", "none": "#ffffff"}


def _esc(value) -> str:
    s = "" if value is None else str(value)
    return (s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
            .replace('"', "&quot;"))


def _f(x: float) -> str:
    return f"{x:.2f}"


def _lerp(ramp: tuple[str, ...], t: float) -> str:
    """Colour at ``t`` in [0, 1], linearly interpolated between the ramp's stops."""
    t = min(max(t, 0.0), 1.0)
    pos = t * (len(ramp) - 1)
    i = min(int(pos), len(ramp) - 2)
    frac = pos - i
    a, b = ramp[i].lstrip("#"), ramp[i + 1].lstrip("#")
    rgb = [round(int(a[k:k + 2], 16) * (1 - frac) + int(b[k:k + 2], 16) * frac) for k in (0, 2, 4)]
    return "#{:02x}{:02x}{:02x}".format(*rgb)


def _normaliser(values: np.ndarray, channel: str, vmin=None, vmax=None):
    """Return ``(to_t, ramp, lo, hi)`` for the channel's colour mapping."""
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return (lambda _v: 0.5), _SEQUENTIAL, 0.0, 1.0
    if channel in SIGNED_CHANNELS:
        span = float(max(abs(np.nanpercentile(finite, 2)), abs(np.nanpercentile(finite, 98)))) or 1.0
        lo, hi = -span, span
        return (lambda v: 0.5 + 0.5 * v / span), _DIVERGING, lo, hi
    lo = float(np.nanpercentile(finite, 2)) if vmin is None else vmin
    hi = float(np.nanpercentile(finite, 98)) if vmax is None else vmax
    rng = (hi - lo) or 1.0
    return (lambda v: (v - lo) / rng), _SEQUENTIAL, lo, hi


def render_surface_map(
    smap,
    channel: str = "h",
    *,
    width: int = 520,
    height: int = 760,
    margin: float = 56.0,
    outline_source: bool = True,
    vmin: float | None = None,
    vmax: float | None = None,
    title: str | None = None,
) -> str:
    """Render a :class:`tcren.surface.SurfaceMap` channel to an SVG string.

    Args:
        smap: the map to draw.
        channel: which channel — ``"h"``, ``"phobic"`` or ``"charge"``.
        width, height, margin: canvas geometry in px. The default is portrait because the groove
            axis (y, peptide N→C) is the long one.
        outline_source: stroke each cell in its source region's colour (peptide vs the two
            helices), so the peptide's footprint is legible on top of the value ramp.
        vmin, vmax: fix the colour range instead of taking the data's 2nd/98th percentile. Pass
            the same values across a set of maps to make the figures directly comparable.
        title: caption; defaults to ``"<structure id> <peptide>"``.

    Returns:
        SVG markup (string).

    Raises:
        KeyError: if ``channel`` is not one of the map's channels.
    """
    from ..surface import SOURCE_NAMES

    grid = smap.channels[channel]
    n_y, n_x = smap.grid
    x0, x1, y0, y1 = smap.extent
    cw = (width - 2 * margin) / n_x
    ch = (height - 2 * margin) / n_y
    to_t, ramp, lo, hi = _normaliser(grid, channel, vmin, vmax)

    cap = title if title is not None else f"{smap.structure_id} {smap.peptide}".strip()
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}" data-channel="{_esc(channel)}" '
        f'data-structure="{_esc(smap.structure_id)}" data-peptide="{_esc(smap.peptide)}" '
        f'data-scale="{_esc(smap.scale)}" data-extent="{_esc(smap.extent)}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        f'<text x="{_f(width / 2)}" y="{_f(margin / 2 + 4)}" text-anchor="middle" '
        f'font-size="14" font-family="sans-serif">{_esc(cap)}</text>',
        '<g class="cells" shape-rendering="crispEdges">',
    ]

    for iy in range(n_y):
        for ix in range(n_x):
            v = grid[iy, ix]
            if not np.isfinite(v):
                continue
            # SVG y grows downward; row 0 is the low-y (peptide N-terminal) end, so flip.
            px = margin + ix * cw
            py = margin + (n_y - 1 - iy) * ch
            src = SOURCE_NAMES[int(smap.source[iy, ix])]
            stroke = (f' stroke="{_OUTLINE.get(src, "#ffffff")}" stroke-width="0.6"'
                      if outline_source and src == "peptide" else "")
            h = smap.channels["h"][iy, ix]
            parts.append(
                f'<rect x="{_f(px)}" y="{_f(py)}" width="{_f(cw + 0.5)}" height="{_f(ch + 0.5)}" '
                f'fill="{_lerp(ramp, to_t(v))}"{stroke} data-ix="{ix}" data-iy="{iy}" '
                f'data-source="{_esc(src)}" data-value="{_f(float(v))}">'
                f'<title>{_esc(src)} {channel}={_f(float(v))} h={_f(float(h))} Å</title></rect>'
            )
    parts.append("</g>")

    # Axes: name the directions rather than printing x/y, matching viz.pymol.CANONICAL_AXES.
    parts.append(
        f'<g font-size="11" font-family="sans-serif" fill="#444">'
        f'<text x="{_f(width / 2)}" y="{_f(height - margin / 3)}" text-anchor="middle">'
        f'groove width ({_f(x0)} to {_f(x1)} Å)</text>'
        f'<text x="{_f(margin / 3)}" y="{_f(height / 2)}" text-anchor="middle" '
        f'transform="rotate(-90 {_f(margin / 3)} {_f(height / 2)})">'
        f'peptide N→C ({_f(y0)} to {_f(y1)} Å)</text></g>'
    )
    parts.append(_legend(ramp, lo, hi, channel, width, height, margin))
    parts.append("</svg>")
    return "".join(parts)


def _legend(ramp, lo, hi, channel, width, height, margin) -> str:
    """Horizontal colour bar under the map."""
    bar_w, bar_h = width - 2 * margin, 9.0
    x, y = margin, height - margin + 12
    stops = "".join(f'<stop offset="{i / (len(ramp) - 1):.3f}" stop-color="{c}"/>'
                    for i, c in enumerate(ramp))
    gid = f"ramp_{channel}"
    return (
        f'<defs><linearGradient id="{gid}" x1="0" x2="1">{stops}</linearGradient></defs>'
        f'<g class="legend" font-size="10" font-family="sans-serif" fill="#444">'
        f'<rect x="{_f(x)}" y="{_f(y)}" width="{_f(bar_w)}" height="{_f(bar_h)}" fill="url(#{gid})"/>'
        f'<text x="{_f(x)}" y="{_f(y + bar_h + 11)}">{_f(lo)}</text>'
        f'<text x="{_f(x + bar_w)}" y="{_f(y + bar_h + 11)}" text-anchor="end">{_f(hi)}</text>'
        f'<text x="{_f(x + bar_w / 2)}" y="{_f(y + bar_h + 11)}" text-anchor="middle">'
        f'{_esc(channel)}</text></g>'
    )
