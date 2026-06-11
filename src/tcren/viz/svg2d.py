"""Direct SVG builder for 2D complementarity maps.

Renders projected interface residues as squares (AA + number), Cα–Cα "chain contacts" as
bold lines, and closest-atom inter-residue contacts as dashed lines. Every element carries
its data as ``data-*`` attributes plus a ``<title>`` tooltip, so the SVG is both a figure
and a queryable, metadata-bearing artifact. Pure string building — no dependencies.
"""

from __future__ import annotations

import polars as pl

from .palette import color_for

_SQUARE = 22.0  # residue square side (px)


def _esc(value) -> str:
    s = "" if value is None else str(value)
    return (s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
            .replace('"', "&quot;"))


def _f(x: float) -> str:
    return f"{x:.2f}"


def _canvas_mapper(uv, width, height, margin):
    us = [p[0] for p in uv]
    vs = [p[1] for p in uv]
    umin, umax, vmin, vmax = min(us), max(us), min(vs), max(vs)
    du, dv = (umax - umin) or 1.0, (vmax - vmin) or 1.0
    sx = (width - 2 * margin) / du
    sy = (height - 2 * margin) / dv

    def to_px(u, v):
        px = margin + (u - umin) * sx
        py = height - margin - (v - vmin) * sy  # flip y (SVG y-down)
        return px, py

    return to_px


def render_complementarity_map(
    markup: pl.DataFrame,
    contacts: pl.DataFrame | None = None,
    ca_contacts: pl.DataFrame | None = None,
    pockets: pl.DataFrame | None = None,
    width: int = 900,
    height: int = 700,
    margin: float = 60.0,
) -> str:
    """Render a complementarity map to an SVG string.

    Args:
        markup: residue markup table (needs ``u``, ``v``; non-null rows are drawn).
        contacts: closest-atom inter-residue contacts (dashed) with ``structure_chain_1/2``,
            ``aa_index_1/2``, ``min_dist``, ``contact_type``.
        ca_contacts: Cα–Cα chain contacts (bold) with ``structure_chain_1/2``,
            ``aa_index_1/2``, ``ca_dist``.
        pockets: optional A–F pocket markers with ``pocket``, ``u``, ``v``.
        width, height, margin: canvas geometry.

    Returns:
        SVG markup (string).
    """
    drawn = markup.filter(pl.col("u").is_not_null() & pl.col("v").is_not_null())
    if drawn.height == 0:
        raise ValueError("no residues with projected (u, v) to render")

    uv = list(zip(drawn["u"].to_list(), drawn["v"].to_list()))
    if pockets is not None and pockets.height:
        uv += list(zip(pockets["u"].to_list(), pockets["v"].to_list()))
    to_px = _canvas_mapper(uv, width, height, margin)

    pos = {}  # (chain, aa_index) -> (px, py)
    for r in drawn.iter_rows(named=True):
        pos[(r["structure_chain"], r["aa_index"])] = to_px(r["u"], r["v"])

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}" font-family="sans-serif">',
        f'<rect width="{width}" height="{height}" fill="white"/>',
    ]

    # Layer 1: dashed inter-residue (closest-atom) contacts.
    if contacts is not None and contacts.height:
        parts.append('<g class="contacts-dashed" stroke="#888" stroke-width="1" '
                     'stroke-dasharray="4 3">')
        for c in contacts.iter_rows(named=True):
            p1 = pos.get((c["structure_chain_1"], c["aa_index_1"]))
            p2 = pos.get((c["structure_chain_2"], c["aa_index_2"]))
            if not p1 or not p2:
                continue
            parts.append(
                f'<line x1="{_f(p1[0])}" y1="{_f(p1[1])}" x2="{_f(p2[0])}" y2="{_f(p2[1])}" '
                f'data-min-dist="{_f(c["min_dist"])}" data-contact-type="{_esc(c["contact_type"])}" '
                f'data-chain-1="{_esc(c["structure_chain_1"])}" data-chain-2="{_esc(c["structure_chain_2"])}" '
                f'data-res-1="{c["aa_index_1"]}" data-res-2="{c["aa_index_2"]}">'
                f'<title>{_esc(c["contact_type"])} {_f(c["min_dist"])} Å '
                f'{_esc(c["structure_chain_1"])}:{c["aa_index_1"]}–'
                f'{_esc(c["structure_chain_2"])}:{c["aa_index_2"]}</title></line>'
            )
        parts.append("</g>")

    # Layer 2: bold Cα–Cα chain contacts.
    if ca_contacts is not None and ca_contacts.height:
        parts.append('<g class="contacts-ca" stroke="#333" stroke-width="3" opacity="0.5">')
        for c in ca_contacts.iter_rows(named=True):
            p1 = pos.get((c["structure_chain_1"], c["aa_index_1"]))
            p2 = pos.get((c["structure_chain_2"], c["aa_index_2"]))
            if not p1 or not p2:
                continue
            parts.append(
                f'<line x1="{_f(p1[0])}" y1="{_f(p1[1])}" x2="{_f(p2[0])}" y2="{_f(p2[1])}" '
                f'data-ca-dist="{_f(c["ca_dist"])}" '
                f'data-chain-1="{_esc(c["structure_chain_1"])}" data-chain-2="{_esc(c["structure_chain_2"])}" '
                f'data-res-1="{c["aa_index_1"]}" data-res-2="{c["aa_index_2"]}">'
                f'<title>Cα {_f(c["ca_dist"])} Å {_esc(c["structure_chain_1"])}:{c["aa_index_1"]}–'
                f'{_esc(c["structure_chain_2"])}:{c["aa_index_2"]}</title></line>'
            )
        parts.append("</g>")

    # Layer 3: residue squares.
    parts.append('<g class="residues">')
    half = _SQUARE / 2
    for r in drawn.iter_rows(named=True):
        px, py = pos[(r["structure_chain"], r["aa_index"])]
        fill = color_for(r["complex_chain"], r["complex_region"])
        label = f'{r["aa"]}{r["residue_index"]} {r["complex_chain"]}/{r["complex_region"] or ""}'
        parts.append(
            f'<g class="residue" data-chain="{_esc(r["structure_chain"])}" '
            f'data-complex-chain="{_esc(r["complex_chain"])}" '
            f'data-region="{_esc(r["complex_region"])}" '
            f'data-residue-index="{r["residue_index"]}" data-aa-index="{r["aa_index"]}" '
            f'data-aa="{_esc(r["aa"])}">'
            f'<title>{_esc(label)}</title>'
            f'<rect x="{_f(px - half)}" y="{_f(py - half)}" width="{_SQUARE}" height="{_SQUARE}" '
            f'rx="3" fill="{fill}" stroke="black" stroke-width="0.6"/>'
            f'<text x="{_f(px)}" y="{_f(py + 1)}" text-anchor="middle" '
            f'dominant-baseline="middle" font-size="11" font-weight="bold">{_esc(r["aa"])}</text>'
            f'<text x="{_f(px)}" y="{_f(py + half + 8)}" text-anchor="middle" '
            f'font-size="7" fill="#444">{r["residue_index"]}</text>'
            f"</g>"
        )
    parts.append("</g>")

    # Layer 4: A–F pocket markers.
    if pockets is not None and pockets.height:
        parts.append('<g class="pockets" fill="none" stroke="#000" stroke-dasharray="2 2">')
        for p in pockets.iter_rows(named=True):
            px, py = to_px(p["u"], p["v"])
            parts.append(
                f'<g class="pocket" data-pocket="{_esc(p["pocket"])}">'
                f'<circle cx="{_f(px)}" cy="{_f(py)}" r="16"/>'
                f'<text x="{_f(px)}" y="{_f(py - 20)}" text-anchor="middle" font-size="12" '
                f'stroke="none" fill="#000">{_esc(p["pocket"])}</text></g>'
            )
        parts.append("</g>")

    parts.append("</svg>")
    return "".join(parts)
