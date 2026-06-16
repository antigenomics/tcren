"""Colorblind-safe palette for complementarity maps (Okabe-Ito), shared by 2D and 3D.

Hue encodes the complex chain (tra/trb/peptide/mhca/mhcb); CDR loops get distinct shades.
"""

from __future__ import annotations

# Okabe-Ito colorblind-safe hues.
_OKABE = {
    "orange": "#E69F00", "skyblue": "#56B4E9", "green": "#009E73", "yellow": "#F0E442",
    "blue": "#0072B2", "vermillion": "#D55E00", "purple": "#CC79A7", "grey": "#999999",
}

CHAIN_COLOR = {
    "tra": _OKABE["blue"],
    "trb": _OKABE["skyblue"],
    "trd": _OKABE["blue"],
    "trg": _OKABE["skyblue"],
    "peptide": _OKABE["vermillion"],
    "mhca": _OKABE["green"],
    "mhcb": _OKABE["yellow"],
    "b2m": _OKABE["grey"],
}

# CDR/region-specific shades (override the chain hue for emphasis).
REGION_COLOR = {
    "cdr1": "#9ecae1", "cdr2": "#4292c6", "cdr3": "#08519c",  # TRA/TRB CDR blues
    "mhc_helix_a1": "#41ab5d", "mhc_helix_a2": "#006d2c", "mhc_helix_b1": "#fec44f",
    "groove_floor": "#d9d9d9", "peptide": _OKABE["vermillion"],
}

DEFAULT_COLOR = _OKABE["grey"]


def color_for(complex_chain: str | None, complex_region: str | None) -> str:
    """Resolve a fill color: region shade if known, else chain hue, else grey."""
    if complex_region and complex_region in REGION_COLOR:
        return REGION_COLOR[complex_region]
    if complex_chain and complex_chain in CHAIN_COLOR:
        return CHAIN_COLOR[complex_chain]
    return DEFAULT_COLOR
