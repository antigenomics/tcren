"""Approximate MHC peptide-binding pocket (A–F) markers along the peptide track.

The six class-I pockets (A–F) line the groove from the peptide N-terminus (A ≈ P1) to its
C-terminus (F ≈ P-Ω). Without explicit pocket geometry we place A–F evenly along the
projected peptide Cα track — a labelling aid for the 2D map and 3D view, not a structural
pocket definition.
"""

from __future__ import annotations

import numpy as np
import polars as pl

_POCKETS = ("A", "B", "C", "D", "E", "F")


def pocket_markers(markup: pl.DataFrame) -> pl.DataFrame:
    """Place A–F markers along the projected peptide (u, v) track.

    Returns columns ``pocket, peptide_pos, u, v`` (empty if the peptide is not projected).
    """
    pep = (
        markup.filter((pl.col("complex_chain") == "peptide") & pl.col("u").is_not_null())
        .sort("aa_index")
    )
    schema = {"pocket": pl.Utf8, "peptide_pos": pl.Float64, "u": pl.Float64, "v": pl.Float64}
    if pep.height < 2:
        return pl.DataFrame(schema=schema)
    u = pep["u"].to_numpy()
    v = pep["v"].to_numpy()
    # Arc-length parameterisation of the peptide track, sampled at 6 even fractions.
    seg = np.r_[0.0, np.cumsum(np.hypot(np.diff(u), np.diff(v)))]
    total = seg[-1] or 1.0
    rows = []
    for k, name in enumerate(_POCKETS):
        frac = k / (len(_POCKETS) - 1)
        target = frac * total
        rows.append(
            {
                "pocket": name,
                "peptide_pos": frac * (len(u) - 1),
                "u": float(np.interp(target, seg, u)),
                "v": float(np.interp(target, seg, v)),
            }
        )
    return pl.DataFrame(rows, schema=schema)
