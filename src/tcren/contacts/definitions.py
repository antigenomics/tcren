"""Flexible multi-threshold contact definition.

Beyond the legacy single 5 Å all-atom contact (the TCRen parity default, ``d1``), this adds
two coarser residue-level layers: ``d2`` over Cβ atoms (Cα for glycine) and ``d3`` over Cα
atoms. The layers nest from tight side-chain proximity to backbone neighbourhood, giving the
2D maps and scoring a tunable contact model without changing the 5 Å default.
"""

from __future__ import annotations

from dataclasses import dataclass

import polars as pl

from ..structure.model import Structure
from .geometry import all_atom_contacts, representative_atom_contacts


@dataclass(frozen=True, slots=True)
class ContactDefinition:
    """Three nested contact thresholds (Å).

    Attributes:
        d1: closest heavy-atom distance (all-atom contact).
        d2: closest Cβ distance (Cα for glycine / missing Cβ).
        d3: closest Cα distance.
    """

    d1: float = 5.0
    d2: float = 8.0
    d3: float = 12.0


TCREN_DEFAULT = ContactDefinition()


def multi_contacts(
    structure: Structure, definition: ContactDefinition = TCREN_DEFAULT
) -> pl.DataFrame:
    """Stacked inter-chain residue contacts across the three layers.

    Returns the union of the ``d1``/``d2``/``d3`` residue-pair tables with a ``layer`` column
    (``"d1"``/``"d2"``/``"d3"``) and the layer's distance. A residue pair can appear in
    several layers; callers filter by ``layer`` as needed.
    """
    layers = {
        "d1": all_atom_contacts(structure, cutoff=definition.d1),
        "d2": representative_atom_contacts(structure, kind="cb", cutoff=definition.d2),
        "d3": representative_atom_contacts(structure, kind="ca", cutoff=definition.d3),
    }
    frames = [df.with_columns(pl.lit(name).alias("layer")) for name, df in layers.items()]
    return pl.concat(frames) if any(f.height for f in frames) else frames[0].with_columns(
        pl.lit("d1").alias("layer")
    )
