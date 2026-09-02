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
from ..structure.model import PEPTIDE_TYPE, RECEPTOR_TYPES
from .geometry import all_atom_contacts, representative_atom_contacts

# Long side chains (Arg, Lys, Trp) put two heavy atoms within 5 A while their Calpha atoms sit far
# apart, so the representative-atom layers are built well past their nominal 8/12 A defaults; the
# per-descriptor thresholds are applied afterwards by filtering.
_REP_BUILD_CUTOFF = 18.0
_KEY = ["chain.id.from", "residue.index.from", "chain.id.to", "residue.index.to"]



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


def _interface_layers(structure: Structure, cutoff: float, partner=(PEPTIDE_TYPE,),
                      receptor=None) -> pl.DataFrame:
    """The d1/d2/d3 layers pivoted onto one row per ``receptor``:``partner`` residue pair.

    Args:
        structure: The chain-typed complex.
        cutoff: The d1 (closest heavy-atom) threshold, Å.
        partner: Chain types on the partner side.
        receptor: Chain types on the receptor side. ``None`` (default) is the TCR/BCR receptor,
            which is what every caller wanted while the only interfaces were TCR:peptide and
            TCR:MHC. Pass the MHC types to make the **groove** the receptor and the peptide the
            partner — the peptide:MHC arm, whose Hamiltonian is a separate field from the
            receptor's and whose partner-side numbering is the plain peptide position.

    Returns a frame keyed by the residue pair with columns ``d1``/``d2``/``d3`` (Angstrom, null where
    that layer does not reach), ``aa.tcr``/``aa.pep`` and the d1 atom names. The receptor side is
    normalised to ``aa.tcr`` regardless of which side the canonical chain ordering put it on — the
    column keeps its name across arms so one site schema serves all three.
    """
    stacked = multi_contacts(
        structure,
        ContactDefinition(d1=cutoff, d2=_REP_BUILD_CUTOFF, d3=_REP_BUILD_CUTOFF),
    )
    ctype = {c.chain_id: c.chain_type for c in structure.chains}
    stacked = stacked.with_columns(
        pl.col("chain.id.from").replace_strict(ctype, default=None).alias("type.from"),
        pl.col("chain.id.to").replace_strict(ctype, default=None).alias("type.to"),
    )
    tcr, pep = list(RECEPTOR_TYPES if receptor is None else receptor), list(partner)
    fwd = pl.col("type.from").is_in(tcr) & pl.col("type.to").is_in(pep)
    rev = pl.col("type.from").is_in(pep) & pl.col("type.to").is_in(tcr)
    stacked = stacked.filter(fwd | rev).with_columns(
        pl.when(fwd).then(pl.col("residue.aa.from")).otherwise(pl.col("residue.aa.to")).alias("aa.tcr"),
        pl.when(fwd).then(pl.col("residue.aa.to")).otherwise(pl.col("residue.aa.from")).alias("aa.pep"),
        # Normalised receptor-side residue key: degree must be counted on the TCR side, because a
        # peptide residue sits *inside* the groove ringed by receptor and has a high degree in every
        # real complex. It is the receptor side-chain reaching too many partners that is the tell.
        pl.when(fwd).then(pl.col("chain.id.from")).otherwise(pl.col("chain.id.to")).alias("key.tcr.chain"),
        pl.when(fwd).then(pl.col("residue.index.from")).otherwise(pl.col("residue.index.to")).alias("key.tcr.res"),
    )
    if stacked.is_empty():
        return stacked.select(*_KEY, "aa.tcr", "aa.pep").with_columns(
            pl.lit(None, dtype=pl.Float64).alias(c) for c in ("d1", "d2", "d3")
        )
    # Identity is taken from the UNION of the three layers, not from d1 alone: a pair present in the
    # Calpha shell but making no contact must still carry its residue types, or every shell
    # descriptor that reads `aa` would collapse to the contact set without saying so.
    ident = stacked.select(*_KEY, "aa.tcr", "aa.pep", "key.tcr.chain", "key.tcr.res").unique(subset=_KEY)
    wide = stacked.filter(pl.col("layer") == "d1").select(*_KEY, pl.col("dist").alias("d1"))
    for layer in ("d2", "d3"):
        part = (stacked.filter(pl.col("layer") == layer)
                .select(*_KEY, pl.col("dist").alias(layer)))
        # A residue pair appears at most once per layer (each keeps its closest atom pair).
        wide = wide.join(part, on=_KEY, how="full", coalesce=True)
    return wide.join(ident, on=_KEY, how="left")
