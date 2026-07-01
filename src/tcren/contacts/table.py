"""Annotate and symmetrise the residue contact table.

Joins per-residue annotations (chain type, region type, region start, amino acid) onto
the raw contacts and mirrors the R ``rbind(contacts, swapped)`` symmetrisation, yielding
the fully annotated, bidirectional contact table the contact map is built from.
"""

from __future__ import annotations

import polars as pl

from ..structure.model import Structure
from .geometry import all_atom_contacts


def residue_annotation(structure: Structure) -> pl.DataFrame:
    """Per-residue annotation table for joining onto contacts.

    Columns: ``chain.id``, ``residue.index``, ``chain.type``, ``chain.supertype``,
    ``region.type``, ``region.start``, ``residue.aa``. ``region.type``/``region.start``
    are null for residues without a region annotation.
    """
    rows = []
    for chain in structure.chains:
        # Map each residue's sequential index to its region (type + region start).
        region_of: dict[int, tuple[str, int]] = {}
        for region in chain.regions:
            for res in region.residues:
                region_of[res.seq_index] = (region.region_type, region.start_seq_index)
        for res in chain.residues:
            region_type, region_start = region_of.get(res.seq_index, (None, None))
            rows.append(
                {
                    "chain.id": chain.chain_id,
                    "residue.index": res.seq_index,
                    "chain.type": chain.chain_type,
                    "chain.supertype": chain.chain_supertype,
                    "region.type": region_type,
                    "region.start": region_start,
                    "residue.aa": res.aa,
                }
            )
    schema = {
        "chain.id": pl.Utf8,
        "residue.index": pl.Int64,
        "chain.type": pl.Utf8,
        "chain.supertype": pl.Utf8,
        "region.type": pl.Utf8,
        "region.start": pl.Int64,
        "residue.aa": pl.Utf8,
    }
    return pl.DataFrame(rows, schema=schema) if rows else pl.DataFrame(schema=schema)


def symmetrize(contacts: pl.DataFrame) -> pl.DataFrame:
    """Return contacts plus their from/to-swapped mirror (R ``rbind`` semantics)."""
    swapped = contacts.rename(
        {
            "chain.id.from": "chain.id.to",
            "chain.id.to": "chain.id.from",
            "residue.index.from": "residue.index.to",
            "residue.index.to": "residue.index.from",
            "residue.aa.from": "residue.aa.to",
            "residue.aa.to": "residue.aa.from",
            "atom.from": "atom.to",
            "atom.to": "atom.from",
        }
    ).select(contacts.columns)
    return pl.concat([contacts, swapped])


def tidy_contacts(
    structure: Structure, cutoff: float = 5.0, count_atoms: bool = False
) -> pl.DataFrame:
    """Symmetrised, fully annotated contact table for a structure.

    Each inter-chain residue contact appears in both directions, with chain type,
    region type, region start and amino acid attached on both the ``from`` and ``to``
    sides — the input to :class:`tcren.contactmap.ContactMap`.

    When ``count_atoms`` is set, the ``n_atom_contacts`` per-residue-pair heavy-atom
    count is carried through (it is symmetric, so it survives the from/to swap
    unchanged). Default ``False`` keeps the table byte-identical to the legacy output.
    """
    contacts = symmetrize(
        all_atom_contacts(structure, cutoff=cutoff, count_atoms=count_atoms)
    )
    ann = residue_annotation(structure)

    from_ann = ann.rename(
        {
            "chain.id": "chain.id.from",
            "residue.index": "residue.index.from",
            "chain.type": "chain.type.from",
            "chain.supertype": "chain.supertype.from",
            "region.type": "region.type.from",
            "region.start": "region.start.from",
            "residue.aa": "residue.aa.from.ann",
        }
    )
    to_ann = ann.rename(
        {
            "chain.id": "chain.id.to",
            "residue.index": "residue.index.to",
            "chain.type": "chain.type.to",
            "chain.supertype": "chain.supertype.to",
            "region.type": "region.type.to",
            "region.start": "region.start.to",
            "residue.aa": "residue.aa.to.ann",
        }
    )
    out = (
        contacts.join(from_ann, on=["chain.id.from", "residue.index.from"], how="left")
        .join(to_ann, on=["chain.id.to", "residue.index.to"], how="left")
        .drop("residue.aa.from.ann", "residue.aa.to.ann")
    )
    return out
