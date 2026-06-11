"""Residue-level contact map and interface partitioning.

A :class:`ContactMap` wraps the annotated, symmetrised contact table and exposes the
three biological interfaces (TCR↔peptide, TCR↔MHC, peptide↔MHC). The TCR↔peptide
interface is the central object for scoring and reproduces the schema of
``data/contact_maps_PDB.csv`` once chains and regions are annotated.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import polars as pl

from .contacts.table import tidy_contacts
from .structure.model import MHC_TYPES, PEPTIDE_TYPE, RECEPTOR_TYPES, Structure

Interface = Literal["tcr_peptide", "tcr_mhc", "peptide_mhc"]


@dataclass(slots=True)
class ContactMap:
    """Annotated, symmetrised residue contacts for one structure."""

    pdb_id: str
    contacts: pl.DataFrame
    peptide_length: int | None = None

    @classmethod
    def from_structure(cls, structure: Structure, cutoff: float = 5.0) -> "ContactMap":
        """Build a contact map from an (annotated) structure."""
        df = tidy_contacts(structure, cutoff=cutoff).with_columns(
            pl.lit(structure.pdb_id).alias("pdb.id")
        )
        peptide_length = next(
            (len(c.residues) for c in structure.chains if c.chain_type == PEPTIDE_TYPE),
            None,
        )
        return cls(pdb_id=structure.pdb_id, contacts=df, peptide_length=peptide_length)

    def _interface(self, from_types: tuple[str, ...], to_types: tuple[str, ...]) -> pl.DataFrame:
        sel = self.contacts.filter(
            pl.col("chain.type.from").is_in(list(from_types))
            & pl.col("chain.type.to").is_in(list(to_types))
        )
        # pos = residue.index - region.start (0-based position within a region).
        return sel.with_columns(
            (pl.col("residue.index.from") - pl.col("region.start.from")).alias("pos.from"),
            (pl.col("residue.index.to") - pl.col("region.start.to")).alias("pos.to"),
        )

    def interface(self, which: Interface) -> pl.DataFrame:
        """Return the contacts of one interface with within-region positions.

        Args:
            which: ``"tcr_peptide"``, ``"tcr_mhc"`` or ``"peptide_mhc"``.

        Returns:
            Filtered contacts with added ``pos.from``/``pos.to`` columns.
        """
        if which == "tcr_peptide":
            return self._interface(RECEPTOR_TYPES, (PEPTIDE_TYPE,))
        if which == "tcr_mhc":
            return self._interface(RECEPTOR_TYPES, MHC_TYPES)
        if which == "peptide_mhc":
            return self._interface((PEPTIDE_TYPE,), MHC_TYPES)
        raise ValueError(f"unknown interface {which!r}")

    def tcr_peptide(self) -> pl.DataFrame:
        """Convenience accessor for the TCR↔peptide interface."""
        return self.interface("tcr_peptide")

    def to_csv(self, path: str | Path) -> None:
        """Write the full annotated contact table to CSV."""
        self.contacts.write_csv(str(path))
