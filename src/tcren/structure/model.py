"""Lightweight structure data model.

These dataclasses wrap the parsed contents of a PDB/mmCIF file in the shape the rest of
the pipeline needs: per-chain residue lists carrying both a 0-based sequential index
(matching the legacy ``mir`` ``residue.index``) and the original author numbering, plus
heavy-atom coordinates for contact computation.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

# Chain type / supertype vocabulary mirrors the legacy general.txt.
TCR_TYPES = ("TRA", "TRB", "TRD", "TRG")
# Some complexes in the dataset use a TCR-mimic antibody (Fab) as the receptor; the
# legacy mir lumped these with the TCR. We treat all antigen-receptor loci uniformly.
BCR_TYPES = ("IGH", "IGK", "IGL")
RECEPTOR_TYPES = TCR_TYPES + BCR_TYPES
PEPTIDE_TYPE = "PEPTIDE"
MHC_TYPES = ("MHCa", "MHCb", "B2M")


@dataclass(frozen=True, slots=True)
class Atom:
    """A single (heavy) atom."""

    name: str
    element: str
    coord: np.ndarray  # float64, shape (3,)


@dataclass(frozen=True, slots=True)
class Residue:
    """A polymer residue.

    Attributes:
        seq_index: 0-based sequential index over the chain's resolved polymer residues
            (the legacy ``residue.index``); independent of author numbering gaps.
        pdb_index: Author residue number (``residue.index.pdb``).
        insertion_code: Author insertion code (``''`` when absent).
        aa: One-letter amino-acid code (``'X'`` for unknown).
        resname: Three-letter residue name (``HIS``, ``MSE`` …).
        atoms: Heavy atoms of the residue.
    """

    seq_index: int
    pdb_index: int
    insertion_code: str
    aa: str
    resname: str
    atoms: tuple[Atom, ...]

    @property
    def ca(self) -> np.ndarray | None:
        """Cα coordinate, or ``None`` if the residue has no Cα atom."""
        for atom in self.atoms:
            if atom.name == "CA":
                return atom.coord
        return None


@dataclass(slots=True)
class Chain:
    """A polymer chain and its annotations."""

    chain_id: str
    residues: list[Residue]
    chain_type: str | None = None  # TRA/TRB/PEPTIDE/MHCa/MHCb/B2M
    chain_supertype: str | None = None  # TRAB/MHCI/MHCII/PEPTIDE
    allele_info: str | None = None
    regions: list["RegionMarkup"] = field(default_factory=list)

    def sequence(self) -> str:
        """One-letter sequence in residue order."""
        return "".join(r.aa for r in self.residues)

    def by_seq_index(self, seq_index: int) -> Residue | None:
        """Return the residue at a given sequential index, or ``None``."""
        # residues are appended in seq_index order, so direct indexing is valid.
        if 0 <= seq_index < len(self.residues):
            r = self.residues[seq_index]
            if r.seq_index == seq_index:
                return r
        for r in self.residues:
            if r.seq_index == seq_index:
                return r
        return None


@dataclass(slots=True)
class RegionMarkup:
    """An annotated region (CDR/FR for TCR, groove regions for MHC)."""

    region_type: str
    start_seq_index: int
    end_seq_index: int
    sequence: str
    residues: list[Residue]


@dataclass(slots=True)
class Structure:
    """A parsed complex: a set of annotated chains."""

    pdb_id: str
    chains: list[Chain]
    complex_species: str | None = None
    cell_type: str | None = None  # "ab" | "gd" | "unknown" (from the TCR constant region)

    def chain(self, chain_id: str) -> Chain:
        """Return the chain with the given id (raises ``KeyError`` if absent)."""
        for c in self.chains:
            if c.chain_id == chain_id:
                return c
        raise KeyError(f"chain {chain_id!r} not in structure {self.pdb_id!r}")

    def by_type(self, *types: str) -> list[Chain]:
        """Return chains whose ``chain_type`` is in ``types``."""
        return [c for c in self.chains if c.chain_type in types]
