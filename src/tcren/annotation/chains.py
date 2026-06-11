"""Chain classification: TRA/TRB (via arda), PEPTIDE, and (provisional) MHC.

Precise MHC sub-typing (MHCa/MHCb/B2M, class I/II) is added in Phase B; here MHC chains
are left with the generic type ``"MHC"`` so the TCR↔peptide scoring path is complete.
"""

from __future__ import annotations

from ..structure.model import RECEPTOR_TYPES, RegionMarkup, Structure
from .arda_adapter import annotate_tcr_chains

_SPECIES = {"human": "Human", "mouse": "Mouse"}


def _tag_peptide(chain) -> None:
    """Mark a chain as PEPTIDE and give it a single full-length PEPTIDE region."""
    chain.chain_type = "PEPTIDE"
    chain.chain_supertype = "PEPTIDE"
    chain.regions = [
        RegionMarkup(
            region_type="PEPTIDE",
            start_seq_index=chain.residues[0].seq_index,
            end_seq_index=chain.residues[-1].seq_index,
            sequence=chain.sequence(),
            residues=list(chain.residues),
        )
    ]


def classify_chains(
    structure: Structure,
    organism: str = "human",
    peptide_max_len: int = 30,
    autodetect_species: bool = True,
) -> None:
    """Classify every chain of ``structure`` in place.

    TRA/TRB are assigned from arda's locus call; the shortest remaining chains (length
    ≤ ``peptide_max_len``) become PEPTIDE; longer remaining chains are tagged ``"MHC"``.

    Args:
        structure: Structure to annotate (mutated in place).
        organism: arda organism (``"human"``/``"mouse"``).
        peptide_max_len: Maximum residue count for a chain to be called PEPTIDE.
        autodetect_species: If arda finds no TCR chain with ``organism``, retry with the
            other supported organism and keep whichever yields TCR chains.
    """
    tcr_ids = annotate_tcr_chains(structure, organism=organism)
    used_organism = organism
    if not tcr_ids and autodetect_species:
        alt = "mouse" if organism == "human" else "human"
        for chain in structure.chains:  # reset before retry
            chain.chain_type = chain.chain_supertype = chain.allele_info = None
            chain.regions = []
        tcr_ids = annotate_tcr_chains(structure, organism=alt)
        if tcr_ids:
            used_organism = alt

    structure.complex_species = _SPECIES.get(used_organism)

    remaining = [c for c in structure.chains if c.chain_type not in RECEPTOR_TYPES]
    for chain in remaining:
        if len(chain.residues) <= peptide_max_len:
            _tag_peptide(chain)
        else:
            chain.chain_type = "MHC"
