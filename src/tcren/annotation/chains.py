"""Chain classification: TRA/TRB (via arda), PEPTIDE, and (provisional) MHC.

Precise MHC sub-typing (MHCa/MHCb/B2M, class I/II) is added in Phase B; here MHC chains
are left with the generic type ``"MHC"`` so the TCR↔peptide scoring path is complete.
"""

from __future__ import annotations

from ..structure.model import RECEPTOR_TYPES, RegionMarkup, Structure
from .arda_adapter import annotate_chains, apply_records, score_records

_SPECIES = {"human": "Human", "mouse": "Mouse"}


def _reset_receptor_chains(structure: Structure) -> None:
    for chain in structure.chains:
        if chain.chain_type in RECEPTOR_TYPES:
            chain.chain_type = chain.chain_supertype = chain.allele_info = None
            chain.regions = []


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


def _records_for(structure, organism, precomputed):
    """Receptor records for ``organism``: re-apply a precomputed map, else call arda."""
    if precomputed is not None and organism in precomputed:
        by_id = precomputed[organism]
        apply_records(structure.chains, by_id)
        return by_id
    return annotate_chains(structure.chains, organism)


def classify_chains(
    structure: Structure,
    organism: str = "human",
    peptide_max_len: int = 30,
    autodetect_species: bool = True,
    precomputed_records: dict[str, dict[str, dict]] | None = None,
) -> None:
    """Classify every chain of ``structure`` in place.

    TRA/TRB are assigned from arda's locus call; the shortest remaining chains (length
    ≤ ``peptide_max_len``) become PEPTIDE; longer remaining chains are tagged ``"MHC"``.

    Args:
        structure: Structure to annotate (mutated in place).
        organism: arda organism (``"human"``/``"mouse"``).
        peptide_max_len: Maximum residue count for a chain to be called PEPTIDE.
        autodetect_species: Annotate against both supported species (human and mouse)
            and keep whichever gives the higher total mmseqs alignment score over the
            receptor chains. TCR/BCR germlines are organism-specific, so the wrong species
            scores measurably lower (e.g. mouse BM3.3 scores ~435 vs ~197 under human);
            this avoids mis-typing a chain under the wrong reference. Ties keep the
            requested ``organism``. Disable to force ``organism`` verbatim.
        precomputed_records: Optional ``{organism: {chain_id: record}}`` of arda records
            for this structure's chains, to reuse instead of calling arda (the batch path
            in :func:`~tcren.paper.helpers.annotate_structure_set` annotates the whole
            dataset in one mmseqs call per organism and injects the per-structure slices).
    """
    used_organism = organism
    req_records = _records_for(structure, organism, precomputed_records)
    _ids, best_score = score_records(structure.chains, req_records)

    # Compare the requested organism against the other supported species and keep the
    # better-scoring annotation. Cached records are re-applied (no extra arda call) when
    # the requested organism wins, so this costs at most two mmseqs invocations.
    other = "mouse" if organism == "human" else "human"
    if autodetect_species and organism in _SPECIES and other in _SPECIES:
        _reset_receptor_chains(structure)
        _alt_ids, alt_score = score_records(
            structure.chains, _records_for(structure, other, precomputed_records)
        )
        if alt_score > best_score:  # strict: ties favour the requested organism
            used_organism = other
        else:  # restore the requested-organism annotation from the cache
            _reset_receptor_chains(structure)
            apply_records(structure.chains, req_records)

    structure.complex_species = _SPECIES.get(used_organism)

    remaining = [c for c in structure.chains if c.chain_type not in RECEPTOR_TYPES]
    for chain in remaining:
        if len(chain.residues) <= peptide_max_len:
            _tag_peptide(chain)
        else:
            chain.chain_type = "MHC"

    # Single-chain pMHC: the peptide can be fused into an MHC chain (engineered constructs,
    # common in class II). If no separate peptide chain was found, split any fused peptide
    # off its MHC chain so the TCR↔peptide interface is recoverable.
    if not any(c.chain_type == "PEPTIDE" for c in structure.chains):
        from ..mhc.linker import split_linked_peptides

        split_linked_peptides(structure)
