"""Project canonical MHC groove regions onto a structure's MHC chains.

Each MHC chain is aligned (global, BLOSUM62) to the canonical chain for its class/role;
the canonical region positions are then mapped through the alignment onto the chain's
residues, producing :class:`~tcren.structure.model.RegionMarkup` entries (``HELIX_A1``,
``HELIX_A2``/``HELIX_B1``, ``GROOVE_FLOOR``) in the same schema as the TCR region markup.
"""

from __future__ import annotations

from functools import lru_cache

from ..structure.model import Chain, RegionMarkup, Structure
from . import domains
from .mapper import MhcCall, apply_mhc_calls, map_mhc


@lru_cache(maxsize=1)
def _aligner():
    from Bio.Align import PairwiseAligner, substitution_matrices

    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    return aligner


def _canonical_to_query(query_seq: str, canonical_seq: str) -> dict[int, int]:
    """Map each canonical position to the aligned query position (best global alignment)."""
    if not query_seq:
        return {}
    alignment = _aligner().align(query_seq, canonical_seq)[0]
    q_blocks, c_blocks = alignment.aligned  # parallel (start, end) block lists
    mapping: dict[int, int] = {}
    for (qs, qe), (cs, ce) in zip(q_blocks, c_blocks):
        for offset in range(qe - qs):
            mapping[cs + offset] = qs + offset
    return mapping


def partition_chain(chain: Chain, mhc_class: str, chain_role: str) -> list[RegionMarkup]:
    """Return groove RegionMarkups for one MHC chain (empty for B2M / unknown roles)."""
    groove = domains.groove_for(mhc_class, chain_role)
    if groove is None:
        return []
    mapping = _canonical_to_query(chain.sequence(), groove["sequence"])

    regions: list[RegionMarkup] = []
    for region_type, canonical_positions in groove["regions"].items():
        residues = [
            chain.residues[mapping[c]]
            for c in canonical_positions
            if c in mapping and mapping[c] < len(chain.residues)
        ]
        if not residues:
            continue
        residues.sort(key=lambda r: r.seq_index)
        regions.append(
            RegionMarkup(
                region_type=region_type,
                start_seq_index=residues[0].seq_index,
                end_seq_index=residues[-1].seq_index,
                sequence="".join(r.aa for r in residues),
                residues=residues,
            )
        )
    return regions


def partition_mhc(structure: Structure, calls: list[MhcCall]) -> None:
    """Assign groove regions to every mapped MHC chain in the structure (in place)."""
    by_id = {c.chain_id: c for c in calls}
    for chain in structure.chains:
        call = by_id.get(chain.chain_id)
        if call is not None:
            chain.regions = partition_chain(chain, call.mhc_class, call.chain_role)


def annotate_mhc(structure: Structure) -> list[MhcCall]:
    """Map and partition the MHC chains of an (already chain-typed) structure.

    Returns the :class:`MhcCall` list and, in place, sets each MHC chain's type
    (MHCa/MHCb/B2M), class supertype, allele and groove regions.
    """
    calls = map_mhc(structure)
    apply_mhc_calls(structure, calls)
    partition_mhc(structure, calls)
    return calls


def annotate_mhc_batch(structures: list[Structure], sensitivity: float = 5.7) -> None:
    """MHC-annotate many (chain-typed) structures with a SINGLE mmseqs search.

    Gathers every candidate MHC chain across all structures, runs one ``easy_search`` (mmseqs
    parallelises internally — no Python threads, no per-structure call), then slices the hits
    back and applies the calls + groove partitioning to each structure in place. This is the
    batched equivalent of calling :func:`annotate_mhc` per structure, for dataset-scale work.
    """
    import tempfile
    from pathlib import Path

    import arda.mmseqs as mmseqs

    from . import reference
    from .mapper import _best_hits, _candidate_chains, calls_from_hits

    flat = [(i, c) for i, s in enumerate(structures) for c in _candidate_chains(s) if c.sequence()]
    best: dict[str, dict] = {}
    if flat:
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            query_fa = tmp / "query.fasta"
            with query_fa.open("w") as fh:
                for i, c in flat:
                    fh.write(f">{i}|{c.chain_id}\n{c.sequence()}\n")
            out_tsv = tmp / "hits.tsv"
            mmseqs.easy_search(query_fa, reference.reference_db(), out_tsv,
                               tmp / "mmseqs_tmp", search_type=1, sensitivity=sensitivity,
                               max_seqs=50)
            best = _best_hits(out_tsv)
    for i, s in enumerate(structures):
        calls = calls_from_hits(_candidate_chains(s), best, key=lambda c, i=i: f"{i}|{c.chain_id}")
        apply_mhc_calls(s, calls)
        partition_mhc(s, calls)
