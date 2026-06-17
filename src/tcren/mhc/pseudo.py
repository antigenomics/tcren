"""MHC pseudosequence (MPS) annotation.

NetMHCpan defines a 34-residue "pseudosequence" per allele — the polymorphic groove positions
that contact the peptide (class I: α1/α2 of MHCa; class II: α1 of MHCa + β1 of MHCb). The
committed ``mhci_pseudo.fa`` / ``mhcii_pseudo.fa`` (see ``scripts/build_pseudo_fasta.py``) hold
the unique pseudosequences.

:func:`annotate_pseudo` adds an ``MPS`` region to a chain-typed + MHC-annotated structure, on
demand. The 34 pseudo positions are **scattered** along the chain (not a contiguous motif), so an
mmseqs/local search can't find them — there is no shared k-mer to seed on. Instead we thread each
candidate 34-mer through the chain with a **fitting alignment** (gaps in the chain are free, the
pseudosequence may not gap), which recovers the positions because NetMHCpan lists them N→C. The
best-scoring pseudosequence is chosen (one hit), and its identically-matched residues are marked —
across MHCa only for class I, split across MHCa+MHCb for class II, never β2m. Scoring all ~4k
pseudosequences this way is ~0.1 s, so no prebuilt index is needed.
"""

from __future__ import annotations

from functools import lru_cache
from importlib import resources
from pathlib import Path

from ..structure.model import Chain, RegionMarkup, Structure

try:  # compiled fitting-alignment hot path (built by scikit-build-core)
    from .. import _align
except ImportError:  # pragma: no cover - pure-Python fallback via Biopython
    _align = None

_PSEUDO_FA = {"MHCI": "mhci_pseudo.fa", "MHCII": "mhcii_pseudo.fa"}


@lru_cache(maxsize=1)
def _pseudo_aligner():
    """Fitting aligner: the pseudosequence is placed in full; chain (target) gaps are free."""
    from Bio.Align import PairwiseAligner, substitution_matrices

    a = PairwiseAligner()
    a.mode = "global"
    a.substitution_matrix = substitution_matrices.load("BLOSUM62")
    a.query_internal_open_gap_score = -11      # the 34-mer should not gap
    a.query_internal_extend_gap_score = -1
    a.target_internal_open_gap_score = 0       # chain gaps are free (fitting)
    a.target_internal_extend_gap_score = 0
    a.target_end_gap_score = 0
    a.query_end_gap_score = 0
    return a


@lru_cache(maxsize=2)
def _pseudo_index(mhc_class: str) -> dict[str, str]:
    """``header-id -> 34-mer`` for the bundled pseudosequence FASTA of a class."""
    fasta = Path(str(resources.files("tcren.data").joinpath(_PSEUDO_FA[mhc_class])))
    index: dict[str, str] = {}
    header = None
    for line in fasta.read_text().splitlines():
        if line.startswith(">"):
            header = line[1:].split()[0]
        elif header is not None:
            index[header] = line.strip()
    return index


@lru_cache(maxsize=2)
def _pseudo_lists(mhc_class: str) -> tuple[list[str], list[str]]:
    """Parallel ``(ids, sequences)`` for the class — the order the C++ ``best_hit`` indexes."""
    index = _pseudo_index(mhc_class)
    return list(index.keys()), list(index.values())


def _aligned_pairs(placed: str, free: str) -> list[tuple[int, int]]:
    """Matched ``(placed_pos, free_pos)`` columns — C++ ``_align`` or the Biopython fallback."""
    if _align is not None:
        return _align.align(placed, free)
    alignment = _pseudo_aligner().align(placed, free)[0]
    return [(ps + k, cs + k) for (ps, pe), (cs, _ce) in zip(*alignment.aligned)
            for k in range(pe - ps)]


def _best_pseudo_hit(query_seq: str, mhc_class: str) -> tuple[str, str] | None:
    """The single best-fitting ``(id, 34-mer)`` for ``query_seq`` (or ``None``)."""
    if not query_seq:
        return None
    ids, seqs = _pseudo_lists(mhc_class)
    if _align is not None:
        best, _score = _align.best_hit(query_seq, seqs)
    else:
        aligner = _pseudo_aligner()
        best = max(range(len(seqs)), key=lambda k: aligner.score(seqs[k], query_seq))
    return ids[best], seqs[best]


def _mark_concatenated(chains: list[Chain], pseudo_seq: str) -> dict[int, list]:
    """Thread the 34-mer through the concatenated chain sequence once; matched residues per chain.

    Aligning to the concatenation (rather than to each chain separately) keeps the α-half of a
    class-II pseudosequence on MHCa and the β-half on MHCb — no spurious cross-chain matches.
    Returns ``{chain_index: [residues]}``.
    """
    seqs = [c.sequence() for c in chains]
    concat = "".join(seqs)
    bounds = []  # (chain_index, start_offset) cumulative
    off = 0
    for ci, s in enumerate(seqs):
        bounds.append((ci, off))
        off += len(s)

    def _locate(pos: int) -> tuple[int, int]:
        ci, start = 0, 0
        for cand_ci, cand_start in bounds:
            if pos >= cand_start:
                ci, start = cand_ci, cand_start
        return ci, pos - start

    out: dict[int, list] = {}
    for p, cpos in _aligned_pairs(pseudo_seq, concat):
        if pseudo_seq[p] == "X" or cpos >= len(concat) or concat[cpos] != pseudo_seq[p]:
            continue
        ci, local = _locate(cpos)
        if local < len(chains[ci].residues):
            out.setdefault(ci, []).append(chains[ci].residues[local])
    return out


def annotate_pseudo(structure: Structure) -> str | None:
    """Add an ``MPS`` region to each groove chain from the best-matching pseudosequence.

    ``structure`` must already be chain-typed + MHC-annotated. Returns the chosen pseudosequence
    id (or ``None`` if there is no MHC). The best hit is selected once over the class groove
    sequence (MHCa for class I; MHCa+MHCb for class II) and its residues marked per chain.
    """
    cls = "MHCII" if any(c.chain_type == "MHCb" for c in structure.chains) else "MHCI"
    order = ("MHCa", "MHCb") if cls == "MHCII" else ("MHCa",)
    chains = [c for t in order for c in structure.chains if c.chain_type == t]
    if not chains:
        return None

    hit = _best_pseudo_hit("".join(c.sequence() for c in chains), cls)
    if hit is None:
        return None
    best_id, pseudo_seq = hit

    per_chain = _mark_concatenated(chains, pseudo_seq)
    for ci, chain in enumerate(chains):
        residues = sorted(per_chain.get(ci, []), key=lambda r: r.seq_index)
        chain.regions = [r for r in chain.regions if r.region_type != "MPS"]
        if residues:
            chain.regions.append(
                RegionMarkup(
                    region_type="MPS",
                    start_seq_index=residues[0].seq_index,
                    end_seq_index=residues[-1].seq_index,
                    sequence="".join(r.aa for r in residues),
                    residues=residues,
                )
            )
    return best_id
