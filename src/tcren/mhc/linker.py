"""Detect and split covalently linked (single-chain) peptides via MHC alignment.

Engineered single-chain pMHC constructs fuse the peptide to an MHC chain through a
flexible (usually Gly/Ser-rich) linker, so the peptide is not a separate chain. Aligning
each chain to the MHC reference reveals this: the MHC domain aligns, leaving an unaligned
terminal segment that — after stripping the linker — is the peptide. This module provides
the alignment check and a splitter that lifts such peptides into their own chain.

No covalently linked peptides occur in the bundled TCR3D / PDB datasets (all conventional,
separate-chain complexes); this is robustness for engineered and predicted structures.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from functools import lru_cache
from importlib import resources

from ..structure.model import Chain, RegionMarkup, Structure

_LINKER_RESIDUES = set("GS")  # flexible linker alphabet


@lru_cache(maxsize=1)
def _mhc_reference_sequences() -> dict[str, str]:
    canon = json.loads(resources.files("tcren.data").joinpath("mhc_canonical.json").read_text())
    return {key: entry["sequence"] for key, entry in canon.items()}


@lru_cache(maxsize=1)
def _aligner():
    from Bio.Align import PairwiseAligner, substitution_matrices

    aligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    return aligner


@dataclass(slots=True)
class MhcAlignmentCheck:
    """Result of aligning a chain onto the MHC reference (a chain-identity check)."""

    chain_id: str
    best_ref: str  # e.g. "MHCI|MHCa"
    score: float
    query_start: int  # first aligned residue index in the chain
    query_end: int  # one past the last aligned residue index
    n_term_extra: int  # unaligned residues before the MHC domain
    c_term_extra: int  # unaligned residues after the MHC domain

    @property
    def is_mhc(self) -> bool:
        return self.score >= 200.0


def check_against_mhc(chain: Chain) -> MhcAlignmentCheck:
    """Align a chain onto the MHC reference and report coverage / terminal extensions."""
    seq = chain.sequence()
    aligner = _aligner()
    best = None
    for key, ref in _mhc_reference_sequences().items():
        aln = aligner.align(seq, ref)[0]
        if best is None or aln.score > best[0]:
            best = (aln.score, key, aln)
    score, key, aln = best
    qs = aln.aligned[0][0][0]
    qe = aln.aligned[0][-1][1]
    return MhcAlignmentCheck(
        chain_id=chain.chain_id,
        best_ref=key,
        score=float(score),
        query_start=qs,
        query_end=qe,
        n_term_extra=qs,
        c_term_extra=len(seq) - qe,
    )


def _strip_linker(residues: list, from_mhc_side: str) -> list:
    """Drop a Gly/Ser-rich linker run adjacent to the MHC domain.

    Args:
        residues: The unaligned terminal residues (in chain order).
        from_mhc_side: ``"n"`` if the MHC domain is C-terminal to this segment (linker is
            at the segment's end), ``"c"`` if the MHC domain is N-terminal (linker first).
    """
    seq = [r.aa for r in residues]
    if from_mhc_side == "n":
        # peptide ... linker | MHC  -> trim trailing linker run
        end = len(seq)
        while end > 0 and seq[end - 1] in _LINKER_RESIDUES:
            end -= 1
        return residues[:end]
    # MHC | linker ... peptide -> trim leading linker run
    start = 0
    while start < len(seq) and seq[start] in _LINKER_RESIDUES:
        start += 1
    return residues[start:]


def detect_linked_peptide(
    chain: Chain, min_len: int = 7, max_len: int = 25
) -> list | None:
    """Return the residues of a peptide fused to an MHC chain, or ``None``.

    Looks for a peptide-length segment (after stripping an adjacent Gly/Ser linker) at the
    N- or C-terminus of a chain whose remainder aligns to the MHC reference.
    """
    check = check_against_mhc(chain)
    if not check.is_mhc:
        return None
    candidates = []
    if check.n_term_extra >= min_len:
        seg = _strip_linker(chain.residues[: check.query_start], from_mhc_side="n")
        candidates.append(seg)
    if check.c_term_extra >= min_len:
        seg = _strip_linker(chain.residues[check.query_end :], from_mhc_side="c")
        candidates.append(seg)
    for seg in candidates:
        if min_len <= len(seg) <= max_len:
            return seg
    return None


def split_linked_peptides(structure: Structure, peptide_chain_id: str = "p") -> list[str]:
    """Split covalently linked peptides off their MHC chains, in place.

    For each chain carrying a fused peptide, the peptide residues are removed and added as
    a new PEPTIDE chain. Returns the list of chain ids that were split (empty if none).
    """
    split_from = []
    new_chains = []
    suffix = 0
    for chain in list(structure.chains):
        if chain.chain_type not in ("MHCa", "MHCb", "MHC", None):
            continue
        seg = detect_linked_peptide(chain)
        if not seg:
            continue
        # Remove the entire terminal extension (peptide + linker) from the MHC chain,
        # leaving the MHC domain; the new peptide chain holds the linker-stripped segment.
        check = check_against_mhc(chain)
        if check.n_term_extra >= check.c_term_extra:
            chain.residues = chain.residues[check.query_start :]
        else:
            chain.residues = chain.residues[: check.query_end]
        new_id = peptide_chain_id if suffix == 0 else f"{peptide_chain_id}{suffix}"
        suffix += 1
        pep = Chain(chain_id=new_id, residues=list(seg), chain_type="PEPTIDE",
                    chain_supertype="PEPTIDE")
        pep.regions = [
            RegionMarkup("PEPTIDE", seg[0].seq_index, seg[-1].seq_index,
                         "".join(r.aa for r in seg), list(seg))
        ]
        new_chains.append(pep)
        split_from.append(chain.chain_id)
    structure.chains.extend(new_chains)
    return split_from
