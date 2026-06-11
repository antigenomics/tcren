"""αβ vs γδ T-cell classification from the TCR constant region (C-gene).

arda annotates the variable V(D)J region but not the constant domain. When a structure
includes an ordered constant domain, aligning each chain to the TCR constant references
(TRAC/TRBC1/TRBC2 → αβ; TRGC/TRDC → γδ) identifies the chain (α/β/γ/δ) unambiguously and
therefore the cell type. This is authoritative for αβ-vs-γδ and independent of the
(occasionally ambiguous, e.g. TRAV/DV) V-gene call. Variable-domain-only chains carry no
constant region and yield no call (cell type ``"unknown"``).
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from importlib import resources

#: Minimum local-alignment score to accept a constant-domain match. The V domain alone
#: scores ~30-45 against any constant; a real constant domain scores in the hundreds.
MIN_CONSTANT_SCORE = 80.0


@dataclass(slots=True)
class ConstantCall:
    """A TCR constant-region identification for one chain."""

    chain_id: str
    gene: str  # TRAC, TRBC1, TRBC2, TRGC1, TRGC2, TRDC
    chain_class: str  # alpha | beta | gamma | delta
    cell_type: str  # ab | gd
    score: float


@lru_cache(maxsize=1)
def _references() -> list[tuple[str, str, str, str]]:
    """Return ``(gene, chain_class, cell_type, sequence)`` constant references."""
    text = resources.files("tcren.data").joinpath("tcr_constant.fasta").read_text()
    out, header, seq = [], None, []
    for line in text.splitlines():
        if line.startswith(">"):
            if header is not None:
                gene, chain_class, cell_type, _species = header.split("|")
                out.append((gene, chain_class, cell_type, "".join(seq)))
            header, seq = line[1:], []
        else:
            seq.append(line.strip())
    if header is not None:
        gene, chain_class, cell_type, _species = header.split("|")
        out.append((gene, chain_class, cell_type, "".join(seq)))
    return out


@lru_cache(maxsize=1)
def _aligner():
    from Bio.Align import PairwiseAligner, substitution_matrices

    aligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    return aligner


def classify_chain_constant(
    sequence: str, min_score: float = MIN_CONSTANT_SCORE
) -> ConstantCall | None:
    """Identify the constant region of a single chain sequence, if one is present."""
    if not sequence:
        return None
    aligner = _aligner()
    best = None
    for gene, chain_class, cell_type, ref_seq in _references():
        score = aligner.score(sequence, ref_seq)
        if best is None or score > best[0]:
            best = (score, gene, chain_class, cell_type)
    if best is None or best[0] < min_score:
        return None
    score, gene, chain_class, cell_type = best
    return ConstantCall(chain_id="", gene=gene, chain_class=chain_class,
                        cell_type=cell_type, score=score)


def classify_constants(structure, min_score: float = MIN_CONSTANT_SCORE) -> list[ConstantCall]:
    """Identify the constant region of every chain that has one."""
    calls = []
    for chain in structure.chains:
        call = classify_chain_constant(chain.sequence(), min_score=min_score)
        if call is not None:
            calls.append(
                ConstantCall(chain.chain_id, call.gene, call.chain_class, call.cell_type, call.score)
            )
    return calls


def cell_type(structure, min_score: float = MIN_CONSTANT_SCORE) -> str:
    """Return ``"ab"``, ``"gd"`` or ``"unknown"`` from the constant regions present.

    γδ wins if any γ/δ constant is found; otherwise αβ if any α/β constant is found;
    otherwise ``"unknown"`` (no ordered constant domain — e.g. variable-only chains).
    """
    types = {c.cell_type for c in classify_constants(structure, min_score=min_score)}
    if "gd" in types:
        return "gd"
    if "ab" in types:
        return "ab"
    return "unknown"
