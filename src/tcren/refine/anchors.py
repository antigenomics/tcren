"""Peptide anchor prediction (which residues bury into the MHC groove).

Class I uses fixed positions counted from both termini (P2 + the C-terminal PΩ); class II has no
fixed register, so a one-pass heuristic slides a 9-mer window and scores each register by a
P1-hydrophobic / P4,P6,P9-not-Pro/Gly rule (a cheap proxy for GibbsCluster/NNAlign register
inference), then reports P1/P4/P6/P9 of the best core. The anchor residues are the groove-facing
("presentation") side; the rest of the peptide is TCR-facing.

This logic is portable and depends only on the standard library. It is lifted from the antigenomics
``mhcmatch``/``seqtree`` decomposition primitives (which themselves carry no tcren dependency at
runtime); :func:`predict_anchors` adds a thin tcren wrapper that prefers the structure's own MHC-class
call (from :mod:`tcren.mhc`) over the length heuristic when a typed structure is available.
"""

from __future__ import annotations

from dataclasses import dataclass

from ..structure.model import PEPTIDE_TYPE, Structure

MASK = "X"
# Pocket-1 of class II strongly favours large hydrophobics (DR/DP and mouse alleles).
_MHC2_P1 = set("FILMVWY")
# Class-I anchors: 1-based P2 and the C-terminus (negative counts from the end, -1 == last).
_MHC1_ANCHORS = (2, -1)


def _resolve(anchors: tuple[int, ...], length: int) -> set[int]:
    """0-based anchor indices for a peptide of ``length`` from 1-based (neg = from C-term) positions."""
    out: set[int] = set()
    for a in anchors:
        idx = (a - 1) if a > 0 else (length + a)
        if 0 <= idx < length:
            out.add(idx)
    return out


def _core_anchor_score(core: str) -> float:
    """Crude register-likelihood of a 9-mer core: P1 hydrophobic dominates; P4/P6/P9 disfavour Pro/Gly."""
    s = 2.0 if core[0] in _MHC2_P1 else 0.0
    for i in (3, 5, 8):  # 0-based P4, P6, P9
        if core[i] not in "PG":
            s += 0.25
    return s


def _mhc2_core_anchors(peptide: str) -> tuple[int, ...]:
    """0-based P1/P4/P6/P9 indices of the best-scoring 9-mer register (one-pass register inference)."""
    if len(peptide) < 9:
        return ()
    best_s = max(range(len(peptide) - 8), key=lambda s: _core_anchor_score(peptide[s : s + 9]))
    return tuple(best_s + j for j in (0, 3, 5, 8))


def infer_class(peptide: str) -> str:
    """Heuristic MHC class from peptide length (≤ 11 → class I, else class II)."""
    return "MHCI" if len(peptide) <= 11 else "MHCII"


def anchor_indices(peptide: str, mhc_class: str) -> tuple[int, ...]:
    """0-based anchor positions: class-I P2/PΩ, class-II core P1/P4/P6/P9. ``mhc_class`` in {MHCI, MHCII}."""
    peptide = peptide.strip().upper()
    if mhc_class == "MHCII":
        return _mhc2_core_anchors(peptide)
    return tuple(sorted(_resolve(_MHC1_ANCHORS, len(peptide))))


@dataclass(frozen=True, slots=True)
class Decomposition:
    """Anchor vs TCR-facing split of a peptide."""

    peptide: str
    mhc_class: str
    anchors: tuple[int, ...]
    tcr_facing: str  # anchor positions masked to 'X'
    presentation: str  # TCR-facing positions masked to 'X'


def decompose(peptide: str, mhc_class: str | None = None) -> Decomposition:
    """Split ``peptide`` into anchor (groove-facing) and TCR-facing parts."""
    peptide = peptide.strip().upper()
    mhc_class = mhc_class or infer_class(peptide)
    anchors = set(anchor_indices(peptide, mhc_class))
    tcr = "".join(MASK if i in anchors else c for i, c in enumerate(peptide))
    present = "".join(c if i in anchors else MASK for i, c in enumerate(peptide))
    return Decomposition(peptide, mhc_class, tuple(sorted(anchors)), tcr, present)


def _structure_mhc_class(structure: Structure) -> str | None:
    """The MHC class a typed structure carries (``chain_supertype`` on its MHC chains), or None."""
    for chain in structure.chains:
        if chain.chain_type in ("MHCa", "MHCb") and chain.chain_supertype in ("MHCI", "MHCII"):
            return chain.chain_supertype
    return None


def predict_anchors(peptide: str, structure: Structure | None = None) -> Decomposition:
    """Predict peptide anchors, preferring the structure's MHC-class call over the length heuristic.

    Pass a chain-typed, MHC-annotated ``structure`` (see :func:`tcren.mhc.annotate_mhc`) to use its
    real class assignment; otherwise the class is inferred from peptide length.
    """
    mhc_class = _structure_mhc_class(structure) if structure is not None else None
    return decompose(peptide, mhc_class)


def native_peptide(structure: Structure) -> str:
    """One-letter sequence of the structure's PEPTIDE chain (raises if absent)."""
    for chain in structure.chains:
        if chain.chain_type == PEPTIDE_TYPE:
            return chain.sequence()
    raise ValueError(f"no peptide chain found in {structure.pdb_id!r}")
