"""Percentile rank of a native peptide's energy against a random pMHC background.

For a given structure, scores the native peptide together with a background of
random (or epitope-sampled) peptides of the same length and reports where the
native score falls in that distribution. Lower TCRen energy = better binder, so a
small ``rank_pct`` means the native peptide scores at least as well as only a small
fraction of the background.
"""

from __future__ import annotations

import random
from collections.abc import Iterable
from pathlib import Path

from .contactmap import ContactMap, Interface
from .potential import Potential
from .scoring import score_peptides

#: The 20 standard amino acids, used to draw uniform-random background peptides.
_AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"


def _read_source(source: str | Path) -> list[str]:
    """Read epitope sequences from a FASTA or plain-text file (one per line)."""
    seqs: list[str] = []
    cur: list[str] = []
    for line in Path(source).read_text().splitlines():
        line = line.strip()
        if not line or line.lower() == "peptide":
            continue
        if line.startswith(">"):
            if cur:
                seqs.append("".join(cur))
                cur = []
            continue
        cur.append(line)
    if cur:
        seqs.append("".join(cur))
    return seqs


def background_peptides(
    length: int,
    n: int = 1000,
    seed: int = 0,
    source: str | None = None,
) -> list[str]:
    """Build a background set of ``n`` peptides of the given length.

    Args:
        length: Required peptide length.
        n: Number of background peptides to return.
        seed: Random seed for reproducibility.
        source: Optional FASTA/text file of epitopes. When given, peptides of the
            requested length are sampled from it (with replacement); each sampled
            sequence is also position-permuted so the background is a shuffled-epitope
            distribution rather than a verbatim copy. When ``None``, peptides are
            drawn uniformly at random over the 20 amino acids.

    Returns:
        A list of ``n`` upper-case peptide strings, each of length ``length``.
    """
    rng = random.Random(seed)
    if source is not None:
        pool = [s for s in _read_source(source) if len(s) == length]
        if not pool:
            raise ValueError(f"no length-{length} epitopes found in {source!r}")
        out = []
        for _ in range(n):
            chars = list(rng.choice(pool))
            rng.shuffle(chars)
            out.append("".join(chars))
        return out
    return ["".join(rng.choice(_AMINO_ACIDS) for _ in range(length)) for _ in range(n)]


def percentile_rank(
    contact_map: ContactMap,
    peptide: str,
    potential: Potential,
    *,
    interface: Interface = "tcr_peptide",
    n_background: int = 1000,
    seed: int = 0,
    tcr_regions: str = "all",
    background: Iterable[str] | None = None,
    contact_weight: str = "residue",
) -> dict:
    """Percentile rank of a peptide's energy against a random pMHC background.

    Scores ``peptide`` together with a background set (supplied or generated) and
    returns the fraction of background peptides whose score is ``<=`` the native
    score. Because lower TCRen energy means a better binder, a smaller ``rank_pct``
    indicates the native peptide is among the strongest binders.

    Args:
        contact_map: The structure's contact map.
        peptide: The native peptide to rank.
        potential: Pairwise potential to score with.
        interface: Which interface to score over (default ``"tcr_peptide"``).
        n_background: Size of the generated background (ignored if ``background`` is
            given).
        seed: Random seed for background generation.
        tcr_regions: Which TCR regions to keep on the TCR side (``"all"`` default,
            ``"cdr"`` or ``"cdr+fr"``); passed through to ``score_peptides``.
        background: Explicit background peptides; when ``None`` a uniform-random
            background of length ``len(peptide)`` is generated.
        contact_weight: ``"residue"`` (default) or ``"atomic"``; passed through to
            ``score_peptides`` (``"atomic"`` needs an ``n_atom_contacts`` contact map).

    Returns:
        Mapping with keys ``peptide``, ``score`` (native energy), ``rank_pct``
        (fraction of background with score ``<=`` native), and ``n_background``.
    """
    if background is None:
        bg = background_peptides(len(peptide), n=n_background, seed=seed)
    else:
        bg = list(background)

    scored = score_peptides(
        contact_map, [peptide] + bg, potential, interface=interface,
        tcr_regions=tcr_regions, contact_weight=contact_weight,
    )
    by_peptide = dict(zip(scored["peptide"], scored["score"]))
    if peptide not in by_peptide:
        raise ValueError(
            f"native peptide {peptide!r} was not scored "
            "(length mismatch with the structure's peptide?)"
        )
    native_score = by_peptide[peptide]

    # Background scores, looked up per generated peptide; entries dropped by the
    # length filter in score_peptides (e.g. a wrong-length supplied background) are
    # excluded from the denominator. A background peptide equal to the native maps
    # to the same score and counts as a tie (score <= native).
    bg_scores = [by_peptide[p] for p in bg if p in by_peptide]
    n_bg = len(bg_scores)
    if n_bg == 0:
        raise ValueError("no background peptides were scored")
    n_le = sum(1 for s in bg_scores if s <= native_score)
    rank_pct = n_le / n_bg

    return {
        "peptide": peptide,
        "score": native_score,
        "rank_pct": rank_pct,
        "n_background": n_bg,
    }
