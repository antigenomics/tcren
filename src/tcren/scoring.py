"""Candidate-peptide scoring by amino-acid substitution.

Ports the second half of ``run_TCRen.R``: for each candidate peptide, substitute its
amino acids at the contacted peptide positions of a structure's contact map and sum the
pairwise potential over all contacts. Lower scores indicate more favourable interactions.
"""

from __future__ import annotations

from collections.abc import Iterable

import numpy as np
import polars as pl

from .contactmap import ContactMap, Interface
from .potential import Potential

# Which side of each interface carries the (substituted) peptide.
_PEPTIDE_SIDE: dict[str, str] = {
    "tcr_peptide": "to",
    "tcr_mhc": "to",  # substitutes the MHC side; peptide is fixed
    "peptide_mhc": "from",
}


def score_peptides(
    contact_map: ContactMap,
    candidates: Iterable[str],
    potential: Potential,
    interface: Interface = "tcr_peptide",
    require_same_length: bool = True,
    substituted_side: str | None = None,
    tcr_regions: str = "all",
) -> pl.DataFrame:
    """Score candidate peptides against a structure's contact map.

    Args:
        contact_map: The structure's contact map.
        candidates: Candidate peptide sequences (one-letter).
        potential: Pairwise potential to score with.
        interface: Which interface to score over (default ``"tcr_peptide"``).
        require_same_length: Only score candidates whose length matches the structure's
            peptide length (mirrors the legacy length join). Ignored when the contact
            map has no recorded peptide length.
        substituted_side: ``"to"`` or ``"from"`` — which contact side the candidate is
            threaded onto. Defaults to the peptide side of ``interface``.
        tcr_regions: which TCR regions to keep on the TCR side (``"all"`` default = no
            filter = legacy behaviour; ``"cdr"`` or ``"cdr+fr"`` to restrict).

    Returns:
        Columns ``complex.id``, ``peptide``, ``potential``, ``score`` sorted by
        ``complex.id`` then ascending ``score``.
    """
    side = substituted_side or _PEPTIDE_SIDE[interface]
    if side not in ("to", "from"):
        raise ValueError(f"substituted_side must be 'to' or 'from', got {side!r}")
    fixed = "from" if side == "to" else "to"

    iface = contact_map.interface(interface, tcr_regions=tcr_regions)
    matrix, index = potential.as_matrix()

    pos = np.asarray(iface[f"pos.{side}"].to_list(), dtype=np.int64)
    fixed_aa = iface[f"residue.aa.{fixed}"].to_list()
    fixed_idx = np.array([index.get(a, -1) for a in fixed_aa], dtype=np.int64)

    candidates = list(candidates)
    rows = []
    for peptide in candidates:
        if require_same_length and contact_map.peptide_length is not None:
            if len(peptide) != contact_map.peptide_length:
                continue
        # Gather the substituted amino acid for each contact from the candidate.
        subst_idx = np.array(
            [index.get(peptide[p], -1) if 0 <= p < len(peptide) else -1 for p in pos],
            dtype=np.int64,
        )
        if side == "to":
            rows_idx, cols_idx = fixed_idx, subst_idx
        else:
            rows_idx, cols_idx = subst_idx, fixed_idx
        valid = (rows_idx >= 0) & (cols_idx >= 0)
        vals = matrix[rows_idx[valid], cols_idx[valid]]
        # Pairs absent from the potential (e.g. Cys on the 'from' axis) are dropped,
        # exactly as the inner join in run_TCRen.R drops unmatched rows.
        score = float(np.nansum(vals))
        rows.append({"complex.id": contact_map.pdb_id, "peptide": peptide, "score": score})

    out = pl.DataFrame(
        rows,
        schema={"complex.id": pl.Utf8, "peptide": pl.Utf8, "score": pl.Float64},
    ).with_columns(pl.lit(potential.name).alias("potential"))
    return out.select("complex.id", "peptide", "potential", "score").sort(
        "complex.id", "score"
    )


def score_structures(
    contact_maps: Iterable[ContactMap],
    candidates: Iterable[str],
    potential: Potential,
    **kwargs,
) -> pl.DataFrame:
    """Score candidates against several structures and stack the results."""
    candidates = list(candidates)
    frames = [score_peptides(cm, candidates, potential, **kwargs) for cm in contact_maps]
    return pl.concat(frames) if frames else pl.DataFrame()
