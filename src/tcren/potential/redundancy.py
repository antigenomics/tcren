"""Non-redundancy clustering for TCRen derivation inputs.

The TCRen potential is derived from a *non-redundant* set of αβ TCR–pMHC structures:
near-duplicate complexes (same/similar CDR3α + CDR3β + peptide) are collapsed to a
single representative so the contact statistics are not dominated by repeated entries.

This is a verbatim lift of the ``nonredundant``/``alphabeta`` helpers from
``notebooks/natcompsci2022/01_nonred_and_derivation.ipynb``: a complete-linkage
hierarchical clustering on the summed Damerau–Levenshtein distance of the
CDR3α+CDR3β+peptide strings, cut at distance ``t`` (default 6.0), keeping one
representative per cluster.
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
import polars as pl
from rapidfuzz.distance import DamerauLevenshtein
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform


def nonredundant_ids(
    markup: pl.DataFrame,
    t: float | None = 6.0,
    fields: Sequence[str] = ("cdr3a", "cdr3b", "peptide"),
    linkage_method: str = "complete",
) -> list[str]:
    """Non-redundant representative ``pdb.id`` for each cluster of similar complexes.

    Clusters structures by the summed Damerau–Levenshtein distance over ``fields``
    (CDR3α + CDR3β + peptide) using hierarchical clustering, then returns one
    representative ``pdb.id`` per cluster (the first in row order).

    Args:
        markup: Per-structure table with a ``pdb.id`` column and the ``fields`` columns.
        t: Distance cutoff for ``fcluster`` (``criterion="distance"``). ``None`` turns
            redundancy filtering **off** and returns every ``pdb.id`` unchanged.
        fields: Sequence columns whose per-pair distances are summed (default
            ``cdr3a``, ``cdr3b``, ``peptide``).
        linkage_method: Linkage method for ``scipy.cluster.hierarchy.linkage``
            (default ``"complete"``).

    Returns:
        Representative ``pdb.id`` values, one per cluster, in first-seen order.
    """
    ids = markup["pdb.id"].to_list()
    if t is None:
        return ids

    mk = markup.with_columns([pl.col(c).fill_null("") for c in fields])
    seqs = {c: mk[c].to_list() for c in fields}
    n = len(ids)
    if n <= 1:
        return ids

    D = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            d = sum(DamerauLevenshtein.distance(seqs[c][i], seqs[c][j]) for c in fields)
            D[i, j] = D[j, i] = d
    cl = fcluster(linkage(squareform(D), method=linkage_method), t=t, criterion="distance")
    seen: set[int] = set()
    return [ids[i] for i in range(n) if not (cl[i] in seen or seen.add(cl[i]))]


def alphabeta_ids(contacts: pl.DataFrame) -> list[str]:
    """``pdb.id`` of complexes whose receptor contacts are exclusively TRA/TRB (αβ TCRs).

    Args:
        contacts: Contact table with ``pdb.id`` and ``chain.type.from`` columns.

    Returns:
        The ``pdb.id`` values whose set of ``chain.type.from`` is a subset of
        ``{"TRA", "TRB"}`` (i.e. no γδ chains).
    """
    cts = contacts.group_by("pdb.id").agg(pl.col("chain.type.from").unique().alias("x"))
    return [
        r["pdb.id"]
        for r in cts.iter_rows(named=True)
        if set(r["x"]) <= {"TRA", "TRB"}
    ]
