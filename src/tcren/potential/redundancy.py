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


def _cluster(
    markup: pl.DataFrame,
    t: float,
    fields: Sequence[str],
    linkage_method: str,
) -> tuple[list[str], np.ndarray]:
    """Cluster structures by summed Damerau–Levenshtein distance over ``fields``.

    Returns ``(ids, labels)`` where ``labels[i]`` is the cluster id of ``ids[i]``.
    For ``n <= 1`` every structure forms its own singleton cluster.
    """
    ids = markup["pdb.id"].to_list()
    n = len(ids)
    if n <= 1:
        return ids, np.arange(n, dtype=int)
    mk = markup.with_columns([pl.col(c).fill_null("") for c in fields])
    seqs = {c: mk[c].to_list() for c in fields}
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            d = sum(DamerauLevenshtein.distance(seqs[c][i], seqs[c][j]) for c in fields)
            D[i, j] = D[j, i] = d
    cl = fcluster(linkage(squareform(D), method=linkage_method), t=t, criterion="distance")
    return ids, cl


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
    if t is None or len(ids) <= 1:
        return ids
    ids, cl = _cluster(markup, t, fields, linkage_method)
    seen: set[int] = set()
    return [ids[i] for i in range(len(ids)) if not (cl[i] in seen or seen.add(cl[i]))]


def cluster_weights(
    markup: pl.DataFrame,
    t: float = 6.0,
    fields: Sequence[str] = ("cdr3a", "cdr3b", "peptide"),
    linkage_method: str = "complete",
) -> dict[str, float]:
    """Inverse-cluster-size (Henikoff-style) weight for each ``pdb.id``.

    Clusters structures exactly as :func:`nonredundant_ids` (same distance, same
    linkage), then assigns every structure the weight ``1 / cluster_size``. A unique
    structure gets weight ``1.0``; each member of a size-``k`` redundancy cluster gets
    ``1/k``, so the cluster contributes a total weight of ``1`` to the derivation.
    Feed the result to :func:`tcren.potential.derive.derive_tcren`'s ``weights`` argument
    to down-weight redundancy while keeping every structure's data.

    Args:
        markup: Per-structure table with a ``pdb.id`` column and the ``fields`` columns.
        t: Distance cutoff for ``fcluster`` (``criterion="distance"``).
        fields: Sequence columns whose per-pair distances are summed.
        linkage_method: Linkage method for ``scipy.cluster.hierarchy.linkage``.

    Returns:
        Mapping ``{pdb.id: 1 / cluster_size}``.
    """
    ids, cl = _cluster(markup, t, fields, linkage_method)
    sizes: dict[int, int] = {}
    for label in cl:
        sizes[label] = sizes.get(label, 0) + 1
    return {ids[i]: 1.0 / sizes[cl[i]] for i in range(len(ids))}


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
