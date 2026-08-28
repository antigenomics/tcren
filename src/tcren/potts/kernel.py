"""Coupling classes: which pairs of sites interact, and the machinery the sampler needs.

Contacts on a real interface are not independent. Three coordinate-free families capture it:

``K(di, dj)``
    **Within-loop**, ``|di| <= 2`` and ``|dj| <= 2`` excluding ``(0, 0)``, both sites in one
    receptor loop and one partner chain. 12 classes after the ``K(d) = K(-d)`` symmetry.
``L(|dj|, same chain?)``
    **Cross-loop**, the two sites in *different* hypervariable loops at partner offset ``|dj|``.
    6 classes. Two CDR loops converging on one partner residue.
``M``
    **Cross-class**, the same receptor residue against both the peptide and the MHC groove.
    1 class, joint models only.

Anything else is asserted uncoupled, which is the model's statement that linkage falls to zero
across loops and beyond two residues.
"""

from __future__ import annotations

import numpy as np
import polars as pl

from .model import CDR_LOOPS, CROSS_DJ, OFFSETS


def _edges_within(q: pl.DataFrame) -> list[np.ndarray]:
    """``(sid_a, sid_b)`` with ``b = a + (di, dj)`` inside one receptor loop and partner chain."""
    key = ["pdb.id", "loop", "pos.rec", "pchain", "pos.par"]
    base = q.select(["sid"] + key)
    out = []
    for di, dj in OFFSETS:
        shifted = base.select(
            pl.col("sid").alias("sid_b"), "pdb.id", "loop",
            (pl.col("pos.rec") - di).alias("pos.rec"), "pchain",
            (pl.col("pos.par") - dj).alias("pos.par"))
        j = base.join(shifted, on=key, how="inner")
        out.append(np.column_stack([j["sid"].to_numpy(), j["sid_b"].to_numpy()]).astype(np.int64))
    return out


def _edges_cross(q: pl.DataFrame) -> list[np.ndarray]:
    """Two sites in DIFFERENT receptor loops, same partner chain, ``|dj| <= 2``.

    ``dj = 0`` is symmetric and is deduped by ``sid_a < sid_b``; ``dj > 0`` already names each
    unordered pair exactly once, because its two members sit at different partner positions.
    """
    key = ["pdb.id", "pchain", "pos.par"]
    base = q.filter(pl.col("loop").is_in(list(CDR_LOOPS))).select(
        ["sid", "loop", "chain.rec"] + key)
    out = []
    for dj in CROSS_DJ:
        shifted = base.select(
            pl.col("sid").alias("sid_b"), pl.col("loop").alias("loop_b"),
            pl.col("chain.rec").alias("chain_b"), "pdb.id", "pchain",
            (pl.col("pos.par") - dj).alias("pos.par"))
        j = base.join(shifted, on=key, how="inner").filter(pl.col("loop") != pl.col("loop_b"))
        if dj == 0:
            j = j.filter(pl.col("sid") < pl.col("sid_b"))
        for same in (True, False):
            k = j.filter((pl.col("chain.rec") == pl.col("chain_b")) == same)
            out.append(np.column_stack([k["sid"].to_numpy(),
                                        k["sid_b"].to_numpy()]).astype(np.int64))
    return out


def _edges_class(q: pl.DataFrame) -> list[np.ndarray]:
    """The same receptor residue against two different partner classes."""
    key = ["pdb.id", "loop", "pos.rec"]
    base = q.select(["sid", "cls"] + key)
    j = (base.join(base.select(pl.col("sid").alias("sid_b"), pl.col("cls").alias("cls_b"), *key),
                   on=key, how="inner")
         .filter((pl.col("cls") != pl.col("cls_b")) & (pl.col("sid") < pl.col("sid_b"))))
    return [np.column_stack([j["sid"].to_numpy(), j["sid_b"].to_numpy()]).astype(np.int64)]


def edges(q: pl.DataFrame, joint: bool = False) -> list[np.ndarray]:
    """One ``(n, 2)`` array of site-index pairs per coupling class, in :func:`kernel_names` order.

    Args:
        q: The frame :func:`tcren.potts.site_codes` returns, carrying ``sid``, ``loop`` and
            ``pchain``.
        joint: Include the cross-class family (both partners in one model).
    """
    out = _edges_within(q) + _edges_cross(q)
    return out + _edges_class(q) if joint else out


def neighbour_counts(edge_lists: list[np.ndarray], sigma: np.ndarray, n_sites: int) -> np.ndarray:
    """``n_k(a)`` — contacting neighbours of each site in each class, both directions of each edge.

    These are the covariates of the pseudolikelihood: the conditional
    ``P(sigma_a = 1 | sigma_{-a})`` is logistic in ``eta_a + sum_k K_k n_k(a)``.
    """
    out = np.zeros((n_sites, len(edge_lists)))
    for k, e in enumerate(edge_lists):
        if not e.size:
            continue
        np.add.at(out[:, k], e[:, 0], sigma[e[:, 1]])
        np.add.at(out[:, k], e[:, 1], sigma[e[:, 0]])
    return out


def bucket_edges(edge_lists: list[np.ndarray], starts: np.ndarray, kv: np.ndarray):
    """Group every non-zero-coefficient edge by structure. Returns ``(E, offsets)``.

    Site indices are contiguous per structure, so one ``searchsorted`` assigns each edge to its
    structure. ``E`` has columns ``(sid_a, sid_b, class)``.
    """
    keep = [(k, e) for k, e in enumerate(edge_lists) if e.size and kv[k] != 0.0]
    E = (np.concatenate([np.column_stack([e, np.full(len(e), k)]) for k, e in keep])
         if keep else np.zeros((0, 3), np.int64))
    g = np.searchsorted(starts, E[:, 0], side="right") - 1
    o = np.argsort(g, kind="stable")
    E, g = E[o], g[o]
    return E, np.searchsorted(g, np.arange(len(starts) + 1))


def colour(n_sites: int, ea: np.ndarray, eb: np.ndarray) -> list[np.ndarray]:
    """Greedy graph colouring of the coupling graph, largest-degree first.

    Sites sharing a colour have no edge between them, so they are conditionally independent given
    the rest and a whole colour class updates in one vectorised Gibbs step. The property is
    *asserted* against the real edge lists wherever it is relied on, never argued.
    """
    adj: list[list[int]] = [[] for _ in range(n_sites)]
    for a, b in zip(ea.tolist(), eb.tolist()):
        adj[a].append(b)
        adj[b].append(a)
    col = np.full(n_sites, -1, np.int64)
    for i in sorted(range(n_sites), key=lambda i: -len(adj[i])):
        used = {col[j] for j in adj[i] if col[j] >= 0}
        c = 0
        while c in used:
            c += 1
        col[i] = c
    return [np.nonzero(col == c)[0] for c in range(col.max() + 1)] if n_sites else []


def coupling_matrix(n_sites: int, ea, eb, ec, kv: np.ndarray) -> np.ndarray:
    """Dense symmetric ``A`` with zero diagonal: ``A[a,b]`` is the coefficient of the class
    linking ``a`` and ``b``."""
    A = np.zeros((n_sites, n_sites))
    np.add.at(A, (ea, eb), kv[ec])
    return A + A.T
