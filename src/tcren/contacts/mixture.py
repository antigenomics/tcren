"""Batched contact likelihood mixture using supplied development effects."""

from __future__ import annotations

from collections.abc import Mapping

import numpy as np
import polars as pl

_ENDPOINTS = ("chain.id.from", "residue.index.from", "chain.id.to", "residue.index.to")


def _mixture_fraction(segment: np.ndarray, effects: np.ndarray, n: int) -> np.ndarray:
    """Run the recorded EM across all structures with one segmented reduction per step."""
    if n == 0:
        return np.empty(0)
    count = np.bincount(segment, minlength=n)
    ratio = np.exp(np.clip(effects, -30, 30))
    weight = np.full(n, 0.5)
    for _ in range(60):
        current = weight[segment]
        responsibility = current * ratio / (current * ratio + (1 - current))
        updated = np.where(count > 0, np.bincount(segment, weights=responsibility, minlength=n)
                           / np.maximum(count, 1), 0.5)
        converged = np.max(np.abs(updated - weight)) < 1e-9
        weight = updated
        if converged:
            break
    return np.where(count > 0, weight, np.nan)


def atom_contact_mixture(
    contacts: pl.DataFrame, effects: Mapping[str, float], ids: list[str],
) -> pl.DataFrame:
    """Estimate the atom mixture after averaging cell effects on each residue pair.

    Supply directed rows from ``native_contact_table`` and development effects
    keyed by interface and atom names, such as ``tcr_peptide/CA|O``. Unsupported
    atoms contribute a zero log ratio and are counted. Structures without
    contacts retain a null mixture estimate.
    """
    if not effects or not np.isfinite(list(effects.values())).all():
        raise ValueError("supply a nonempty finite development cell-effect reference")
    models = pl.DataFrame({"pdb.id": ids}, schema={"pdb.id": pl.String}).with_row_index("segment")
    if len(set(ids)) != len(ids) or contacts.join(models, on="pdb.id", how="anti").height:
        raise ValueError("declare unique structure IDs including every contact")
    cells = pl.concat_str(
        "interface", pl.concat_str("atom.from", "atom.to", separator="|"), separator="/")
    edges = (contacts.with_columns(
        cells.replace_strict(dict(effects), default=None, return_dtype=pl.Float64).alias("effect"),
    ).group_by("pdb.id", *_ENDPOINTS).agg(
        pl.col("effect").fill_null(0).mean(), pl.col("effect").is_null().sum().alias("n_unknown"),
    ).join(models, on="pdb.id", how="inner"))
    counts = edges.group_by("pdb.id").agg(
        pl.len().alias("n_residue_pairs"), pl.col("n_unknown").sum().alias("n_atom_pairs_unknown"),
    )
    fraction = _mixture_fraction(edges["segment"].to_numpy(), edges["effect"].to_numpy(), len(ids))
    return (models.join(counts, on="pdb.id", how="left", maintain_order="left")
            .with_columns(pl.col("n_residue_pairs", "n_atom_pairs_unknown").fill_null(0),
                          pl.Series("pi_atom", fraction, nan_to_null=True))
            .drop("segment").rename({"pdb.id": "id"}))
