"""Dataset-level analyses for the TCRen contact statistics and potentials.

Helpers for the analysis notebook / benchmarks, following the TCRen manuscript logic:
potential heatmaps and comparisons, the distribution of TCR↔peptide contacts per
structure and per region, and how contacts distribute over peptide / CDR3 positions as a
function of peptide / CDR3 length. They take the manuscript contact + summary tables as
explicit paths (the committed oracle lives under ``tests/assets/oracle/``).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import polars as pl

from .potential import Potential

# chain.type.from -> the summary CDR3 column whose length applies
_CDR3_LEN_COL = {"TRA": "cdr3a", "TRB": "cdr3b"}


def load_interface_contacts(
    contact_maps: str | Path, summary: str | Path
) -> pl.DataFrame:
    """Load and enrich the manuscript TCR↔peptide contact table.

    Adds, per contact: ``peptide_pos`` (0-based peptide position = ``residue.index.to``),
    ``peptide_len``, ``cdr3_len`` (CDR3α/β length for the contacting TCR chain),
    ``cdr3_rel_pos`` (residue index relative to the chain's first contacting CDR3 residue —
    a relative position, since the committed table carries no region start), and the
    ``nonred`` flag.
    """
    c = pl.read_csv(contact_maps)
    s = pl.read_csv(summary).select(
        "pdb.id", "peptide", "cdr3a", "cdr3b", "nonred"
    ).with_columns(
        pl.col("peptide").str.len_chars().alias("peptide_len"),
        pl.col("cdr3a").str.len_chars().alias("cdr3a_len"),
        pl.col("cdr3b").str.len_chars().alias("cdr3b_len"),
    )
    df = c.join(s, on="pdb.id", how="left").with_columns(
        pl.col("residue.index.to").alias("peptide_pos"),
        pl.when(pl.col("chain.type.from") == "TRA")
        .then(pl.col("cdr3a_len"))
        .otherwise(pl.col("cdr3b_len"))
        .alias("cdr3_len"),
    )
    # Relative CDR3 position: index minus the first CDR3-contacting index per (pdb, chain).
    cdr3 = df.filter(pl.col("region.type.from") == "CDR3").with_columns(
        (pl.col("residue.index.from")
         - pl.col("residue.index.from").min().over(["pdb.id", "chain.id.from"]))
        .alias("cdr3_rel_pos")
    )
    return df.join(
        cdr3.select("pdb.id", "chain.id.from", "residue.index.from", "cdr3_rel_pos"),
        on=["pdb.id", "chain.id.from", "residue.index.from"],
        how="left",
    )


def contacts_per_structure(df: pl.DataFrame, nonred_only: bool = True) -> pl.DataFrame:
    """Number of TCR↔peptide contacts per structure (with TRA/TRB split)."""
    if nonred_only:
        df = df.filter(pl.col("nonred"))
    return (
        df.group_by("pdb.id", "chain.type.from")
        .agg(pl.len().alias("n_contacts"))
        .sort("pdb.id", "chain.type.from")
    )


def region_contact_counts(df: pl.DataFrame, nonred_only: bool = True) -> pl.DataFrame:
    """Total contacts by TCR region (CDR1/2/3, FR) and chain (TRA/TRB)."""
    if nonred_only:
        df = df.filter(pl.col("nonred"))
    return (
        df.group_by("chain.type.from", "region.type.from")
        .agg(pl.len().alias("n_contacts"))
        .sort("chain.type.from", "region.type.from")
    )


def position_distribution(
    df: pl.DataFrame, side: str = "peptide", nonred_only: bool = True
) -> pl.DataFrame:
    """Contact counts by position, stratified by chain/molecule length.

    Args:
        df: enriched contacts (see :func:`load_interface_contacts`).
        side: ``"peptide"`` (peptide position vs peptide length) or ``"cdr3a"`` / ``"cdr3b"``
            (relative CDR3 position vs CDR3 length, for that TCR chain).
        nonred_only: restrict to non-redundant structures.

    Returns:
        Long counts: ``length, position, n_contacts``.
    """
    if nonred_only:
        df = df.filter(pl.col("nonred"))
    if side == "peptide":
        sub = df.select(pl.col("peptide_len").alias("length"),
                        pl.col("peptide_pos").alias("position"))
    elif side in ("cdr3a", "cdr3b"):
        chain = "TRA" if side == "cdr3a" else "TRB"
        sub = (
            df.filter((pl.col("region.type.from") == "CDR3")
                      & (pl.col("chain.type.from") == chain))
            .select(pl.col("cdr3_len").alias("length"),
                    pl.col("cdr3_rel_pos").alias("position"))
        )
    else:
        raise ValueError(f"side must be 'peptide'|'cdr3a'|'cdr3b', got {side!r}")
    return (
        sub.drop_nulls()
        .group_by("length", "position")
        .agg(pl.len().alias("n_contacts"))
        .sort("length", "position")
    )


def potential_long(potential: Potential) -> pl.DataFrame:
    """Heatmap-ready long form of a potential: ``residue.aa.from, residue.aa.to, value``."""
    return potential.matrix.select("residue.aa.from", "residue.aa.to", "value")


def compare_potentials(a: Potential, b: Potential) -> pl.DataFrame:
    """Join two potentials on their amino-acid pairs and add their difference.

    Returns ``residue.aa.from, residue.aa.to, value_a, value_b, diff`` over shared pairs.
    """
    left = potential_long(a).rename({"value": "value_a"})
    right = potential_long(b).rename({"value": "value_b"})
    return (
        left.join(right, on=["residue.aa.from", "residue.aa.to"], how="inner")
        .with_columns((pl.col("value_a") - pl.col("value_b")).alias("diff"))
    )


def potential_matrix(potential: Potential) -> tuple[np.ndarray, list[str], list[str]]:
    """Dense matrix + ``(from_labels, to_labels)`` for plotting a potential heatmap."""
    long = potential_long(potential)
    froms = sorted(long["residue.aa.from"].unique().to_list())
    tos = sorted(long["residue.aa.to"].unique().to_list())
    fi = {a: i for i, a in enumerate(froms)}
    ti = {a: i for i, a in enumerate(tos)}
    m = np.full((len(froms), len(tos)), np.nan)
    for row in long.iter_rows(named=True):
        m[fi[row["residue.aa.from"]], ti[row["residue.aa.to"]]] = row["value"]
    return m, froms, tos
