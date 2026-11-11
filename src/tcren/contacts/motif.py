"""Frozen all-six-CDR/peptide mock enrichment for production screening."""

import gzip
import json
from importlib.resources import files

import polars as pl

from .motif_mixture import empirical_survival, fisher_summary
from .residue_surprise import residue_counts, score_residues


def bundled_enrichment(residues: pl.DataFrame) -> pl.DataFrame:
    """Score residue-count rows against frozen histograms and all 2,056 mock graphs.

    Input is :func:`residue_counts` output. Whole mock graphs use the recorded
    leave-one-model-out Fisher values. Query graphs use the complete reference.
    Returns one row per ``id`` with enrichment, Fisher statistic and support.
    """
    scopes = (residues.filter(pl.col("chain") != "PEPTIDE")
              .select("id", "chain", "region").unique().group_by("id").len())
    if scopes.height != residues["id"].n_unique() or (scopes["len"] != 6).any():
        raise ValueError("Every graph requires all six CDR loops")
    data = files("tcren.data")
    reference = json.loads(data.joinpath("motif_residue_reference.json").read_bytes())
    controls = pl.read_csv(gzip.decompress(data.joinpath("motif_mock_controls.tsv.gz").read_bytes()),
                           separator="\t")
    summary = fisher_summary(score_residues(residues, reference, "cycles"))
    p = empirical_survival(summary["fisher_statistic"], controls["fisher_statistic"],
                           controls["weight"])
    return summary.with_columns(pl.Series("mock_fisher_p", p),
                                pl.lit(controls.height).alias("null_models")).with_columns(
        pl.when(pl.col("n_local_tests") > 0).then(1 - pl.col("mock_fisher_p"))
        .otherwise(0.0).alias("enrichment"))


def contact_enrichment(atoms: pl.DataFrame, contacts: pl.DataFrame) -> pl.DataFrame:
    """Frozen enrichment from complete atom inventory and physical contact tables."""
    return bundled_enrichment(residue_counts(atoms, contacts))
