"""Contact-mechanism diagnostics from frozen parameters; no score or model changes."""
from __future__ import annotations

import polars as pl

from .plausibility import marginal_compatibilities
from ..contact_types import TYPES_V2, UNTYPED

EDGE = ["pdb.id", "chain.id.from", "residue.index.from", "chain.id.to", "residue.index.to"]
CHEMICAL = ["is_"+t for t in TYPES_V2 if t not in UNTYPED]


def geometry_diagnostics(observations: pl.DataFrame, resource: dict) -> pl.DataFrame:
    """One CDR:peptide edge, all-atom chemistry, and exact minimum-tail lineage.

    Chemical flags are the existing v2 heavy-atom heuristics; they do not certify
    physical correctness. The minimum tail is not a naturalness probability.
    Orientation winners concern the residue frame, not the arbitrarily tied atom.
    """
    obs = observations.with_row_index("row")
    swap = pl.col("role.from") != "tcr"
    obs = obs.with_columns([
        pl.when(swap).then(pl.col(field+".to")).otherwise(pl.col(field+".from")).alias(field+".from")
        for field in ("chain.id", "residue.index", "chem")] + [
        pl.when(swap).then(pl.col(field+".from")).otherwise(pl.col(field+".to")).alias(field+".to")
        for field in ("chain.id", "residue.index", "chem")])
    tails = marginal_compatibilities(observations, resource).join(
        obs.select("row", *EDGE, "chem.from", "chem.to", "dist"), on="row", validate="m:1")
    tails = tails.with_columns(pl.col("key").str.split("|").list.first().alias("family"))
    supported = tails.filter(pl.col("compatibility").is_finite())
    winner = supported.sort("compatibility", "key", "row").unique(EDGE, keep="first", maintain_order=True)
    winner = winner.select(*EDGE, *[pl.col(c).alias("winner_"+c) for c in
        ("family", "key", "pool", "x", "n", "mean", "m2", "parameter_source", "compatibility", "chem.from", "chem.to", "dist")])
    counts = supported.group_by(EDGE).agg(
        pl.len().alias("supported_marginal_rows"),
        pl.struct("key", "x").n_unique().alias("n_distinct_marginals"),
        *[pl.col("compatibility").filter(pl.col("family")==f).min().alias(f+"_compatibility_min")
          for f in ("distance", "omega", "theta", "phi")])
    base = obs.group_by(EDGE).agg(pl.len().alias("n_atom_pairs"),
        pl.col("dist").n_unique().alias("n_distinct_distances"),
        pl.col("dist").min().alias("min_distance_A"), pl.col("dist").max().alias("max_distance_A"),
        pl.any_horizontal(CHEMICAL).any().alias("chemical_supported"),
        *[pl.col("is_"+t).any() for t in TYPES_V2])
    nearest = obs.sort("dist", "row").unique(EDGE, keep="first").select(*EDGE,
        pl.col("contact.type").alias("closest_type"),
        pl.any_horizontal(CHEMICAL).alias("closest_chemical_supported"))
    result = base.join(nearest, on=EDGE, validate="1:1").join(winner, on=EDGE, how="left", validate="1:1").join(counts, on=EDGE, how="left", validate="1:1")
    return result.with_columns((1-pl.col("winner_compatibility")).alias("audit_unnaturalness"),
        pl.when(pl.col("winner_family")=="distance").then(pl.lit("atom_pair")).otherwise(pl.lit("residue_frame")).alias("winner_scope"),
        pl.coalesce([pl.when(pl.col("is_"+t)).then(pl.lit(t)) for t in TYPES_V2]).alias("union_type"))
