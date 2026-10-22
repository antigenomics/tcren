"""Label-free observed contact counts and selected-atom availability marginals."""
from hashlib import sha256

import numpy as np
import polars as pl
from scipy.spatial import cKDTree

from .observations import _HELICES, atom_inventory

_BINDER_FIELDS = ["binder.chain.type", "binder.region.type", "binder.residue.aa", "binder.atom.name"]
_TARGET_FIELDS = ["target.chain.type", "target.region.type", "target.residue.aa", "target.atom.name",
                  "target.peptide.pos", "target.residue.index.pdb", "target.insertion_code"]
_DICT_SCHEMA = {"feature": pl.String, "binder_feature": pl.String, "target_feature": pl.String,
                "level": pl.String, "view": pl.String,
                **{key: pl.String for key in _BINDER_FIELDS + _TARGET_FIELDS},
                "target.peptide.pos": pl.Int64, "target.residue.index.pdb": pl.Int64}
_PAIR_SCHEMA = {"pdb.id": pl.String, "binder.atom.source_index": pl.Int64,
                "target.atom.source_index": pl.Int64, "dist": pl.Float64}
_COUNT_SCHEMA = {"pdb.id": pl.String, "feature": pl.String, "C": pl.Int64}
_BINDER_SCHEMA = {"pdb.id": pl.String, "binder_feature": pl.String, "N_binder": pl.Int64}
_TARGET_SCHEMA = {"pdb.id": pl.String, "target_feature": pl.String, "N_target": pl.Int64}
_OPPORTUNITY_SCHEMA = {**_COUNT_SCHEMA, "O": pl.Int64}


def _endpoint_features(atoms, side, level, view):
    """Vectorized endpoint keys; hash only distinct keys, never physical pairs."""
    fields = _BINDER_FIELDS if side == "binder" else _TARGET_FIELDS
    table = atoms.select(
        "pdb.id", "atom.source_index", "chain.id", "residue.index",
        pl.col("chain.type").alias(f"{side}.chain.type"),
        (pl.when(pl.col("chain.type") == "PEPTIDE").then(pl.lit("PEPTIDE"))
         .otherwise(pl.col("region.type").fill_null("UNASSIGNED"))).alias(f"{side}.region.type"),
        pl.col("residue.aa").alias(f"{side}.residue.aa"),
        (pl.col("atom.name") if level == "atom" else pl.lit(None, dtype=pl.String))
        .alias(f"{side}.atom.name"))
    if side == "target":
        # Chemistry pools sites within the declared case; site keys never pool positions.
        table = table.with_columns(
            atoms["peptide.pos"].alias("target.peptide.pos") if view == "site"
            else pl.lit(None, dtype=pl.Int64).alias("target.peptide.pos"),
            (pl.when(pl.Series(atoms["chain.type"]) != "PEPTIDE")
             .then(pl.Series(atoms["residue.index.pdb"])).otherwise(None)
             if view == "site" else pl.lit(None, dtype=pl.Int64)).alias("target.residue.index.pdb"),
            (pl.when(pl.Series(atoms["chain.type"]) != "PEPTIDE")
             .then(pl.Series(atoms["insertion_code"])).otherwise(None)
             if view == "site" else pl.lit(None, dtype=pl.String)).alias("target.insertion_code"))
    table = table.with_columns(pl.struct(fields).struct.json_encode().alias("_key"))
    keys = table.select("_key").unique().with_columns(
        pl.col("_key").map_elements(lambda key: sha256(f"{side}:{level}:{view}:{key}".encode()).hexdigest(),
                                    return_dtype=pl.String).alias(f"{side}_feature"))
    return table.join(keys, on="_key", how="left", validate="m:1").drop("_key")


def _feature_id():
    return pl.concat_str("binder_feature", "target_feature", separator=":").alias("feature")


def features_for_structure(structure, cutoff=5.0):
    """Return inventory, physical pairs, motifs, two marginals and feature dictionary.

    All selected TCR atoms (including framework/unassigned regions) are binders;
    peptide and annotated MHC helices are targets. Full inventory keeps other
    components as evidence. One sparse neighbour search supplies atom pairs at an
    inclusive cutoff in Angstroms. Residue contacts collapse physical residue pairs.
    Availability counts all selected endpoints, including uncontacted endpoints.

    Fixed ``site`` and ``chemistry`` views are each exported at ``atom`` and
    ``residue`` levels. Peptide site is its one-based observed ordinal; MHC site is
    chain type, annotated region, author residue number and insertion code. Callers
    must keep cases/alleles separate. No labels, fitted coefficients or scores enter.
    """
    if not np.isfinite(cutoff) or cutoff <= 0:
        raise ValueError("cutoff must be finite and positive")
    inventory = atom_inventory(structure)
    for chain in structure.by_type("MHCa", "MHCb"):
        if not any(region.region_type in _HELICES for region in chain.regions):
            raise ValueError("MHC helix annotation required for contact-class audit")
        if inventory.filter((pl.col("chain.id") == chain.chain_id) &
                            (pl.col("contact_role") == "mhc")).is_empty():
            raise ValueError("observed MHC helix heavy atoms required for contact-class audit")
    binders = inventory.filter(pl.col("role") == "tcr")
    targets = inventory.filter((pl.col("role") == "peptide") |
                               ((pl.col("role") == "mhc") & (pl.col("contact_role") == "mhc")))
    pairs = pl.DataFrame(schema=_PAIR_SCHEMA)
    if binders.height and targets.height:
        neighbors = cKDTree(binders.select("x", "y", "z").to_numpy()).sparse_distance_matrix(
            cKDTree(targets.select("x", "y", "z").to_numpy()), cutoff, output_type="coo_matrix")
        pairs = pl.DataFrame({
            "pdb.id": [structure.pdb_id] * neighbors.nnz,
            "binder.atom.source_index": binders["atom.source_index"].to_numpy()[neighbors.row],
            "target.atom.source_index": targets["atom.source_index"].to_numpy()[neighbors.col],
            "dist": neighbors.data}, schema=_PAIR_SCHEMA).sort(
                "binder.atom.source_index", "target.atom.source_index")
    if pairs.select("binder.atom.source_index", "target.atom.source_index").is_duplicated().any():
        raise ValueError("duplicate physical atom pair")
    counts, binder_counts, target_counts, dictionaries = [], [], [], []
    for level in ("atom", "residue"):
        for view in ("site", "chemistry"):
            b = _endpoint_features(binders, "binder", level, view)
            t = _endpoint_features(targets, "target", level, view)
            for endpoints, side, output in ((b, "binder", binder_counts), (t, "target", target_counts)):
                available = endpoints if level == "atom" else endpoints.unique(
                    subset=["pdb.id", "chain.id", "residue.index"])
                output.append(available.group_by("pdb.id", f"{side}_feature").agg(
                    pl.len().cast(pl.Int64).alias(f"N_{side}")))
            contacts = pairs.join(b.select(
                pl.col("atom.source_index").alias("binder.atom.source_index"), "binder_feature",
                pl.col("chain.id").alias("binder.chain.id"),
                pl.col("residue.index").alias("binder.residue.index")),
                on="binder.atom.source_index", how="left", validate="m:1").join(t.select(
                    pl.col("atom.source_index").alias("target.atom.source_index"), "target_feature",
                    pl.col("chain.id").alias("target.chain.id"),
                    pl.col("residue.index").alias("target.residue.index")),
                    on="target.atom.source_index", how="left", validate="m:1")
            if contacts["binder_feature"].null_count() or contacts["target_feature"].null_count():
                raise ValueError("physical atom indices failed inventory join")
            if level == "residue":
                contacts = contacts.unique(subset=["pdb.id", "binder.chain.id", "binder.residue.index",
                                                   "target.chain.id", "target.residue.index"])
            motifs = contacts.group_by("pdb.id", "binder_feature", "target_feature").agg(
                pl.len().cast(pl.Int64).alias("C")).with_columns(_feature_id())
            counts.append(motifs.select(*_COUNT_SCHEMA))
            dictionaries.append(motifs.select("feature", "binder_feature", "target_feature").unique()
                .join(b.select("binder_feature", *_BINDER_FIELDS).unique(), on="binder_feature", validate="m:1")
                .join(t.select("target_feature", *_TARGET_FIELDS).unique(), on="target_feature", validate="m:1")
                .with_columns(pl.lit(level).alias("level"), pl.lit(view).alias("view"))
                .select(*_DICT_SCHEMA))
    result = (inventory, pairs,
              pl.concat(counts).sort("pdb.id", "feature"),
              pl.concat(binder_counts).sort("pdb.id", "binder_feature"),
              pl.concat(target_counts).sort("pdb.id", "target_feature"),
              pl.concat(dictionaries).sort("feature"))
    opportunity_rows(*result[2:], [structure.pdb_id])
    return result


def opportunity_rows(motif_counts, binder_marginals, target_marginals, feature_dictionary, model_ids):
    """Positive opportunities for observed feature union, including zero-contact models.

    Join required dictionary endpoints onto compact availability marginals; no
    model×feature or binder×target universe is constructed. O=0 combinations are
    absent explicitly. Counts with missing support, duplicates or C>O fail.
    """
    ids = pl.DataFrame({"pdb.id": list(model_ids)}, schema={"pdb.id": pl.String})
    if ids["pdb.id"].null_count() or ids.is_duplicated().any():
        raise ValueError("duplicate or null model identity")
    for table, keys in ((motif_counts, ["pdb.id", "feature"]),
                        (binder_marginals, ["pdb.id", "binder_feature"]),
                        (target_marginals, ["pdb.id", "target_feature"]),
                        (feature_dictionary, ["feature"])):
        if table.select(keys).is_duplicated().any():
            raise ValueError("duplicate motif, marginal or feature identity")
        if table.select(pl.any_horizontal(pl.col(keys).is_null()).any()).item():
            raise ValueError("null motif, marginal or feature identity")
    for table, column in ((motif_counts, "C"), (binder_marginals, "N_binder"), (target_marginals, "N_target")):
        if table[column].null_count() or table.filter(pl.col(column) <= 0).height:
            raise ValueError("observed counts and available endpoint counts must be positive")
        if table.schema[column] != pl.Int64:
            raise ValueError("counts must be Int64 integers")
    if feature_dictionary.select(pl.any_horizontal(
            pl.col("binder_feature", "target_feature").is_null()).any()).item():
        raise ValueError("null dictionary endpoint identity")
    counts = motif_counts.join(ids, on="pdb.id", how="inner", validate="m:1")
    binders = binder_marginals.join(ids, on="pdb.id", how="inner", validate="m:1")
    targets = target_marginals.join(ids, on="pdb.id", how="inner", validate="m:1")
    rows = (feature_dictionary.select("feature", "binder_feature", "target_feature")
        .join(binders, on="binder_feature", how="inner", validate="m:m")
        .join(targets, on=["pdb.id", "target_feature"], how="inner", validate="m:1"))
    if rows.filter(pl.col("N_binder") > (pl.lit(np.iinfo(np.int64).max) // pl.col("N_target"))).height:
        raise ValueError("integer overflow in opportunity product")
    rows = (rows.with_columns((pl.col("N_binder") * pl.col("N_target")).alias("O"))
        .join(counts, on=["pdb.id", "feature"], how="left", validate="1:1")
        .with_columns(pl.col("C").fill_null(0)).select(*_OPPORTUNITY_SCHEMA))
    if rows.filter((pl.col("O") <= 0) | (pl.col("C") > pl.col("O"))).height:
        raise ValueError("contact counts require 0 <= C <= O with positive integer O")
    if counts.join(rows.select("pdb.id", "feature"), on=["pdb.id", "feature"], how="anti").height:
        raise ValueError("contact C has no positive opportunity O or dictionary feature")
    return rows.sort("pdb.id", "feature")
