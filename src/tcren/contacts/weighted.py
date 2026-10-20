"""Multiplicity-conditioned reference ranks for frozen weighted contact energies.

The target calibration is native CDR:peptide. External antibody:antigen edges
supply exact-opportunity strata absent in native observations. These are ranks
of reference compatibility, not probabilities of physical or binding validity.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import polars as pl
import typer

from .plausibility import marginal_compatibilities

VERSION = "contact-rank-calibration-1"
STRATA = ["n_distance", "n_omega", "n_theta", "n_phi"]
KEYS = ["pdb.id", "chain.id.from", "residue.index.from", "chain.id.to", "residue.index.to"]


def multiplicity_edges(observations, resource):
    """One edge; distinct supported angular terms never multiply with atom rows."""
    obs = observations.with_row_index("row")
    if "source_resid.from" in obs.columns:
        keys = ["pdb.id", "source_chain.from", "source_resid.from", "source_chain.to", "source_resid.to"]
    else:
        keys = KEYS
    tails = marginal_compatibilities(observations, resource).join(obs.select("row", *keys), on="row", validate="m:1")
    tails = tails.filter(pl.col("compatibility").is_finite()).with_columns(
        pl.col("key").str.split("|").list.first().alias("family"))
    counts = tails.group_by(keys).agg(pl.col("compatibility").min().alias("minimum_tail"),
        pl.col("row").filter(pl.col("family")=="distance").n_unique().alias("n_distance"),
        *[pl.struct("key", "x").filter(pl.col("family")==family).n_unique().alias("n_"+family)
          for family in ("omega", "theta", "phi")])
    base = obs.group_by(keys).agg(pl.len().alias("n_atom_pairs"),
        pl.col("dist").min().alias("min_distance_A"), pl.col("dist").max().alias("max_distance_A"))
    result = base.join(counts, on=keys, how="left", validate="1:1").with_columns(pl.col(STRATA).fill_null(0))
    if keys == KEYS:
        # Native and benchmark observations may store the physical TCR endpoint second.
        roles = obs.group_by(keys).agg(pl.col("role.from").first())
        result = result.join(roles, on=keys, validate="1:1")
        swap = pl.col("role.from") != "tcr"
        result = result.with_columns(*[
            pl.when(swap).then(pl.col(field+".to")).otherwise(pl.col(field+".from")).alias(field+".from")
            for field in ("chain.id", "residue.index")], *[
            pl.when(swap).then(pl.col(field+".from")).otherwise(pl.col(field+".to")).alias(field+".to")
            for field in ("chain.id", "residue.index")]).drop("role.from")
    return result


def fit_calibration(external, native, resource):
    """Freeze sorted empirical minimum-tail distributions without benchmark input."""
    native = native.filter(pl.col("interface")=="tcr_peptide")
    records, checks = [], []
    for source, observations in (("native_tcr_peptide", native), ("external_antibody_antigen", external)):
        edges = multiplicity_edges(observations, resource)
        for stratum, part in edges.filter(pl.col("minimum_tail").is_finite()).group_by(STRATA):
            record = dict(source=source, **dict(zip(STRATA, stratum)), n_edges=part.height,
                n_structures=part["pdb.id"].n_unique(), minimum_tails=sorted(part["minimum_tail"].to_list()))
            records.append(record)
        checks.append(dict(source=source, n_atom_pairs=observations.height, n_edges=edges.height,
            n_supported_edges=edges["minimum_tail"].count(), n_structures=edges["pdb.id"].n_unique()))
    return dict(version=VERSION, strata=STRATA, records=records, reference_counts=checks,
        rank="w=(1+count(reference minimum_tail <= query minimum_tail))/(N+1); ties retained conservatively",
        source_policy="native CDR:peptide exact stratum; external antibody:antigen exact stratum only if native absent; otherwise neutral 1 with explicit uncalibrated flag",
        angular_opportunities="unique(key,transformed_value) per physical edge; repeated angle copies across atom rows removed",
        interpretation="In-reference conditional empirical ranks; correlated atoms and complexes, fitted reference marginals; no conformal coverage, naturalness, occupancy or binding probability claim",
        energy="sum(w_r*e_r) on CDR peptide contacts; unchanged framework, lookup matrix and missing-entry omission")


def calibrated_weights(edges, calibration):
    """Vectorized searchsorted within coarse opportunity groups; no label access."""
    if calibration.get("version") != VERSION or calibration.get("strata") != STRATA:
        raise ValueError("unsupported contact rank calibration")
    lookup = {(r["source"], *[r[c] for c in STRATA]): r for r in calibration["records"]}
    parts = []
    for key, part in edges.with_row_index("_order").group_by(STRATA, maintain_order=True):
        record = lookup.get(("native_tcr_peptide", *key)) or lookup.get(("external_antibody_antigen", *key))
        s = part["minimum_tail"].to_numpy()
        finite = np.isfinite(s)
        weight = np.ones(part.height)
        if record:
            weight[finite] = (1+np.searchsorted(record["minimum_tails"], s[finite], side="right"))/(record["n_edges"]+1)
        parts.append(part.with_columns(pl.Series("contact_weight", weight),
            pl.Series("calibration_source", [record["source"] if record and valid else "unavailable_geometry" if not valid else "unseen_multiplicity" for valid in finite]),
            pl.lit(record["n_edges"] if record else 0).alias("calibration_n_edges"),
            pl.lit(record["n_structures"] if record else 0).alias("calibration_n_structures")))
    if not parts:
        return edges.with_columns(pl.lit(1.).alias("contact_weight"), pl.lit("unseen_multiplicity").alias("calibration_source"),
            pl.lit(0).alias("calibration_n_edges"), pl.lit(0).alias("calibration_n_structures"))
    return pl.concat(parts).sort("_order").drop("_order")


def apply_weighted(table, weights, summary):
    """Multiply original finite entries once; framework is explicitly unweighted."""
    result = table.join(weights, on=KEYS, how="left", validate="1:1")
    cdr = pl.col("geometry_status") != "unassessed_framework"
    result = result.with_columns(pl.when(cdr).then(pl.col("contact_weight").fill_null(1.)).otherwise(1.).alias("contact_weight"),
        pl.when(cdr).then(pl.col("calibration_source").fill_null("unavailable_geometry")).otherwise(pl.lit("unassessed_framework")).alias("calibration_source"))
    if result["contact_weight"].null_count():
        raise ValueError("CDR geometry edge missing from weight join")
    result = result.with_columns((pl.col("original_contribution")*pl.col("contact_weight")).alias("weighted_contribution"))
    result = result.with_columns((pl.col("weighted_contribution")-pl.col("original_contribution")).alias("weighted_increment"))
    summary = dict(summary, weighted=result["weighted_contribution"].sum(),
        uncalibrated_cdr_contacts=result.filter(pl.col("calibration_source").is_in(["unseen_multiplicity", "unavailable_geometry"])).height)
    return result, summary


def calibration_command(
    external: Path = typer.Option(..., "--external", exists=True, dir_okay=False),
    native: Path = typer.Option(..., "--native", exists=True, dir_okay=False),
    resource: Path = typer.Option(..., "--resource", exists=True, dir_okay=False),
    out: Path = typer.Option(..., "-o", "--out"),
):
    """Freeze multiplicity-conditioned compatibility ranks from reference structures."""
    from ..exports import _sha256
    if out.exists():
        raise FileExistsError(out)
    keys = ["source_chain.from", "source_resid.from", "source_chain.to", "source_resid.to"]
    cols = ["pdb.id", "chem.from", "chem.to", "dist", "omega_rad", "theta.from_rad", "theta.to_rad", "phi.from_rad", "phi.to_rad"]
    ex = pl.read_csv(external, separator="\t", columns=keys+cols, schema_overrides=dict.fromkeys(keys, pl.String), infer_schema_length=10000)
    calibration = fit_calibration(ex, pl.read_parquet(native), json.loads(resource.read_text()))
    calibration["sources"] = {str(p): _sha256(p) for p in (external, native, resource)}
    calibration["marginal_resource_sha256"] = _sha256(resource)
    calibration["producer_sha256"] = _sha256(Path(__file__))
    out.write_text(json.dumps(calibration, indent=2, allow_nan=False)+"\n")
    print(json.dumps(dict(resource_sha256=_sha256(out), reference_counts=calibration["reference_counts"], n_strata=len(calibration["records"]))))
