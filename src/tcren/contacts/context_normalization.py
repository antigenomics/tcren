"""Frozen native opportunity-conditioned reference partition for geometry penalties.

This normalizes a dimensionless factor against empirical native contacts; it is
not a probability of validity, a binding free energy or a SASA-matched potential.
"""
from __future__ import annotations

import json
from pathlib import Path

import polars as pl
import typer

from .plausibility import VERSION as PRIOR_VERSION
from .weighted import KEYS, STRATA, multiplicity_edges

VERSION = "contact-context-partition-1"


def fit_partition(edges):
    """Freeze log(mean(exp(-u))) per exact opportunity stratum; no labels."""
    if edges.select(KEYS).is_duplicated().any():
        raise ValueError("partition requires unique directed physical edges")
    supported = edges.filter(pl.col("minimum_tail").is_finite())
    if supported.is_empty():
        raise ValueError("native reference has no supported edges")
    if supported.filter(~pl.col("minimum_tail").is_between(0., 1.)).height:
        raise ValueError("reference tails must be in [0,1]")
    result = supported.with_columns((pl.col("minimum_tail")-1).exp().alias("factor")).group_by(STRATA).agg(
        pl.col("factor").mean().log().alias("logZ"), pl.len().alias("n_edges"),
        pl.col("pdb.id").n_unique().alias("n_structures"))
    checked = supported.join(result.select(*STRATA, "logZ"), on=STRATA, validate="m:1").group_by(STRATA).agg(
        (pl.col("minimum_tail")-1-pl.col("logZ")).exp().mean().alias("mean_reference_normalized_factor"))
    return result.join(checked, on=STRATA, validate="1:1").sort(STRATA)


def fit_normalizer(native, resource, resource_sha256, mhc_class="I"):
    """Estimate only from explicit native CDR→peptide contacts in one MHC class."""
    if resource.get("version") != PRIOR_VERSION:
        raise ValueError("normalizer requires directed v2 marginal resource")
    if mhc_class not in {"I", "II"}:
        raise ValueError("MHC class must be I or II")
    required = {"mhc_class", "interface", "role.from", "role.to", "region.type.from", "region.type.to"}
    if not required <= set(native.columns):
        raise ValueError("native reference requires explicit class, interface, roles and residue regions")
    cdr = pl.when(pl.col("role.from") == "tcr").then(pl.col("region.type.from")).otherwise(pl.col("region.type.to"))
    selected = native.filter((pl.col("mhc_class") == mhc_class) & (pl.col("interface") == "tcr_peptide") &
                             cdr.is_in(["CDR1", "CDR2", "CDR3"]).fill_null(False))
    if selected.is_empty():
        raise ValueError("no native CDR-peptide observations in declared MHC class")
    edges = multiplicity_edges(selected, resource)
    partition = fit_partition(edges)
    return dict(version=VERSION, strata=STRATA, records=partition.to_dicts(), mhc_class=mhc_class,
        marginal_resource_version=resource["version"], marginal_resource_sha256=resource_sha256,
        reference_counts=dict(n_atom_pairs=selected.height, n_edges=edges.height,
            n_supported_edges=edges["minimum_tail"].is_finite().sum(), n_structures=selected["pdb.id"].n_unique()),
        reference_context="native CDR-to-peptide only; exact MHC class and opportunity stratum; no source fallback",
        definition="u=1-minimum_tail; logZ=log(mean_native exp(-u)); g=u+logZ; E_normalized=E_original+sum_CDR g",
        dependencies="Empirically estimated native reference partition and v2 geometry marginals; no benchmark labels or fitted score coefficients",
        missing_policy="unseen strata neutral g=0; missing actual CDR geometry unavailable; framework g=0",
        interpretation="Relative dimensionless reference likelihood factor exp(-g); not contact validity, binding probability or binding free energy")


def validate_normalizer(normalizer, resource_sha256):
    if normalizer.get("version") != VERSION or normalizer.get("strata") != STRATA:
        raise ValueError("unsupported contact context normalizer")
    if normalizer.get("marginal_resource_version") != PRIOR_VERSION:
        raise ValueError("normalizer requires directed v2 marginal resource")
    if normalizer.get("marginal_resource_sha256") != resource_sha256:
        raise ValueError("normalizer belongs to a different marginal resource")
    if normalizer.get("mhc_class") not in {"I", "II"}:
        raise ValueError("normalizer requires explicit MHC class")
    records = pl.DataFrame(normalizer["records"], infer_schema_length=None)
    if records.is_empty() or records.select(STRATA).is_duplicated().any():
        raise ValueError("normalizer requires unique nonempty strata")
    if records.filter((~pl.col("logZ").is_finite() | ~pl.col("logZ").is_between(-1., 0.) |
                      (pl.col("n_edges") <= 0) | (pl.col("n_structures") <= 0)).fill_null(True)).height:
        raise ValueError("finite logZ in [-1,0] and positive reference counts required")


def normalized_edges(edges, normalizer):
    """One dimensionless increment per directed edge, retaining unavailable status."""
    if edges.select(KEYS).is_duplicated().any():
        raise ValueError("normalization requires unique directed physical edges")
    reference = pl.DataFrame(normalizer["records"], infer_schema_length=None).rename({
        "n_edges": "normalizer_n_edges", "n_structures": "normalizer_n_structures"})
    joined = edges.join(reference, on=STRATA, how="left", validate="m:1")
    finite = pl.col("minimum_tail").is_finite().fill_null(False)
    if joined.filter(finite & ~pl.col("minimum_tail").is_between(0., 1.)).height:
        raise ValueError("query tails must be in [0,1]")
    return joined.with_columns(
        pl.when(~finite).then(None).when(pl.col("logZ").is_null()).then(0.)
        .otherwise(1-pl.col("minimum_tail")+pl.col("logZ")).alias("normalized_increment"),
        pl.when(~finite).then(pl.lit("unavailable_geometry"))
        .when(pl.col("logZ").is_null()).then(pl.lit("unseen_multiplicity"))
        .otherwise(pl.lit("native_tcr_peptide")).alias("normalizer_status"))


def apply_normalized(table, increments, summary):
    """Preserve original matrix lookup and framework, summing each CDR edge once."""
    if table.select(KEYS).is_duplicated().any() or increments.select(KEYS).is_duplicated().any():
        raise ValueError("normalization requires unique directed physical edges")
    cols = [*KEYS, "normalized_increment", "normalizer_status", "logZ", "normalizer_n_edges", "normalizer_n_structures"]
    result = table.join(increments.select(cols), on=KEYS, how="left", validate="1:1")
    cdr = pl.col("geometry_status") != "unassessed_framework"
    result = result.with_columns(
        pl.when(cdr).then(pl.col("normalized_increment")).otherwise(0.).alias("normalized_increment"),
        pl.when(cdr).then(pl.col("normalizer_status").fill_null("unavailable_geometry"))
        .otherwise(pl.lit("unassessed_framework")).alias("normalizer_status"))
    result = result.with_columns((pl.col("original_contribution")+pl.col("normalized_increment")).alias("normalized_contribution"))
    missing = result.filter(cdr & (pl.col("normalizer_status") == "unavailable_geometry")).height
    return result, dict(summary, normalized=None if missing else result["normalized_contribution"].sum(),
        normalized_missing_geometry=missing,
        normalized_unseen_stratum=result.filter(pl.col("normalizer_status") == "unseen_multiplicity").height)


def normalizer_command(
    native: Path = typer.Option(..., "--native", exists=True, dir_okay=False),
    resource: Path = typer.Option(..., "--resource", exists=True, dir_okay=False),
    mhc_class: str = typer.Option(..., "--mhc-class"),
    out: Path = typer.Option(..., "-o", "--out"),
    reference_checks: Path | None = typer.Option(None, "--reference-checks"),
):
    """Freeze native opportunity-conditioned dimensionless geometry partitions."""
    from ..exports import _sha256
    if out.exists():
        raise FileExistsError(out)
    if reference_checks is not None and reference_checks.exists():
        raise FileExistsError(reference_checks)
    model = fit_normalizer(pl.read_parquet(native), json.loads(resource.read_text()), _sha256(resource), mhc_class)
    model["sources"] = {str(p): _sha256(p) for p in (native, resource)}
    model["producer_sha256"] = _sha256(Path(__file__))
    validate_normalizer(model, _sha256(resource))
    out.write_text(json.dumps(model, indent=2, allow_nan=False)+"\n")
    if reference_checks is not None:
        pl.DataFrame(model["records"]).write_csv(reference_checks, separator="\t")
    print(json.dumps(dict(normalizer_sha256=_sha256(out), reference_counts=model["reference_counts"], n_strata=len(model["records"]))))
