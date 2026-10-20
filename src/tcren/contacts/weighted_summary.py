"""Table-only signed accounting for frozen multiplicity-weighted contacts."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import polars as pl
import typer

from .mechanism import EDGE
from .mechanism_summary import _mean_profiles, _read, rank_changes
from ..associations import GROUPS, _unique
from ..runtime import single_threaded

TERMS = ["original_positive_magnitude", "original_favourable_magnitude",
         "weighted_positive_magnitude", "weighted_favourable_magnitude",
         "positive_magnitude_lost", "favourable_magnitude_lost", "weighted_increment"]
COUNTS = ["contact_weight", "uncalibrated_cdr_contacts", "external_fallback_contacts",
          "framework_contacts", "known_potential_contacts", "n_atom_pairs", "cdr_contacts",
          "native_calibrated_contacts", "cdr_weight_equivalent", "cdr_downweighted_equivalent"]


def weighted_terms(edges):
    """Preserve positive/favourable signs; missing matrix entries stay unavailable."""
    e, w = pl.col("original_contribution"), pl.col("contact_weight")
    cdr = pl.col("geometry_status") != "unassessed_framework"
    return edges.with_columns(
        cdr.alias("cdr_contacts"),
        (pl.col("calibration_source")=="native_tcr_peptide").alias("native_calibrated_contacts"),
        pl.when(cdr).then(w).otherwise(0.).alias("cdr_weight_equivalent"),
        pl.when(cdr).then(1-w).otherwise(0.).alias("cdr_downweighted_equivalent"),
        e.clip(lower_bound=0).alias("original_positive_magnitude"),
        (-e).clip(lower_bound=0).alias("original_favourable_magnitude"),
        (w*e.clip(lower_bound=0)).alias("weighted_positive_magnitude"),
        (w*(-e).clip(lower_bound=0)).alias("weighted_favourable_magnitude"),
        ((1-w)*e.clip(lower_bound=0)).alias("positive_magnitude_lost"),
        ((1-w)*(-e).clip(lower_bound=0)).alias("favourable_magnitude_lost"),
        pl.col("calibration_source").is_in(["unseen_multiplicity", "unavailable_geometry"]).alias("uncalibrated_cdr_contacts"),
        (pl.col("calibration_source")=="external_antibody_antigen").alias("external_fallback_contacts"),
        (pl.col("geometry_status")=="unassessed_framework").alias("framework_contacts"),
        e.is_not_null().alias("known_potential_contacts"))


def model_summary(scores, labels, edges):
    """Sum all physical edges, retaining score rows with zero contacts."""
    sums = edges.group_by("pdb.id", "potential").agg(
        *[pl.col(c).sum().alias(c) for c in [*TERMS, *COUNTS]],
        pl.len().alias("summarized_contacts"),
        pl.col("contact_weight").mean().alias("mean_contact_weight"),
        pl.col("contact_weight").filter(pl.col("cdr_contacts")).mean().alias("mean_cdr_contact_weight"))
    result = scores.drop([c for c in [*TERMS, *COUNTS] if c in scores.columns]).join(labels.drop("mhc_class"), on="id", validate="m:1").join(
        sums.rename({"pdb.id": "id"}), on=["id", "potential"], how="left", validate="1:1")
    return result.with_columns(pl.col([*TERMS, *COUNTS, "summarized_contacts"]).fill_null(0))


def _exact(left, right, keys, columns, name):
    _unique(left, keys, name)
    _unique(right, keys, name+" accepted")
    if not left.select(*keys, *columns).sort(keys).equals(right.select(*keys, *columns).sort(keys)):
        raise ValueError("unchanged accepted "+name+" required")


@single_threaded
def summarize_weighted(run: Path, accepted: Path, mechanism: Path, out: Path):
    """Verify immutable populations/baselines, then explain weighted signed losses."""
    from ..exports import _sha256
    if out.exists():
        raise FileExistsError(out)
    labels, scores = _read(run/"labels.tsv"), _read(run/"scoring/scores.tsv")
    edges = _read(run/"scoring/contacts.tsv.gz")
    metrics = _read(run/"comparison/metrics.tsv")
    old_labels, old_scores = _read(accepted/"labels.tsv"), _read(accepted/"scoring/scores.tsv")
    _exact(labels, old_labels, ["id"], [*GROUPS, "y", "iptm"], "labels")
    _exact(scores, old_scores, ["id", "potential"], ["original", "n_contacts", "missing_potential_entries", "mhc_class"], "scores")
    _exact(edges, _read(accepted/"scoring/contacts.tsv.gz"), [*EDGE, "potential"], ["original_contribution"], "contact graph and matrix entries")
    ledger, old_ledger = _read(run/"scoring/accounting.tsv"), _read(accepted/"scoring/accounting.tsv")
    _exact(ledger, old_ledger, ["id"], ["sha256", "status"], "raw sources")
    if set(ledger["id"]) != set(labels["id"]) or ledger["sha256"].null_count() or set(ledger["status"]) != {"ok"}:
        raise ValueError("complete successful raw-source identities required")
    if set(labels["y"])-{0, 1} or labels["y"].null_count() or not labels["iptm"].is_finite().all() or labels.select(pl.any_horizontal(pl.col(GROUPS).is_null()).any()).item():
        raise ValueError("finite binary labels and ipTM required")
    if not scores.select(pl.all_horizontal(pl.col(["original", "weighted"]).is_finite()).all()).item():
        raise ValueError("finite original and weighted scores required")
    for _, part in scores.group_by("potential"):
        if set(part["id"]) != set(labels["id"]):
            raise ValueError("every potential must retain every accepted model")
    membership = _read(run/"comparison/membership.tsv")
    _exact(membership, scores, ["id", "potential"], [], "comparison membership")
    if not membership["used"].all():
        raise ValueError("all accepted comparison rows must be used")
    if edges["contact_weight"].null_count() or not edges["contact_weight"].is_finite().all() or edges.filter(~pl.col("contact_weight").is_between(0, 1)).height:
        raise ValueError("finite contact weights in [0,1] required")
    expected = edges["original_contribution"]*edges["contact_weight"]
    if not np.allclose(expected.to_numpy(), edges["weighted_contribution"].to_numpy(), rtol=0, atol=1e-12, equal_nan=True):
        raise ValueError("edge weighted identity failed")
    expected_increment = edges["weighted_contribution"]-edges["original_contribution"]
    if not np.allclose(expected_increment.to_numpy(), edges["weighted_increment"].to_numpy(), rtol=0, atol=1e-12, equal_nan=True):
        raise ValueError("edge increment identity failed")
    multiplicity = edges.group_by("potential", "n_distance", "n_omega", "n_theta", "n_phi", "calibration_source").agg(
        pl.len().alias("n_edges"), pl.col("pdb.id").n_unique().alias("n_models"),
        pl.col("calibration_n_edges").first().alias("n_reference_edges"),
        pl.col("calibration_n_structures").first().alias("n_reference_structures"),
        pl.col("contact_weight").mean().alias("mean_contact_weight"),
        pl.col("contact_weight").median().alias("median_contact_weight"),
        pl.col("minimum_tail").mean().alias("mean_minimum_tail"))
    joined = weighted_terms(edges)
    models = model_summary(scores, labels, joined)
    identities = {
        "original": pl.col("original_positive_magnitude")-pl.col("original_favourable_magnitude"),
        "weighted": pl.col("weighted_positive_magnitude")-pl.col("weighted_favourable_magnitude"),
        "weighted_increment": pl.col("favourable_magnitude_lost")-pl.col("positive_magnitude_lost"),
        "cdr_contacts": pl.col("cdr_weight_equivalent")+pl.col("cdr_downweighted_equivalent"),
    }
    errors = {c: models.select((pl.col(c)-expr).abs().max()).item() for c, expr in identities.items()}
    errors["calibration_support"] = models.select((pl.col("cdr_contacts")-pl.col("native_calibrated_contacts")-pl.col("external_fallback_contacts")-pl.col("uncalibrated_cdr_contacts")).abs().max()).item()
    errors["score_increment"] = models.select((pl.col("weighted")-pl.col("original")-pl.col("weighted_increment")).abs().max()).item()
    if any(error is None or error>1e-10 for error in errors.values()) or models.filter(pl.col("n_contacts")!=pl.col("summarized_contacts")).height:
        raise ValueError("signed model decomposition or contact accounting failed")

    previous = _read(mechanism/"models.tsv")
    _exact(previous, scores, ["id", "potential"], ["original", "n_contacts"], "previous mechanism population")
    chemistry = [c for c in ["chemical_real", "chemical_proximity", "closest_misses_chemical"] if c in previous.columns]
    models = models.join(previous.select("id", "potential", *[pl.col(c).alias("previous_"+c) for c in chemistry]),
        on=["id", "potential"], validate="1:1")
    # The existing final geometry carries chemical labels only; old geometry weights are excluded.
    geometry = pl.read_parquet(mechanism/"contacts.parquet", columns=["id", *EDGE[1:], "potential", "union_type", "chemical_real", "chemical_proximity"]).rename({"id": "pdb.id"})
    _exact(geometry, edges, [*EDGE, "potential"], [], "previous chemical graph")
    joined = joined.join(geometry, on=[*EDGE, "potential"], validate="1:1").rename({"pdb.id": "id"}).join(labels, on="id", validate="m:1")
    values = [*TERMS, *COUNTS, "chemical_real", "chemical_proximity"]
    denominators = models.group_by(*GROUPS, "potential", "y").len(name="n_models")
    motifs = joined.group_by(*GROUPS, "potential", "y", "chain.type.from", "region.type.from", "pos.to", "residue.aa.from", "residue.aa.to").agg(
        pl.len().alias("n_edges"), pl.col("id").n_unique().alias("n_models_with_motif"),
        *[pl.col(c).sum().alias(c) for c in values]).join(denominators, on=[*GROUPS, "potential", "y"], validate="m:1")
    motifs = motifs.with_columns(*[(pl.col(c)/pl.col("n_models")).alias("per_model_"+c) for c in ["n_edges", *values]])
    motif_keys = [*GROUPS, "potential", "chain.type.from", "region.type.from", "pos.to", "residue.aa.from", "residue.aa.to"]
    contrast_values = ["per_model_"+c for c in ["weighted_increment", "favourable_magnitude_lost", "positive_magnitude_lost"]]
    positive = motifs.filter(pl.col("y")==1).select(*motif_keys, *[pl.col(c).alias("positive_"+c) for c in contrast_values])
    negative = motifs.filter(pl.col("y")==0).select(*motif_keys, *[pl.col(c).alias("negative_"+c) for c in contrast_values])
    contrasts = positive.join(negative, on=motif_keys, how="full", coalesce=True, validate="1:1").fill_null(0)
    support = denominators.group_by(*GROUPS, "potential").agg(
        (pl.col("y")==1).any().alias("has_positive_class"), (pl.col("y")==0).any().alias("has_negative_class"))
    contrasts = contrasts.join(support, on=[*GROUPS, "potential"], validate="m:1").with_columns(*[
        pl.when(pl.col("has_"+side+"_class")).then(pl.col(side+"_"+c)).otherwise(None).alias(side+"_"+c)
        for side in ["positive", "negative"] for c in contrast_values])
    contrasts = contrasts.with_columns(*[(pl.col("positive_"+c)-pl.col("negative_"+c)).alias("positive_minus_negative_"+c) for c in contrast_values])
    classification = joined.group_by(*GROUPS, "potential", "y", "union_type", "geometry_status", "calibration_source").agg(
        pl.len().alias("n_edges"), pl.col("id").n_unique().alias("n_models_with_edges"), *[pl.col(c).sum().alias(c) for c in values])
    rank_rows, pairs = [], []
    for key, part in models.sort("id").group_by(*GROUPS, "potential", maintain_order=True):
        context = dict(zip([*GROUPS, "potential"], key))
        y, old = part["y"].to_numpy(), part["original"].to_numpy()
        ledger = metrics.filter(pl.all_horizontal([pl.col(c)==v for c, v in context.items()]) & (pl.col("metric")=="roc_auc"))
        for arm in ["weighted", "iptm"]:
            new = -part[arm].to_numpy() if arm=="iptm" else part[arm].to_numpy()
            result = rank_changes(y, old, new)
            row = ledger.filter(pl.col("arm")==arm)
            if row.height != 1:
                raise ValueError("one recorded AUC per case/potential/arm required")
            target, base = row["value"].item(), row["original_value"].item()
            expected = None if target is None or base is None else target-base
            actual = result["delta_auc"]
            if not ((actual is None and expected is None) or (actual is not None and expected is not None and np.isclose(actual, expected, rtol=0, atol=2e-15))):
                raise ValueError("rank accounting differs from recorded AUC delta")
            n = result["n_pairs"]
            if n and (not np.isclose(result["old_credit"]/n, base, rtol=0, atol=2e-15) or not np.isclose(result["new_credit"]/n, target, rtol=0, atol=2e-15)):
                raise ValueError("rank accounting differs from recorded absolute AUC")
            for kind, (a, b) in result.pop("examples").items():
                pairs.append(dict(**context, arm=arm, change=kind, positive_id=part["id"][a], negative_id=part["id"][b],
                    original_positive=old[a], original_negative=old[b], new_positive_energy=new[a], new_negative_energy=new[b],
                    **{side+"_"+c: part[c][i] for side, i in [("positive", a), ("negative", b)] for c in TERMS}))
            rank_rows.append(dict(**context, reference_arm="original", arm=arm, **result, recorded_delta_auc=expected))
    out.mkdir(parents=True, exist_ok=False)
    multiplicity.write_csv(out/"multiplicity_profiles.tsv", separator="\t")
    models.write_csv(out/"models.tsv", separator="\t")
    _mean_profiles(models, ["n_contacts", *TERMS, *COUNTS, "mean_contact_weight", "mean_cdr_contact_weight", "original", "weighted", "iptm", *["previous_"+c for c in chemistry]]).write_csv(out/"label_profiles.tsv", separator="\t")
    motifs.write_csv(out/"motifs.tsv", separator="\t")
    contrasts.write_csv(out/"motif_contrasts.tsv", separator="\t")
    classification.write_csv(out/"classification.tsv", separator="\t")
    pl.DataFrame(rank_rows).write_csv(out/"rank_changes.tsv", separator="\t")
    pair_schema = {**{c: pl.String for c in [*GROUPS, "potential", "arm", "change", "positive_id", "negative_id"]},
                   **{c: pl.Float64 for c in ["original_positive", "original_negative", "new_positive_energy", "new_negative_energy", *[side+"_"+c for side in ["positive", "negative"] for c in TERMS]]}}
    pl.DataFrame(pairs, schema=pair_schema).write_csv(out/"case_pairs.tsv", separator="\t")
    (out/"VERIFICATION.json").write_text(json.dumps(dict(complete=True, rows=labels.height, model_potential_rows=models.height,
        cases=models.select(GROUPS).unique().height, potentials=sorted(scores["potential"].unique().to_list()),
        signed_identity_max_abs_errors=errors, rank_accounting_rows=len(rank_rows),
        unchanged_raw_sources=True, unchanged_original=True, unchanged_membership=True,
        no_fitting=True, no_new_accuracy_metrics=True,
        sources={str(p): _sha256(p) for p in [run/"labels.tsv", run/"scoring/scores.tsv", run/"scoring/contacts.tsv.gz", run/"comparison/metrics.tsv", accepted/"scoring/accounting.tsv", mechanism/"models.tsv", mechanism/"contacts.parquet"]}), indent=2)+"\n")


def weighted_summary_command(
    run: Path = typer.Option(..., "--run", exists=True, file_okay=False),
    accepted: Path = typer.Option(..., "--accepted", exists=True, file_okay=False),
    mechanism: Path = typer.Option(..., "--mechanism", exists=True, file_okay=False),
    out: Path = typer.Option(..., "-o", "--out"),
):
    """Explain signed weighted contacts and exact matched AUC rank changes."""
    summarize_weighted(run, accepted, mechanism, out)
