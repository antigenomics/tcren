"""Exploratory, epitope-held-out proxies for unlabelled cohort confidence.

These diagnostics never replace the default cohort-mean score. Frozen score
coefficients are reused; confidence proxy fits never read binding outcomes.
"""

import argparse
import gzip
from importlib.resources import files
import json
from pathlib import Path
import time

import numpy as np
import polars as pl
import typer

from .motif_optimize import COHORT, apply_confidence_score, complementary_score
from .profile import digest


def allele_group(values, fields):
    """Retain one or two colon-separated HLA fields; unknown alleles stay null."""
    if fields not in (1, 2):
        raise ValueError("MHC resolution must be one or two fields")
    pattern = (
        r"^(?:HLA-)?([A-Z]+[0-9]*\*\d{2,3}" + (r":\d{2,3}" if fields == 2 else "") + r")(?:[:A-Z].*)?$"
    )
    return pl.Series(values).str.extract(pattern, 1).to_numpy()


def fit_allele_reference(rows):
    """Equal-cohort confidence means; full reference, distinct from held-out tests."""
    cohort = rows.group_by(COHORT).agg(
        pl.col("mhc").n_unique().alias("alleles"),
        pl.col("mhc").first(),
        pl.col("plddt").mean().alias("mean_plddt"),
        pl.col("iptm_scored").mean().alias("mean_iptm"),
        pl.len().alias("n_models"),
    )
    if (cohort["alleles"] != 1).any():
        raise ValueError("Each confidence-reference cohort must have one MHC allele")
    parts = []
    for fields in (0, 1, 2):
        keys = ["*"] * cohort.height if fields == 0 else allele_group(cohort["mhc"], fields)
        part = cohort.with_columns(pl.Series("allele", keys)).filter(pl.col("allele").is_not_null())
        parts.append(
            part.group_by("allele")
            .agg(
                pl.col("mean_plddt", "mean_iptm").mean(),
                pl.col("n_models").sum(),
                pl.len().alias("n_cohorts"),
                pl.col("epitope").n_unique().alias("n_epitopes"),
            )
            .with_columns(pl.lit(fields).alias("fields"))
        )
    return (
        pl.concat(parts)
        .select(
            "fields", "allele", "mean_plddt", "mean_iptm", "n_models", "n_cohorts", "n_epitopes"
        )
        .sort("fields", "allele")
    )


def allele_confidence(mhc, fields=2, reference=None):
    """Explicit allele lookup with a flagged global fallback for unseen alleles.

    pLDDT is on the 0–100 scale; ipTM is on [0,1]. The bundled reference is a
    full fit of unlabelled confidence means, not independent validation.
    """
    if reference is None:
        reference = pl.read_csv(
            gzip.decompress(files("tcren.data").joinpath("confidence_alleles.tsv.gz").read_bytes()),
            separator="\t",
        )
    required = [
        "fields",
        "allele",
        "mean_plddt",
        "mean_iptm",
        "n_models",
        "n_cohorts",
        "n_epitopes",
    ]
    if any(name not in reference.columns for name in required) or any(
        reference[name].null_count() for name in required
    ):
        raise ValueError("Complete allele reference columns are required")
    if reference.select("fields", "allele").unique().height != reference.height:
        raise ValueError("Duplicate allele reference entries")
    numeric = reference.select(
        "mean_plddt", "mean_iptm", "n_models", "n_cohorts", "n_epitopes"
    ).to_numpy()
    if (
        not np.isfinite(numeric).all()
        or np.any((numeric[:, 0] < 0) | (numeric[:, 0] > 100))
        or np.any((numeric[:, 1] < 0) | (numeric[:, 1] > 1))
        or np.any(numeric[:, 2:] <= 0)
        or np.any(numeric[:, 2:] != np.floor(numeric[:, 2:]))
    ):
        raise ValueError("Finite bounded confidence and positive integer support required")
    global_rows = reference.filter((pl.col("fields") == 0) & (pl.col("allele") == "*"))
    if global_rows.height != 1:
        raise ValueError("Exactly one global confidence fallback is required")
    columns = ["mean_plddt", "mean_iptm", "n_models", "n_cohorts", "n_epitopes"]
    selected = reference.filter(pl.col("fields") == fields).select("allele", *columns)
    out = pl.DataFrame({"allele": allele_group(mhc, fields)}).join(
        selected, on="allele", how="left", validate="m:1", maintain_order="left"
    )
    fallback = global_rows.row(0, named=True)
    out = out.with_columns(pl.col("n_cohorts").is_null().alias("confidence_proxy_fallback"))
    out = out.with_columns([pl.col(c).fill_null(fallback[c]) for c in columns])
    if (
        not out["mean_plddt"].is_between(0, 100).all()
        or not out["mean_iptm"].is_between(0, 1).all()
    ):
        raise ValueError("Reference confidence must be finite and bounded")
    return out


def nearest_confidence(plddt, groups=None, model=None):
    """Choose the nearest frozen cohort-pLDDT representative (midpoint ties lower).

    Optional explicit groups use their observed mean before choosing a bin.
    Support is the count of observed input confidences, never reference support.
    Empty-confidence groups abstain. No model is fitted or smoothed.
    """
    if model is None:
        model = json.loads(
            files("tcren.data").joinpath("motif_confidence_reference.json").read_bytes()
        )
    reference = 100 * np.asarray(model["confidence_cdf"], dtype=float)
    if (
        reference.ndim != 1
        or not len(reference)
        or not np.isfinite(reference).all()
        or np.any(np.diff(reference) <= 0)
        or np.any((reference < 0) | (reference > 100))
    ):
        raise ValueError("Distinct sorted bounded confidence representatives required")
    value = pl.Series("observed_mean_plddt", plddt, dtype=pl.Float64).fill_nan(None)
    if value.drop_nulls().is_infinite().any() or not value.drop_nulls().is_between(0, 100).all():
        raise ValueError("Generator plddt must be finite or missing and lie in [0,100]")
    table = pl.DataFrame(value).with_columns(
        pl.col("observed_mean_plddt").is_not_null().cast(pl.Int64).alias("context_n")
    )
    if groups is not None:
        group = pl.Series("_group", groups)
        if len(group) != len(value) or group.null_count():
            raise ValueError("One nonmissing explicit group per confidence is required")
        table = (
            table.with_columns(group)
            .with_columns(
                pl.col("observed_mean_plddt").mean().over("_group"),
                pl.col("context_n").sum().over("_group"),
            )
            .drop("_group")
        )
    observed = table["observed_mean_plddt"].to_numpy()
    index = np.searchsorted((reference[:-1] + reference[1:]) / 2, observed, side="left")
    return table.with_columns(
        pl.Series("mean_plddt", reference[index]),
        pl.Series("representative_bin", index + 1),
        pl.Series("distance_plddt", np.abs(observed - reference[index])),
        (pl.col("context_n") == 0).alias("missing_context"),
    ).with_columns(
        [
            pl.when(pl.col("missing_context")).then(None).otherwise(pl.col(name)).alias(name)
            for name in ("mean_plddt", "representative_bin", "distance_plddt")
        ]
    )


def score_command(
    features: Path = typer.Option(..., exists=True),
    output: Path = typer.Option(..., "--output", "-o"),
    fields: int = typer.Option(2, min=1, max=2),
    reference: Path | None = typer.Option(None, exists=True),
    context: str = "allele",
    group_column: str | None = None,
):
    """Apply an explicit allele or nearest-plddt context to score components.

    Requires iptm, betti and enrichment, plus mhc for allele or plddt for
    nearest-plddt. --group-column explicitly averages nearest-plddt context.
    Ordinary screening retains its cohort mean; no fit is performed here.
    """
    if output.exists():
        raise ValueError("Output exists")
    rows = (
        pl.read_parquet(features)
        if features.suffix == ".parquet"
        else pl.read_csv(features, separator="\t")
    )
    if context not in ("allele", "nearest-plddt"):
        raise ValueError("--context must be allele or nearest-plddt")
    model = json.loads(files("tcren.data").joinpath("motif_confidence_reference.json").read_bytes())
    if context == "allele":
        if group_column is not None:
            raise ValueError("--group-column applies only to nearest-plddt context")
        ref = None if reference is None else pl.read_csv(reference, separator="\t")
        context_rows = allele_confidence(rows["mhc"], fields, ref)
        context_name = "allele_proxy_full_reference"
    else:
        if reference is not None:
            raise ValueError("--reference applies only to the allele lookup")
        context_rows = nearest_confidence(
            rows["plddt"], None if group_column is None else rows[group_column], model
        )
        context_name = (
            "nearest_frozen_plddt" if group_column is None else "nearest_group_mean_plddt"
        )
    features_only = rows.select("iptm", "betti", "enrichment").with_columns(
        pl.Series("plddt", context_rows["mean_plddt"])
    )
    # Each row receives the explicitly selected context; the score formula stays frozen.
    scored = apply_confidence_score(features_only, {**model, "gate": "complex-plddt"})
    if "missing_context" in context_rows.columns:
        scored = scored.with_columns(
            [
                pl.when(context_rows["missing_context"])
                .then(None)
                .otherwise(pl.col(name))
                .alias(name)
                for name in ("score", "gate_q", "w_I", "w_B", "w_M")
            ]
        )
    output_columns = [
        "score",
        "gate_q",
        "w_I",
        "w_B",
        "w_M",
        "missing_iptm",
        "missing_betti",
        "missing_enrichment",
    ]
    rows.with_columns(
        *context_rows.rename({c: "proxy_" + c for c in context_rows.columns}).get_columns(),
        *scored.select(output_columns).get_columns(),
        pl.lit(context_name).alias("confidence_context"),
    ).write_csv(output, separator="\t")


def proxy_predictions(rows):
    """Predict pLDDT/100 and ipTM cohort means, excluding the test epitope.

    Cohorts have equal training mass. Missing predictors use training means and
    explicit indicators. Fixed ipTM bins have width 0.1; unseen bins/alleles use the
    training cohort mean. Neither labels nor held-out cohort means enter fitting.
    """
    if rows["id"].n_unique() != rows.height:
        raise ValueError("Unique model IDs required")
    rows = rows.sort(*COHORT, "id")
    targets = rows.select(
        (pl.col("plddt").mean().over(COHORT) / 100).alias("mean_plddt"),
        pl.col("iptm_scored").mean().over(COHORT).alias("mean_iptm"),
    ).to_numpy()
    if not np.isfinite(targets).all():
        raise ValueError("Complete cohort confidence targets required")
    counts = rows.select(pl.len().over(COHORT))["len"].to_numpy()
    weight = 1 / counts
    epitope = rows["epitope"].to_numpy()
    raw = rows.select(
        "iptm_scored", "plddt", "ranking_confidence", "betti_scored", "enrichment_scored"
    ).to_numpy()
    raw[:, 1] /= 100
    methods = ["linear_confidence", "linear_structure", "iptm_decile", "mhc_field1", "mhc_field2"]
    predictions = {name: np.empty_like(targets) for name in methods}
    support = {name: np.zeros(rows.height, dtype=int) for name in methods[2:]}
    keys = {"iptm_decile": np.minimum((raw[:, 0] * 10).astype(int), 9)}
    keys.update({f"mhc_field{n}": allele_group(rows["mhc"], n) for n in (1, 2)})
    for held in np.unique(epitope):
        test = epitope == held
        train = ~test
        if not train.any():
            raise ValueError("At least two epitopes are required")
        w = weight[train]
        fallback = np.average(targets[train], axis=0, weights=w)
        for name, width in [("linear_confidence", 3), ("linear_structure", 5)]:
            x = raw[:, :width]
            missing = ~np.isfinite(x)
            means = np.array(
                [
                    np.average(x[train & ~missing[:, j], j], weights=weight[train & ~missing[:, j]])
                    if (train & ~missing[:, j]).any()
                    else 0.0
                    for j in range(width)
                ]
            )
            design = np.column_stack([np.ones(len(x)), np.where(missing, means, x), missing])
            root_w = np.sqrt(w)[:, None]
            beta = np.linalg.lstsq(design[train] * root_w, targets[train] * root_w, rcond=None)[0]
            predictions[name][test] = np.clip(design[test] @ beta, 0, 1)
        for name, key in keys.items():
            predictions[name][test] = fallback
            for value in set(key[test]):
                matching = (
                    train & (key == value)
                    if value is not None
                    else np.zeros(len(train), dtype=bool)
                )
                selected = test & (key == value)
                if matching.any():
                    predictions[name][selected] = np.average(
                        targets[matching], axis=0, weights=weight[matching]
                    )
                    support[name][selected] = np.unique(epitope[matching]).size
    predictions["single_model"] = raw[:, :2][:, ::-1]
    # Small sets are deterministic disjoint ID-sorted batches, not cherry-picked subsets.
    for size in (3, 5):
        group = rows.with_columns((pl.int_range(pl.len()).over(COHORT) // size).alias("batch"))
        predictions[f"small_set_{size}"] = group.select(
            (pl.col("plddt").mean().over([*COHORT, "batch"]) / 100),
            pl.col("iptm_scored").mean().over([*COHORT, "batch"]),
        ).to_numpy()
    predictions["cohort_default"] = targets.copy()
    return rows, targets, predictions, support


def evaluate(rows, resource):
    """Return flat individual predictions and cohort errors/score metrics."""
    from sklearn.metrics import average_precision_score, roc_auc_score

    rows, targets, predictions, support = proxy_predictions(rows)
    reference = np.round(resource["confidence_cdf"], 12)
    x = rows.select("iptm_scored", "betti_scored", "enrichment_scored").to_numpy()
    p = resource["parameters"]
    expected_q = np.searchsorted(reference, np.round(targets[:, 0], 12), side="right") / len(
        reference
    )
    records, details = [], []
    for method, predicted in predictions.items():
        q = np.searchsorted(reference, np.round(predicted[:, 0], 12), side="right") / len(reference)
        score = complementary_score(x, q, [p["b"], p["m"]])
        detail = rows.select("id", *COHORT, "y").with_columns(
            pl.lit(method).alias("method"),
            pl.Series("predicted_mean_plddt", predicted[:, 0] * 100),
            pl.Series("predicted_mean_iptm", predicted[:, 1]),
            pl.Series("score", score),
            pl.Series("gate_q", q),
            pl.Series("gate_error", np.abs(q - expected_q)),
            pl.Series("plddt_error", np.abs(predicted[:, 0] - targets[:, 0]) * 100),
            pl.Series("iptm_error", np.abs(predicted[:, 1] - targets[:, 1])),
            pl.Series("support_epitopes", support.get(method, np.full(rows.height, -1))),
        )
        details.append(detail)
        for key, part in detail.group_by(COHORT):
            called = part.filter(pl.col("y").is_in([0, 1]))
            record = dict(
                zip(COHORT, key),
                method=method,
                n=part.height,
                n_called=called.height,
                binders=int(called["y"].sum()),
                nonbinders=called.height - int(called["y"].sum()),
                plddt_MAE=part["plddt_error"].mean(),
                iptm_MAE=part["iptm_error"].mean(),
                gate_MAE=part["gate_error"].mean(),
                fallback_n=int((part["support_epitopes"] == 0).sum()),
            )
            record.update(
                AP=float(average_precision_score(called["y"], called["score"])),
                AUROC=float(roc_auc_score(called["y"], called["score"])),
            )
            records.append(record)
    return pl.concat(details), pl.DataFrame(records).sort(*COHORT, "method")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--features", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, nargs="+", required=True)
    parser.add_argument("--resource", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    start = time.monotonic()
    meta = pl.concat(
        [
            pl.read_csv(p, separator="\t", infer_schema_length=10000).select(
                "id", "mhc", "ranking_confidence"
            )
            for p in args.metadata
        ]
    )
    rows = pl.read_parquet(args.features).join(meta, on="id", how="left", validate="1:1")
    if rows["mhc"].null_count():
        raise ValueError("Missing model metadata")
    details, metrics = evaluate(rows, json.loads(args.resource.read_text()))
    args.output.mkdir(parents=True, exist_ok=False)
    details.write_parquet(args.output / "predictions.parquet")
    metrics.write_csv(args.output / "cohort_metrics.tsv", separator="\t")
    metrics.group_by("dataset", "stratum", "method").agg(
        pl.len().alias("cohorts"),
        pl.col("n").sum(),
        pl.col("n_called").sum(),
        pl.col("AP", "AUROC", "plddt_MAE", "iptm_MAE", "gate_MAE").mean(),
        pl.col("fallback_n").sum(),
    ).sort("dataset", "stratum", "method").write_csv(args.output / "summary.tsv", separator="\t")
    (args.output / "SOURCES.json").write_text(
        json.dumps(
            {
                "inputs": {
                    str(p): digest(p) for p in [args.features, args.resource, *args.metadata]
                },
                "models": rows.height,
                "seconds": time.monotonic() - start,
                "interpretation": "Exploratory leave-epitope-out confidence regressions; all 27 cohorts supply unlabelled training confidence, no binding labels in fits; frozen accepted score; small sets sorted by ID; no score-default change",
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
