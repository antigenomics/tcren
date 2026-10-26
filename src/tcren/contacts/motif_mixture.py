"""Fisher contact-motif enrichment and a three-proportion linear score.

S = w_I*ipTM + w_B*recognition_beta1/E + w_M*mock_enrichment, w>=0, sum(w)=1.
No arbitrary powers, peptide predictors, regression features or new dependencies.
"""

from pathlib import Path
import json
import time

import numpy as np
import polars as pl
import typer
from scipy.special import gammaincc

from .residue_surprise import CELL, KEY, score_residues
from .profile import digest, connected_groups
from .atom_formula import ap_columns
from ..runtime import parallel_map, resolve_threads

FEATURES = ["iptm", "betti", "enrichment"]


def fisher_summary(rows):
    """Raw Fisher statistic and independent-test diagnostic; absent support abstains."""
    p = rows["tail_probability"].drop_nulls().to_numpy()
    if not np.isfinite(p).all() or np.any((p <= 0) | (p > 1)):
        raise ValueError("Local tail probabilities must be in (0,1]")
    out = rows.group_by("id").agg(
        (-2 * pl.col("tail_probability").log().fill_null(0).sum()).alias("fisher_statistic"),
        pl.col("tail_probability").count().alias("n_local_tests"),
        pl.len().alias("n_residues"),
    )
    independent = gammaincc(out["n_local_tests"].to_numpy(), out["fisher_statistic"].to_numpy() / 2)
    return out.with_columns(pl.Series("fisher_independent_p", independent)).sort("id")


def empirical_survival(query, reference, weights):
    """Inclusive plus-one empirical upper tail, with explicit positive model weights."""
    q, r, w = map(lambda a: np.asarray(a, dtype=float), (query, reference, weights))
    if (
        not np.isfinite(q).all()
        or not np.isfinite(r).all()
        or not np.isfinite(w).all()
        or len(r) != len(w)
        or not len(r)
        or np.any(w <= 0)
    ):
        raise ValueError("Finite reference statistics and positive weights required")
    order = np.argsort(r, kind="stable")
    r, w = r[order], w[order]
    tail = np.r_[np.cumsum(w[::-1])[::-1], 0.0]
    return np.clip((1 + tail[np.searchsorted(r, q, side="left")]) / (1 + w.sum()), 0, 1)


def fisher_command(
    features: Path = typer.Option(..., exists=True),
    reference: Path = typer.Option(..., exists=True),
    controls: Path = typer.Option(..., exists=True),
    output: Path = typer.Option(..., "--output", "-o"),
    statistic: str = typer.Option("cycles"),
    joint_residue: bool = False,
    positive_reference: Path | None = typer.Option(None, exists=True),
    condition_epitope: bool = False,
):
    """Calibrate Fisher sums against whole control graphs; retain every eligible control."""
    output.mkdir(exist_ok=False)
    resource = json.loads(reference.read_text())
    if resource.get("reference_label") != 0 or resource.get("input_sha256") != digest(controls):
        raise ValueError("Control inventory must match the frozen negative reference")
    selected = pl.read_csv(reference.with_suffix(".models.tsv"), separator="\t").select(
        "id", "weight", *(["PEPTIDE"] if condition_epitope else [])
    )
    raw = pl.read_parquet(controls).join(selected, on="id", how="inner", validate="m:1")
    if raw["id"].n_unique() != selected.height:
        raise ValueError("Missing eligible control graphs")
    detail = score_residues(raw, resource, statistic, joint_residue=joint_residue)
    # Whole-model leave-one-out ranks for the empirical null. This is vectorized
    # subtraction of own sufficient counts, not a fit or loop per query model.
    own = (
        raw.filter("complete")
        .group_by("id", *CELL, statistic)
        .agg(pl.col("weight").sum().alias("own_mass"))
        .sort("id", *CELL, statistic)
        .with_columns(
            pl.col("own_mass")
            .sum()
            .over(["id", "chain", "region"] if joint_residue else ["id", *CELL])
            .alias("own_total"),
            pl.col("own_mass").reverse().cum_sum().reverse().over(["id", *CELL]).alias("own_tail"),
        )
    )
    own = own.with_columns(pl.col("available").cast(pl.Int64), pl.col(statistic).cast(pl.Int64))
    detail = detail.join(own, on=["id", *CELL, statistic], how="left", validate="m:1")
    detail = detail.with_columns(
        pl.when(pl.col("complete") & (pl.col("reference_mass") > pl.col("own_total")))
        .then(
            pl.when(pl.col(statistic) == 0)
            .then(1.0)
            .otherwise(
                (
                    (1 + pl.col("tail_count") - pl.col("own_tail"))
                    / (1 + pl.col("reference_mass") - pl.col("own_total"))
                ).clip(0, 1)
            )
        )
        .otherwise(None)
        .alias("tail_probability")
    )
    positive = (
        json.loads(positive_reference.read_text()) if isinstance(positive_reference, Path) else None
    )
    if positive is not None and positive.get("reference_label") != 1:
        raise ValueError("Declared source-positive reference required")

    def summarize(base, local):
        if positive is None:
            return local, fisher_summary(local)
        from .residue_surprise import KEY

        bound = score_residues(base, positive, statistic, joint_residue=joint_residue).select(
            *KEY, pl.col("tail_probability").alias("positive_tail_probability")
        )
        local = local.join(bound, on=KEY, validate="1:1", maintain_order="left").with_columns(
            (pl.col("positive_tail_probability") / pl.col("tail_probability"))
            .log(2)
            .alias("evidence_bits")
        )
        summary = local.group_by("id").agg(
            (2 * np.log(2) * pl.col("evidence_bits").sum()).alias("fisher_statistic"),
            pl.col("evidence_bits").count().alias("n_local_tests"),
            pl.len().alias("n_residues"),
        )
        return local, summary.with_columns(
            pl.lit(None, dtype=pl.Float64).alias("fisher_independent_p")
        )

    detail, null = summarize(raw, detail)
    null = null.join(selected, on="id", validate="1:1")
    if null.height != selected.height:
        raise ValueError("Control row conservation failed")
    null.write_parquet(output / "control_fisher.parquet")
    queries = pl.read_parquet(features)
    for table in [raw, queries]:
        scopes = (
            table.filter(pl.col("chain") != "PEPTIDE")
            .select("id", "chain", "region")
            .unique()
            .group_by("id")
            .len()
        )
        if scopes.height != table["id"].n_unique() or (scopes["len"] != 6).any():
            raise ValueError("Every graph requires all six CDR loops")
    query = score_residues(queries, resource, statistic, joint_residue=joint_residue)
    query, summary = summarize(queries, query)
    query.write_parquet(output / "residues.parquet")
    supported = null  # A control without supported local tests has F=0; retain it.
    null_n = np.full(summary.height, supported.height, dtype=int)
    specific = np.zeros(summary.height, dtype=bool)
    if condition_epitope:
        identities = pl.read_csv(features.with_suffix(".identities.tsv"), separator="\t").select(
            "id", "PEPTIDE"
        )
        summary = summary.join(identities, on="id", validate="1:1", maintain_order="left")
        if summary["PEPTIDE"].null_count():
            raise ValueError("Missing query epitope identity")
        p = np.empty(summary.height)
        for peptide in summary["PEPTIDE"].unique():
            mask = (summary["PEPTIDE"] == peptide).to_numpy()
            pool = supported.filter(pl.col("PEPTIDE") == peptide)
            if pool.is_empty():
                pool = supported
            else:
                specific[mask] = True
            null_n[mask] = pool.height
            p[mask] = empirical_survival(
                summary.filter(pl.Series(mask))["fisher_statistic"],
                pool["fisher_statistic"],
                pool["weight"],
            )
    else:
        p = empirical_survival(
            summary["fisher_statistic"], supported["fisher_statistic"], supported["weight"]
        )
    summary = summary.with_columns(
        pl.Series("mock_fisher_p", p),
        pl.Series("null_models", null_n),
        pl.Series("epitope_specific_null", specific),
    ).with_columns(
        pl.when(pl.col("n_local_tests") > 0)
        .then(1 - pl.col("mock_fisher_p"))
        .otherwise(0.0)
        .alias("enrichment")
    )
    summary.write_csv(output / "scores.tsv", separator="\t")
    (output / "reference.json").write_text(
        json.dumps(
            dict(
                schema="tcren-motif-fisher-1",
                statistic=statistic,
                condition_epitope=condition_epitope,
                positive_reference_sha256=digest(positive_reference)
                if positive is not None
                else None,
                evidence="Signed sum of log bound/mock tail-frequency ratios"
                if positive is not None
                else "Fisher negative-tail surprise",
                local_tail="Joint residue-type/contact event within chain/region"
                if joint_residue
                else "Contact tail conditional on residue type/atom count",
                residue_reference_sha256=digest(reference),
                controls_sha256=digest(controls),
                control_ids=supported["id"].to_list(),
                fisher=supported["fisher_statistic"].to_list(),
                weights=supported["weight"].to_list(),
                formula=(
                    "F=2 sum ln(p_bound_tail/p_mock_tail); "
                    if positive is not None
                    else "F=-2 sum ln(p_residue); "
                )
                + "enrichment=1-(1+sum w*I(F_control>=F))/(1+sum w)",
                null="Whole control graphs with whole-model leave-one-out local ranks; empirical calibration, not independent chi-square",
                interpretation="Dependent/discrete motif tails; diagnostic reference enrichment, not calibrated binding probability",
            ),
            indent=2,
        )
        + "\n"
    )


def simplex_weights(resolution=100):
    """A two-dimensional simplex, at one-percent coefficient precision by default."""
    if resolution < 1 or int(resolution) != resolution:
        raise ValueError("Positive integer coefficient resolution required")
    return (
        np.array(
            [
                (i, j, resolution - i - j)
                for i in range(resolution + 1)
                for j in range(resolution - i + 1)
            ],
            dtype=float,
        )
        / resolution
    )


def linear_score(x, weights):
    """The entire final model, with identical arithmetic for batch and single fits."""
    x, w = np.asarray(x, dtype=float), np.asarray(weights, dtype=float)
    if (
        w.ndim not in [1, 2]
        or w.shape[-1] != 3
        or not np.isfinite(w).all()
        or np.any(w < 0)
        or not np.allclose(w.sum(axis=-1), 1, atol=1e-12, rtol=0)
    ):
        raise ValueError("Three nonnegative proportions must sum to one")
    if x.shape[-1] != 3 or not np.isfinite(x).all() or np.any((x < 0) | (x > 1)):
        raise ValueError("Three finite component scores in [0,1] required")
    if w.ndim == 1:
        return x[..., 0] * w[0] + x[..., 1] * w[1] + x[..., 2] * w[2]
    return x[:, 0, None] * w[:, 0] + x[:, 1, None] * w[:, 1] + x[:, 2, None] * w[:, 2]


def _weight_chunk(task):
    x, y, cohort, w = task
    scores = linear_score(x, w)
    return np.mean(
        [ap_columns(y[cohort == c], scores[cohort == c]) for c in np.unique(cohort)], axis=0
    )


def fit_weights(x, y, cohort, threads=0, weights=None):
    """Fit only the three proportions by equal-cohort macro AP; no fitted powers."""
    x, y, cohort = np.asarray(x), np.asarray(y), np.asarray(cohort)
    weights = simplex_weights() if weights is None else weights
    budget = resolve_threads(threads)
    tasks = [
        (x, y, cohort, w) for w in np.array_split(weights, min(budget, len(weights))) if len(w)
    ]
    ap = np.concatenate(parallel_map(_weight_chunk, tasks, budget))
    tied = np.flatnonzero(ap >= ap.max() - 1e-12)
    best = min(tied, key=lambda i: (np.count_nonzero(weights[i]), -weights[i, 0], weights[i, 2]))
    return weights[best], ap


def _impute(x, train):
    means = np.nanmean(x[train], axis=0)
    if not np.isfinite(means).all():
        raise ValueError("Training lacks an observed component")
    return np.where(np.isfinite(x), x, means), means


def _partition_fit(task):
    """One coarse fold/full-fit job; all native work stays serial inside it."""
    fold, x, y, cohort, a, b = task
    imputed, mean = _impute(x, a)
    w, ap = fit_weights(imputed[a], y[a], cohort[a], threads=1)
    return (
        fold,
        b,
        linear_score(imputed[b], w),
        imputed[b],
        w,
        mean,
        len(a),
        float(ap.max()),
        ap if fold == -1 else None,
    )


def _metrics(rows):
    from sklearn.metrics import average_precision_score, roc_auc_score

    records = []
    for key, part in rows.filter(pl.col("y").is_in([0, 1])).group_by(
        "dataset", "stratum", "epitope"
    ):
        for population, subset in [
            ("all_labelled", part),
            ("observed_components", part.filter("components_observed")),
        ]:
            if subset.is_empty() or subset["y"].n_unique() != 2:
                continue
            y = subset["y"].to_numpy()
            for arm in [
                "iptm_scored",
                "betti_scored",
                "enrichment_scored",
                "prepared_scored",
                "mixture_fullfit",
                "mixture_oof",
            ]:
                if subset[arm].null_count() or not subset[arm].is_finite().all():
                    continue
                records.append(
                    dict(
                        zip(["dataset", "stratum", "epitope"], key),
                        population=population,
                        arm=arm,
                        n=subset.height,
                        binders=int(y.sum()),
                        nonbinders=int(len(y) - y.sum()),
                        AP=float(average_precision_score(y, subset[arm])),
                        AUROC=float(roc_auc_score(y, subset[arm])),
                    )
                )
    return pl.DataFrame(records).sort("population", "dataset", "stratum", "epitope", "arm")


def diagnose_command(
    predictions: Path = typer.Option(..., exists=True),
    folds: Path = typer.Option(..., exists=True),
    fold_weights: Path = typer.Option(..., exists=True),
    resource: Path = typer.Option(..., exists=True),
    output: Path = typer.Option(..., "--output", "-o"),
):
    """Separate within-fold ranking from between-fold shifts of the same formula."""
    from sklearn.metrics import average_precision_score, roc_auc_score

    output.mkdir(exist_ok=False)
    rows = pl.read_parquet(predictions).join(
        pl.read_parquet(folds).select("id", "fold"), on="id", how="inner", validate="1:1"
    )
    fits = json.loads(fold_weights.read_text())
    model = json.loads(resource.read_text())
    records, shifts = [], []
    for epitope, part in rows.group_by("epitope"):
        # The same feature population under each weight vector isolates a
        # coefficient-induced shift from differences in fold composition.
        common = part.select(FEATURES).to_numpy()
        common = np.where(np.isfinite(common), common, model["imputation"])
        anchor = linear_score(common, model["weights"])
        for fit in fits:
            fold = fit["fold"]
            s = part.filter(pl.col("fold") == fold)
            if s["y"].n_unique() != 2:
                raise ValueError("Within-fold AP needs both labels")
            for arm in ["prepared_scored", "mixture_oof", "mixture_fullfit"]:
                records.append(
                    dict(
                        epitope=epitope[0],
                        fold=fold,
                        arm=arm,
                        n=s.height,
                        binders=int(s["y"].sum()),
                        nonbinders=int(s.height - s["y"].sum()),
                        AP=float(average_precision_score(s["y"], s[arm])),
                        AUROC=float(roc_auc_score(s["y"], s[arm])),
                    )
                )
            delta = linear_score(common, fit["weights"]) - anchor
            shifts.append(
                dict(
                    epitope=epitope[0],
                    fold=fold,
                    n=part.height,
                    w_I=fit["weights"][0],
                    w_B=fit["weights"][1],
                    w_M=fit["weights"][2],
                    common_mean_shift=float(delta.mean()),
                    common_shift_sd=float(delta.std()),
                )
            )
    metrics = pl.DataFrame(records).sort("epitope", "fold", "arm")
    metrics.write_csv(output / "within_fold.tsv", separator="\t")
    metrics.group_by("epitope", "arm").agg(
        pl.len().alias("folds"), pl.col("n").sum(), pl.col("AP", "AUROC").mean()
    ).sort("epitope", "arm").write_csv(output / "fold_mean.tsv", separator="\t")
    pl.DataFrame(shifts).sort("epitope", "fold").write_csv(
        output / "scale_shifts.tsv", separator="\t"
    )


def fit_command(
    baseline: Path = typer.Option(..., exists=True),
    topology: list[Path] = typer.Option(..., exists=True),
    motifs: Path = typer.Option(..., exists=True),
    residues: Path | None = typer.Option(None, exists=True),
    folds: Path | None = typer.Option(None, exists=True),
    output: Path = typer.Option(..., "--output", "-o"),
    threads: int = typer.Option(0, "--threads", "-t"),
):
    """Fit one global ipTM/complex-Betti/motif weight triplet; folds are optional diagnostics."""
    output.mkdir(exist_ok=False)
    start = time.monotonic()
    budget = resolve_threads(threads)
    rows = pl.read_parquet(baseline).sort("id")
    g = pl.concat([pl.read_csv(p, separator="\t") for p in topology]).select(
        pl.col("pdb.id").alias("id"),
        pl.col("complex_cyclo_frac").alias("betti"),
        "receptor_contacts_available",
    )
    m = pl.read_csv(motifs, separator="\t").select("id", "enrichment")
    if set(rows["id"]) != set(g["id"]) or set(rows["id"]) != set(m["id"]):
        raise ValueError("Component model membership differs")
    rows = rows.join(g, on="id", validate="1:1", maintain_order="left").join(
        m, on="id", validate="1:1", maintain_order="left"
    )
    x = rows.select(FEATURES).to_numpy().astype(float)
    if np.any(np.isfinite(x) & ((x < 0) | (x > 1))):
        raise ValueError("Components must be bounded in [0,1]")
    eligible = (rows["dataset"] == "tcrvdb") & rows["y"].is_in([0, 1]).fill_null(False)
    indices = np.flatnonzero(eligible.to_numpy())
    train_rows = rows[indices]
    if folds is not None:
        if residues is None:
            raise ValueError("--residues is required with --folds")
        d = (
            pl.read_parquet(residues)
            .sort(KEY)
            .filter(pl.col("region") == "CDR3")
            .group_by("id", "chain")
            .agg(pl.col("aa").str.join("").alias("seq"))
            .pivot(on="chain", index="id", values="seq")
            .rename({"TRA": "cdr3_alpha", "TRB": "cdr3_beta"})
        )
        rows = rows.join(d, on="id", validate="1:1", maintain_order="left")
        train_rows = rows[indices]
        saved = pl.read_parquet(folds).select("id", "fold")
        if set(saved["id"]) != set(train_rows["id"]):
            raise ValueError("Recorded fold population changed; no resplit is allowed")
        train_rows = train_rows.join(saved, on="id", validate="1:1", maintain_order="left")
        group = connected_groups(train_rows["cdr3_alpha"], train_rows["cdr3_beta"])
        if (
            train_rows.with_columns(pl.Series("_g", group))
            .group_by("_g")
            .agg(pl.col("fold").n_unique())["fold"]
            .max()
            != 1
        ):
            raise ValueError("Recorded folds split an exact receptor chain component")
    xt = x[indices]
    y = train_rows["y"].to_numpy()
    cohort = train_rows["epitope"].to_numpy()
    oof = np.full(len(rows), np.nan)
    scored = np.full_like(x, np.nan)
    fits = []
    if folds is None:
        values, full_mean = _impute(xt, np.arange(len(xt)))
        full_w, full_ap = fit_weights(values, y, cohort, threads=budget)
        scored[indices] = values
    else:
        tasks = []
        for fold in sorted(train_rows["fold"].unique()):
            a = np.flatnonzero(train_rows["fold"].to_numpy() != fold)
            b = np.flatnonzero(train_rows["fold"].to_numpy() == fold)
            tasks.append((int(fold), xt, y, cohort, a, b))
        tasks.append((-1, xt, y, cohort, np.arange(len(xt)), np.array([], dtype=int)))
        results = parallel_map(_partition_fit, tasks, budget)
        for fold, b, pred, values, w, mean, n, objective, ap in results:
            if fold == -1:
                full_w, full_mean, full_ap = w, mean, ap
                continue
            oof[indices[b]] = pred
            scored[indices[b]] = values
            fits.append(
                dict(
                    fold=fold,
                    weights=w.tolist(),
                    imputation=mean.tolist(),
                    train_n=n,
                    test_n=len(b),
                    training_macro_AP=objective,
                )
            )
            print("fold", fold, "weights", w.tolist(), flush=True)
    w, mean, ap = full_w, full_mean, full_ap
    imputed = np.where(np.isfinite(x), x, mean)
    final = linear_score(imputed, w)
    others = np.flatnonzero(~eligible.to_numpy())
    scored[others] = imputed[others]
    conf = np.where(rows["stratum"].to_numpy() == "template_absent", 0.90, 0.55)
    prepared = conf * scored[:, 0] + (1 - conf) * scored[:, 1]
    old = rows["prepared"].to_numpy()
    finite = np.isfinite(old) & np.isfinite(prepared)
    np.testing.assert_allclose(old[finite], prepared[finite], atol=1e-12, rtol=0)
    rows = rows.with_columns(
        pl.Series("mixture_fullfit", final),
        pl.Series("mixture_oof", oof).fill_nan(None),
        pl.Series("components_observed", np.isfinite(x).all(1)),
        *[pl.Series(f + "_missing", ~np.isfinite(x[:, j])) for j, f in enumerate(FEATURES)],
        *[pl.Series(f + "_scored", scored[:, j]) for j, f in enumerate(FEATURES)],
        pl.Series("prepared_scored", prepared),
    )
    rows.write_parquet(output / "predictions.parquet")
    metrics = _metrics(rows)
    metrics.write_csv(output / "metrics.tsv", separator="\t")
    metrics.group_by("population", "dataset", "stratum", "arm").agg(
        pl.len().alias("cohorts"), pl.col("n").sum(), pl.col("AP", "AUROC").mean()
    ).sort("population", "dataset", "stratum", "arm").write_csv(
        output / "summary.tsv", separator="\t"
    )
    resource = dict(
        schema="tcren-motif-mixture-1",
        features=FEATURES,
        weights=w.tolist(),
        imputation=mean.tolist(),
        formula="S=w_I*ipTM+w_B*(recognition beta1/E)+w_M*(1-p_mock_Fisher); weights>=0, sum=1",
        objective="Equal-cohort macro AP across five labelled TCRvdb epitopes",
        coefficient_resolution=0.01,
        weight_estimation="Direct global training macro AP maximizer",
        training_macro_AP=float(ap.max()),
        training_n=len(indices),
        workers=min(budget, 6 if folds is not None else len(ap)),
        seconds=time.monotonic() - start,
        inputs={
            str(p): digest(p)
            for p in [baseline, *topology, motifs, residues, folds]
            if p is not None
        },
        assessment=(
            "Recorded five-fold exact-chain-grouped OOF diagnostic; "
            if folds is not None
            else "No folds; "
        )
        + "global fit is in-sample; VDJdb is a transferred readout, not untouched external validation",
    )
    (output / "weights.json").write_text(json.dumps(resource, indent=2) + "\n")
    (output / "fold_weights.json").write_text(json.dumps(fits, indent=2) + "\n")
    pl.DataFrame(simplex_weights(), schema=FEATURES).with_columns(
        pl.Series("training_macro_AP", ap)
    ).write_csv(output / "full_fit_objective.tsv", separator="\t")
    if folds is not None:
        diagnose_command(
            predictions=output / "predictions.parquet",
            folds=folds,
            fold_weights=output / "fold_weights.json",
            resource=output / "weights.json",
            output=output / "fold-diagnosis",
        )


def _cohort_mask(table, group):
    mask = table["dataset"] == group["dataset"]
    if group["stratum"] is not None:
        mask = mask & (table["stratum"] == group["stratum"])
    return mask.to_numpy()


def _apply_resource(table, model):
    x = table.select(FEATURES).to_numpy().astype(float)
    missing = ~np.isfinite(x)
    values, score = np.full_like(x, np.nan), np.full(len(x), np.nan)
    groups = model.get("groups", [model])
    assigned = np.zeros(len(x), dtype=int)
    for group in groups:
        mask = _cohort_mask(table, group) if "groups" in model else np.ones(len(x), dtype=bool)
        values[mask] = np.where(missing[mask], group["imputation"], x[mask])
        score[mask] = linear_score(values[mask], group["weights"])
        assigned[mask] += 1
    if np.any(assigned != 1):
        raise ValueError("Every row must match exactly one fitted cohort")
    score[missing.all(1)] = np.nan
    return values, score, missing


def fit_cohorts_command(
    features: Path = typer.Option(..., exists=True),
    output: Path = typer.Option(..., "--output", "-o"),
    threads: int = typer.Option(0, "--threads", "-t"),
):
    """Fit three coefficients separately for TCRvdb and both VDJdb template strata; no folds."""
    output.mkdir(exist_ok=False)
    start = time.monotonic()
    rows = pl.read_parquet(features).sort("id")
    groups = [
        dict(name="tcrvdb", dataset="tcrvdb", stratum=None),
        dict(name="vdjdb_templated", dataset="vdjdb", stratum="template_present"),
        dict(name="vdjdb_nontemplated", dataset="vdjdb", stratum="template_absent"),
    ]
    objectives = []
    for group in groups:
        train = rows.filter(pl.Series(_cohort_mask(rows, group)) & pl.col("y").is_in([0, 1]))
        if train.is_empty() or train["y"].n_unique() != 2:
            raise ValueError("Each fitted cohort needs both labels")
        x, mean = _impute(train.select(FEATURES).to_numpy(), np.arange(train.height))
        w, ap = fit_weights(x, train["y"].to_numpy(), train["epitope"].to_numpy(), threads=threads)
        group.update(
            weights=w.tolist(),
            imputation=mean.tolist(),
            training_n=train.height,
            epitopes=train["epitope"].n_unique(),
            training_macro_AP=float(ap.max()),
        )
        objectives.append(
            pl.DataFrame(simplex_weights(), schema=FEATURES).with_columns(
                pl.lit(group["name"]).alias("cohort"), pl.Series("training_macro_AP", ap)
            )
        )
    model = dict(
        schema="tcren-cohort-motif-mixture-1",
        features=FEATURES,
        groups=groups,
        formula="S=w_I*ipTM+w_B*(complex beta1/E)+w_M*enrichment; weights>=0, sum=1",
        objective="Equal-epitope macro AP independently in each declared cohort",
        coefficient_resolution=0.01,
        input_sha256=digest(features),
        assessment="Full-fit in-sample accuracy in all three cohorts; no folds or score-scale transform",
        seconds=time.monotonic() - start,
    )
    values, score, missing = _apply_resource(rows, model)
    conf = np.where(rows["stratum"].to_numpy() == "template_absent", 0.90, 0.55)
    rows = rows.with_columns(
        pl.Series("mixture_fullfit", score),
        pl.lit(None, dtype=pl.Float64).alias("mixture_oof"),
        pl.Series("components_observed", ~missing.any(1)),
        *[pl.Series(f + "_scored", values[:, j]) for j, f in enumerate(FEATURES)],
        *[pl.Series(f + "_missing", missing[:, j]) for j, f in enumerate(FEATURES)],
        pl.Series("prepared_scored", conf * values[:, 0] + (1 - conf) * values[:, 1]),
    )
    rows.write_parquet(output / "predictions.parquet")
    metrics = _metrics(rows)
    metrics.write_csv(output / "metrics.tsv", separator="\t")
    metrics.group_by("population", "dataset", "stratum", "arm").agg(
        pl.len().alias("cohorts"), pl.col("n").sum(), pl.col("AP", "AUROC").mean()
    ).sort("population", "dataset", "stratum", "arm").write_csv(
        output / "summary.tsv", separator="\t"
    )
    pl.concat(objectives).write_csv(output / "objectives.tsv", separator="\t")
    (output / "weights.json").write_text(json.dumps(model, indent=2) + "\n")


def score_command(
    features: Path = typer.Option(..., exists=True),
    resource: Path = typer.Option(..., exists=True),
    output: Path = typer.Option(..., "--output", "-o"),
):
    """Apply one frozen three-proportion formula; no score-time fitting."""
    if output.exists():
        raise ValueError("Output already exists")
    model = json.loads(resource.read_text())
    if model.get("schema") not in ["tcren-motif-mixture-1", "tcren-cohort-motif-mixture-1"]:
        raise ValueError("Unknown mixture resource")
    table = (
        pl.read_parquet(features)
        if features.suffix == ".parquet"
        else pl.read_csv(features, separator="\t")
    )
    _, score, missing = _apply_resource(table, model)
    table = table.with_columns(
        pl.Series("score", score),
        *[pl.Series(f + "_missing", missing[:, j]) for j, f in enumerate(FEATURES)],
    )
    table.write_csv(output, separator="\t")
