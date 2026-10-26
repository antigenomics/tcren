"""Exploratory, grouped contact-profile evaluation; no default score changes.

Input tables are fresh CLI atom-contact/inventory exports and an accepted baseline.
All fitted transforms and mixture choices are confined to outer training folds.
"""

from pathlib import Path
import gzip
import json
import time
import warnings

import numpy as np
import polars as pl
import typer
from scipy.special import expit
from scipy.stats import mannwhitneyu
from threadpoolctl import threadpool_limits

SEED = 20261025
CS = (0.01, 0.1, 1.0)
WEIGHTS = (0.0, 0.25, 0.5, 1.0, 2.0, 4.0)
PRIORITY = ("TTDPSFLGRY", "NLVPMVATV")
KEY = ["id", "chain_id", "residue_index"]
BACKBONE = ["N", "CA", "C", "O", "OXT"]


def digest(path):
    from ..exports import _sha256

    return _sha256(Path(path))


def save_flat(table, path):
    """Columnar plus compressed, flat TSV deliverables, with no list cells."""
    table.write_parquet(path.with_suffix(".parquet"))
    with gzip.open(path.with_suffix(".tsv.gz"), "wt") as stream:
        table.write_csv(stream, separator="\t")


def connected_groups(alpha, beta):
    """Transitive exact-chain components, shared across epitope cohorts."""
    parent = list(range(len(alpha)))

    def root(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    seen = {}
    for i, pair in enumerate(zip(alpha, beta)):
        for chain, seq in enumerate(pair):
            if not seq:
                raise ValueError("Missing CDR3 sequence")
            key = chain, seq
            if key in seen:
                a, b = root(i), root(seen[key])
                parent[max(a, b)] = min(a, b)
            seen[key] = i
    return np.array([root(i) for i in range(len(alpha))])


def separation(y, x):
    """Two-sided tie-corrected asymptotic MW; positive effect favors binders."""
    y, x = np.asarray(y), np.asarray(x, dtype=float)
    keep = np.isfinite(x)
    a, b = x[keep & (y == 1)], x[keep & (y == 0)]
    result = dict(n_binder=len(a), n_nonbinder=len(b))
    if not len(a) or not len(b):
        return {
            **result,
            **dict.fromkeys(
                [
                    "p_value",
                    "rank_biserial",
                    "mean_binder",
                    "mean_nonbinder",
                    "mean_difference",
                    "median_binder",
                    "median_nonbinder",
                    "median_difference",
                ]
            ),
        }
    test = mannwhitneyu(a, b, alternative="two-sided", method="asymptotic")
    return dict(
        **result,
        p_value=float(test.pvalue),
        rank_biserial=float(2 * test.statistic / (len(a) * len(b)) - 1),
        mean_binder=float(a.mean()),
        mean_nonbinder=float(b.mean()),
        mean_difference=float(a.mean() - b.mean()),
        median_binder=float(np.median(a)),
        median_nonbinder=float(np.median(b)),
        median_difference=float(np.median(a) - np.median(b)),
    )


def ordered_matrix(models, wide, names):
    """Align every feature family explicitly before dropping model IDs."""
    frame = models.select("id").join(
        wide.select("id", *names), on="id", how="left", validate="1:1", maintain_order="left"
    )
    if not frame["id"].equals(models["id"]):
        raise ValueError("Feature row alignment changed")
    return frame.select(names).to_numpy().astype(float)


def contact_states(prof):
    """Separate endpoint masses from partial occupation; gaps remain absent."""
    return (
        prof.with_columns(
            (pl.col("contacted_atoms") > 0).cast(pl.Float64).alias("engaged"),
            pl.when((pl.col("occupation") > 0) & (pl.col("occupation") < 1))
            .then(pl.col("occupation"))
            .otherwise(0.0)
            .alias("partial_occupation"),
        )
        .unpivot(
            index=["id", "feature"],
            on=["contacted_atoms", "engaged", "fully_contacted", "partial_occupation"],
            variable_name="measure",
            value_name="value",
        )
        .with_columns(pl.concat_str("feature", "measure", separator=":").alias("feature"))
    )


def prepare(contacts_path, atoms_path, baseline_path, out, *, include_unlabelled=False):
    labels = pl.read_parquet(baseline_path)
    if not include_unlabelled:
        labels = labels.filter(pl.col("y").is_in([0, 1]))
    labels = labels.sort("id")
    if labels["id"].n_unique() != labels.height:
        raise ValueError("Baseline IDs must be unique")
    ids = labels.select("id")
    rename = {
        "pdb.id": "id",
        "chain.id": "chain_id",
        "chain.type": "chain",
        "residue.index": "residue_index",
        "residue.aa": "aa",
        "region.type": "region",
        "peptide.pos": "peptide_position",
        "atom.name": "atom",
    }
    atoms = (
        pl.scan_parquet(atoms_path)
        .rename(rename)
        .filter(pl.col("chain").is_in(["TRA", "TRB", "PEPTIDE"]))
        .select(*KEY, "chain", "aa", "region", "peptide_position", "atom")
        .collect()
        .join(ids, on="id", how="semi")
    )
    residues = atoms.select(*KEY, "chain", "aa", "region", "peptide_position").unique().sort(KEY)
    seq = (
        residues.group_by("id", "chain", maintain_order=True)
        .agg(pl.col("aa").str.join("").alias("sequence"))
        .pivot(on="chain", index="id", values="sequence")
        .rename(
            {"TRA": "tcr_alpha_sequence", "TRB": "tcr_beta_sequence", "PEPTIDE": "epitope_sequence"}
        )
    )
    cdrs = (
        residues.filter(pl.col("region") == "CDR3")
        .group_by("id", "chain", maintain_order=True)
        .agg(pl.col("aa").str.join("").alias("sequence"))
        .pivot(on="chain", index="id", values="sequence")
        .rename({"TRA": "cdr3_alpha", "TRB": "cdr3_beta"})
    )
    models = (
        labels.join(seq, on="id", validate="1:1")
        .join(cdrs, on="id", validate="1:1")
        .with_columns(
            pl.col("epitope").alias("cohort"),
            pl.when(pl.col("y") == 1)
            .then(pl.lit("binder"))
            .when(pl.col("y") == 0)
            .then(pl.lit("nonbinder"))
            .otherwise(pl.lit("ambiguous"))
            .alias("binder_status"),
        )
    )
    if (
        models.height != labels.height
        or not (models["epitope"] == models["epitope_sequence"]).all()
    ):
        raise ValueError("Raw sequence/model membership mismatch")
    models = models.sort("id")
    models = models.with_columns(
        pl.Series("group", connected_groups(models["cdr3_alpha"], models["cdr3_beta"]))
    )
    meta = models.select(
        "id",
        "cohort",
        "binder_status",
        "y",
        "tcr_alpha_sequence",
        "tcr_beta_sequence",
        "cdr3_alpha",
        "cdr3_beta",
        "epitope_sequence",
        "group",
    )
    contacts = (
        pl.read_parquet(contacts_path).rename({"pdb.id": "id"}).join(ids, on="id", how="semi")
    )
    if (
        set(contacts["role.from"]) != {"tcr"}
        or set(contacts["role.to"]) != {"peptide"}
        or not (contacts["dist"] <= 5).all()
    ):
        raise ValueError("Expected oriented TCR-peptide contacts at 5 A")
    cols = {
        "chain.id.from": "tcr_chain_id",
        "chain.type.from": "chain",
        "region.type.from": "region",
        "residue.index.from": "tcr_residue_index",
        "residue.index.pdb.from": "tcr_author_position",
        "residue.aa.from": "tcr_aa",
        "chain.id.to": "peptide_chain_id",
        "residue.index.to": "peptide_residue_index",
        "residue.aa.to": "peptide_aa",
        "atom.from": "tcr_atom",
        "atom.to": "peptide_atom",
        "dist": "distance_A",
    }
    edges = contacts.select("id", *cols).rename(cols)
    peptide_pos = residues.filter(pl.col("chain") == "PEPTIDE").select(
        "id",
        pl.col("chain_id").alias("peptide_chain_id"),
        pl.col("residue_index").alias("peptide_residue_index"),
        "peptide_position",
    )
    edges = edges.join(
        peptide_pos, on=["id", "peptide_chain_id", "peptide_residue_index"], validate="m:1"
    )
    save_flat(edges.join(meta, on="id", validate="m:1"), out / "atom_contacts")
    edgekey = [c for c in edges.columns if c not in ["tcr_atom", "peptide_atom", "distance_A"]]
    redges = edges.group_by(edgekey).agg(
        pl.len().alias("atom_pairs"),
        pl.col("tcr_atom").n_unique().alias("tcr_contacted_atoms"),
        pl.col("peptide_atom").n_unique().alias("peptide_contacted_atoms"),
        pl.col("distance_A").min().alias("minimum_distance_A"),
    )
    save_flat(redges.join(meta, on="id", validate="m:1"), out / "residue_contacts")
    inventory = atoms.group_by(KEY).agg(
        pl.col("atom").n_unique().alias("available_atoms"),
        pl.col("atom")
        .filter(pl.col("atom").is_in(BACKBONE))
        .n_unique()
        .alias("available_backbone_atoms"),
    )
    tcr = residues.filter(
        pl.col("chain").is_in(["TRA", "TRB"]) & pl.col("region").is_in(["CDR1", "CDR2", "CDR3"])
    )
    tcr = tcr.with_columns(
        (
            pl.col("residue_index")
            - pl.col("residue_index").min().over("id", "chain", "region")
            + 1
        ).alias("position"),
        pl.len().over("id", "chain", "region").alias("length"),
    )
    maxlens = tcr.group_by("chain", "region").agg(pl.col("length").max().alias("max_length"))
    tcr = tcr.join(maxlens, on=["chain", "region"], validate="m:1").with_columns(
        (
            pl.col("position")
            + pl.when(pl.col("position") > 3)
            .then(pl.col("max_length") - pl.col("length"))
            .otherwise(0)
        ).alias("aligned_position")
    )
    tc = (
        edges.group_by("id", "tcr_chain_id", "tcr_residue_index")
        .agg(
            pl.col("tcr_atom").n_unique().alias("contacted_atoms"),
            pl.col("tcr_atom")
            .filter(pl.col("tcr_atom").is_in(BACKBONE))
            .n_unique()
            .alias("contacted_backbone_atoms"),
            pl.col("peptide_position").n_unique().alias("degree"),
        )
        .rename({"tcr_chain_id": "chain_id", "tcr_residue_index": "residue_index"})
    )
    tcr = (
        tcr.join(tc, on=KEY, how="left", validate="1:1")
        .join(inventory, on=KEY, validate="1:1")
        .with_columns(pl.lit("tcr").alias("side"))
    )
    pep = (
        residues.filter(pl.col("chain") == "PEPTIDE")
        .drop("chain", "region")
        .join(pl.DataFrame({"chain": ["TRA", "TRB"]}), how="cross")
        .with_columns(
            pl.lit("PEPTIDE").alias("region"),
            pl.col("peptide_position").alias("position"),
            pl.col("peptide_position").alias("aligned_position"),
            pl.lit("peptide").alias("side"),
        )
    )
    pc = (
        edges.group_by("id", "peptide_chain_id", "peptide_residue_index", "chain")
        .agg(
            pl.col("peptide_atom").n_unique().alias("contacted_atoms"),
            pl.col("peptide_atom")
            .filter(pl.col("peptide_atom").is_in(BACKBONE))
            .n_unique()
            .alias("contacted_backbone_atoms"),
            pl.col("tcr_residue_index").n_unique().alias("degree"),
        )
        .rename({"peptide_chain_id": "chain_id", "peptide_residue_index": "residue_index"})
    )
    pep = pep.join(pc, on=KEY + ["chain"], how="left", validate="1:1").join(
        inventory, on=KEY, validate="m:1"
    )
    detail = (
        pl.concat([tcr, pep], how="diagonal_relaxed")
        .with_columns(pl.col("contacted_atoms", "contacted_backbone_atoms", "degree").fill_null(0))
        .with_columns(
            (pl.col("contacted_atoms") / pl.col("available_atoms")).alias("occupation"),
            (pl.col("contacted_atoms") == pl.col("available_atoms"))
            .cast(pl.Int64)
            .alias("fully_contacted"),
            (pl.col("contacted_atoms") - pl.col("contacted_backbone_atoms")).alias(
                "contacted_sidechain_atoms"
            ),
        )
        .join(meta, on="id", validate="m:1")
    )
    if not detail["occupation"].is_between(0, 1).all():
        raise ValueError("Contact atoms exceed observed inventory")
    save_flat(detail.sort("id", "side", "chain", "region", "position"), out / "residues")
    # CDR3 position profiles retain presence separately: gap != zero contact.
    prof = detail.filter((pl.col("side") == "peptide") | (pl.col("region") == "CDR3")).with_columns(
        pl.concat_str("side", "chain", "aligned_position", separator=":").alias("feature")
    )
    long = prof.select(
        "id",
        "cohort",
        "y",
        "feature",
        "aa",
        "contacted_atoms",
        "occupation",
        "fully_contacted",
        "degree",
    ).unpivot(
        index=["id", "cohort", "y", "feature", "aa"],
        on=["contacted_atoms", "occupation", "fully_contacted", "degree"],
        variable_name="measure",
        value_name="value",
    )
    summary = (
        detail.group_by("id", "side", "chain", "region")
        .agg(
            pl.col(
                "contacted_atoms",
                "contacted_backbone_atoms",
                "contacted_sidechain_atoms",
                "fully_contacted",
                "degree",
            ).sum(),
            pl.col("occupation").mean(),
            (pl.col("contacted_atoms") > 0).sum().alias("involved_residues"),
        )
        .with_columns(pl.concat_str("side", "chain", "region", separator=":").alias("feature"))
        .unpivot(
            index=["id", "feature"],
            on=[
                "contacted_atoms",
                "contacted_backbone_atoms",
                "contacted_sidechain_atoms",
                "fully_contacted",
                "degree",
                "occupation",
                "involved_residues",
            ],
            variable_name="measure",
            value_name="value",
        )
        .join(models.select("id", "cohort", "y"), on="id", validate="m:1")
    )
    statrows = []
    for scope, data in [("position", long), ("aggregate", summary)]:
        for keys, part in data.group_by("cohort", "feature", "measure"):
            statrows.append(
                dict(
                    scope=scope,
                    cohort=keys[0],
                    feature=keys[1],
                    measure=keys[2],
                    **separation(part["y"], part["value"]),
                )
            )
    pl.DataFrame(statrows).sort("cohort", "scope", "feature", "measure").write_csv(
        out / "separation_first.tsv", separator="\t"
    )
    # Fixed formulas are shown before any fitted readout.
    save_flat(summary, out / "aggregate_features")
    matrices = {}
    for measure in ["contacted_atoms", "occupation"]:
        wide = prof.select("id", "feature", measure).pivot(on="feature", index="id", values=measure)
        for side in ["tcr", "peptide", "joint"] if measure == "contacted_atoms" else ["joint"]:
            names = sorted(
                c
                for c in wide.columns
                if c != "id" and (side == "joint" or c.startswith(side + ":"))
            )
            matrices[side if measure == "contacted_atoms" else "occupation"] = (
                ordered_matrix(models, wide, names),
                names,
            )
    wide = (
        summary.with_columns(pl.concat_str("feature", "measure", separator=":").alias("feature"))
        .select("id", "feature", "value")
        .pivot(on="feature", index="id", values="value")
    )
    names = sorted(c for c in wide.columns if c != "id")
    matrices["commitment"] = (ordered_matrix(models, wide, names), names)
    x, names = matrices["tcr"]
    for arm, chain in [("alpha", "TRA"), ("beta", "TRB")]:
        indices = [j for j, name in enumerate(names) if name.startswith(f"tcr:{chain}:")]
        matrices[arm] = (x[:, indices], [names[j] for j in indices])
    wide = (
        contact_states(prof)
        .select("id", "feature", "value")
        .pivot(on="feature", index="id", values="value")
    )
    names = sorted(c for c in wide.columns if c != "id")
    matrices["states"] = (ordered_matrix(models, wide, names), names)
    models.write_parquet(out / "models.parquet")
    return models, matrices


def splits(models, n, seed):
    from sklearn.model_selection import StratifiedGroupKFold

    y = models["y"].to_numpy()
    groups = models["group"].to_numpy()
    strata = np.array([f"{e}:{v}" for e, v in zip(models["cohort"], y)])
    folds = list(
        StratifiedGroupKFold(n_splits=n, shuffle=True, random_state=seed).split(
            np.zeros(len(y)), strata, groups
        )
    )
    for train, test in folds:
        if set(groups[train]) & set(groups[test]):
            raise ValueError("Group leakage")
        for ep in set(models["cohort"]):
            for idx in [train, test]:
                if set(y[idx][models["cohort"].to_numpy()[idx] == ep]) != {0, 1}:
                    raise ValueError(f"Class missing in grouped fold: {ep}")
    return folds


def fit_predict(x, y, train, test, c):
    """Train-only imputation/scaling; explicit missing-position indicator."""
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        means = np.nanmean(x[train], axis=0)
    means = np.nan_to_num(means, nan=0.0)

    def design(idx):
        return np.concatenate(
            [np.where(np.isfinite(x[idx]), x[idx], means), ~np.isfinite(x[idx])], axis=1
        )

    scaler = StandardScaler().fit(design(train))
    a = scaler.transform(design(train))
    b = scaler.transform(design(test))
    model = LogisticRegression(C=c, solver="lbfgs", max_iter=1000, random_state=SEED).fit(
        a, y[train]
    )
    if model.n_iter_.max() >= 1000:
        raise RuntimeError("Logistic convergence limit reached")
    contribution = b * model.coef_[0]
    score = expit(model.intercept_[0] + contribution.sum(axis=1))
    np.testing.assert_allclose(score, model.predict_proba(b)[:, 1], atol=1e-12)
    params = dict(
        C=c,
        imputation=means.tolist(),
        mean=scaler.mean_.tolist(),
        scale=scaler.scale_.tolist(),
        coef=model.coef_[0].tolist(),
        intercept=float(model.intercept_[0]),
    )
    return score, contribution, params


def imputed_baselines(models, train, test):
    a = models.select("iptm", "prepared").to_numpy().astype(float)
    means = np.nanmean(a[train], axis=0)
    if not np.isfinite(means).all():
        raise ValueError("No baseline observations in training fold")
    return np.where(np.isfinite(a[test]), a[test], means)


def target_ap(models, values):
    from sklearn.metrics import average_precision_score

    return float(
        np.mean(
            [
                average_precision_score(
                    models.filter(pl.col("cohort") == e)["y"],
                    values[models["cohort"].to_numpy() == e],
                )
                for e in PRIORITY
            ]
        )
    )


def outer_fold(task):
    models, matrices, fold, train, test = task
    with threadpool_limits(limits=1):
        return _outer_fold(models, matrices, fold, train, test)


def _outer_fold(models, matrices, fold, train, test):
    started = time.monotonic()
    dev = models[train]
    y = models["y"].to_numpy()
    ep = models["cohort"].to_numpy()
    inner = splits(dev, 3, SEED + 100 + fold)
    base_inner = np.full((len(train), 2), np.nan)
    for it, iv in inner:
        base_inner[iv] = imputed_baselines(models, train[it], train[iv])
    base_test = imputed_baselines(models, train, test)
    scores = {"iptm": base_test[:, 0], "prepared": base_test[:, 1]}
    candidates = {}
    params = {}
    contributions = {}
    search = []
    choices = {}
    for arm, (x, names) in matrices.items():
        cands = []
        for c in CS:
            pred = np.full(len(train), np.nan)
            for it, iv in inner:
                for e in np.unique(ep):
                    a = train[it][ep[train[it]] == e]
                    b = train[iv][ep[train[iv]] == e]
                    pred[iv[ep[train[iv]] == e]] = fit_predict(x, y, a, b, c)[0]
            ap = target_ap(dev, pred)
            cands.append((ap, c, pred))
            search.append(
                dict(
                    fold=fold,
                    stage="regularization",
                    arm=arm,
                    C=c,
                    w1=None,
                    w2=None,
                    inner_target_AP=ap,
                )
            )
        _, c, inner_pred = max(cands, key=lambda z: (z[0], -z[1]))
        candidates[arm] = inner_pred
        choices[arm] = c
        result = np.full(len(test), np.nan)
        parts = np.zeros((len(test), 2 * len(names)))
        parameters = []
        for e in np.unique(ep):
            a = train[ep[train] == e]
            b = test[ep[test] == e]
            pred, contr, par = fit_predict(x, y, a, b, c)
            result[ep[test] == e] = pred
            parts[ep[test] == e] = contr
            parameters.append(
                dict(cohort=e, features=names + ["missing:" + n for n in names], **par)
            )
        scores[arm] = result
        contributions[arm] = parts
        params[arm] = parameters
    best = None
    best_base = None
    for arm, pred in candidates.items():
        for w1 in WEIGHTS:
            for w2 in WEIGHTS:
                mix = (base_inner[:, 0] + w1 * base_inner[:, 1] + w2 * pred) / (1 + w1 + w2)
                ap = target_ap(dev, mix)
                key = (ap, -w2, -w1, arm)
                search.append(
                    dict(
                        fold=fold,
                        stage="mixture",
                        arm=arm,
                        C=choices[arm],
                        w1=w1,
                        w2=w2,
                        inner_target_AP=ap,
                    )
                )
                if best is None or key > best[0]:
                    best = (key, arm, w1, w2)
                if w2 == 0 and (best_base is None or key > best_base[0]):
                    best_base = (key, arm, w1, w2)
    _, arm, w1, w2 = best
    scores["combined"] = (base_test[:, 0] + w1 * base_test[:, 1] + w2 * scores[arm]) / (1 + w1 + w2)
    scores["no_contact"] = (base_test[:, 0] + best_base[2] * base_test[:, 1]) / (1 + best_base[2])
    scores["selected_contact"] = scores[arm]
    rows = (
        models[test]
        .select("id", "cohort", "y", "group", "matched")
        .with_columns(pl.lit(fold).alias("fold"), *[pl.Series(k, v) for k, v in scores.items()])
    )
    names = matrices[arm][1]
    names = names + ["missing:" + n for n in names]
    cont = (
        pl.DataFrame(contributions[arm], schema=names)
        .with_columns(pl.Series("id", models[test]["id"]))
        .unpivot(index="id", variable_name="feature", value_name="logit_contribution")
        .with_columns(
            pl.lit(fold).alias("fold"),
            pl.lit(arm).alias("arm"),
            pl.lit(w2 / (1 + w1 + w2)).alias("contact_mixture_weight"),
        )
    )
    selected = dict(
        fold=fold,
        arm=arm,
        C=choices[arm],
        w1=w1,
        w2=w2,
        no_contact_w1=best_base[2],
        inner_target_AP=best[0][0],
        seconds=time.monotonic() - started,
        train_n=len(train),
        test_n=len(test),
    )
    # Component assignments, including inner folds, are retained for independent audit.
    inner_assign = []
    for j, (it, iv) in enumerate(inner):
        inner_assign.extend(
            dict(outer_fold=fold, inner_fold=j, id=s, role="validation") for s in dev[iv]["id"]
        )
    return rows, cont, selected, search, params, inner_assign


def evaluate(models, matrices, out, threads):
    from sklearn.metrics import average_precision_score, roc_auc_score

    from tcren.runtime import parallel_map

    folds = splits(models, 5, SEED)
    tasks = [(models, matrices, i, a, b) for i, (a, b) in enumerate(folds)]
    # Five coarse outer-fold jobs; native pools are one-thread per worker.
    results = parallel_map(outer_fold, tasks, threads)
    rows = pl.concat([r[0] for r in results]).sort("id")
    rows.write_parquet(out / "oof_predictions.parquet")
    rows.write_csv(out / "oof_predictions.tsv", separator="\t")
    if rows["id"].n_unique() != models.height or rows.height != models.height:
        raise ValueError("OOF row loss")
    names = ["iptm", "prepared", "no_contact", *matrices, "selected_contact", "combined"]
    records = []
    for cohort, p in rows.group_by("cohort"):
        for population, q in [("all_labelled", p), ("observed_baselines", p.filter("matched"))]:
            for arm in names:
                values = q[arm].to_numpy()
                y = q["y"].to_numpy()
                if not np.isfinite(values).all() or not ((values >= 0) & (values <= 1)).all():
                    raise ValueError("Invalid OOF score")
                records.append(
                    dict(
                        cohort=cohort[0],
                        population=population,
                        arm=arm,
                        n=q.height,
                        AP=float(average_precision_score(y, values)),
                        AUROC=float(roc_auc_score(y, values)),
                        **separation(y, values),
                    )
                )
    metrics = pl.DataFrame(records)
    for baseline in ["iptm", "prepared", "no_contact"]:
        ref = metrics.filter(pl.col("arm") == baseline).select(
            "cohort", "population", pl.col("AP").alias("reference_AP")
        )
        metrics = (
            metrics.join(ref, on=["cohort", "population"], validate="m:1")
            .with_columns((pl.col("AP") - pl.col("reference_AP")).alias("delta_AP_" + baseline))
            .drop("reference_AP")
        )
    metrics.sort("population", "cohort", "arm").write_csv(out / "metrics.tsv", separator="\t")
    cont = pl.concat([r[1] for r in results]).join(
        models.select("id", "cohort", "y", "cdr3_alpha", "cdr3_beta", "epitope_sequence"),
        on="id",
        validate="m:1",
    )
    save_flat(cont, out / "oof_contributions")
    summary = cont.group_by("cohort", "y", "arm", "feature").agg(
        pl.len().alias("n"),
        pl.col("logit_contribution").mean().alias("mean_logit_contribution"),
        pl.col("logit_contribution").median().alias("median_logit_contribution"),
        pl.col("logit_contribution").abs().mean().alias("mean_absolute_logit_contribution"),
        pl.col("contact_mixture_weight").mean(),
    )
    summary.sort(
        "cohort", "y", "mean_absolute_logit_contribution", descending=[False, False, True]
    ).write_csv(out / "contribution_summary.tsv", separator="\t")
    pl.DataFrame([s for r in results for s in r[3]]).write_csv(
        out / "inner_search.tsv", separator="\t"
    )
    pl.DataFrame([s for r in results for s in r[5]]).write_csv(
        out / "inner_folds.tsv", separator="\t"
    )
    (out / "selected.json").write_text(json.dumps([r[2] for r in results], indent=2) + "\n")
    (out / "fold_models.json").write_text(json.dumps([r[4] for r in results]) + "\n")
    return metrics


def command(
    contacts: Path = typer.Option(..., exists=True),
    atoms: Path = typer.Option(..., exists=True),
    baseline: Path = typer.Option(..., exists=True),
    output: Path = typer.Option(..., "--output", "-o"),
    threads: int = typer.Option(0, "--threads", "-t"),
):
    """Evaluate simple residue profiles with nested exact-CDR3-grouped folds."""
    from tcren.runtime import resolve_threads

    budget = resolve_threads(threads)
    output.mkdir(parents=True, exist_ok=False)
    start = time.monotonic()
    models, matrices = prepare(contacts, atoms, baseline, output)
    print(
        f"First separation readout ready: {models.height} models; "
        + str(output / "separation_first.tsv"),
        flush=True,
    )
    metrics = evaluate(models, matrices, output, budget)
    provenance = dict(
        seed=SEED,
        outer_folds=5,
        inner_folds=3,
        workers=min(budget, 5),
        C=CS,
        weights=WEIGHTS,
        priority=PRIORITY,
        models=models.height,
        groups=models["group"].n_unique(),
        seconds=time.monotonic() - start,
        matrix_shapes={k: list(v[0].shape) for k, v in matrices.items()},
        inputs={str(p): digest(p) for p in [contacts, atoms, baseline]},
        interpretation="Internal exploratory OOF after prior panel inspection; no independent external validation; no default change",
        p_values="Two-sided asymptotic tie-corrected Mann-Whitney; model-level descriptive P-values; no multiplicity correction; shared sequence groups can be dependent",
    )
    (output / "PROVENANCE.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print(metrics.filter(pl.col("arm") == "combined"), flush=True)
