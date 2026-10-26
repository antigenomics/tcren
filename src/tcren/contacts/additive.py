"""Small, explicit additive contact search; no change to production scores.

CLI consumes fresh atom, inventory and recognition-topology exports. Exact-chain
components stay together in five outer folds. Selection maximizes equal-epitope
mean AP; p-values are reported afterwards, never used to select coefficients.
"""

from pathlib import Path
import json

import numpy as np
import polars as pl
import typer
from sklearn.metrics import average_precision_score

from .profile import connected_groups, digest, separation, splits, SEED
from ..runtime import parallel_map, resolve_threads

FEATURES = (
    "max_pairs",
    "max_atoms",
    "peptide_atoms",
    "pairs_per_peptide_atom",
    "tcr_atoms_per_peptide_atom",
)
WEIGHTS = np.array([0.0] + [m * 10.0**e for e in range(-5, 2) for m in (1, 2, 5)])


def oriented_contacts(contacts):
    """Select peptide contacts and orient every paired endpoint field TCR→peptide."""
    c = contacts.filter(pl.col("interface") == "tcr_peptide")
    forward = (pl.col("role.from") == "tcr") & (pl.col("role.to") == "peptide")
    reverse = (pl.col("role.from") == "peptide") & (pl.col("role.to") == "tcr")
    if (not c.select((forward | reverse).fill_null(False).all()).item()
            or c["dist"].null_count() or not c["dist"].is_between(0, 5).all()):
        raise ValueError("Expected TCR/peptide heavy-atom contacts at 5 Angstrom")
    if c.select(reverse.any()).item():
        if "chain.type.to" not in c.columns:
            raise ValueError("Reversed contacts require paired chain.type endpoint columns")
        fields = [name[:-5] for name in c.columns
                  if name.endswith(".from") and name[:-5] + ".to" in c.columns]
        c = c.with_columns(*[
            pl.when(reverse).then(pl.col(field + "." + other))
            .otherwise(pl.col(field + "." + side)).alias(field + "." + side)
            for field in fields for side, other in [("from", "to"), ("to", "from")]])
    return c


def contact_counts(contacts, ids):
    """Count physical pairs and unique endpoint atoms across either CDR chain."""
    c = oriented_contacts(contacts)
    keys = ["pdb.id"] + [
        f"{k}.{s}" for s in ("from", "to") for k in ("chain.id", "residue.index", "atom")
    ]
    if c.unique(keys).height != c.height:
        raise ValueError("Duplicate physical atom pairs")
    grouped = c.group_by("pdb.id", "chain.type.from").agg(
        pl.len().alias("pairs"),
        pl.struct("chain.id.from", "residue.index.from", "atom.from").n_unique().alias("atoms"),
    )
    maxima = grouped.group_by("pdb.id").agg(
        pl.col("pairs").max().alias("max_pairs"), pl.col("atoms").max().alias("max_atoms")
    )
    peptide = c.group_by("pdb.id").agg(
        pl.struct("chain.id.to", "residue.index.to", "atom.to").n_unique().alias("peptide_atoms")
    )
    result = (
        ids.join(maxima, on="pdb.id", how="left", validate="1:1")
        .join(peptide, on="pdb.id", how="left", validate="1:1")
        .fill_null(0)
    )
    return result.with_columns(
        (pl.col("max_pairs") / pl.col("peptide_atoms").clip(1)).alias("pairs_per_peptide_atom"),
        (pl.col("max_atoms") / pl.col("peptide_atoms").clip(1)).alias("tcr_atoms_per_peptide_atom"),
    )


def average_precision_matrix(y, scores):
    """Exact threshold AP for many candidate columns, including tied scores."""
    y = np.asarray(y)
    if set(y) != {0, 1} or not np.isfinite(scores).all():
        raise ValueError("AP needs both classes and finite scores")
    order = np.argsort(-scores, axis=0, kind="stable")
    ranked = np.take_along_axis(scores, order, axis=0)
    hits = y[order].cumsum(axis=0)
    end = np.vstack([ranked[:-1] != ranked[1:], np.ones((1, scores.shape[1]), bool)])
    previous = np.maximum.accumulate(np.where(end, hits, 0), axis=0)
    previous = np.vstack([np.zeros((1, scores.shape[1])), previous[:-1]])
    return (
        np.sum(
            np.where(end, (hits - previous) * hits / np.arange(1, len(y) + 1)[:, None], 0), axis=0
        )
        / y.sum()
    )


def candidate_spec():
    return [
        (name, mode, float(a), float(b))
        for name in FEATURES
        for mode in ("shared", "separate")
        for a in WEIGHTS
        for b in ([a] if mode == "shared" else WEIGHTS)
    ]


def scores_for(models, specs):
    x = models["iptm"].to_numpy()[:, None]
    g = models["betti"].to_numpy()[:, None]
    features = models.select(FEATURES).to_numpy()
    return (
        x
        + g * np.array([s[2] for s in specs])
        + features[:, [FEATURES.index(s[0]) for s in specs]] * np.array([s[3] for s in specs])
    )


def select(models, scores, specs, train):
    ap = np.mean(
        [
            average_precision_matrix(models["y"].to_numpy()[take], scores[take])
            for e in sorted(models["cohort"].unique())
            for take in [train[models["cohort"].to_numpy()[train] == e]]
        ],
        axis=0,
    )
    groups = {"best": np.arange(len(specs))}
    for feature in FEATURES:
        for mode in ("shared", "separate"):
            groups[feature + ":" + mode] = np.array(
                [j for j, s in enumerate(specs) if s[:2] == (feature, mode)]
            )
    # Deterministic ties prefer smaller contact and topology coefficients.
    chosen = {
        name: min(indices, key=lambda j: (-ap[j], specs[j][3], specs[j][2], j))
        for name, indices in groups.items()
    }
    return chosen, ap


def training_percentile(train, query):
    """Midrank transform fitted without held-out labels, to compare fold scales."""
    train = np.sort(train)
    return (np.searchsorted(train, query, "left") + np.searchsorted(train, query, "right")) / (
        2 * len(train)
    )


def fold_job(args):
    models, specs, fold, train, test = args
    scores = scores_for(models, specs)
    chosen, ap = select(models, scores, specs, train)
    epi = models["cohort"].to_numpy()
    rows, selected = [], []
    arms = {"iptm": models["iptm"].to_numpy(), "prepared": models["prepared"].to_numpy()}
    arms.update({arm: scores[:, j] for arm, j in chosen.items()})
    for arm, j in chosen.items():
        feature, mode, a, b = specs[j]
        selected.append(
            dict(
                fold=fold,
                arm=arm,
                feature=feature,
                mode=mode,
                a=a,
                b=b,
                training_macro_AP=float(ap[j]),
            )
        )
    for arm, values in arms.items():
        for e in sorted(set(epi)):
            tr, te = train[epi[train] == e], test[epi[test] == e]
            rows.append(
                models[te]
                .select("id", "cohort", "y")
                .with_columns(
                    pl.lit(fold).alias("fold"),
                    pl.lit(arm).alias("arm"),
                    pl.Series("raw_score", values[te]),
                    pl.Series("score", training_percentile(values[tr], values[te])),
                )
            )
    return pl.concat(rows), selected


def metrics(rows):
    result = []
    for (e, arm), part in rows.group_by("cohort", "arm"):
        y, x = part["y"].to_numpy(), part["score"].to_numpy()
        result.append(
            dict(cohort=e, arm=arm, AP=float(average_precision_score(y, x)), **separation(y, x))
        )
    return pl.DataFrame(result).sort("cohort", "arm")


def command(
    contacts: Path = typer.Option(...),
    atoms: Path = typer.Option(...),
    topology: Path = typer.Option(...),
    baseline: Path = typer.Option(...),
    output: Path = typer.Option(..., "-o"),
    threads: int = typer.Option(0, "-t"),
):
    """Evaluate raw and peptide-normalized maxima with shared/separate weights."""
    output.mkdir(parents=True, exist_ok=False)
    budget = resolve_threads(threads)
    b = pl.read_parquet(baseline).filter(pl.col("dataset") == "tcrvdb").sort("id")
    counts = contact_counts(
        pl.read_parquet(contacts), b.select(pl.col("id").alias("pdb.id"))
    ).rename({"pdb.id": "id"})
    g = (
        pl.read_csv(topology, separator="\t", schema_overrides={"pdb.id": pl.String})
        .filter(pl.col("graph") == "tp_pm_tm")
        .select(
            pl.col("pdb.id").alias("id"),
            pl.when(pl.col("available")).then("cyclo_frac").otherwise(None).alias("betti"),
        )
    )
    inv = pl.read_parquet(atoms)
    cdr = (
        inv.filter(pl.col("chain.type").is_in(["TRA", "TRB"]) & (pl.col("region.type") == "CDR3"))
        .select("pdb.id", "chain.type", "residue.index", "residue.aa")
        .unique()
        .sort("pdb.id", "chain.type", "residue.index")
        .group_by("pdb.id", "chain.type", maintain_order=True)
        .agg(pl.col("residue.aa").str.join("").alias("seq"))
        .pivot(on="chain.type", index="pdb.id", values="seq")
        .rename({"pdb.id": "id"})
    )
    data = (
        b.join(counts, on="id", validate="1:1")
        .join(g, on="id", validate="1:1")
        .join(cdr, on="id", validate="1:1")
    )
    if data.height != b.height:
        raise ValueError("Missing models in fresh exports")
    data = data.with_columns(
        pl.col("epitope").alias("cohort"),
        pl.Series("group", connected_groups(data["TRA"], data["TRB"])),
    )
    data.write_parquet(output / "all_models.parquet")
    # Keep every called row with the tested variables, not an unrelated historical mask.
    models = data.filter(
        pl.col("y").is_in([0, 1])
        & pl.all_horizontal(pl.col("iptm", "prepared", "betti").is_finite())
    ).sort("id")
    if not np.isfinite(models.select("iptm", "prepared", "betti", *FEATURES).to_numpy()).all():
        raise ValueError("Missing active scores in the established matched population")
    specs = candidate_spec()
    folds = splits(models, 5, SEED)
    results = parallel_map(
        fold_job, [(models, specs, i, a, b) for i, (a, b) in enumerate(folds)], budget
    )
    oof = pl.concat([r[0] for r in results])
    if oof.unique(["id", "arm"]).height != oof.height:
        raise ValueError("Duplicate held-out predictions")
    oof.write_parquet(output / "oof_predictions.parquet")
    metrics(oof).write_csv(output / "oof_metrics.tsv", separator="\t")
    pl.DataFrame([s for r in results for s in r[1]]).write_csv(
        output / "fold_selection.tsv", separator="\t"
    )
    scores = scores_for(models, specs)
    chosen, ap = select(models, scores, specs, np.arange(models.height))
    pl.DataFrame(
        [
            dict(feature=s[0], mode=s[1], a=s[2], b=s[3], macro_AP=float(v))
            for s, v in zip(specs, ap)
        ]
    ).write_csv(output / "full_search.tsv", separator="\t")
    full = []
    for arm, values in [
        ("iptm", models["iptm"].to_numpy()),
        ("prepared", models["prepared"].to_numpy()),
    ] + [(arm, scores[:, j]) for arm, j in chosen.items()]:
        full.append(
            models.select("id", "cohort", "y").with_columns(
                pl.lit(arm).alias("arm"), pl.Series("score", values)
            )
        )
    full = pl.concat(full)
    full.write_parquet(output / "full_fit_predictions.parquet")
    metrics(full).write_csv(output / "full_fit_metrics.tsv", separator="\t")
    best = specs[chosen["best"]]
    (output / "SELECTION.json").write_text(
        json.dumps(
            dict(
                feature=best[0],
                mode=best[1],
                a=best[2],
                b=best[3],
                formula="iptm + a * recognition_beta1_over_E + b * feature",
                selection="five-epitope macro AP; p-values not selected",
                inputs={str(p): digest(p) for p in [contacts, atoms, topology, baseline]},
                models=models.height,
                all_models=data.height,
                called=int(data["y"].is_in([0, 1]).sum()),
                candidates=len(specs),
                seed=SEED,
                contact_scope="CDR1/2/3 heavy atoms at 5 Angstrom; union over TRA/TRB for peptide atoms",
                oof="Exact-chain grouped; within-epitope training midrank transform for all arms",
                caveat="Exploratory after prior panel inspection; full-fit p-values are post-selection descriptive",
            ),
            indent=2,
        )
        + "\n"
    )
    print(json.dumps(dict(best=best, models=models.height, candidates=len(specs))), flush=True)


if __name__ == "__main__":
    typer.run(command)
