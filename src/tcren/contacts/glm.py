"""Requested in-sample, unpenalized binomial GLMs for contact totals.

The historical prepared score is named betti here. No folds, regularization,
epitope indicators, class weights or production default changes are introduced.
"""

from pathlib import Path
import json
import warnings

import numpy as np
import polars as pl
from scipy.special import expit
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
import typer

from .additive import contact_counts, metrics, oriented_contacts
from .profile import digest

COUNTS = ["peptide_atoms", "alpha_atoms", "beta_atoms"]


def atom_totals(contacts, ids):
    """Unique endpoint atoms and unique residue edges across TRA/TRB at 5 Å."""
    c = oriented_contacts(contacts)
    result = contact_counts(c, ids).select("pdb.id", "peptide_atoms")
    totals = c.group_by("pdb.id").agg(
        *[
            pl.struct("chain.id.from", "residue.index.from", "atom.from")
            .filter(pl.col("chain.type.from") == chain)
            .n_unique()
            .alias(name)
            for chain, name in [("TRA", "alpha_atoms"), ("TRB", "beta_atoms")]
        ],
        pl.struct("chain.id.from", "residue.index.from", "chain.id.to", "residue.index.to")
        .n_unique()
        .alias("tcr_peptide_edges"),
    )
    return result.join(totals, on="pdb.id", how="left", validate="1:1").fill_null(0)


def fit_glm(x, y):
    """Unpenalized maximum-likelihood logistic fit, coefficients in input units."""
    x, y = np.asarray(x, float), np.asarray(y)
    if not np.isfinite(x).all() or set(y) != {0, 1}:
        raise ValueError("GLM needs finite predictors and both classes")
    scaler = StandardScaler().fit(x)
    z = scaler.transform(x)
    with warnings.catch_warnings():
        warnings.simplefilter("error", ConvergenceWarning)
        model = LogisticRegression(C=np.inf, solver="lbfgs", tol=1e-12, max_iter=10000).fit(z, y)
    coef = model.coef_[0] / scaler.scale_
    intercept = float(model.intercept_[0] - scaler.mean_ @ coef)
    probability = expit(intercept + x @ coef)
    gradient = np.r_[np.mean(probability - y), z.T @ (probability - y) / len(y)]
    if np.max(np.abs(gradient)) > 1e-6:
        raise ValueError("Unpenalized GLM did not reach the likelihood optimum")
    np.testing.assert_allclose(probability, model.predict_proba(z)[:, 1], atol=1e-12)
    return probability, dict(
        intercept=intercept,
        coefficients=coef.tolist(),
        iterations=int(model.n_iter_[0]),
        max_abs_gradient=float(abs(gradient).max()),
    )


def command(
    contacts: Path = typer.Option(...),
    topology: Path = typer.Option(...),
    baseline: Path = typer.Option(...),
    output: Path = typer.Option(..., "-o"),
):
    """Fit raw-count and edge-normalized GLMs to all called rows, in-sample."""
    output.mkdir(parents=True, exist_ok=False)
    b = (
        pl.read_parquet(baseline)
        .filter(pl.col("dataset") == "tcrvdb")
        .rename({"prepared": "betti"})
    )
    counts = atom_totals(pl.read_parquet(contacts), b.select(pl.col("id").alias("pdb.id"))).rename(
        {"pdb.id": "id"}
    )
    g = (
        pl.read_csv(topology, separator="\t", schema_overrides={"pdb.id": pl.String})
        .filter(pl.col("graph") == "tp_pm_tm")
        .select(pl.col("pdb.id").alias("id"), "cyclo_frac", "available")
    )
    data = (
        b.join(counts, on="id", validate="1:1")
        .join(g, on="id", validate="1:1")
        .with_columns(pl.col("epitope").alias("cohort"))
    )
    if data.height != b.height:
        raise ValueError("Fresh export population differs from the input inventory")
    rebuilt = 0.55 * data["iptm"].to_numpy() + 0.45 * data["cyclo_frac"].to_numpy()
    np.testing.assert_allclose(rebuilt, data["betti"].to_numpy(), atol=1e-14, equal_nan=True)
    data = data.with_columns(pl.Series("betti", rebuilt))
    data.write_parquet(output / "all_models.parquet")
    models = data.filter(
        pl.col("y").is_in([0, 1]) & pl.all_horizontal(pl.col("iptm", "betti").is_finite())
    ).sort("id")
    y = models["y"].to_numpy()
    names = ["iptm", "betti", *COUNTS]
    parameters, predictions = {}, []
    for arm in ["raw_glm", "normalized_glm"]:
        x = models.select(names).to_numpy().astype(float)
        if arm == "normalized_glm":
            # A zero-contact model has zero atoms; define its three densities as zero.
            x[:, 2:] /= np.maximum(models["tcr_peptide_edges"].to_numpy(), 1)[:, None]
        probability, fit = fit_glm(x, y)
        parameters[arm] = dict(features=names, **fit)
        predictions.append(
            models.select("id", "cohort", "y").with_columns(
                pl.lit(arm).alias("arm"), pl.Series("score", probability)
            )
        )
    for arm in ["iptm", "betti"]:
        predictions.append(
            models.select("id", "cohort", "y").with_columns(
                pl.lit(arm).alias("arm"), pl.Series("score", models[arm])
            )
        )
    predictions = pl.concat(predictions)
    predictions.write_parquet(output / "predictions.parquet")
    metrics(predictions).write_csv(output / "metrics.tsv", separator="\t")
    models.write_parquet(output / "models.parquet")
    (output / "FIT.json").write_text(
        json.dumps(
            dict(
                models=models.height,
                positive=int(y.sum()),
                negative=int(len(y) - y.sum()),
                fits=parameters,
                evaluation="in-sample, intentionally overfit",
                fit_scope="one pooled classifier per representation; no epitope effects",
                betti="previously called prepared: 0.55 iptm + 0.45 recognition beta1/E",
                normalization="each atom count divided by total unique TCR-peptide residue edges; zero/zero assigned zero",
                inputs={str(p): digest(p) for p in [contacts, topology, baseline]},
            ),
            indent=2,
        )
        + "\n"
    )
    print(json.dumps(parameters), flush=True)


if __name__ == "__main__":
    typer.run(command)
