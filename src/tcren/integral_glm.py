"""Small frozen logistic readout trained on an explicitly separate real/mock pool.

The probability scale describes the balanced training task, not binding affinity
or calibrated experimental binding probability. Evaluation panels never fit it.
"""
from __future__ import annotations

import json
from pathlib import Path
import warnings

import numpy as np
import polars as pl
from scipy.special import expit
import typer

from .integral import _potential_hash, _read, _sha, _validated

FORMAT = "integral-glm-1"
PAIRWISE_FORMAT = "integral-pairwise-glm-1"
INPUTS = ("iptm", "energy", "geometry", "topology")


def _matrix(table):
    values, cls = _validated(table)
    if (cls != 0).any():
        raise ValueError("integral-glm-1 supports only MHC class I")
    return np.column_stack([values[k] for k in INPUTS])


def fit_integral_glm(table, labels, *, excluded_epitopes, excluded_receptors, provenance, pairwise=False):
    """Fit nine or fewer coefficients; labels require explicit epitope/receptor IDs.

    Caller verifies current positive assignments, mock construction, coordinate
    identity and disjoint evaluation provenance. Exclusions are mandatory inputs.
    Each epitope/class cell has equal total weight, normalized to sample count.
    """
    from sklearn import __version__ as sklearn_version
    from sklearn.exceptions import ConvergenceWarning
    from sklearn.linear_model import LogisticRegression
    from threadpoolctl import threadpool_limits

    x = _matrix(table)
    columns = ["complex.id", "epitope", "receptor_id", "y"]
    lab = labels.select(columns)
    if lab.select(pl.any_horizontal(pl.all().is_null()).any()).item():
        raise ValueError("training labels and identities must be nonmissing")
    for col in columns[:3]:
        if lab[col].dtype != pl.String or (lab[col] == "").any():
            raise ValueError(f"nonempty string {col} required")
    if lab["complex.id"].is_duplicated().any() or set(lab["complex.id"]) != set(table["complex.id"]):
        raise ValueError("training feature/label IDs must match one to one")
    if not excluded_epitopes or not excluded_receptors:
        raise ValueError("explicit nonempty evaluation exclusions required")
    if set(lab["epitope"]) & set(excluded_epitopes) or set(lab["receptor_id"]) & set(excluded_receptors):
        raise ValueError("evaluation epitope or receptor occurs in training")
    if not set(lab["y"]).issubset({0, 1}):
        raise ValueError("binary training labels required")
    lab = table.select("complex.id").join(lab, on="complex.id", validate="1:1", maintain_order="left")
    cells = lab.group_by("epitope", "y").len()
    if cells.group_by("epitope").len().filter(pl.col("len") != 2).height or not lab.height:
        raise ValueError("every training epitope requires both real and mock models")
    if not np.isfinite(x).any(axis=1).all() or not np.isfinite(x).any(axis=0).all():
        raise ValueError("training requires an observation in every row and every component")
    digest = provenance.get("input_sha256", "")
    if not isinstance(digest, str) or len(digest) != 64 or any(c not in "0123456789abcdef" for c in digest):
        raise ValueError("provenance requires input_sha256")
    counts = lab.join(cells, on=["epitope", "y"], validate="m:1", maintain_order="left")["len"].to_numpy()
    weights = len(lab) / (len(cells) * counts)
    if pairwise:
        sizes = lab.group_by("epitope").len()
        n_epitope = lab.join(sizes, on="epitope", validate="m:1", maintain_order="left")["len"].to_numpy()
        weights = len(lab) * np.sqrt(n_epitope) / (2 * counts * np.sqrt(sizes["len"].to_numpy()).sum())
    missing = ~np.isfinite(x)
    means = np.array([np.average(x[~missing[:, j], j], weights=weights[~missing[:, j]]) for j in range(4)])
    scales = np.array([np.sqrt(np.average((x[~missing[:, j], j] - means[j])**2,
        weights=weights[~missing[:, j]])) for j in range(4)])
    scales[scales == 0] = 1
    design = np.column_stack([np.where(missing, 0, (x - means) / scales), missing.astype(float)])
    target = lab["y"].to_numpy()
    n_pairs = None
    if pairwise:
        design, target, weights, n_pairs = _pairwise_design(design, target, lab["epitope"].to_numpy())
    model = LogisticRegression(C=1., solver="lbfgs", tol=1e-8, max_iter=1000,
                               random_state=0, fit_intercept=not pairwise)
    with warnings.catch_warnings(), threadpool_limits(limits=1):
        warnings.simplefilter("error", ConvergenceWarning)
        model.fit(design, target, sample_weight=weights)
    return dict(format=PAIRWISE_FORMAT if pairwise else FORMAT, mhc_class="I", potential_sha256=_potential_hash(),
        inputs=list(INPUTS), means=means.tolist(), scales=scales.tolist(),
        coefficients=model.coef_[0].tolist(), intercept=float(model.intercept_[0]),
        n_training=len(lab), n_epitopes=lab["epitope"].n_unique(), n_pairs=n_pairs,
        n_missing=missing.sum(axis=0).tolist(), C=1.,
        sklearn_version=sklearn_version, numpy_version=np.__version__,
        objective=("within-epitope pairwise logistic, square-root model-count epitope weights, sum weights equals n_training"
                   if pairwise else "equal epitope/class cell weights, sum weights equals n_training"),
        provenance=dict(provenance))


def _pairwise_design(design, target, epitopes):
    groups = [(np.flatnonzero((epitopes == e) & (target == 1)),
               np.flatnonzero((epitopes == e) & (target == 0))) for e in np.unique(epitopes)]
    n_pairs = sum(len(pos) * len(neg) for pos, neg in groups)
    # ponytail: materialize at most one million pairs; stream the objective for larger cohorts.
    if n_pairs > 1_000_000:
        raise ValueError("pairwise GLM exceeds one million pairs; no subsampling is performed")
    paired = np.empty((2 * n_pairs, design.shape[1]))
    labels, weights = np.empty(2 * n_pairs, dtype=int), np.empty(2 * n_pairs)
    support = np.sqrt([len(pos) + len(neg) for pos, neg in groups])
    start = 0
    for group_index, (pos, neg) in enumerate(groups):
        delta = (design[pos, None, :] - design[neg][None, :, :]).reshape(-1, design.shape[1])
        n = len(delta)
        paired[start:start+n], paired[start+n:start+2*n] = delta, -delta
        labels[start:start+n], labels[start+n:start+2*n] = 1, 0
        weights[start:start+2*n] = len(design) * support[group_index] / (2 * n * support.sum())
        start += 2 * n
    return paired, labels, weights, n_pairs


def integral_glm_scores(table, model):
    """Apply frozen training transformations and expose additive logit contributions."""
    if model.get("format") not in {FORMAT, PAIRWISE_FORMAT} or model.get("mhc_class") != "I" or model.get("inputs") != list(INPUTS):
        raise ValueError("unsupported integral GLM format, class or inputs")
    if model.get("potential_sha256") != _potential_hash():
        raise ValueError("GLM potential differs from bundled TCRen2")
    means = np.asarray(model.get("means", []), dtype=float)
    scales = np.asarray(model.get("scales", []), dtype=float)
    beta = np.asarray(model.get("coefficients", []), dtype=float)
    intercept = np.asarray(model.get("intercept"), dtype=float)
    if (means.shape != (4,) or scales.shape != (4,) or beta.shape != (8,) or intercept.shape != ()
            or not all(np.isfinite(a).all() for a in (means, scales, beta, intercept)) or (scales <= 0).any()):
        raise ValueError("invalid finite GLM parameters")
    if model["format"] == PAIRWISE_FORMAT and intercept != 0:
        raise ValueError("pairwise GLM has no intercept")
    x = _matrix(table)
    missing = ~np.isfinite(x)
    contribution = np.where(missing, 0, (x - means) / scales) * beta[:4]
    missing_contribution = missing @ beta[4:]
    logit = intercept + contribution.sum(axis=1) + missing_contribution
    unavailable = pl.Series(missing.all(axis=1))
    return pl.DataFrame({"complex.id": table["complex.id"],
        **{"glm_" + key + "_contribution": contribution[:, j] for j, key in enumerate(INPUTS)},
        **{"missing_" + key: missing[:, j] for j, key in enumerate(INPUTS)},
        "glm_missing_contribution": missing_contribution,
        "n_components": (~missing).sum(axis=1),
        "integral_glm_logit": pl.Series(logit).set(unavailable, None),
        "integral_glm_score": pl.Series(expit(logit)).set(unavailable, None)})


def fit_command(
    features: Path = typer.Option(..., "--features", exists=True, dir_okay=False),
    labels: Path = typer.Option(..., "--labels", exists=True, dir_okay=False),
    excluded_epitopes: Path = typer.Option(..., "--excluded-epitopes", exists=True, dir_okay=False),
    excluded_receptors: Path = typer.Option(..., "--excluded-receptors", exists=True, dir_okay=False),
    output: Path = typer.Option(..., "-o", "--output"),
    pairwise: bool = typer.Option(False, "--pairwise", help="fit within-epitope ordering without an intercept"),
):
    """Freeze a separate-training GLM; exclusion files contain epitope/receptor_id columns."""
    if output.exists():
        raise FileExistsError(output)
    provenance = dict(input_sha256=_sha(features), labels_sha256=_sha(labels),
        excluded_epitopes_sha256=_sha(excluded_epitopes), excluded_receptors_sha256=_sha(excluded_receptors))
    model = fit_integral_glm(_read(features), _read(labels),
        excluded_epitopes=set(_read(excluded_epitopes)["epitope"]),
        excluded_receptors=set(_read(excluded_receptors)["receptor_id"]), provenance=provenance, pairwise=pairwise)
    with output.open("x") as stream:
        stream.write(json.dumps(model, indent=2, allow_nan=False) + "\n")
    typer.echo(f"Frozen GLM from {model['n_training']} models in {model['n_epitopes']} training epitopes")


def predict_command(
    features: Path = typer.Option(..., "--features", exists=True, dir_okay=False),
    model: Path = typer.Option(..., "--model", exists=True, dir_okay=False),
    output: Path = typer.Option(..., "-o", "--output"),
):
    """Apply a frozen GLM without fitting, cohort normalization or new thresholds."""
    from .provenance import stamp
    resource = json.loads(model.read_text())
    result = integral_glm_scores(_read(features), resource)
    result.write_csv(output, separator="\t")
    stamp(output, command="integral-predict", extra=dict(input_sha256=_sha(features),
        model_sha256=_sha(model), format=resource["format"], potential_sha256=_potential_hash()))
    typer.echo(f"Wrote {result.height} frozen GLM scores")
