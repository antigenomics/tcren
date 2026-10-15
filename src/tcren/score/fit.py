"""Fitting the frozen hold-out model, and the reason this file is shipped rather than kept aside.

``P_native`` was withdrawn from this project because its coefficients were frozen against a
training set that no longer existed, which made it the one part of the package a reader could not
reproduce. The repair is not to stop fitting -- a two-class covariance is what reads a
variance break, and a variance break is what separates binders on the hardest stratum -- but to
ship the fitter, the manifest and the frozen output together, so that

    tcren fit-holdout --features <table> --manifest <csv> -o holdout_model.npz

implements its historical recipe from named public inputs. The shipped numerical arrays remain
frozen through the 2026-10-14 definition migration; refits under corrected descriptor semantics
write a separate output and require matched comparison before promotion.

What is estimated, and on what: a mean and a covariance per class over the transformed descriptor
coordinates, on out-of-panel structures only, weighted 1/n_epitope so that one deeply sampled
epitope does not set the shape of the binder manifold. ipTM enters as one further coordinate of
the **binder** Gaussian alone, because every hold-out positive carries it and one whole negative
arm does not -- a two-class joint over it would learn a property of the deposit rather than of the
interface.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from .. import __version__
from ..provenance import registry_digest
from . import CONF_COORD, _logit
from .model import Joint, conditional_variance
from .transform import Transformer, working_set


def epitope_weights(epitopes) -> np.ndarray:
    """1/n per epitope, normalised: even coverage without discarding a structure.

    Subsampling to an even epitope count reaches the same first and second moments and throws rows
    away; the whole construction here is a covariance, so thinning it is the one thing not to do.
    """
    e = np.asarray(["?" if x is None else str(x) for x in epitopes])
    _, inv, cnt = np.unique(e, return_inverse=True, return_counts=True)
    w = 1.0 / cnt[inv]
    return w / w.sum()


def fit_holdout(features, manifest, *, out: Path | None = None) -> dict:
    """Fit the frozen model. ``features`` and ``manifest`` are polars frames keyed on ``pdb.id``.

    ``manifest`` needs the id column, ``y`` (1 binder / 0 non-binder) and ``epitope``; an ``iptm``
    column is used for the confidence coordinate when over 200 finite binder observations exist;
    otherwise the confidence residual is unavailable. Keys must be unique and non-null. Fitting
    a feature subset of the manifest is supported; unmatched-row counts are recorded in metadata.
    """
    import polars as pl

    cols = working_set(receptor_task=False)          # the wider set; the receptor read is a marginal
    key = next((c for c in ("complex.id", "pdb.id") if c in features.columns and c in manifest.columns), None)
    if key is None:
        raise KeyError("features and manifest must share a 'complex.id' or 'pdb.id' column")
    for name, frame in (("features", features), ("manifest", manifest)):
        if frame[key].null_count() or frame[key].n_unique() != frame.height:
            raise ValueError(f"{name} {key!r} keys must be unique and non-null")
    labels = manifest["y"]
    if (labels.null_count() or not (labels.dtype.is_numeric() or labels.dtype == pl.Boolean)
            or not np.isin(labels.to_numpy(), [0, 1]).all()):
        raise ValueError("manifest y must contain only binary labels 0 and 1")
    t = manifest.join(features, on=key, how="inner", validate="1:1")
    n_joined = t.height
    missing = [c for c in cols if c not in t.columns]
    if missing:
        raise KeyError(f"the feature table is missing {len(missing)} descriptors, "
                       f"first few: {missing[:5]}")
    X = t.select(cols).to_numpy().astype(float)
    ok = np.isfinite(X).all(1)
    t, X = t.filter(pl.Series(ok)), X[ok]
    if t.height < 200:
        raise ValueError(f"only {t.height} complete-case rows; this is not a hold-out.")

    tr = Transformer(names=cols).fit(X)
    Z = tr.transform(X)
    names = tr.out_names()
    y = t["y"].to_numpy().astype(int)

    Xc = {c: Z[y == c] for c in (0, 1)}
    w = {c: epitope_weights(t.filter(pl.col("y") == c)["epitope"]) for c in (0, 1)}
    if min(len(Xc[0]), len(Xc[1])) < 50:
        raise ValueError(f"class sizes {len(Xc[0])} / {len(Xc[1])}: too few to fit a covariance.")
    j = Joint(names=names).fit(Xc, w)

    # ipTM as one more coordinate of the binder Gaussian: its mean, its covariance with every
    # descriptor coordinate, and its own variance. Weighted with the same epitope weights.
    conf_mu, conf_cov, conf_var = 0.0, np.zeros(len(names)), 0.0
    confidence_available = False
    if "iptm" in t.columns:
        ip = _logit(t.filter(pl.col("y") == 1)["iptm"].to_numpy().astype(float))
        good = np.isfinite(ip)
        if good.sum() > 200:
            confidence_available = True
            ww = w[1][good]
            ww = ww / ww.sum()
            conf_mu = float(ww @ ip[good])
            Zc = Xc[1][good] - j.mu[1]
            conf_cov = ((ip[good] - conf_mu) * ww) @ Zc / (1.0 - (ww ** 2).sum())
            conf_var = float(((ip[good] - conf_mu) ** 2 * ww).sum() / (1.0 - (ww ** 2).sum()))

    conf_cond_var = (conditional_variance(j.cov[1], conf_cov, conf_var)
                     if confidence_available else None)
    meta = {
        "descriptors": list(cols),
        "coordinates": list(names),
        "receptor_coordinates": [n for n in names
                                 if (n[:-4] if n.endswith(("_cos", "_sin")) else n)
                                 in set(working_set(receptor_task=True))],
        "lam": tr.lam, "loc": tr.loc, "scale": tr.scale,
        "prior": [j.prior[0], j.prior[1]], "alpha": [j.alpha[0], j.alpha[1]],
        "conf_mu": conf_mu, "conf_var": conf_var, "conf_coord": CONF_COORD,
        "confidence_available": confidence_available,
        "conf_conditional_variance": conf_cond_var,
        "transform_training": "all joined complete rows, both classes",
        "binder_moments_training": "label-selected positive rows",
        "catalogue_digest": registry_digest(),
        "n_pos": int((y == 1).sum()), "n_neg": int((y == 0).sum()),
        "n_epitopes": int(t["epitope"].n_unique()),
        "n_unmatched_features": features.height - n_joined,
        "n_unmatched_manifest": manifest.height - n_joined,
        "n_incomplete": n_joined - t.height,
        "tcren_version": __version__,
    }
    if out is not None:
        out = Path(out)
        out.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(out, meta=json.dumps(meta), mu0=j.mu[0], mu1=j.mu[1],
                            cov0=j.cov[0], cov1=j.cov[1], conf_cov=conf_cov)
    return meta
