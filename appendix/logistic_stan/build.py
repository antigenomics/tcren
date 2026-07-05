#!/usr/bin/env python3
"""Fit + freeze the distribution-aware Bayesian logistic recognizer (PyMC); emit appendix figures.

Real (Oriented2026, label 1) vs Shuffled2026 10x decoys (label 0), full fresh tcren geometry + TCRen
energetics + MHC-class label. Each feature is entered by its natural-family canonical form
(tcren.recognition.encode_features: dock_torsion -> cos/sin, chain_balance -> logit, counts/continuous linear),
standardised, then a Bayesian logistic ``y ~ Bernoulli(sigmoid(a + Z b))`` is fit with PyMC (weakly-informative
Normal priors; a regularized-horseshoe variant is available via --prior horseshoe). 5-fold CV gives the honest
ROC/PR; the full-data posterior mean is frozen into src/tcren/data/shuffle_logistic.json.gz.

Run (in the PyMC venv):  results/eda/.venv/bin/python appendix/logistic_stan/build.py
2026-07-06
"""
from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pymc as pm
import arviz as az
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import (average_precision_score, balanced_accuracy_score, f1_score,
                             matthews_corrcoef, precision_recall_curve, roc_auc_score, roc_curve)

HERE = Path(__file__).resolve().parent
FIG = HERE / "figures"; FIG.mkdir(exist_ok=True)
DATA = HERE.parent / "shuffle_bn" / "data"          # reuse the same feature CSVs
META = HERE.parents[1] / "data" / "orient_metadata.json"
MODEL_OUT = HERE.parents[1] / "src" / "tcren" / "data" / "shuffle_logistic.json.gz"
DROP = {"key", "species", "mhc_class", "y", "mc"}

# load tcren.recognition standalone (numpy/json/gzip/math only — no tcren install needed in this venv)
_spec = importlib.util.spec_from_file_location(
    "tcren_recognition", HERE.parents[1] / "src" / "tcren" / "recognition.py")
_rec = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(_rec)
encode_features = _rec.encode_features
BayesianLogisticRecognizer = _rec.BayesianLogisticRecognizer


def load():
    meta = {m["pdb.id"]: m for m in json.load(open(META))}
    real = pd.read_csv(DATA / "real_features.csv"); real["y"] = 1
    real["mhc_class_bin"] = (real["mhc_class"] == "MHCII").astype(float)
    shuf = pd.read_csv(DATA / "shuffled_features.csv"); shuf["y"] = 0
    shuf["mhc_class_bin"] = shuf["key"].map(
        lambda k: 1.0 if meta.get(k.split("__")[0], {}).get("mhc.class") == "MHCII" else 0.0)
    feat = [c for c in real.columns if c in shuf.columns and c not in DROP]
    df = pd.concat([real[feat + ["y"]], shuf[feat + ["y"]]], ignore_index=True)
    return df, feat


def fit_pymc(Z, y, prior="normal", draws=1000, tune=1000, seed=0):
    p = Z.shape[1]
    with pm.Model():
        a = pm.Normal("a", 0.0, 2.5)
        if prior == "horseshoe":
            tau = pm.HalfCauchy("tau", 1.0)
            lam = pm.HalfCauchy("lam", 1.0, shape=p)
            b = pm.Normal("b", 0.0, tau * lam, shape=p)
        else:
            b = pm.Normal("b", 0.0, 1.0, shape=p)
        pm.Bernoulli("obs", logit_p=a + pm.math.dot(Z, b), observed=y)
        idata = pm.sample(draws, tune=tune, chains=2, cores=1, target_accept=0.9,
                          random_seed=seed, progressbar=False)
    return idata


def coefs(idata):
    return (float(idata.posterior["a"].mean()),
            idata.posterior["b"].mean(("chain", "draw")).values)


def standardize(Ztr):
    mean = np.nanmean(Ztr, 0); sd = np.nanstd(Ztr, 0) + 1e-9
    return mean, sd


def prep(Z, mean, sd):
    Z = np.where(np.isfinite(Z), Z, mean[None, :])
    return (Z - mean) / sd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prior", choices=["normal", "horseshoe"], default="normal")
    ap.add_argument("--cv-draws", type=int, default=600)
    args = ap.parse_args()

    df, feat = load()
    Zraw, enc = encode_features(df[feat].to_numpy(float), feat)
    y = df["y"].to_numpy()
    print(f"real={int(y.sum())} shuffled={int((1-y).sum())}  raw-feats={len(feat)}  encoded={len(enc)} "
          f"(prior={args.prior})")

    # 5-fold CV (standardise on train, refit PyMC per fold)
    oof = np.full(len(y), np.nan)
    for k, (tr, te) in enumerate(StratifiedKFold(5, shuffle=True, random_state=0).split(Zraw, y)):
        mean, sd = standardize(Zraw[tr])
        idata = fit_pymc(prep(Zraw[tr], mean, sd), y[tr], args.prior, draws=args.cv_draws, tune=args.cv_draws)
        a, b = coefs(idata)
        oof[te] = 1.0 / (1.0 + np.exp(-(a + prep(Zraw[te], mean, sd) @ b)))
        print(f"  fold {k+1}/5 done")

    auc = roc_auc_score(y, oof); ap_ = average_precision_score(y, oof); base = y.mean()
    fpr, tpr, thr = roc_curve(y, oof); t = thr[int(np.argmax(tpr - fpr))]
    yhat = (oof >= t).astype(int)
    bacc = balanced_accuracy_score(y, yhat); f1 = f1_score(y, yhat); mcc = matthews_corrcoef(y, yhat)
    prec, rec, _ = precision_recall_curve(y, oof)
    np.savetxt(FIG / "roc.dat", np.c_[fpr, tpr], header="fpr tpr", comments="# ")
    np.savetxt(FIG / "pr.dat", np.c_[rec, prec], header="recall precision", comments="# ")
    (FIG / "base.txt").write_text(f"{base:.4f}")
    print(f"5-fold CV: ROC-AUC={auc:.3f}  PR-AUC={ap_:.3f} (base {base:.3f})  bal-acc={bacc:.3f}  "
          f"F1={f1:.3f}  MCC={mcc:.3f}   (BN 0.865 / plain-logit 0.870)")
    (HERE / "metrics.tex").write_text(
        "\\begin{tabular}{lr}\n\\toprule\nmetric & value \\\\\n\\midrule\n"
        f"ROC--AUC & {auc:.3f} \\\\\nPR--AUC (baseline {base:.3f}) & {ap_:.3f} \\\\\n"
        f"balanced accuracy & {bacc:.3f} \\\\\n$F_1$ & {f1:.3f} \\\\\nMCC & {mcc:.3f} \\\\\n"
        f"\\bottomrule\n\\end{{tabular}}\n")

    # full-data fit -> freeze recognizer + posterior forest
    mean, sd = standardize(Zraw)
    idata = fit_pymc(prep(Zraw, mean, sd), y, args.prior, draws=1000, tune=1000)
    ndiv = int(idata.sample_stats["diverging"].sum())
    try:
        rhat = float(az.rhat(idata, var_names=["a", "b"]).to_array().max())
    except Exception:
        rhat = float("nan")
    print(f"convergence: max R-hat={rhat:.3f}  divergences={ndiv}")
    a, b = coefs(idata)
    MODEL_OUT.parent.mkdir(parents=True, exist_ok=True)
    BayesianLogisticRecognizer(feat, enc, mean, sd, a, b, prior=args.prior).save(MODEL_OUT)
    print(f"saved recognizer -> {MODEL_OUT}")

    bpost = idata.posterior["b"]                         # (chain, draw, p)
    bmean = bpost.mean(("chain", "draw")).values
    blo = bpost.quantile(0.03, dim=("chain", "draw")).values
    bhi = bpost.quantile(0.97, dim=("chain", "draw")).values
    order = np.argsort(bmean)
    with open(FIG / "forest.dat", "w") as fh:
        fh.write("# idx mean lo hi name\n")
        for rank, j in enumerate(order):
            fh.write(f"{rank} {bmean[j]:.4f} {blo[j]:.4f} {bhi[j]:.4f} {enc[j]}\n")
    (HERE / "forest_labels.tex").write_text(
        "\\def\\forestlabels{" + "; ".join(
            f"{r}: \\texttt{{{enc[j].replace('_', chr(92)+'_')}}}" for r, j in enumerate(order)) + "}\n")
    print("wrote forest.dat + roc/pr.dat + metrics.tex")


if __name__ == "__main__":
    main()
