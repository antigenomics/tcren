#!/usr/bin/env python3
"""Train + evaluate the real-vs-shuffled Gaussian BN classifier; emit figures data + store the model.

Reads data/real_features.csv (Oriented2026, label 1) and data/shuffled_features.csv (Shuffled2026 10x decoys,
label 0), both with the full fresh tcren geometry + TCRen energetics feature set + the MHC-class label. Runs a
5-fold stratified CV of tcren.recognition.GaussianBNClassifier, writes ROC/PR/marginal .dat for gnuplot and a
balanced-metrics LaTeX table, renders the learned BN with graphviz, and saves the final model (fit on all data)
to src/tcren/data/shuffle_bn.json.gz.

Run:  python appendix/shuffle_bn/build.py     (from the repo root, with tcren installed)
2026-07-06
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import (average_precision_score, balanced_accuracy_score, f1_score,
                             matthews_corrcoef, precision_recall_curve, roc_auc_score, roc_curve)

from tcren.recognition import GaussianBNClassifier

HERE = Path(__file__).resolve().parent
FIG = HERE / "figures"
FIG.mkdir(exist_ok=True)
# Repo-relative: this was an absolute path into `tcren-ms`, the repo's name before it was renamed
# to `tcren`, so it had stopped resolving on every machine including the author's.
META = HERE.parents[1] / "data" / "orient_metadata.json"
MODEL_OUT = HERE.parents[1] / "src" / "tcren" / "data" / "shuffle_bn.json.gz"
DROP = {"key", "species", "mhc_class", "y", "mc"}


def load():
    meta = {m["pdb.id"]: m for m in json.load(open(META))}
    real = pd.read_csv(HERE / "data" / "real_features.csv"); real["y"] = 1
    shuf = pd.read_csv(HERE / "data" / "shuffled_features.csv"); shuf["y"] = 0
    real["mc"] = (real["mhc_class"] == "MHCII").astype(int)
    shuf["mc"] = shuf["key"].map(lambda k: 1 if meta.get(k.split("__")[0], {}).get("mhc.class") == "MHCII" else 0)
    feat = [c for c in real.columns if c in shuf.columns and c not in DROP]
    df = pd.concat([real[feat + ["y", "mc"]], shuf[feat + ["y", "mc"]]], ignore_index=True)
    return df, feat


def main():
    df, feat = load()
    X = df[feat].to_numpy(float); y = df["y"].to_numpy(); m = df["mc"].to_numpy()
    print(f"real={int(y.sum())}  shuffled={int((1-y).sum())}  features={len(feat)}")

    # 5-fold stratified CV
    oof = np.full(len(y), np.nan)
    for tr, te in StratifiedKFold(5, shuffle=True, random_state=0).split(X, y):
        clf = GaussianBNClassifier(feat, max_parents=3).fit(X[tr], y[tr], m[tr])
        oof[te] = clf.predict_proba(X[te], m[te])[:, 1]

    auc = roc_auc_score(y, oof); ap = average_precision_score(y, oof); base = y.mean()
    fpr, tpr, thr = roc_curve(y, oof)
    jopt = int(np.argmax(tpr - fpr)); t = thr[jopt]
    yhat = (oof >= t).astype(int)
    bacc = balanced_accuracy_score(y, yhat); f1 = f1_score(y, yhat); mcc = matthews_corrcoef(y, yhat)
    prec, rec, _ = precision_recall_curve(y, oof)

    np.savetxt(FIG / "roc.dat", np.c_[fpr, tpr], header="fpr tpr", comments="# ")
    np.savetxt(FIG / "pr.dat", np.c_[rec, prec], header="recall precision", comments="# ")
    (FIG / "base.txt").write_text(f"{base:.4f}")
    print(f"ROC-AUC={auc:.3f}  PR-AUC={ap:.3f} (base {base:.3f})  bal-acc={bacc:.3f}  F1={f1:.3f}  MCC={mcc:.3f}")

    # balanced-metrics LaTeX table
    (HERE / "metrics.tex").write_text(
        "\\begin{tabular}{lr}\n\\toprule\nmetric & value \\\\\n\\midrule\n"
        f"ROC--AUC & {auc:.3f} \\\\\nPR--AUC (baseline {base:.3f}) & {ap:.3f} \\\\\n"
        f"balanced accuracy & {bacc:.3f} \\\\\n$F_1$ & {f1:.3f} \\\\\nMCC & {mcc:.3f} \\\\\n"
        f"\\bottomrule\n\\end{{tabular}}\n")

    # final model on all data -> store in tcren
    clf = GaussianBNClassifier(feat, max_parents=3).fit(X, y, m)
    MODEL_OUT.parent.mkdir(parents=True, exist_ok=True)
    clf.save(MODEL_OUT)
    print(f"saved model -> {MODEL_OUT}")
    (FIG / "bn.dot").write_text(clf.to_dot(coef_threshold=0.15))
    try:
        subprocess.run(["dot", "-Tpdf", str(FIG / "bn.dot"), "-o", str(FIG / "bn_dag.pdf")], check=True)
    except Exception as e:
        print("dot render failed:", e)

    # marginals: 6 most class-separating features, real vs shuffled normalized histograms
    dsep = sorted(feat, key=lambda c: -abs(
        (df.loc[df.y == 0, c].mean() - df.loc[df.y == 1, c].mean()) /
        (df[c].std() + 1e-9)))[:6]
    with open(FIG / "marginals.dat", "w") as fh:
        for c in dsep:
            a = df.loc[df.y == 1, c].dropna().to_numpy(); b = df.loc[df.y == 0, c].dropna().to_numpy()
            lo, hi = np.percentile(np.r_[a, b], [1, 99])
            bins = np.linspace(lo, hi, 25)
            ha, _ = np.histogram(a, bins, density=True); hb, _ = np.histogram(b, bins, density=True)
            ctr = 0.5 * (bins[:-1] + bins[1:])
            fh.write(f'# {c}\n')
            for x, va, vb in zip(ctr, ha, hb):
                fh.write(f"{x:.4g} {va:.4g} {vb:.4g}\n")
            fh.write("\n\n")                                     # gnuplot index separator
    (HERE / "marginal_labels.tex").write_text(
        "\\def\\marginallabels{" + ", ".join(f"({i+1})~\\texttt{{{c.replace('_','\\_')}}}" for i, c in enumerate(dsep)) + "}\n")
    print("marginals:", dsep)


if __name__ == "__main__":
    main()
