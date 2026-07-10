#!/usr/bin/env python
"""Can C-beta ride on the Frenet frame, and does Ramachandran add anything?

The Frenet descriptor (kappa, tau) sees only the C-alpha trace. It knows the loop's shape but
not which way the side chains point -- and which way they point is what decides whether a
residue faces the peptide. This script asks three questions of the crystal set:

  1. Is the idealised virtual C-beta good enough to stand in for glycine's missing one?
     Gate: RMSD against observed C-beta on non-Gly residues must be < 0.15 A.

  2. Does the side-chain direction in the local Frenet frame -- (polar, azimuth) of the
     Ca->Cb unit vector, both rigid-motion invariant -- predict peptide contact? Compared
     against the obvious baseline, position along the loop.

  3. Does Ramachandran (phi, psi) carry information the C-alpha (kappa, tau) does not?
     Gate: if a k-NN regressor recovers (phi, psi) from (kappa, tau) with circular R^2 >= 0.8,
     Ramachandran is redundant and does not get built into anything.

Residues are the unit of measurement but junctions are the unit of independence, so every
interval is bootstrapped over junctions.

Usage:
    python scripts/cb_contacts.py [--contact-cutoff 4.5] [--boot 1000]

2026-07-10
"""
from __future__ import annotations

import argparse
import glob
import gzip
import os
import random
import sys
import warnings

import numpy as np

warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
sys.path.insert(0, os.path.dirname(__file__))

from _harvest import _markup  # noqa: E402

from tcren.loops import cb_orientation, find_junctions, frenet, is_omega_loop, ramachandran, virtual_cb  # noqa: E402

AA = "ACDEFGHIKLMNPQRSTVWY"


def collect(pattern, markup_csv):
    """Per junction: backbone arrays, observed/virtual CB, and a peptide-contact label."""
    from Bio.PDB import PDBParser
    from Bio.PDB.Polypeptide import index_to_one, three_to_index

    mk = _markup(markup_csv)
    parser = PDBParser(QUIET=True)
    out = []
    for path in sorted(glob.glob(pattern)):
        pdb = os.path.basename(path).split(".")[0].lower()
        row = mk.get(pdb)
        if row is None or not row.get("peptide"):
            continue
        try:
            with gzip.open(path, "rt") as fh:
                model = parser.get_structure(pdb, fh)[0]
        except Exception:
            continue

        # the peptide chain: the short one whose sequence matches the curated epitope
        pep_atoms = None
        for chain in model:
            res = [r for r in chain if r.id[0] == " "]
            if not 7 <= len(res) <= 20:
                continue
            try:
                seq = "".join(index_to_one(three_to_index(r.get_resname())) for r in res)
            except Exception:
                continue
            if seq == row["peptide"]:
                pep_atoms = np.array([a.get_coord() for r in res for a in r
                                      if a.element != "H"], dtype=float)
                break
        if pep_atoms is None or len(pep_atoms) == 0:
            continue

        for chain in model:
            res = [r for r in chain if r.id[0] == " " and {"N", "CA", "C"} <= {a.get_id() for a in r}]
            if not 90 <= len(res) <= 130:
                continue
            try:
                seq = "".join(index_to_one(three_to_index(r.get_resname())) for r in res)
            except Exception:
                continue
            ca = np.array([r["CA"].get_coord() for r in res], dtype=float)
            for j in find_junctions(seq, ca):
                if not (set(j.seq) <= set(AA) and is_omega_loop(j.ca, relax_length=True)):
                    continue
                sub = res[j.cys:j.fw + 1]
                if j.cdr3 not in (row.get("cdr3a"), row.get("cdr3b")):
                    continue
                ct = "TRA" if j.cdr3 == row.get("cdr3a") else "TRB"
                nn = np.array([r["N"].get_coord() for r in sub], dtype=float)
                cc = np.array([r["C"].get_coord() for r in sub], dtype=float)
                obs_cb = np.array([r["CB"].get_coord() if "CB" in r else [np.nan] * 3
                                   for r in sub], dtype=float)
                # min heavy-atom distance from each junction residue to the peptide
                dmin = np.array([
                    min(np.linalg.norm(np.array(a.get_coord()) - pep_atoms, axis=1).min()
                        for a in r if a.element != "H")
                    for r in sub], dtype=float)
                out.append({"pdb": pdb, "chain_type": ct, "seq": j.seq, "ca": j.ca,
                            "n": nn, "c": cc, "obs_cb": obs_cb, "dmin": dmin})
    return out


def auc(labels, scores):
    """ROC-AUC with ties at 0.5 credit; no sklearn dependency."""
    order = sorted(range(len(scores)), key=lambda i: scores[i])
    pos = sum(labels)
    neg = len(labels) - pos
    if pos == 0 or neg == 0:
        return float("nan")
    rank, i = {}, 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and scores[order[j + 1]] == scores[order[i]]:
            j += 1
        r = (i + j) / 2 + 1
        for k in range(i, j + 1):
            rank[order[k]] = r
        i = j + 1
    s = sum(rank[i] for i in range(len(labels)) if labels[i])
    return (s - pos * (pos + 1) / 2) / (pos * neg)


def boot_auc(groups, boot, rng):
    """Bootstrap over junctions: resample whole loops, not residues."""
    vals = []
    for _ in range(boot):
        draw = [groups[rng.randrange(len(groups))] for _ in range(len(groups))]
        lab = [l for g in draw for l in g[0]]
        sc = [s for g in draw for s in g[1]]
        a = auc(lab, sc)
        if not np.isnan(a):
            vals.append(a)
    vals.sort()
    return vals[int(0.025 * len(vals))], vals[int(0.975 * len(vals))]


def circular_r2(true_deg, pred_deg):
    """1 - E[1 - cos(err)] / E[1 - cos(true - circular mean)]. 1.0 = perfect, 0 = no better than the mean."""
    err = np.radians(true_deg - pred_deg)
    resid = np.mean(1 - np.cos(err))
    mu = np.arctan2(np.mean(np.sin(np.radians(true_deg))), np.mean(np.cos(np.radians(true_deg))))
    total = np.mean(1 - np.cos(np.radians(true_deg) - mu))
    return 1 - resid / total if total > 0 else float("nan")


def knn_circular(x_train, y_train, x_test, k=15):
    """k-NN circular-mean regression; y in degrees."""
    pred = []
    yr = np.radians(y_train)
    for row in x_test:
        d = np.linalg.norm(x_train - row, axis=1)
        idx = np.argpartition(d, k)[:k]
        pred.append(np.degrees(np.arctan2(np.sin(yr[idx]).mean(), np.cos(yr[idx]).mean())))
    return np.array(pred)


def _residue_table(loops, cutoff):
    """One row per interior junction residue: features, contact label, junction id."""
    from tcren.loops import frenet as _fr

    rows, labels, gid = [], [], []
    for g, lp in enumerate(loops):
        if len(lp["ca"]) < 6:
            continue
        cb = np.where(np.isnan(lp["obs_cb"]), virtual_cb(lp["n"], lp["ca"], lp["c"]), lp["obs_cb"])
        polar, azim = cb_orientation(lp["ca"], cb)
        kap, tau = _fr(lp["ca"])
        phi, psi = ramachandran(lp["n"], lp["ca"], lp["c"])
        m = min(len(polar), len(kap), len(phi), len(tau))
        L = len(lp["ca"])
        for i in range(m):
            relpos = abs(2 * (i + 1) - (L - 1)) / (L - 1)
            v = [relpos, kap[i], abs(tau[i]), polar[i], abs(azim[i]),
                 np.cos(np.radians(azim[i])), np.cos(np.radians(phi[i])), np.sin(np.radians(phi[i]))]
            if any(np.isnan(v)):
                continue
            rows.append(v)
            labels.append(int(lp["dmin"][i + 1] < cutoff))
            gid.append(g)
    return np.array(rows), np.array(labels), np.array(gid)


def incremental(loops, args):
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import roc_auc_score
    from sklearn.model_selection import GroupKFold
    from sklearn.preprocessing import StandardScaler

    X, y, g = _residue_table(loops, args.contact_cutoff)
    cols = {"relpos": [0], "+ shape (kappa,|tau|)": [0, 1, 2],
            "+ Cbeta (polar,|azim|,cos azim)": [0, 1, 2, 3, 4, 5],
            "+ Ramachandran phi": [0, 1, 2, 3, 4, 5, 6, 7]}
    print(f"  residues {len(y):,}, junctions {len(set(g))}, contacting {100*y.mean():.1f}%\n")
    print(f"  {'feature set':<34}{'AUC':>8}{'gain':>8}")
    prev = None
    for name, idx in cols.items():
        aucs = []
        for tr, te in GroupKFold(n_splits=5).split(X, y, g):
            sc = StandardScaler().fit(X[tr][:, idx])
            clf = LogisticRegression(max_iter=2000).fit(sc.transform(X[tr][:, idx]), y[tr])
            aucs.append(roc_auc_score(y[te], clf.predict_proba(sc.transform(X[te][:, idx]))[:, 1]))
        a = float(np.mean(aucs))
        gain = "" if prev is None else f"{a-prev:+.3f}"
        print(f"  {name:<34}{a:>8.3f}{gain:>8}")
        prev = a


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--structures", default="data/Canonical2026/*.pdb.gz")
    ap.add_argument("--markup", default="notebooks/natcompsci2022/results_new/markup_2026.csv")
    ap.add_argument("--contact-cutoff", type=float, default=4.5)
    ap.add_argument("--boot", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    rng = random.Random(args.seed)

    loops = collect(args.structures, args.markup)
    print(f"{len(loops)} junctions with a resolved peptide chain "
          f"({sum(1 for l in loops if l['chain_type']=='TRA')} TRA, "
          f"{sum(1 for l in loops if l['chain_type']=='TRB')} TRB)")

    # ---------------------------------------------------------------- 1. virtual CB gate
    err, n_gly = [], 0
    for lp in loops:
        vcb = virtual_cb(lp["n"], lp["ca"], lp["c"])
        ok = ~np.isnan(lp["obs_cb"]).any(1)
        n_gly += (~ok).sum()
        err.extend(np.linalg.norm(vcb[ok] - lp["obs_cb"][ok], axis=1))
    rmsd = float(np.sqrt(np.mean(np.square(err))))
    print(f"\n=== 1. virtual C-beta gate ===")
    print(f"  non-Gly residues checked: {len(err):,}   glycines needing the virtual CB: {n_gly:,}")
    print(f"  RMSD(virtual, observed) = {rmsd:.3f} A   (gate: < 0.15 A)")
    if rmsd >= 0.15:
        print("  !! construction rejected; do not use it")
        return 1
    print("  passed -- glycine can be given a side-chain direction.")

    # ---------------------------------------------------------------- 2. contact prediction
    print(f"\n=== 2. does the Frenet-frame C-beta direction predict peptide contact? ===")
    print(f"  contact := min heavy-atom distance to the peptide < {args.contact_cutoff} A")
    feats: dict[str, list] = {"polar": [], "azimuth": [], "|azimuth|": [],
                              "kappa": [], "|tau|": [], "relpos": []}
    labels: list[int] = []
    groups: dict[str, list] = {k: [] for k in feats}
    for lp in loops:
        if len(lp["ca"]) < 5:
            continue
        cb = np.where(np.isnan(lp["obs_cb"]), virtual_cb(lp["n"], lp["ca"], lp["c"]), lp["obs_cb"])
        polar, azim = cb_orientation(lp["ca"], cb)
        kap, tau = frenet(lp["ca"])
        m = len(polar)                       # interior residues 1 .. len-2
        lab = (lp["dmin"][1:-1] < args.contact_cutoff).astype(int).tolist()
        col = {
            "polar": polar.tolist(),
            "azimuth": azim.tolist(),
            "|azimuth|": np.abs(azim).tolist(),
            "kappa": kap[:m].tolist(),
            "|tau|": np.abs(np.concatenate([tau, [np.nan]]))[:m].tolist(),
            "relpos": [abs(2 * (i + 1) - (len(lp["ca"]) - 1)) / (len(lp["ca"]) - 1) for i in range(m)],
        }
        keep = [i for i in range(m) if not any(np.isnan(col[k][i]) for k in col)]
        if not keep:
            continue
        labels.extend(lab[i] for i in keep)
        for k in feats:
            feats[k].extend(col[k][i] for i in keep)
            groups[k].append(([lab[i] for i in keep], [col[k][i] for i in keep]))
    print(f"  residues: {len(labels):,}   contacting: {sum(labels):,} "
          f"({100*sum(labels)/len(labels):.1f}%)\n")
    print(f"  {'feature':<12}{'ROC-AUC':>10}{'95% CI (junction bootstrap)':>32}")
    for k in feats:
        a = auc(labels, feats[k])
        a = max(a, 1 - a)                    # direction-agnostic: report discriminability
        lo, hi = boot_auc(groups[k], args.boot, rng)
        lo, hi = (lo, hi) if a == auc(labels, feats[k]) else (1 - hi, 1 - lo)
        print(f"  {k:<12}{a:>10.3f}{f'[{lo:.3f}, {hi:.3f}]':>32}")
    print("\n  relpos (distance from the loop apex) is the baseline any shape feature must beat.")

    # ---------------------------------------------------------------- 2b. incremental value
    print("\n=== 2b. does any of it add ON TOP of position? (5-fold logistic, grouped by junction) ===")
    print("  A marginal AUC is not evidence: every one of these features correlates with position.")
    incremental(loops, args)

    # ---------------------------------------------------------------- 3. Ramachandran gate
    print("\n=== 3. Ramachandran gate: is (phi, psi) recoverable from (kappa, tau)? ===")
    X, Yp, Ys = [], [], []
    for lp in loops:
        if len(lp["ca"]) < 6:
            continue
        kap, tau = frenet(lp["ca"])
        phi, psi = ramachandran(lp["n"], lp["ca"], lp["c"])
        m = min(len(kap) - 1, len(tau), len(phi))
        for i in range(m):
            v = [kap[i], tau[i], kap[i + 1]]
            if any(np.isnan(v)) or np.isnan(phi[i]) or np.isnan(psi[i]):
                continue
            X.append(v)
            Yp.append(phi[i])
            Ys.append(psi[i])
    X, Yp, Ys = np.array(X), np.array(Yp), np.array(Ys)
    idx = np.arange(len(X))
    np.random.default_rng(args.seed).shuffle(idx)
    folds = np.array_split(idx, 5)
    print(f"  residues: {len(X):,}, 5-fold k-NN (k=15) on (kappa_i, tau_i, kappa_i+1)")
    for name, Y in (("phi", Yp), ("psi", Ys)):
        preds = np.zeros(len(X))
        for f in folds:
            tr = np.setdiff1d(idx, f)
            preds[f] = knn_circular(X[tr], Y[tr], X[f])
        r2 = circular_r2(Y, preds)
        verdict = "REDUNDANT" if r2 >= 0.8 else "adds information"
        print(f"    circular R^2({name} | kappa, tau) = {r2:.3f}   -> {verdict}")
    print("\n  Gate: R^2 >= 0.80 on BOTH means Ramachandran is not built.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
