#!/usr/bin/env python3
"""Validate the native _geom binder descriptors against the manuscript TCRvdb features + AF.

Resolves the "618 structures" question and guards against its two hazards:
  * all 618 rows are labeled (padj present) — no NaN-padj mislabeling;
  * 4 TCR-pMHC complexes are duplicated as two structural models each (identical name/seq/padj,
    different TCR_hash + PDB) -> we DEDUP to 614 unique complexes (by `name`) so a complex's two
    models cannot straddle CV folds.

The TCRvdb model PDBs use canonical chains A,B = TCR, C = peptide, D = MHC (per the manuscript
geom_*/relax_worker scripts), so features are extracted by chain letter — no arda annotation needed.

Checks:
  1. native _geom {pm_cov_ntcr, chain_balance, n_hbond} vs the manuscript columns
     {network_pm_cov_ntcr, geometry_chain_balance, geometry_n_hbond} — Pearson r (reproduction gate).
  2. per-epitope + pooled marginal-over-AF: AUC of AF confidence vs native geometry vs native+AF
     (5-fold CV logistic), on the deduped 614. A term earns its place only if it beats/adds to AF.

Usage: python scripts/binder_validate.py [--data <tcrvdb dir>] [--limit N]
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import numpy as np

from tcren.structure import parse_structure

MS_DEFAULT = Path("/Users/mikesh/vcs/manuscripts/2026-tcren2/data/tcrvdb")
AF_COLS = ["iptm", "ranking_confidence", "tcr-pmhc_iptm", "ptm", "plddt"]
POLAR = {"N", "O"}
# Biopython ShrakeRupley element radii (Å), to match the manuscript geometry_dSASA_interface.
BONDI = {"H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "P": 1.80}


def native_dsasa(structure):
    """Interface ΔSASA = (SASA(TCR alone)+SASA(pMHC alone)) − (their SASA in the full complex).

    One full-complex SASA pass (split into TCR A,B and pMHC C,D) + two "alone" passes.
    """
    from tcren import _geom
    xyz, rad, tag = [], [], []  # tag 0 = TCR(A,B), 1 = pMHC(C,D)
    for c in structure.chains:
        role = 0 if c.chain_id in ("A", "B") else (1 if c.chain_id in ("C", "D") else -1)
        if role < 0:
            continue
        for r in c.residues:
            for a in r.atoms:
                if a.element == "H":
                    continue
                xyz.append(a.coord)
                rad.append(BONDI.get(a.element, 1.70))
                tag.append(role)
    xyz = np.asarray(xyz, float).reshape(-1, 3)
    rad = np.asarray(rad, float)
    tag = np.asarray(tag)
    if len(xyz) == 0 or (tag == 0).sum() == 0 or (tag == 1).sum() == 0:
        return float("nan")
    full = _geom.shrake_rupley(xyz, rad, 1.4, 100)
    tcr_full, pmhc_full = full[tag == 0].sum(), full[tag == 1].sum()
    tcr_alone = _geom.shrake_rupley(xyz[tag == 0], rad[tag == 0], 1.4, 100).sum()
    pmhc_alone = _geom.shrake_rupley(xyz[tag == 1], rad[tag == 1], 1.4, 100).sum()
    return float((tcr_alone + pmhc_alone) - (tcr_full + pmhc_full))


def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return math.nan


def auc(scores, labels):
    pos = [s for s, l in zip(scores, labels) if l and not math.isnan(s)]
    neg = [s for s, l in zip(scores, labels) if not l and not math.isnan(s)]
    if not pos or not neg:
        return float("nan")
    return sum((p > n) + 0.5 * (p == n) for p in pos for n in neg) / (len(pos) * len(neg))


def _chain_atoms(structure, chain_id):
    """(heavy-atom xyz (N,3), per-atom residue seq_index (N,), polar-atom xyz (M,3)) for one chain."""
    xyz, res, polar = [], [], []
    ch = next((c for c in structure.chains if c.chain_id == chain_id), None)
    if ch is None:
        return np.zeros((0, 3)), np.zeros(0, np.int32), np.zeros((0, 3))
    for r in ch.residues:
        for a in r.atoms:
            if a.element == "H":
                continue
            xyz.append(a.coord)
            res.append(r.seq_index)
            if a.element in POLAR:
                polar.append(a.coord)
    return (np.asarray(xyz, float).reshape(-1, 3),
            np.asarray(res, np.int32),
            np.asarray(polar, float).reshape(-1, 3))


def native_features(pdb_path):
    from tcren import _geom

    s = parse_structure(pdb_path)
    a_xyz, a_res, a_pol = _chain_atoms(s, "A")
    b_xyz, b_res, b_pol = _chain_atoms(s, "B")
    c_xyz, _c_res, c_pol = _chain_atoms(s, "C")
    d_xyz, _d_res, _d_pol = _chain_atoms(s, "D")
    if len(c_xyz) == 0 or (len(a_xyz) == 0 and len(b_xyz) == 0):
        return None
    cd = _geom.contact_descriptors(a_xyz, a_res, b_xyz, b_res, c_xyz, d_xyz, 5.0, 4.5)
    tcr_pol = np.vstack([a_pol, b_pol]) if len(a_pol) or len(b_pol) else np.zeros((0, 3))
    n_hbond = _geom.interface_hbonds(tcr_pol, c_pol, 3.5) if len(tcr_pol) and len(c_pol) else 0
    return {"pm_cov_ntcr": cd["pm_cov_ntcr"], "chain_balance": cd["chain_balance"],
            "n_hbond": n_hbond, "dSASA": native_dsasa(s)}


def cvp(X, y):
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import StratifiedKFold, cross_val_predict
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler
    return cross_val_predict(make_pipeline(StandardScaler(), LogisticRegression(max_iter=3000)),
                             X, y, cv=StratifiedKFold(5, shuffle=True, random_state=0),
                             method="predict_proba")[:, 1]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data", type=Path, default=MS_DEFAULT)
    ap.add_argument("--limit", type=int, default=None)
    args = ap.parse_args()

    scored = list(csv.DictReader(open(args.data / "tcrvdb_scored.csv")))
    phys = {r["TCR_hash"]: r for r in csv.DictReader(open(args.data / "tcrvdb_physics_features.csv"))}
    clust = {(r["epitope_aa"], r["cdr3aa"]) for r in csv.DictReader(open(args.data / "tcrnet_clusters.csv"))}

    # DEDUP the 4 ensemble-duplicate complexes: keep the first model per `name`.
    seen, rows = set(), []
    for r in scored:
        if r["name"] in seen:
            continue
        seen.add(r["name"])
        rows.append(r)
    print(f"scored rows={len(scored)}  unique complexes (deduped by name)={len(rows)}")
    if args.limit:
        rows = rows[: args.limit]

    nat, keep = [], []
    for i, r in enumerate(rows):
        f = native_features(args.data / (r["TCR_hash"] + ".pdb"))
        if f is None:
            continue
        r["_bind"] = fnum(r["padj"]) < 1e-5
        r["_ep"] = r["epitope_aa"]
        tcrnet = (r["_ep"], r.get("cdr3_alpha_aa")) in clust or (r["_ep"], r.get("cdr3_beta_aa")) in clust
        # denoised (manuscript primary): keep tcrnet-consistent rows (binder&in-cluster or nonbinder&out).
        r["_den"] = (r["_bind"] and tcrnet) or (not r["_bind"] and not tcrnet)
        ph = phys.get(r["TCR_hash"], {})
        r["_ms"] = {c: fnum(ph.get(c)) for c in
                    ("network_pm_cov_ntcr", "geometry_chain_balance", "geometry_n_hbond",
                     "geometry_dSASA_interface")}
        # CDR1/2-vs-CDR3a TCRen potential term (a tcren-native quantity; native wiring pending).
        r["_pp"] = fnum(ph.get("potential_PP_combo_z12_minus_z3a"))
        r["_af"] = [fnum(r.get(c)) for c in AF_COLS]
        nat.append(f)
        keep.append(r)
        if (i + 1) % 100 == 0:
            print(f"  ... {i + 1}/{len(rows)}")

    print(f"\nfeaturized {len(keep)} complexes")

    # 1) native vs manuscript feature reproduction
    pairs = [("pm_cov_ntcr", "network_pm_cov_ntcr"), ("chain_balance", "geometry_chain_balance"),
             ("n_hbond", "geometry_n_hbond"), ("dSASA", "geometry_dSASA_interface")]
    print("\n=== native _geom vs manuscript feature (Pearson r) ===")
    for nk, mk in pairs:
        x = np.array([f[nk] for f in nat], float)
        y = np.array([r["_ms"][mk] for r in keep], float)
        ok = ~(np.isnan(x) | np.isnan(y))
        r = np.corrcoef(x[ok], y[ok])[0, 1] if ok.sum() > 2 else float("nan")
        print(f"  {nk:14s} vs {mk:26s} r={r:.3f}  (n={ok.sum()})")

    # 2) marginal-over-AF (per-epitope + pooled): native geometry (3) and geometry+potential (4).
    y = np.array([r["_bind"] for r in keep], int)
    Xg = np.array([[f["pm_cov_ntcr"], f["chain_balance"], f["n_hbond"]] for f in nat], float)
    pp = np.array([[r["_pp"]] for r in keep], float)
    pp = np.where(np.isnan(pp), np.nanmedian(pp), pp)
    ds = np.array([[f["dSASA"]] for f in nat], float)
    ds = np.where(np.isnan(ds), np.nanmedian(ds), ds)
    X4 = np.hstack([Xg, pp])       # geom3 + CDR1/2-vs-CDR3a potential
    X5 = np.hstack([Xg, pp, ds])   # + interface ΔSASA
    AF = np.array([r["_af"] for r in keep], float)
    AF = np.where(np.isnan(AF), np.nanmedian(AF, axis=0), AF)
    den = np.array([r["_den"] for r in keep], bool)
    eps = sorted({r["_ep"] for r in keep})
    print("\n=== marginal over AF: AUC (5-fold CV logistic); native features ===")
    print("  scope    set      n  bind    AF   struct5  s5+AF")
    scopes = [("POOLED", np.ones(len(keep), bool))] + \
             [(e[:5], np.array([r["_ep"] == e for r in keep])) for e in eps]
    for scope, base in scopes:
        for tag, idx in [("RAW", base), ("DEN", base & den)]:
            yi = y[idx]
            if yi.sum() < 3 or (len(yi) - yi.sum()) < 3:
                print(f"  {scope:8s} {tag:3s}  too few, skip"); continue
            a_af = auc(list(cvp(AF[idx], yi)), list(yi.astype(bool)))
            a5 = auc(list(cvp(X5[idx], yi)), list(yi.astype(bool)))
            a5af = auc(list(cvp(np.hstack([X5[idx], AF[idx]]), yi)), list(yi.astype(bool)))
            star = "  **" if a5af > a_af else ""
            print(f"  {scope:8s} {tag:3s} {idx.sum():4d} {yi.sum():4d}  {a_af:.3f}  {a5:.3f}   {a5af:.3f}{star}")


if __name__ == "__main__":
    main()
