#!/usr/bin/env python3
"""Does peptide conformational stability explain where an additive contact model fails?

    "Supplementary Figure 4 in our JCI paper shows quite a substantial intra-peptide interaction
     between P3 and P6 that stabilises the central peptide bulge recognised by the TCR. It made me
     wonder whether this sort of internal peptide stabilisation might explain why 4C6 appears to be
     less well described by an additive contact model than the other systems. Poor binders could
     perhaps still make many contacts but fail to stabilise the productive peptide conformation."
        -- Sewell, 2026-08, manuscripts/2026-tcren/suggestions/sewell.txt

The CPL set gives the design for free: ~160 modelled complexes for the BEST peptides and ~160 for
the WORST, for each of seven clones. For every structure compute

  * the additive contact energy (the model whose failures we are explaining), and
  * peptide conformational stability from flexible-backbone MC (`tcren.dynamics`), with the
    intra-peptide term ON and OFF at the same seed, so the comparison is paired,

then ask, per clone, whether stability separates best from worst where the contact energy does not.

    python3 scripts/sewell_stability.py            # full set, ~28 min on 14 cores
    python3 scripts/sewell_stability.py 6          # 6 structures per group, a smoke test
    MC_STEPS=8000 python3 scripts/sewell_stability.py

Result (2102 structures, 2026-08-17): stability beats the contact energy in 4/4 clones where the
contact model fails and 0/3 where it works. See README and CHANGELOG for the numbers.
"""
from __future__ import annotations

import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import polars as pl

CPL = Path(os.environ.get("CPL_DIR", "data/cpl/pdb_cpl"))
STEPS = int(os.environ.get("MC_STEPS", 4000))
TEMP, SIGMA = 20.0, 10.0
LIMIT = int(sys.argv[1]) if len(sys.argv) > 1 else 0


def one(task):
    """One structure -> contact energy + stability with the intra term on and off."""
    from tcren.annotation import classify_chains
    from tcren.contactmap import ContactMap
    from tcren.dynamics import peptide_stability
    from tcren.mhc import annotate_mhc
    from tcren.pipeline import _interface_energy
    from tcren.potential import mj, tcren
    from tcren.structure import parse_structure

    path, clone, is_best = task
    try:
        s = parse_structure(path, pdb_id=path.stem)
        classify_chains(s, organism="human")
        annotate_mhc(s)
        cm = ContactMap.from_structure(s, peptide_internal=True)

        row = {"clone": clone, "is_best": is_best, "pdb": path.stem,
               "peptide": "".join(r.aa for c in s.chains if c.chain_type == "PEPTIDE"
                                  for r in c.residues)}
        row["Phi_tcr_pep"] = float(_interface_energy(cm.interface("tcr_peptide"), tcren()))
        row["Phi_pep_mhc"] = float(_interface_energy(cm.interface("peptide_mhc"), mj()))
        # The intra-peptide CONTACT energy: the additive model's own view of the same interactions
        # whose dynamical effect we are measuring.
        from tcren.scoring import intra_peptide_energy
        row["Phi_intra"] = float(intra_peptide_energy(cm, mj()))
        row["n_intra"] = int(cm.peptide_internal.height) if cm.peptide_internal is not None else 0

        for w, tag in ((1.0, "intra1"), (0.0, "intra0")):
            st = peptide_stability(s, intra_weight=w, n_steps=STEPS, temperature=TEMP,
                                   sigma_deg=SIGMA, seed=0)
            row[f"rmsf_{tag}"] = st.rmsf
            row[f"drift_{tag}"] = st.drift
        row["delta_rmsf"] = row["rmsf_intra0"] - row["rmsf_intra1"]
        return row
    except Exception as exc:                                            # noqa: BLE001 - survey
        return {"clone": clone, "is_best": is_best, "pdb": path.stem,
                "error": f"{type(exc).__name__}: {str(exc)[:70]}"}


def auc(score, label):
    """ROC-AUC via the rank-sum identity; higher score = positive class."""
    from scipy.stats import rankdata

    score, label = np.asarray(score, float), np.asarray(label, bool)
    ok = np.isfinite(score)
    score, label = score[ok], label[ok]
    n1, n0 = label.sum(), (~label).sum()
    if n1 == 0 or n0 == 0:
        return float("nan")
    r = rankdata(score)
    return float((r[label].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


if __name__ == "__main__":
    tasks = []
    for d in sorted(CPL.iterdir()):
        if not d.is_dir() or "_" not in d.name:
            continue
        clone, kind = d.name.rsplit("_", 1)
        pdbs = sorted(d.glob("*.pdb"))
        if LIMIT:
            pdbs = pdbs[:LIMIT]
        tasks += [(p, clone, kind == "best") for p in pdbs]
    print(f"{len(tasks)} structures, {STEPS} MC steps each x2 weights", flush=True)

    t0 = time.perf_counter()
    n_workers = max((os.cpu_count() or 2) - 2, 1)
    with ProcessPoolExecutor(max_workers=n_workers) as ex:
        rows = list(ex.map(one, tasks, chunksize=8))
    print(f"done in {time.perf_counter() - t0:.0f}s on {n_workers} workers", flush=True)

    d = pl.DataFrame(rows, infer_schema_length=None)
    out_path = Path(os.environ.get("OUT", "scratch/sewell_stability.parquet"))
    out_path.parent.mkdir(parents=True, exist_ok=True)
    d.write_parquet(out_path)
    bad = d.filter(pl.col("error").is_not_null()) if "error" in d.columns else d.head(0)
    print(f"{d.height} rows, {bad.height} failed")
    if bad.height:
        print(bad.select("pdb", "error").head(5))
    d = d.filter(pl.col("rmsf_intra1").is_not_null()) if "error" in d.columns else d

    pl.Config.set_tbl_rows(30)
    pl.Config.set_tbl_width_chars(170)

    print("\n=== per clone: AUC(best vs worst) — does each signal separate them? ===")
    print("    contact = -Phi_tcr_pep (lower energy = better binder)")
    print("    stable  = -rmsf      (less motion  = better binder)")
    out = []
    for clone in sorted(set(d["clone"])):
        sub = d.filter(pl.col("clone") == clone)
        lab = sub["is_best"].to_numpy()
        out.append({
            "clone": clone, "n_best": int(lab.sum()), "n_worst": int((~lab).sum()),
            "contact": auc(-sub["Phi_tcr_pep"].to_numpy(), lab),
            "pep_mhc": auc(-sub["Phi_pep_mhc"].to_numpy(), lab),
            "stable": auc(-sub["rmsf_intra1"].to_numpy(), lab),
            "stable_nointra": auc(-sub["rmsf_intra0"].to_numpy(), lab),
            "delta_rmsf": auc(sub["delta_rmsf"].to_numpy(), lab),
            "Phi_intra": auc(-sub["Phi_intra"].to_numpy(), lab),
        })
    res = pl.DataFrame(out).with_columns(
        pl.col(c).round(3) for c in ("contact", "pep_mhc", "stable", "stable_nointra",
                                     "delta_rmsf", "Phi_intra"))
    print(res)

    lab = d["is_best"].to_numpy()
    print("\npooled: contact %.3f | stable %.3f | delta_rmsf %.3f"
          % (auc(-d["Phi_tcr_pep"].to_numpy(), lab), auc(-d["rmsf_intra1"].to_numpy(), lab),
             auc(d["delta_rmsf"].to_numpy(), lab)))

    print("\n=== the hypothesis: does stability carry signal where the contact model fails? ===")
    r = res.with_columns((pl.col("stable") - 0.5).abs().alias("stable_eff"),
                         (pl.col("contact") - 0.5).abs().alias("contact_eff"))
    from scipy.stats import spearmanr
    rho, p = spearmanr(r["contact"].to_numpy(), r["stable"].to_numpy())
    print(f"  Spearman(contact AUC, stability AUC) over {r.height} clones = {rho:+.3f} (p={p:.3f})")
    print("  (a negative correlation is the prediction: stability helps where contacts do not)")

    print("\n=== does the intra-peptide term steady the peptide at all? ===")
    for grp, name in ((True, "best"), (False, "worst")):
        v = d.filter(pl.col("is_best") == grp)["delta_rmsf"].to_numpy()
        v = v[np.isfinite(v)]
        print(f"  {name:5s} n={len(v):4d}  mean delta_rmsf = {v.mean():+.4f} A "
              f"(SE {v.std(ddof=1)/np.sqrt(len(v)):.4f})")
    from scipy.stats import mannwhitneyu
    a = d.filter(pl.col("is_best"))["delta_rmsf"].to_numpy()
    b = d.filter(~pl.col("is_best"))["delta_rmsf"].to_numpy()
    a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    u, p = mannwhitneyu(a, b, alternative="two-sided")
    print(f"  best vs worst: p = {p:.3g}, AUC = {u / (len(a) * len(b)):.3f}")
