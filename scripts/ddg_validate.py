#!/usr/bin/env python3
"""ΔΔG-recovery validation: does native interface energy + sampling recover binding affinity?

The ATLAS set (data/oracle/ddg_direct/ddg_direct_deltas.csv) is 178 TCR-pMHC complexes each with an
experimental binding ``ddG``. Baselines on that set: the dead single-structure TCR:peptide delta
``d_tcr_pep`` gives R²≈0.004; the physical/sampled ``delta_sum`` gives R²≈0.641 (the target). The
scientific claim is that no single-structure scorer recovers affinity — sampling does.

Here we measure, over the 178 complexes, the R² vs experimental ddG of:
  * e_native  — tcren._relax.interface_energy on the crystal pose (a single-structure DOPE scorer);
  * e_relax   — the same after rigid-body DOPE relaxation (tcren.refine.refine_peptide);
  * gap       — e_native − e_relax (the relaxation gap, cf. relax_worker.py abc_gap).
The rotamer repack + flexible-backbone relax (native _relax, later phases) are expected to close the
gap toward 0.641; this establishes where the cheap (rigid) sampling already lands vs the floor.

Usage: python scripts/ddg_validate.py [--limit N] [--relax]   (--relax adds the rigid-relax pass)
"""

from __future__ import annotations

import argparse
import csv
import glob
import math
import time
from pathlib import Path

import numpy as np

from tcren.annotation import classify_chains
from tcren.mhc import annotate_mhc
from tcren.refine import refine_peptide
from tcren.refine.interface import interface_energy
from tcren.structure import parse_structure

MS = Path("/Users/mikesh/vcs/manuscripts/2026-tcren2")
ATLAS = MS / "data/oracle/ddg_direct/ddg_direct_deltas.csv"
SEARCH = ["/Users/mikesh/vcs/code/tcren-ms/data/**", str(MS / "data/**")]


def _find_pdb(name):
    for pat in SEARCH:
        hits = glob.glob(f"{pat}/{name}.pdb.gz", recursive=True) + glob.glob(f"{pat}/{name}.pdb", recursive=True)
        if hits:
            return hits[0]
    return None


def _r2(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 3:
        return float("nan")
    c = np.corrcoef(x[ok], y[ok])[0, 1]
    return c * c


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--relax", action="store_true", help="also compute the rigid-relax pass (slower)")
    ap.add_argument("--out", type=Path, default=Path("scratch/ddg_validate.csv"))
    args = ap.parse_args()

    rows = list(csv.DictReader(open(ATLAS)))
    if args.limit:
        rows = rows[: args.limit]
    out, t0, done = [], time.perf_counter(), 0
    for r in rows:
        pdb = _find_pdb(r["structure_name"])
        if pdb is None:
            continue
        try:
            s = parse_structure(pdb, pdb_id=r["structure_name"])
            classify_chains(s, organism="human")
            annotate_mhc(s)
            e_nat = interface_energy(s)
            e_rel = interface_energy(refine_peptide(s, seed=0)[0]) if args.relax else float("nan")
        except Exception as exc:
            print(f"  skip {r['structure_name']}: {exc}")
            continue
        out.append({"structure": r["structure_name"], "ddG": float(r["ddG"]),
                    "d_tcr_pep": float(r["d_tcr_pep"]), "delta_sum": float(r["delta_sum"]),
                    "e_native": e_nat, "e_relax": e_rel, "gap": e_nat - e_rel})
        done += 1
        if done % 40 == 0:
            print(f"  ... {done}/{len(rows)}  {time.perf_counter()-t0:.0f}s")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(out[0]))
        w.writeheader(); w.writerows(out)

    ddg = [o["ddG"] for o in out]
    print(f"\nn = {len(out)} complexes | wall {time.perf_counter()-t0:.0f}s")
    print("=== R² vs experimental ddG (target delta_sum=0.641, floor d_tcr_pep=0.004) ===")
    print(f"  d_tcr_pep (dead floor)   R² = {_r2([o['d_tcr_pep'] for o in out], ddg):.3f}")
    print(f"  delta_sum (target)       R² = {_r2([o['delta_sum'] for o in out], ddg):.3f}")
    print(f"  e_native (single-struct) R² = {_r2([o['e_native'] for o in out], ddg):.3f}")
    if args.relax:
        print(f"  e_relax  (rigid sampling) R² = {_r2([o['e_relax'] for o in out], ddg):.3f}")
        print(f"  gap      (relax gap)      R² = {_r2([o['gap'] for o in out], ddg):.3f}")
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
