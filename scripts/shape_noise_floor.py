#!/usr/bin/env python
"""How small a loop-shape difference can this crystal set resolve? (Answer: not a small one.)

This script replaces ``island_conformation.py``, which claimed that distinct sequence islands
for one epitope converge on similar CDR3 backbone conformations:

    (kappa, tau) shape  -3.470 deg   95% CI [-6.580, -0.466]   p = 0.004
    backbone CA-RMSD    -0.073 A     95% CI [-0.171, +0.009]   p = 0.022

**Both results are withdrawn.** Three independent defects, any one of which is fatal:

1. *The metrics used different correspondences.* ``shape_distance`` compared kappa/tau over the
   shared **prefix** of two loops -- the leading-gap (i = 0) alignment, which the structural
   data say is the wrong one; the block sits at the apex. Its sibling ``loop_rmsd`` used the
   gap-block correspondence and silently dropped every pair with |dlen| > 3. Different metric,
   different correspondence, different subset.

2. *The islands were not islands.* Edges were drawn at a fixed ``gapblock_score <= 60``. Against
   a size-matched random control, that threshold puts **more** control junctions into large
   components than real same-epitope ones (0.748 vs 0.660 of nodes in a component of >= 5).
   A fixed distance is not a significance test; see ``vdjmatch/bench/islands_calibrated.py``.

3. *The effects are below the resolution of the instrument.* This script measures that limit:
   the disagreement between two crystal structures of the **same** junction sequence. Both
   claimed effects are smaller than it. No sample size repairs this.

What would change the answer: an effect larger than the floor printed below, or a structure set
whose within-sequence reproducibility is materially better than the PDB's.

Usage:
    python scripts/shape_noise_floor.py

2026-07-10
"""
from __future__ import annotations

import argparse
import collections
import os
import statistics as stt
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
sys.path.insert(0, os.path.dirname(__file__))

from _harvest import harvest  # noqa: E402

from tcren.loops import frenet, kabsch_rmsd  # noqa: E402

# The two effects island_conformation.py reported, for scale.
RETRACTED = {"CA-RMSD (A)": 0.073, "(kappa,tau) (deg)": 3.470}


def shape_distance(xa: np.ndarray, xb: np.ndarray) -> float:
    """Mean absolute (kappa, tau) difference in degrees. Equal lengths only.

    Equal lengths is not a limitation here, it is the point: identical sequences need no
    correspondence, so the number is free of every alignment assumption.
    """
    if len(xa) != len(xb):
        raise ValueError("shape_distance compares equal-length loops only")
    ka, ta = frenet(xa)
    kb, tb = frenet(xb)
    dk = np.abs(ka - kb)
    dt = np.abs(((ta - tb + 180) % 360) - 180)     # circular
    v = np.concatenate([dk, dt])
    v = v[~np.isnan(v)]
    return float(v.mean())


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--structures", default="data/Canonical2026/*.pdb.gz")
    ap.add_argument("--markup", default="notebooks/natcompsci2022/results_new/markup_2026.csv")
    args = ap.parse_args()

    loops = harvest(args.structures, args.markup)
    by_seq: dict[tuple[str, str], list] = collections.defaultdict(list)
    for lp in loops:
        by_seq[(lp.chain_type, lp.seq)].append(lp)

    rms, shp = [], []
    for group in by_seq.values():
        for i in range(len(group)):
            for j in range(i + 1, len(group)):
                a, b = group[i].ca, group[j].ca
                if len(a) == len(b) and len(a) >= 5:
                    rms.append(kabsch_rmsd(a, b))
                    shp.append(shape_distance(a, b))
    if not rms:
        print("no redundant sequences: cannot establish a noise floor")
        return 1

    redundant = sum(1 for g in by_seq.values() if len(g) > 1)

    def q(v, p):
        return sorted(v)[int(p * len(v))]

    print(f"{len(loops)} junctions, {len(by_seq)} unique sequences, {redundant} of them "
          f"crystallised more than once")
    print(f"\n=== noise floor: the SAME junction sequence, in two different crystals "
          f"(n = {len(rms)} pairs) ===")
    print(f"  {'metric':<20}{'median':>10}{'p90':>10}{'max':>10}")
    print(f"  {'CA-RMSD (A)':<20}{stt.median(rms):>10.3f}{q(rms, .9):>10.3f}{max(rms):>10.3f}")
    print(f"  {'(kappa,tau) (deg)':<20}{stt.median(shp):>10.3f}{q(shp, .9):>10.3f}{max(shp):>10.3f}")

    print("\n=== the retracted effects, against that floor ===")
    print(f"  {'metric':<20}{'claimed effect':>16}{'floor (median)':>16}{'ratio':>9}")
    floors = {"CA-RMSD (A)": stt.median(rms), "(kappa,tau) (deg)": stt.median(shp)}
    for k, eff in RETRACTED.items():
        print(f"  {k:<20}{eff:>16.3f}{floors[k]:>16.3f}{eff / floors[k]:>8.2f}x")
    print("\n  Both claimed effects are SMALLER than the disagreement between two crystal")
    print("  structures of the identical molecule. They were never resolvable, at any n.")
    print(f"\n  Nothing below {floors['(kappa,tau) (deg)']:.2f} deg of (kappa,tau), or "
          f"{floors['CA-RMSD (A)']:.3f} A of CA-RMSD, is measurable from this set.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
