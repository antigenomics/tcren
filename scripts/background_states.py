"""What the three outcomes actually look like on the crystals. 2026-09-08.

Verification that the collision branch is doing work rather than decorating the code: over
Native2026, per interface position, how many candidates collide, how many are clean non-contacts and
how many contact, plus the distribution of |A_i| -- the size of the admissible background. A scan
where |A_i| = 20 everywhere would mean the branch never fires; one where it is small at buried
positions and large at exposed ones is the result.

    python scripts/background_states.py -n 60
"""
from __future__ import annotations

import argparse
import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from tcren.annotation.batch import iter_annotated_set  # noqa: E402
from tcren.energetics import background as bg  # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--structures", default=str(Path.home() / "hf/tcren_structures/Native2026"))
    ap.add_argument("-n", "--limit", type=int, default=60)
    ap.add_argument("--side", default="to", choices=("to", "from"))
    ap.add_argument("--interface", default="tcr_peptide")
    args = ap.parse_args()

    sizes, natives, per_pos, times = [], [], [], []
    coll = cont = none = 0
    by_aa: Counter = Counter()
    seen_aa: Counter = Counter()
    n_struct = 0
    for s in iter_annotated_set(args.structures):
        if n_struct >= args.limit:
            break
        t0 = time.perf_counter()
        try:
            r, pos, _par = bg.scan_structure(s, args.interface, args.side)
        except Exception:
            continue
        times.append(time.perf_counter() - t0)
        n_struct += 1
        adm = r.admissible[r.pair_i]
        coll += int((~adm).sum())
        none += int((adm & (r.p <= 0)).sum())
        cont += int((adm & (r.p > 0)).sum())
        sizes.append(r.admissible.sum(1))
        per_pos.append(np.array([len(pos)]))
        for i, (_c, res) in enumerate(pos):
            if res.aa in bg.AA:
                seen_aa[res.aa] += 1
                if r.admissible[i, bg.AA.index(res.aa)]:
                    natives.append(1)
                else:
                    natives.append(0)
                    by_aa[res.aa] += 1

    a = np.concatenate(sizes)
    tot = coll + none + cont
    t = np.asarray(times)
    print(f"# background scan states -- {args.interface}, side={args.side}, "
          f"step={bg.DEFAULT_STEP:.0f} deg, clash_tol={bg.DEFAULT_CLASH_TOL} A")
    print(f"{n_struct} structures, {a.size} scanned positions, {tot} (pair, candidate) cells\n")
    print(f"  collision   {coll:8d}  {coll / tot:7.3%}")
    print(f"  no contact  {none:8d}  {none / tot:7.3%}")
    print(f"  contact     {cont:8d}  {cont / tot:7.3%}\n")
    print(f"|A_i|: mean {a.mean():.2f}  median {np.median(a):.0f}  "
          f"p05 {np.percentile(a, 5):.0f}  p95 {np.percentile(a, 95):.0f}  "
          f"frac at 20 {np.mean(a == 20):.3f}  frac <= 10 {np.mean(a <= 10):.3f}")
    print(f"native admissible: {np.mean(natives):.4f} over {len(natives)} deposited residues")
    if by_aa:
        print("  rejections by type (count / deposited): "
              + ", ".join(f"{k} {v}/{seen_aa[k]}" for k, v in by_aa.most_common(8)))
    print(f"wall time per structure: mean {t.mean() * 1000:.1f} ms, "
          f"median {np.median(t) * 1000:.1f} ms, max {t.max() * 1000:.1f} ms")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
