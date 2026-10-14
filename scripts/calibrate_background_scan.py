"""Pick the chi grid and the clash tolerance against one invariant. 2026-09-08.

`background.scan_structure` calls a candidate INADMISSIBLE when every rotamer of it collides. That
is only meaningful if the grid can reproduce a conformation the residue actually adopts -- otherwise
"cannot fit" is confounded with "my grid is too coarse". The invariant that separates them needs no
new data: **a residue the crystal deposited must be admissible at its own position.** The crystal is
the proof that it fits.

So this sweeps (step, clash_tol) and reports, per setting, the fraction of deposited residues their
own scan calls admissible, the mean |A_i|, and the wall time. The setting to ship is the cheapest
one that reaches ~1.0 on the invariant; anything below that is reporting grid coarseness as
chemistry.

    python scripts/calibrate_background_scan.py --structures ~/hf/tcren_structures/Native2026 -n 40
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from tcren.annotation.batch import iter_annotated_set  # noqa: E402
from tcren.energetics import background as bg  # noqa: E402

SETTINGS = [(120.0, 0.4), (120.0, 0.6), (120.0, 0.8),
            (60.0, 0.4), (60.0, 0.6), (60.0, 0.8),
            (30.0, 0.4), (30.0, 0.6)]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--structures", default=str(Path.home() / "hf/tcren_structures/Native2026"))
    ap.add_argument("-n", "--limit", type=int, default=40)
    ap.add_argument("--side", default="to", choices=("to", "from"))
    ap.add_argument("--interface", default="tcr_peptide")
    args = ap.parse_args()

    # One batched annotation pass over the set, never one call per structure (CLAUDE.md 0-mmseqs).
    structures = []
    for s in iter_annotated_set(args.structures):
        structures.append(s)
        if len(structures) >= args.limit:
            break
    print(f"{len(structures)} structures annotated", file=sys.stderr)

    print(f"{'step':>5s} {'tol':>4s} {'native adm':>11s} {'n_res':>7s} {'mean |A|':>9s} "
          f"{'|A|=20':>7s} {'ms/struct':>10s}")
    for step, tol in SETTINGS:
        nat, sizes, t0 = [], [], time.perf_counter()
        for s in structures:
            try:
                r, pos, _par = bg.scan_structure(s, args.interface, args.side,
                                                 step=step, clash_tol=tol)
            except Exception:
                continue
            sizes.append(r.admissible.sum(1))
            nat.extend(r.admissible[i, bg.AA.index(res.aa)]
                       for i, (_c, res) in enumerate(pos) if res.aa in bg.AA)
        el = (time.perf_counter() - t0) / max(1, len(structures)) * 1000
        a = np.concatenate(sizes)
        print(f"{step:5.0f} {tol:4.1f} {np.mean(nat):11.4f} {len(nat):7d} {a.mean():9.2f} "
              f"{np.mean(a == 20):7.3f} {el:10.1f}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
