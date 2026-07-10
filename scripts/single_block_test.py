#!/usr/bin/env python
"""Is one gap block enough? And where does the block go when lengths differ?

The single-gap-block aligner restricts alignment to exactly one contiguous indel. Whether that
restriction costs anything cannot be answered with ``structural_block_position``, which is an
argmin over the single-block family and so cannot express any other answer. This script uses
``loops.structural_align`` -- iterative superposition plus an *unrestricted* affine DP, free to
open any number of gap blocks -- as a model-independent oracle, and asks:

  G2  what fraction of true structural correspondences ARE a single contiguous block?
  G1  what does the restriction cost, in CA-RMSD, when it is wrong?
  P1  when it is one block, where does the block sit, in Cys-anchor coordinates? Does it match
      the sequence centre, or the germline untemplated span?

Everything is resampled by JUNCTION, never by pair: 372 crystal junctions collapse to 199
unique sequences and generate ~7,000 pairs, so pair-level intervals are anticonservative by
several fold.

Usage:
    python scripts/single_block_test.py [--d-max 4] [--per-cell 400] [--boot 2000]

2026-07-10
"""
from __future__ import annotations

import argparse
import collections
import itertools
import os
import random
import statistics as stt
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
sys.path.insert(0, os.path.dirname(__file__))

from _harvest import collapse, crystal_noise_floor, harvest  # noqa: E402

from tcren.loops import (  # noqa: E402
    gap_runs, is_single_block, structural_align, structural_block_position,
)

# Germline untemplated-span centre, measured model-free from VDJdb V/J calls against the IMGT
# germline residue strings (n = 55,078 TRA / 113,354 TRB unique (cdr3, v, j)).
# Templated prefix p, templated suffix s, both medians.
GERMLINE = {"TRA": (3, 8), "TRB": (4, 5)}


def block_start_predictions(chain_type, m, d):
    """Where each rule says the block of length ``d`` starts in a length-``m`` junction."""
    p, s = GERMLINE[chain_type]
    core_lo, core_hi = p, m - s          # the untemplated span [p, m - s)
    centre = (m - d) // 2                # central_prior's argmin
    germ = max(0, min(m - d, (core_lo + core_hi - d) // 2))   # centred in the untemplated span
    return {"central": centre, "germline": germ}


def cluster_bootstrap(idx_pairs, value, n_loops, boot, rng):
    """Resample JUNCTIONS with replacement and rebuild the induced pair set."""
    keyed = {(a, b): v for a, b, v in idx_pairs}
    out = []
    for _ in range(boot):
        draw = [rng.randrange(n_loops) for _ in range(n_loops)]
        vals = []
        for x, y in itertools.combinations(range(len(draw)), 2):
            a, b = draw[x], draw[y]
            if a == b:
                continue
            v = keyed.get((a, b)) or keyed.get((b, a))
            if v is not None:
                vals.append(v)
        if vals:
            out.append(value(vals))
    out.sort()
    if not out:
        return float("nan"), float("nan")
    return out[int(0.025 * len(out))], out[int(0.975 * len(out))]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--structures", default="data/Canonical2026/*.pdb.gz")
    ap.add_argument("--markup", default="notebooks/natcompsci2022/results_new/markup_2026.csv")
    ap.add_argument("--d-max", type=int, default=4)
    ap.add_argument("--per-cell", type=int, default=400, help="pairs sampled per (chain, d)")
    ap.add_argument("--boot", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    rng = random.Random(args.seed)

    t0 = time.perf_counter()
    loops = harvest(args.structures, args.markup)
    reps, mult = collapse(loops)
    noise = crystal_noise_floor(loops)
    print(f"{len(loops)} typed junctions -> {len(reps)} unique sequences "
          f"({len(loops)/len(reps):.2f}x crystal redundancy), {time.perf_counter()-t0:.0f}s")
    if noise:
        print(f"crystal noise floor (same sequence, different crystal): n={len(noise)}, "
              f"median CA-RMSD {stt.median(noise):.3f} A, p90 {sorted(noise)[int(.9*len(noise))]:.3f} A")
        print("  -> no structural effect smaller than this is resolvable, whatever the sample size.\n")

    by_chain = collections.defaultdict(list)
    for lp in reps:
        by_chain[lp.chain_type].append(lp)

    results = collections.defaultdict(list)   # (chain, d) -> [(ia, ib, rec)]
    t0 = time.perf_counter()
    for ct, group in sorted(by_chain.items()):
        cand = collections.defaultdict(list)
        for i, j in itertools.combinations(range(len(group)), 2):
            a, b = group[i], group[j]
            d = abs(len(a.seq) - len(b.seq))
            if 1 <= d <= args.d_max and min(len(a.seq), len(b.seq)) >= 8:
                cand[d].append((i, j))
        for d in sorted(cand):
            pairs = cand[d]
            if len(pairs) > args.per_cell:
                pairs = random.Random(f"{args.seed}:{ct}:{d}").sample(pairs, args.per_cell)
            for i, j in pairs:
                a, b = group[i], group[j]
                # orient so `a` is the longer loop; the block lives in it
                if len(a.seq) < len(b.seq):
                    i, j, a, b = j, i, b, a
                _, rmsd_free, ops = structural_align(a.ca, b.ca)
                blk, rmsd_blk, _ = structural_block_position(a.ca, b.ca)
                runs = gap_runs(ops)
                single = is_single_block(ops)
                start = runs[0][1] if single and runs else None
                results[(ct, d)].append((i, j, {
                    "single": single, "n_runs": len(runs), "start": start,
                    "rmsd_free": rmsd_free, "rmsd_block": rmsd_blk, "blk": blk,
                    "m": len(a.seq), "d": d, "chain": ct,
                }))
    print(f"aligned {sum(len(v) for v in results.values()):,} pairs in {time.perf_counter()-t0:.0f}s\n")

    # ---------------------------------------------------------------- G2
    print("=== G2: is one gap block enough? (model-independent oracle) ===")
    print(f"  {'chain':<6}{'d':>3}{'pairs':>8}{'single-block':>15}{'95% CI (junction bootstrap)':>32}")
    for (ct, d), recs in sorted(results.items()):
        frac = sum(r["single"] for _, _, r in recs) / len(recs)
        lo, hi = cluster_bootstrap(
            [(a, b, r["single"]) for a, b, r in recs],
            lambda vs: sum(vs) / len(vs), len(by_chain[ct]), args.boot, rng)
        print(f"  {ct:<6}{d:>3}{len(recs):>8}{100*frac:>14.1f}%{f'[{100*lo:.1f}, {100*hi:.1f}]':>32}")

    # ---------------------------------------------------------------- G1
    print("\n=== G1: what does the single-block restriction cost in CA-RMSD? ===")
    print(f"  {'chain':<6}{'d':>3}{'free (A)':>11}{'1-block (A)':>13}{'excess (A)':>12}"
          f"{'excess when NOT 1-block':>26}")
    for (ct, d), recs in sorted(results.items()):
        free = stt.median(r["rmsd_free"] for _, _, r in recs)
        blk = stt.median(r["rmsd_block"] for _, _, r in recs)
        multi = [r for _, _, r in recs if not r["single"]]
        exc_m = stt.median(r["rmsd_block"] - r["rmsd_free"] for r in multi) if multi else float("nan")
        print(f"  {ct:<6}{d:>3}{free:>11.3f}{blk:>13.3f}{blk-free:>12.3f}{exc_m:>26.3f}")

    # ---------------------------------------------------------------- seed robustness
    print("\n=== robustness: does the anchor-seeded DP agree with an exhaustive, unseeded search? ===")
    print("  structural_block_position superposes EVERY single-block layout and takes the argmin.")
    print("  It has no seed and no iteration. If the free DP's block matches it, the anchor seed")
    print("  is not steering the answer.")
    print(f"  {'chain':<6}{'d':>3}{'n(1-block)':>12}{'exact agree':>14}{'within 1':>11}")
    for (ct, d), recs in sorted(results.items()):
        sb = [r for _, _, r in recs if r["single"] and r["start"] is not None]
        if not sb:
            continue
        ex = sum(r["start"] == r["blk"] for r in sb) / len(sb)
        w1 = sum(abs(r["start"] - r["blk"]) <= 1 for r in sb) / len(sb)
        print(f"  {ct:<6}{d:>3}{len(sb):>12}{100*ex:>13.1f}%{100*w1:>10.1f}%")

    # ---------------------------------------------------------------- P1
    print("\n=== P1: when it IS one block, where does the block start? "
          "(Cys-anchor residue offset) ===")
    print(f"  {'chain':<6}{'d':>3}{'n':>6}{'observed':>11}{'central':>10}{'germline':>10}"
          f"{'|obs-cen|':>12}{'|obs-germ|':>12}")
    agree = collections.Counter()
    for (ct, d), recs in sorted(results.items()):
        sb = [r for _, _, r in recs if r["single"] and r["start"] is not None]
        if len(sb) < 10:
            continue
        obs = stt.median(r["start"] for r in sb)
        cen = stt.median(block_start_predictions(ct, r["m"], d)["central"] for r in sb)
        ger = stt.median(block_start_predictions(ct, r["m"], d)["germline"] for r in sb)
        e_cen = stt.median(abs(r["start"] - block_start_predictions(ct, r["m"], d)["central"]) for r in sb)
        e_ger = stt.median(abs(r["start"] - block_start_predictions(ct, r["m"], d)["germline"]) for r in sb)
        agree[(ct, "central")] += sum(r["start"] == block_start_predictions(ct, r["m"], d)["central"] for r in sb)
        agree[(ct, "germline")] += sum(r["start"] == block_start_predictions(ct, r["m"], d)["germline"] for r in sb)
        agree[(ct, "n")] += len(sb)
        print(f"  {ct:<6}{d:>3}{len(sb):>6}{obs:>11.1f}{cen:>10.1f}{ger:>10.1f}"
              f"{e_cen:>12.2f}{e_ger:>12.2f}")

    print("\n  exact-hit rate of each rule on the single-block pairs:")
    for ct in sorted(by_chain):
        n = agree[(ct, "n")]
        if not n:
            continue
        print(f"    {ct}: central {100*agree[(ct,'central')]/n:5.1f}%   "
              f"germline {100*agree[(ct,'germline')]/n:5.1f}%   (n = {n})")
    print("\n  P2 (the powered test) lives in vdjmatch: same question, n = 51,852 TRA sequences,")
    print("  no structures needed. This crystal set can only corroborate it.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
