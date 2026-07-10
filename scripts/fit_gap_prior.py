#!/usr/bin/env python
"""Fit the gap-position prior for single-gap-block junction alignment.

A sequence score alone cannot decide *where* the gap in a length-different CDR3 pair goes:
against structure, minimum-BLOSUM62 agrees with the structurally correct block position
about as often as picking at random. seqtree.gapblock therefore takes a positional prior.
This script measures that prior from backbone geometry.

Method. Harvest every ``C ... [FW]GXG`` junction from a structure set, keep the ones that
satisfy all three Leszczynski-Rose omega-loop conditions (tcren.loops.is_omega_loop), then
for every pair of junctions from the same chain type differing by 1..d_max residues, find
the block position whose induced residue correspondence superposes best (loop-local Kabsch).
Pairs are NOT restricted to a shared epitope: the prior is a statement about loop geometry,
not specificity, so every pair of anchored loops is informative -- which is what makes the
fit set large enough to trust.

Outputs data/gap_prior.tsv: for each (block length d, normalised position bin) the empirical
probability that the structurally-best block sits there, plus the integer penalty
``lam * |i - L/2|`` that best reproduces it.

Usage:
    python scripts/fit_gap_prior.py [--structures data/Canonical2026] [--d-max 3] [--out data/gap_prior.tsv]

2026-07-10
"""
from __future__ import annotations

import argparse
import collections
import csv
import glob
import gzip
import itertools
import math
import os
import random
import statistics as stt
import sys
import warnings

import numpy as np

warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from tcren.loops import (  # noqa: E402
    find_junctions, is_omega_loop, kabsch_rmsd, block_layouts,
    structural_block_position, omega_stats,
)

# BLOSUM62 Gram penalty, parsed from seqtree's own table so the two stay in lockstep.
SEQTREE_INC = "/Users/mikesh/vcs/code/seqtree/src/blosum62.inc"
AA = "ACDEFGHIKLMNPQRSTVWY"


def blosum_gram():
    import re
    order = "ARNDCQEGHILKMFPSTWYVBZX*"
    txt = open(SEQTREE_INC).read()
    mark = "kBlosum62[24 * 24] = {"
    body = txt[txt.index(mark) + len(mark):]
    nums = [int(x) for x in re.findall(r"-?\d+", re.sub(r"//.*", "", body))][:576]
    sim = {(order[i], order[j]): nums[i * 24 + j] for i in range(24) for j in range(24)}
    return {(a, b): max(0, sim[a, a] + sim[b, b] - 2 * sim[a, b]) for a in AA for b in AA}


PEN = blosum_gram()
SCALE = int(stt.median([PEN[a, b] for a in AA for b in AA if a != b]))


def load_junctions(pattern: str) -> list[tuple[str, str, np.ndarray]]:
    """(chain_type_hint, junction_seq, CA coords) for every omega-loop junction on disk."""
    from Bio.PDB import PDBParser
    from Bio.PDB.Polypeptide import index_to_one, three_to_index

    parser = PDBParser(QUIET=True)
    out = []
    for path in sorted(glob.glob(pattern)):
        try:
            with gzip.open(path, "rt") as fh:
                model = parser.get_structure(os.path.basename(path), fh)[0]
        except Exception:
            continue
        for chain in model:
            res = [r for r in chain if r.id[0] == " " and "CA" in r]
            if not 90 <= len(res) <= 130:      # a TCR/Ig variable domain
                continue
            try:
                seq = "".join(index_to_one(three_to_index(r.get_resname())) for r in res)
            except Exception:
                continue
            ca = np.array([r["CA"].get_coord() for r in res], dtype=float)
            for j in find_junctions(seq, ca):
                if set(j.seq) <= set(AA) and is_omega_loop(j.ca, relax_length=True):
                    out.append((chain.id, j.seq, j.ca))
    return out


def sequence_block_position(q: str, r: str) -> int:
    """argmin BLOSUM62-Gram block position, no prior. The thing we are testing."""
    best, best_i = None, 0
    for i, pairs in enumerate(block_layouts(len(q), len(r))):
        s = sum(PEN[q[x], r[y]] for x, y in pairs)
        if best is None or s < best:
            best, best_i = s, i
    return best_i


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--structures", default="data/Canonical2026/*.pdb.gz")
    ap.add_argument("--d-max", type=int, default=3)
    ap.add_argument("--out", default="data/gap_prior.tsv")
    ap.add_argument("--max-pairs", type=int, default=6000)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    rng = random.Random(args.seed)

    print(f"BLOSUM62 Gram scale (median mismatch) = {SCALE}")
    juncs = load_junctions(args.structures)
    print(f"harvested {len(juncs):,} omega-loop junctions from {args.structures}")
    if not juncs:
        print("no junctions found", file=sys.stderr)
        return 1
    lens = [len(s) for _, s, _ in juncs]
    print(f"  junction length: {min(lens)}-{max(lens)}, median {stt.median(lens):.0f}")

    by_chain = collections.defaultdict(list)
    for cid, seq, ca in juncs:
        by_chain[cid].append((seq, ca))

    # Every pair of same-chain-type junctions differing by 1..d_max residues.
    pairs = []
    for cid, group in by_chain.items():
        for (sa, xa), (sb, xb) in itertools.combinations(group, 2):
            d = abs(len(sa) - len(sb))
            if 1 <= d <= args.d_max and sa != sb and min(len(sa), len(sb)) >= 8:
                pairs.append((sa, xa, sb, xb, d))
    rng.shuffle(pairs)
    pairs = pairs[:args.max_pairs]
    print(f"  usable pairs (1 <= d <= {args.d_max}, same chain id): {len(pairs):,}\n")

    # --- structurally-best block position per pair -----------------------------------
    rows = []       # (d, L, i_struct, i_seq, rmsd_struct, rmsd_seq, rmsd_random)
    for sa, xa, sb, xb, d in pairs:
        if len(sa) < len(sb):
            sa, xa, sb, xb = sb, xb, sa, xa
        i_str, rm_str, all_rm = structural_block_position(xa, xb)
        i_seq = sequence_block_position(sa, sb)
        i_rnd = rng.randrange(len(all_rm))
        rows.append((d, min(len(sa), len(sb)), i_str, i_seq, rm_str, all_rm[i_seq], all_rm[i_rnd]))

    print("=== where does the structurally-correct block sit? ===")
    norm = [i / L for _, L, i, _, _, _, _ in rows]
    print(f"  normalised position i/L: mean {stt.mean(norm):.3f}  median {stt.median(norm):.3f}"
          f"  sd {stt.stdev(norm):.3f}   (n = {len(norm):,})")
    hist = collections.Counter(round(v * 10) / 10 for v in norm)
    for k in sorted(hist):
        bar = "#" * max(1, int(60 * hist[k] / max(hist.values())))
        print(f"    {k:.1f} {bar} {hist[k]}")

    print("\n=== can a sequence score find it? ===")
    exact = sum(1 for r in rows if r[2] == r[3]) / len(rows)
    near = sum(1 for r in rows if abs(r[2] - r[3]) <= 1) / len(rows)
    print(f"  BLOSUM62-Gram picks the structural block exactly : {100*exact:5.1f}%")
    print(f"  ... within one position                          : {100*near:5.1f}%")
    print(f"  median CA-RMSD  oracle {stt.median([r[4] for r in rows]):.3f} A"
          f"   sequence {stt.median([r[5] for r in rows]):.3f} A"
          f"   random {stt.median([r[6] for r in rows]):.3f} A")

    # --- fit lam for the central prior, and compare to the empirical prior ------------
    print("\n=== fitting the prior ===")
    print(f"  {'chooser':<40}{'exact':>8}{'within 1':>10}{'median RMSD':>14}")

    def evaluate(name, choose):
        ex = wi = 0
        rms = []
        for (sa, xa, sb, xb, d), row in zip(pairs, rows):
            if len(sa) < len(sb):
                sa, xa, sb, xb = sb, xb, sa, xa
            L = min(len(sa), len(sb))
            i = choose(sa, sb, L, d)
            _, _, all_rm = structural_block_position(xa, xb)
            ex += i == row[2]
            wi += abs(i - row[2]) <= 1
            rms.append(all_rm[i])
        print(f"  {name:<40}{100*ex/len(rows):7.1f}%{100*wi/len(rows):9.1f}%{stt.median(rms):13.3f} A")
        return ex / len(rows), stt.median(rms)

    def seq_only(q, r, L, d):
        return sequence_block_position(q, r)

    def centre_only(q, r, L, d):
        return min(range(L + 1), key=lambda i: (abs(2 * i - L), i))

    def with_prior(lam):
        def choose(q, r, L, d):
            best, best_i = None, 0
            for i, pairs_ in enumerate(block_layouts(len(q), len(r))):
                s = sum(PEN[q[x], r[y]] for x, y in pairs_) + (lam * abs(2 * i - L)) // 2
                if best is None or s < best:
                    best, best_i = s, i
            return best_i
        return choose

    evaluate("BLOSUM62 only (no prior)", seq_only)
    evaluate("central block, sequence ignored", centre_only)
    best_lam, best_ex = None, -1.0
    for mult in (0.5, 1.0, 1.5, 2.0, 3.0, 5.0):
        lam = int(mult * SCALE)
        ex, _ = evaluate(f"BLOSUM62 + central prior (lam={lam} = {mult}*scale)", with_prior(lam))
        if ex > best_ex:
            best_ex, best_lam = ex, lam
    print(f"  {'ORACLE (min structural RMSD)':<40}{100.0:7.1f}%{100.0:9.1f}%"
          f"{stt.median([r[4] for r in rows]):13.3f} A")
    print(f"\n  best lam = {best_lam}  ({best_lam/SCALE:.1f} * matrix scale)")

    # --- emit the empirical prior ------------------------------------------------------
    counts = collections.defaultdict(collections.Counter)
    for d, L, i_str, *_ in rows:
        counts[d][round(10 * i_str / L)] += 1
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["# fitted by scripts/fit_gap_prior.py from", args.structures])
        w.writerow(["# n_pairs", len(rows), "blosum62_scale", SCALE, "best_lam", best_lam])
        w.writerow(["d", "rel_pos_decile", "n", "p"])
        for d in sorted(counts):
            tot = sum(counts[d].values())
            for b in sorted(counts[d]):
                w.writerow([d, b / 10, counts[d][b], f"{counts[d][b]/tot:.5f}"])
    print(f"\nwrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
