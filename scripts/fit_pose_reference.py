#!/usr/bin/env python3
"""Fit the bundled crystal reference for the pose-consistency descriptors. 2026-08-23

`tcren.pose.pose_consistency` returns raw descriptors; to score a *single* user structure they must
be standardized against the natural interface manifold, exactly as `tcren.cohort.q_score` is
standardized against `q_native_reference()`. This writes that reference.

It writes a SECOND file, `src/tcren/data/pose_native_reference.csv`, and never touches
`q_native_reference.csv`: regenerating that one would move every published `Q`.

Scope matches the potential derivation's hard rule -- alpha-beta TCR : peptide-MHC only, i.e. both
CDR3 loops and a 20-letter peptide present. A pMHC-only or single-chain file has no TCR:peptide
interface and would enter the reference at the origin.

TWO references are fitted, because they answer different questions:

  native  Native2026 crystals. C then reads *departure from the crystal manifold*, i.e. provenance
          -- "was this solved or generated?".
  af      vdjdb_binder_benchmark AlphaFold models, label-blind (both classes). C then reads
          *departure from the typical generated pose*, i.e. which AlphaFold models are outliers
          among their own kind. This is the reference to use when the input is itself generated,
          which is the deployment case: a crystal reference measures the constant AlphaFold offset
          shared by every model and swamps the within-cohort signal.

    python scripts/fit_pose_reference.py                        # native (default)
    python scripts/fit_pose_reference.py \
        --struct-dir ../../projects/2026-tcren2-code/data-large/vdjdb_binder_bm \
        --out src/tcren/data/pose_af_reference.csv
"""
from __future__ import annotations

import argparse
import glob
import sys
import time
from pathlib import Path

import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import tcren  # noqa: E402
from tcren.annotation import classify_chains  # noqa: E402
from tcren.mhc import annotate_mhc_batch  # noqa: E402
from tcren.paper.helpers import annotate_batch  # noqa: E402
from tcren.pose import POSE_FEATURES, pose_consistency  # noqa: E402

AA20 = set("ACDEFGHIKLMNPQRSTVWY")
CHUNK = 128  # one batched mmseqs call per chunk; the whole set in one call is a large index build


def _in_scope(s) -> bool:
    """The derivation's hard rule: both CDR3 loops resolved and a 20-letter peptide."""
    def cdr3(ctype):
        ch = next((c for c in s.chains if c.chain_type == ctype), None)
        return ch is not None and any(r.region_type == "CDR3" for r in (ch.regions or []))

    pep = next((c for c in s.chains if c.chain_type == "PEPTIDE"), None)
    if pep is None or not cdr3("TRA") or not cdr3("TRB"):
        return False
    return bool(pep.sequence()) and set(pep.sequence()) <= AA20


def annotate(structures):
    """Batched TCR then MHC annotation -- the two-pass path; a process pool here deadlocks."""
    recs = annotate_batch(structures, organisms=("human", "mouse"))
    for i, s in enumerate(structures):
        try:
            classify_chains(s, organism="human", autodetect_species=True,
                            precomputed_records=recs[i])
        except Exception as exc:  # a chain arda cannot call is not a reason to lose the batch
            print(f"  ! classify {s.pdb_id}: {type(exc).__name__}: {exc}", file=sys.stderr)
    annotate_mhc_batch(structures)  # SEPARATE pass: classify_chains types TRA/TRB/PEPTIDE only


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--struct-dir", default="data/Native2026")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--out", default="src/tcren/data/pose_native_reference.csv")
    args = ap.parse_args()

    paths = sorted(glob.glob(f"{args.struct_dir}/**/*.pdb*", recursive=True)
                   + glob.glob(f"{args.struct_dir}/**/*.cif*", recursive=True))
    if args.limit:
        paths = paths[: args.limit]
    print(f"{len(paths)} structures in {args.struct_dir}")

    rows, dropped = [], {"parse": 0, "scope": 0, "score": 0}
    t0 = time.time()
    for start in range(0, len(paths), CHUNK):
        batch = []
        for p in paths[start : start + CHUNK]:
            try:
                batch.append(tcren.import_structure(p))
            except Exception:
                dropped["parse"] += 1
        annotate(batch)
        for s in batch:
            if not _in_scope(s):
                dropped["scope"] += 1
                continue
            try:
                d = pose_consistency(s)
            except Exception as exc:
                print(f"  ! score {s.pdb_id}: {type(exc).__name__}: {exc}", file=sys.stderr)
                dropped["score"] += 1
                continue
            rows.append({"pdb.id": s.pdb_id, **d})
        print(f"  {min(start + CHUNK, len(paths))}/{len(paths)}  kept={len(rows)}  "
              f"{time.time() - t0:.0f}s", flush=True)

    df = pl.DataFrame(rows)
    # A descriptor that is nan for this structure cannot contribute a mean/sd; drop rows that are
    # nan on any reference column rather than letting nanmean quietly change the n per column.
    before = df.height
    df = df.filter(pl.all_horizontal([pl.col(c).is_finite() for c in POSE_FEATURES]))
    print(f"\nkept {df.height} of {len(paths)} "
          f"(parse {dropped['parse']}, out-of-scope {dropped['scope']}, "
          f"score {dropped['score']}, non-finite {before - df.height})")
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    df.select("pdb.id", *POSE_FEATURES, "n_contacts", "n_cb_close").write_csv(args.out)
    print(f"wrote {args.out}")
    for c in POSE_FEATURES:
        v = df[c]
        print(f"  {c:24s} median={v.median(): .4f}  mean={v.mean(): .4f}  sd={v.std(): .4f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
