#!/usr/bin/env python3
"""Fit a JOINT reference carrying the Q descriptors and the pose descriptors together. 2026-08-23

`q_native_reference.csv` and `pose_native_reference.csv` are fitted on different structure sets
(369 and 359 rows), so an *extended* Q -- one whitened over both descriptor families at once --
cannot be standardized against either: the whitening needs a covariance, and a covariance needs the
descriptors measured on the SAME structures.

This writes one table per manifold with every column an extended score can draw on: the 4 Q
descriptors, the energies, all 16 pose descriptors and the 4 placement descriptors.

    python scripts/fit_joint_reference.py --struct-dir data/Native2026 \\
        --out src/tcren/data/q_pose_native_reference.csv

Neither `q_native_reference.csv` nor `pose_native_reference.csv` is touched, so no published score
moves.
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
from tcren.cohort import Q_FEATURES_GEOM  # noqa: E402
from tcren.mhc import annotate_mhc_batch  # noqa: E402
from tcren.orient import tcr_placement  # noqa: E402
from tcren.paper.helpers import annotate_batch  # noqa: E402
from tcren.pose import POSE_FEATURES, pose_consistency  # noqa: E402
from tcren.recognition import recognition_table  # noqa: E402

AA20 = set("ACDEFGHIKLMNPQRSTVWY")
CHUNK = 128
ENERGY = ["F_tcr_pep", "F_tcr_mhc", "dF_tcr_pep", "F_cdr12", "F_cdr3a"]
PLACE = ["tcr_height", "shift_u", "shift_w", "offset"]


def _in_scope(s) -> bool:
    def cdr3(ctype):
        ch = next((c for c in s.chains if c.chain_type == ctype), None)
        return ch is not None and any(r.region_type == "CDR3" for r in (ch.regions or []))

    pep = next((c for c in s.chains if c.chain_type == "PEPTIDE"), None)
    if pep is None or not cdr3("TRA") or not cdr3("TRB"):
        return False
    return bool(pep.sequence()) and set(pep.sequence()) <= AA20


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--struct-dir", default="data/Native2026")
    ap.add_argument("--out", default="src/tcren/data/q_pose_native_reference.csv")
    ap.add_argument("--limit", type=int, default=0)
    args = ap.parse_args()

    paths = sorted(glob.glob(f"{args.struct_dir}/**/*.pdb*", recursive=True)
                   + glob.glob(f"{args.struct_dir}/**/*.cif*", recursive=True))
    if args.limit:
        paths = paths[: args.limit]
    print(f"{len(paths)} structures")

    t0 = time.time()
    items = [(Path(p).name.split(".pdb")[0].split(".cif")[0], p) for p in paths]
    # threads=1 deliberately: recognition_table's process pool forks after mmseqs/BLAS have spawned
    # threads and dies with BrokenProcessPool.
    rec = [r for r in recognition_table(items, full=True, scores=True, threads=1, chunk=CHUNK)
           if "error" not in r]
    print(f"recognition: {len(rec)} rows in {time.time() - t0:.0f}s")

    pose_rows = []
    for start in range(0, len(paths), CHUNK):
        batch, ids = [], []
        for p in paths[start : start + CHUNK]:
            try:
                batch.append(tcren.import_structure(p))
                ids.append(Path(p).name.split(".pdb")[0].split(".cif")[0])
            except Exception:
                pass
        recs = annotate_batch(batch, organisms=("human", "mouse"))
        for i, s in enumerate(batch):
            try:
                classify_chains(s, organism="human", autodetect_species=True,
                                precomputed_records=recs[i])
            except Exception:
                pass
        annotate_mhc_batch(batch)
        for sid, s in zip(ids, batch):
            if not _in_scope(s):
                continue
            row = {"complex.id": sid}
            try:
                row.update(pose_consistency(s))
            except Exception:
                continue
            try:
                p = tcr_placement(s)
                row.update({"tcr_height": p.height, "shift_u": p.shift_u,
                            "shift_w": p.shift_w, "offset": p.offset})
            except Exception:
                pass
            pose_rows.append(row)
        print(f"  pose {min(start + CHUNK, len(paths))}/{len(paths)} {time.time() - t0:.0f}s",
              flush=True)

    df = (pl.DataFrame(rec, infer_schema_length=None)
          .join(pl.DataFrame(pose_rows, infer_schema_length=None), on="complex.id", how="inner"))
    keep = ["complex.id", *Q_FEATURES_GEOM, *[c for c in ENERGY if c in df.columns],
            *POSE_FEATURES, *[c for c in PLACE if c in df.columns]]
    df = df.select(keep)
    before = df.height
    core = [*Q_FEATURES_GEOM, *POSE_FEATURES]
    df = df.filter(pl.all_horizontal([pl.col(c).is_finite() for c in core]))
    print(f"kept {df.height} of {len(paths)} (dropped {before - df.height} non-finite)")
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    df.write_csv(args.out)
    print(f"wrote {args.out} ({df.height} x {df.width})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
