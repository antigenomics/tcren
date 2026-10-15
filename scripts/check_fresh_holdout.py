"""Compare a full fresh reproduction's hold-out with the shipped frozen object.

The reproduction driver owns freshness: it rebuilds the harmonization parquet from raw
structures before invoking this gate. Refit only into a temporary directory, never the library.
The original full-scale contract is exact array equality; this gate does not relax it.
"""
from __future__ import annotations

import argparse
import json
import tempfile
from pathlib import Path

import numpy as np
import polars as pl

from tcren.score import MODEL_FILE, holdout_manifest
from tcren.score.fit import fit_holdout


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--features", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()
    manifest = holdout_manifest()
    report = {"n_manifest": manifest.height, "contract": "exact frozen numerical identity"}
    try:
        with np.load(MODEL_FILE, allow_pickle=False) as z:
            ref = {k: z[k].copy() for k in z.files}
        meta = json.loads(str(ref["meta"]))
        scan = pl.scan_parquet(args.features)
        key = next(k for k in ("id", "pdb.id", "complex.id") if k in scan.collect_schema())
        ids = scan.select(key).unique().collect()[key]
        missing = manifest.filter(~pl.col("pdb.id").is_in(ids.implode()))
        report["n_missing"] = missing.height
        if missing.height:
            raise ValueError(f"{missing.height} missing hold-out structures; full reproduction required")
        features = (scan.filter(pl.col(key).is_in(manifest["pdb.id"].implode()))
                    .select(key, *meta["descriptors"]).unique().collect()
                    .rename({key: "pdb.id"}).sort("pdb.id"))
        if features["pdb.id"].n_unique() != features.height:
            raise ValueError("conflicting descriptor rows for the same hold-out structure")
        report["n_features"] = features.height
        with tempfile.TemporaryDirectory(prefix="tcren-fresh-fit-") as td:
            out = Path(td) / "model.npz"
            got_meta = fit_holdout(features, manifest, out=out)
            with np.load(out, allow_pickle=False) as z:
                got = {k: z[k].copy() for k in z.files if k != "meta"}
        report["arrays"] = {k: {"equal": bool(np.array_equal(got[k], ref[k])),
                                  "changed_cells": int(np.count_nonzero(got[k] != ref[k])),
                                  "max_abs_delta": float(np.max(np.abs(got[k] - ref[k])))}
                             for k in got}
        fields = ("descriptors", "coordinates", "receptor_coordinates", "prior", "alpha",
                  "n_pos", "n_neg", "n_epitopes", "conf_mu", "conf_var", "lam", "loc", "scale")
        report["metadata_equal"] = {k: got_meta[k] == meta[k] for k in fields}
        report["n_pos"] = got_meta["n_pos"]
        report["n_neg"] = got_meta["n_neg"]
        report["n_incomplete"] = got_meta["n_incomplete"]
        if (not all(d["equal"] for d in report["arrays"].values())
                or not all(report["metadata_equal"].values())):
            raise ValueError("fresh full hold-out differs from frozen parameters; preserve baseline and review")
        report["passed"] = True
    except (ValueError, KeyError, StopIteration) as exc:
        report.update(passed=False, error=str(exc))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2) + "\n")
    if not report["passed"]:
        p.exit(1, report["error"] + "\n")
    print(f"full hold-out: {manifest.height} structures; frozen numerical identity verified")


if __name__ == "__main__":
    main()
