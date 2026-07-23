"""Populate scanning_angle / pitch_angle in data_patch/patch_metadata.tsv via STCRpy.

2026-07-19

A separate step because STCRpy needs a working ANARCI (`anarci-mhc` + built germline models, see
README). It reads the already-written aligned PDBs in data_patch/ — no tcren required — and fills
the two angle columns in place. Structures whose angles cannot be computed keep blank angles.

    python 21_angles.py
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

os.environ["PATH"] = "/opt/homebrew/bin:" + os.environ.get("PATH", "")  # hmmscan / muscle for ANARCI
sys.path.insert(0, str(Path(__file__).parent))
import reformat as rf  # noqa: E402

DATA_PATCH = Path(__file__).parent / "data_patch"


def main() -> None:
    import pandas as pd
    meta = DATA_PATCH / "patch_metadata.tsv"
    df = pd.read_csv(meta, sep="\t", dtype=str, keep_default_na=False)
    pdb_dir = DATA_PATCH / "pdb_files_native"

    ok = 0
    for i, row in df.iterrows():
        pid = row["meta.structure.id"].lower()
        apdb = pdb_dir / f"aligned_{pid}.pdb"
        if not apdb.exists():
            print(f"  {pid}: no aligned PDB")
            continue
        a = rf.angles(apdb)
        if a:
            df.at[i, "scanning_angle"] = str(a["scanning_angle"])
            df.at[i, "pitch_angle"] = str(a["pitch_angle"])
            ok += 1
            print(f"  {pid}: scanning={a['scanning_angle']} pitch={a['pitch_angle']}")
        else:
            print(f"  {pid}: angles unavailable")

    df.to_csv(meta, sep="\t", index=False)
    print(f"\n{ok}/{len(df)} rows got angles -> {meta}")


if __name__ == "__main__":
    main()
