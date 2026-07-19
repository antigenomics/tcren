"""Backfill metadata rows for native structures that have files but no metadata row.

2026-07-19

Some shipped natives (e.g. the A6 TCR 3D3V) have a PDB/coords/contacts file but no row in
``vdjdb_structures_metadata.tsv.gz``. This adds the missing rows — **no new files**. The row's
``tcr_pmhc_hash`` is taken from the structure's own existing files (the `<hash>_<pdbid>_aa_contacts`
name) so it links to them; the identity fields come from the VDJdb record carrying that hash (the
metadata itself, or the annotated TSV). ``num_contacts`` is recomputed with tcren when the source
PDB is in Native2026. Rows are appended to ``data_patch/patch_metadata.tsv`` (run after 20_process).

Structures whose hash/fields cannot be recovered are skipped and logged.
"""

from __future__ import annotations

import csv
import gzip
import os
import re
import sys
import tarfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import reformat as rf  # noqa: E402

HERE = Path(__file__).parent
DATA_PATCH = HERE / "data_patch"
HF_CLONE = Path(os.environ.get("TCREN_HF_CLONE", os.path.expanduser("~/hf/vdjdb_structure_models")))


def _strip_aligned(stem: str) -> str:
    while stem.startswith("aligned_"):
        stem = stem[len("aligned_"):]
    return stem.lower()


def _present_natives() -> set:
    with tarfile.open(HF_CLONE / "data" / "pdb_files_native.tgz") as t:
        return {_strip_aligned(os.path.basename(n)[:-4]) for n in t.getnames() if n.endswith(".pdb")}


def _meta_native_ids() -> set:
    ids = set()
    with gzip.open(HF_CLONE / "vdjdb_structures_metadata.tsv.gz", "rt") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row.get("is_native") == "True":
                sid = (row.get("meta.structure.id") or "").strip().lower()
                if sid:
                    ids.add(sid)
    return ids


def _file_hashes() -> dict:
    h = {}
    with tarfile.open(HF_CLONE / "data" / "contacts_aa.tgz") as t:
        for n in t.getnames():
            m = re.match(r"^([0-9a-f]{64})_([0-9a-zA-Z]{4})_aa_contacts\.tsv$", os.path.basename(n))
            if m:
                h[m.group(2).lower()] = m.group(1)
    return h


def _index_by_hash(rows, key) -> dict:
    d = {}
    for row in rows:
        h = (row.get(key) or "").strip()
        if h and h not in d:
            d[h] = row
    return d


def main() -> None:
    file_hash = _file_hashes()
    missing = sorted(_present_natives() - _meta_native_ids())
    with gzip.open(HF_CLONE / "vdjdb_structures_metadata.tsv.gz", "rt") as fh:
        by_meta = _index_by_hash(list(csv.DictReader(fh, delimiter="\t")), "tcr_pmhc_hash")
    with open(rf.ANNOTATED_TSV) as fh:
        by_ann = _index_by_hash(list(csv.DictReader(fh, delimiter="\t")), "TCR_hash")
    native2026 = {f[:-4].lower() for f in os.listdir(rf.NATIVE2026_DIR) if f.endswith(".pdb")}

    rows, skipped = [], []
    for pid in missing:
        if pid not in native2026:  # native ⟺ in Native2026; anything else is not a native entry
            skipped.append((pid, "not in Native2026 (cannot be native)"))
            continue
        h = file_hash.get(pid)
        rec = (by_meta.get(h) or by_ann.get(h)) if h else None
        if not h or rec is None:
            skipped.append((pid, "no recoverable hash/fields"))
            continue
        row = {c: "" for c in rf.PUBLISHED_COLUMNS}
        for c in rf.VDJDB_COLUMNS:  # present in both metadata and annotated schemas
            row[c] = rec.get(c, "") or ""
        row["idx"] = rec.get("idx", "") or ""
        row["meta.structure.id"] = pid
        row["tcr_pmhc_hash"] = h
        row["is_native"] = "True"
        if pid in native2026:
            try:
                oriented, _, _ = rf.annotate_and_orient(pid)
                row["num_contacts"] = float(rf.native_num_contacts(oriented))
            except Exception as e:  # noqa: BLE001
                print(f"  ! {pid}: num_contacts recompute failed ({type(e).__name__})")
        rows.append(row)
        print(f"  + {pid}  hash={h[:10]}… src={'metadata' if h in by_meta else 'annotated'} "
              f"num_contacts={row['num_contacts']}")

    import pandas as pd
    out = DATA_PATCH / "patch_metadata.tsv"
    if rows:
        df = pd.DataFrame(rows)[rf.PUBLISHED_COLUMNS].astype(str)
        if out.exists():
            prev = pd.read_csv(out, sep="\t", dtype=str, keep_default_na=False)
            df = pd.concat([prev, df]).drop_duplicates("meta.structure.id", keep="last")
        df.to_csv(out, sep="\t", index=False)

    print(f"\nbackfilled {len(rows)} metadata rows (no new files); skipped {len(skipped)}")
    for pid, msg in skipped:
        print(f"   SKIP {pid}: {msg}")


if __name__ == "__main__":
    main()
