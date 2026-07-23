"""Complete 'orphan' natives — files but no metadata AND no recoverable VDJdb hash.

2026-07-19

25_backfill_metadata.py handles orphans whose hash is recoverable from VDJdb. This handles the
rest by annotating their *shipped aligned* PDB directly: those that are a complete αβ TCR:pMHC get
a tcren-annotated metadata row + contacts + skeleton maps (pdb/coords already exist, so are not
regenerated). Structures that are pMHC-only / MHC-only are reported as not-a-TCR:pMHC and left out.

    python 26_complete_orphans.py
"""

from __future__ import annotations

import csv
import gzip
import os
import re
import shutil
import sys
import tarfile
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import reformat as rf  # noqa: E402
from importlib import import_module

_bf = import_module("25_backfill_metadata")

HERE = Path(__file__).parent
DATA_PATCH = HERE / "data_patch"
HF_CLONE = Path(os.environ.get("TCREN_HF_CLONE", os.path.expanduser("~/hf/vdjdb_structure_models")))


def _extract_aligned(pdbid: str, dest: Path) -> Path | None:
    with tarfile.open(HF_CLONE / "data" / "pdb_files_native.tgz") as t:
        for m in t.getmembers():
            b = os.path.basename(m.name)
            if not b.endswith(".pdb"):
                continue
            if _bf._strip_aligned(b[:-4]) == pdbid:
                m.name = b
                t.extract(m, dest)
                return dest / b
    return None


def main() -> None:
    file_hash = _bf._file_hashes()
    missing = sorted(_bf._present_natives() - _bf._meta_native_ids())
    with gzip.open(HF_CLONE / "vdjdb_structures_metadata.tsv.gz", "rt") as fh:
        by_meta = _bf._index_by_hash(list(csv.DictReader(fh, delimiter="\t")), "tcr_pmhc_hash")
    with open(rf.ANNOTATED_TSV) as fh:
        by_ann = _bf._index_by_hash(list(csv.DictReader(fh, delimiter="\t")), "TCR_hash")
    native2026 = {f[:-4].lower() for f in os.listdir(rf.NATIVE2026_DIR) if f.endswith(".pdb")}
    orphans = [p for p in missing
               if not (file_hash.get(p) and (by_meta.get(file_hash[p]) or by_ann.get(file_hash[p])))]

    rows, not_tcr = [], []
    for pid in orphans:
        if pid not in native2026:  # native ⟺ in Native2026
            not_tcr.append((pid, "not in Native2026 (cannot be native)"))
            continue
        tmp = Path(tempfile.mkdtemp())
        try:
            apdb = _extract_aligned(pid, tmp)
            if apdb is None:
                not_tcr.append((pid, "no aligned PDB in dataset"))
                continue
            # Metadata-only: pdb/coords already ship. num_contacts from the asymmetric unit is
            # unreliable for these (the complex forms under crystal symmetry) → leave it blank.
            row = rf.process_prealigned(pid, apdb, tmp, make_map=False, do_angles=True)
            if str(row["num_contacts"]) in ("0", "0.0"):
                row["num_contacts"] = ""
            rows.append(row)
            print(f"  + {pid}  metadata-only  hash={row['tcr_pmhc_hash'][:12]}… "
                  f"epitope={row['antigen.epitope']} cdr3b={row['cdr3.beta']}")
        except Exception as e:  # noqa: BLE001
            not_tcr.append((pid, str(e)))
            print(f"  - {pid}  {e}")
        finally:
            shutil.rmtree(tmp, ignore_errors=True)

    if rows:
        import pandas as pd
        df = pd.DataFrame(rows)[rf.PUBLISHED_COLUMNS].astype(str)
        out = DATA_PATCH / "patch_metadata.tsv"
        if out.exists():
            prev = pd.read_csv(out, sep="\t", dtype=str, keep_default_na=False)
            df = pd.concat([prev, df]).drop_duplicates("meta.structure.id", keep="last")
        df.to_csv(out, sep="\t", index=False)

    print(f"\ncompleted {len(rows)} orphan(s); {len(not_tcr)} not a TCR:pMHC (excluded):")
    for pid, msg in not_tcr:
        print(f"   - {pid}: {msg}")


if __name__ == "__main__":
    main()
