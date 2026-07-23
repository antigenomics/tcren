"""Process the missing Native2026 structures into a patch under ``scripts/legacy/data_patch/``.

2026-07-19

Reads ``data_dump/worklist.tsv`` (from 00_bootstrap.py) and, per pdbid, runs the reformat
pipeline (tcren canonical orientation + legacy coord/plot/contact format) and drops the files
into archive-shaped subdirs, accumulating one metadata row per structure:

    data_patch/pdb_files_native/aligned_<pdbid>.pdb
    data_patch/coordinates_aa/<pdbid>_aa_coordinates.tsv
    data_patch/contacts_aa/<hash>_<pdbid>_aa_contacts.tsv
    data_patch/complementarity_maps[_simplified]/<hash>[_simplified].svg   (only if hash is new)
    data_patch/patch_metadata.tsv

Everything under data_patch/ is gitignored. Usage:
    python 20_process.py [--limit N] [--only PDBID[,PDBID...]] [--no-angles]
"""

from __future__ import annotations

import argparse
import shutil
import sys
import tempfile
import traceback
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import reformat as rf  # noqa: E402

HERE = Path(__file__).parent
DATA_DUMP = HERE / "data_dump"
DATA_PATCH = HERE / "data_patch"

ARCHIVE_DIR = {
    "pdb": "pdb_files_native", "coords": "coordinates_aa", "contacts": "contacts_aa",
    "map": "complementarity_maps", "map_simplified": "complementarity_maps_simplified",
}


def _dest(name: str) -> str:
    if name.endswith("_simplified.svg"):
        return "map_simplified"
    if name.endswith(".svg"):
        return "map"
    if name.endswith("_aa_coordinates.tsv"):
        return "coords"
    if name.endswith("_aa_contacts.tsv"):
        return "contacts"
    if name.startswith("aligned_") and name.endswith(".pdb"):
        return "pdb"
    return ""


def _distribute(tmp: Path, present_map_hashes: set, hash_: str) -> None:
    for f in tmp.iterdir():
        kind = _dest(f.name)
        if not kind:
            continue
        if kind in ("map", "map_simplified") and hash_ in present_map_hashes:
            continue  # a predicted model already provides this hash's map — no duplicate
        out = DATA_PATCH / ARCHIVE_DIR[kind]
        out.mkdir(parents=True, exist_ok=True)
        shutil.copy2(f, out / f.name)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--only", type=str, default=None, help="comma-separated pdbids")
    ap.add_argument("--no-angles", action="store_true")
    args = ap.parse_args()

    worklist = [ln.split("\t") for ln in (DATA_DUMP / "worklist.tsv").read_text().splitlines()[1:]]
    if args.only:
        keep = {x.strip().lower() for x in args.only.split(",")}
        worklist = [w for w in worklist if w[0] in keep]
    if args.limit:
        worklist = worklist[: args.limit]

    vdjdb = rf.load_vdjdb_index()
    pmh = set((DATA_DUMP / "present_map_hashes.txt").read_text().split())
    DATA_PATCH.mkdir(exist_ok=True)

    rows, ok, failed = [], 0, []
    for pdbid, status, _hash in worklist:
        rec = vdjdb.get(pdbid) if status == "joinable" else None
        tmp = Path(tempfile.mkdtemp())
        try:
            row = rf.process(pdbid, tmp, vdjdb_rec=rec, make_map=True, do_angles=not args.no_angles)
            _distribute(tmp, pmh, row["tcr_pmhc_hash"])
            rows.append(row)
            ok += 1
            print(f"  OK   {pdbid:6} [{status:8}] hash={row['tcr_pmhc_hash'][:12]}… "
                  f"num_contacts={row['num_contacts']} angles={'y' if row['scanning_angle']!='' else 'n'}")
        except Exception as e:  # noqa: BLE001 — log & continue; one bad structure must not stop the run
            failed.append((pdbid, f"{type(e).__name__}: {e}"))
            print(f"  FAIL {pdbid:6} [{status:8}] {type(e).__name__}: {e}")
            traceback.print_exc()
        finally:
            shutil.rmtree(tmp, ignore_errors=True)

    if rows:
        import pandas as pd
        df = pd.DataFrame(rows)[rf.PUBLISHED_COLUMNS]
        out = DATA_PATCH / "patch_metadata.tsv"
        # append if it exists (idempotent-ish across chunked runs), else write with header
        if out.exists():  # key on the structure, NOT the hash (distinct crystals can share a hash)
            prev = pd.read_csv(out, sep="\t", dtype=str, keep_default_na=False)
            df = pd.concat([prev, df.astype(str)]).drop_duplicates("meta.structure.id", keep="last")
        df.to_csv(out, sep="\t", index=False)
        print(f"\nwrote {len(df)} metadata rows -> {out}")

    print(f"\ndone: {ok} ok, {len(failed)} failed")
    for pid, msg in failed:
        print(f"   FAIL {pid}: {msg}")


if __name__ == "__main__":
    main()
