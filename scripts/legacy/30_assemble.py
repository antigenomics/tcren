"""Stage data_patch/ into the local HF clone ~/hf/vdjdb_structure_models.

2026-07-19

Repacks the 5 affected ``data/*.tgz`` archives (existing members + new native files), appends the
patch metadata rows to ``vdjdb_structures_metadata.tsv.gz``, and copies the legacy PCA model into
the dataset so the reproduction is self-contained. Dry-run by default.

    python 30_assemble.py            # report deltas only, no changes
    python 30_assemble.py --apply    # repack + append + copy asset + git commit (NO push)
    python 30_assemble.py --apply --push   # also `git push` (publishing is the user's call)
"""

from __future__ import annotations

import argparse
import gzip
import os
import shutil
import subprocess
import sys
import tarfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import reformat as rf  # noqa: E402

HERE = Path(__file__).parent
DATA_PATCH = HERE / "data_patch"
HF_CLONE = Path(os.environ.get("TCREN_HF_CLONE", os.path.expanduser("~/hf/vdjdb_structure_models")))
METADATA = HF_CLONE / "vdjdb_structures_metadata.tsv.gz"
# patch subdir -> archive filename under data/
ARCHIVES = {
    "pdb_files_native": "pdb_files_native.tgz",
    "coordinates_aa": "coordinates_aa.tgz",
    "contacts_aa": "contacts_aa.tgz",
    "complementarity_maps": "complementarity_maps.tgz",
    "complementarity_maps_simplified": "complementarity_maps_simplified.tgz",
}


def _strip_aligned(stem: str) -> str:
    while stem.startswith("aligned_"):
        stem = stem[len("aligned_"):]
    return stem.lower()


def _native_file_pdbid(name: str):
    """pdbid a native-archive member belongs to (pdb / coords / contacts), else None (e.g. maps)."""
    if name.endswith(".pdb"):
        return _strip_aligned(name[:-4])
    if name.endswith("_aa_coordinates.tsv"):
        return name[: -len("_aa_coordinates.tsv")].lower()
    if name.endswith("_aa_contacts.tsv"):
        parts = name[: -len("_aa_contacts.tsv")].split("_")
        return parts[-1].lower() if len(parts) >= 2 else None
    return None


def repack(subdir: str, archive: str, apply: bool, remove_pdbids: set) -> tuple[int, int, int]:
    tgz = HF_CLONE / "data" / archive
    with tarfile.open(tgz) as t:
        names = [os.path.basename(n) for n in t.getnames()]
    existing = set(names)
    patch_dir = DATA_PATCH / subdir
    patch_files = sorted(patch_dir.iterdir()) if patch_dir.exists() else []
    new = [f for f in patch_files if f.name not in existing]
    dup = [f for f in patch_files if f.name in existing]
    drop = {n for n in names if _native_file_pdbid(n) in remove_pdbids}  # bogus-native files
    if apply and (new or drop):
        tmp = tgz.with_name(tgz.name + ".new")
        with tarfile.open(tgz) as told, tarfile.open(tmp, "w:gz") as tnew:
            for m in told.getmembers():
                if os.path.basename(m.name) in drop:
                    continue
                tnew.addfile(m, told.extractfile(m) if m.isfile() else None)
            for f in new:
                tnew.add(f, arcname=f.name)
        shutil.move(str(tmp), str(tgz))
    return len(new), len(dup), len(drop)


def append_metadata(apply: bool, native2026: set) -> tuple[int, int, int, int]:
    import pandas as pd
    patch = pd.read_csv(DATA_PATCH / "patch_metadata.tsv", sep="\t", dtype=str, keep_default_na=False)
    with gzip.open(METADATA, "rt") as fh:
        existing = pd.read_csv(fh, sep="\t", dtype=str, keep_default_na=False)
    if list(patch.columns) != list(existing.columns):
        raise SystemExit(f"column mismatch:\n patch={list(patch.columns)}\n meta ={list(existing.columns)}")
    combined = pd.concat([existing, patch], ignore_index=True)
    # Enforce is_native ⟺ in Native2026: drop bogus native rows (structure not in Native2026).
    bogus = (combined["is_native"] == "True") & \
            (~combined["meta.structure.id"].str.lower().isin(native2026))
    n_bogus = int(bogus.sum())
    combined = combined[~bogus]
    if apply:
        with gzip.open(METADATA, "wt") as fh:
            combined.to_csv(fh, sep="\t", index=False)
    return len(existing), len(patch), len(combined), n_bogus


def copy_pca_asset(apply: bool) -> str:
    dest = HF_CLONE / "legacy_assets" / rf.PCA_PATH.name
    if apply:
        dest.parent.mkdir(exist_ok=True)
        shutil.copy2(rf.PCA_PATH, dest)
    return str(dest.relative_to(HF_CLONE))


def _present_native_pdbids() -> set:
    with tarfile.open(HF_CLONE / "data" / "pdb_files_native.tgz") as t:
        return {_strip_aligned(os.path.basename(n)[:-4]) for n in t.getnames() if n.endswith(".pdb")}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--apply", action="store_true")
    ap.add_argument("--push", action="store_true")
    args = ap.parse_args()
    apply = args.apply

    if apply:
        subprocess.run(["git", "-C", str(HF_CLONE), "pull", "--ff-only"], check=True)

    native2026 = {f[:-4].lower() for f in os.listdir(rf.NATIVE2026_DIR) if f.endswith(".pdb")}
    remove_pdbids = _present_native_pdbids() - native2026  # present natives not in Native2026 = bogus

    print(f"{'APPLY' if apply else 'DRY-RUN'} — clone: {HF_CLONE}")
    print(f"enforce is_native ⟺ Native2026: {len(remove_pdbids)} bogus native(s) removed "
          f"-> {sorted(remove_pdbids)}\n")
    total_new = total_drop = 0
    for subdir, archive in ARCHIVES.items():
        n, d, r = repack(subdir, archive, apply, remove_pdbids)
        total_new += n
        total_drop += r
        print(f"  {archive:38} +{n:3} new  -{r:2} bogus  ({d} already present)")
    e, p, c, nb = append_metadata(apply, native2026)
    print(f"\n  metadata rows: {e} existing + {p} patch - {nb} bogus native = {c}")
    asset = copy_pca_asset(apply)
    print(f"  PCA asset -> {asset}")

    if not apply:
        print("\n(dry-run: nothing changed — re-run with --apply)")
        return

    subprocess.run(["git", "-C", str(HF_CLONE), "add", "-A"], check=True)
    msg = (f"Add {p} Native2026 structures (+{total_new} files) + PCA asset; "
           f"remove {nb} bogus native rows / {total_drop} files (not in Native2026)")
    subprocess.run(["git", "-C", str(HF_CLONE), "commit", "-m", msg], check=True)
    print(f"\ncommitted: {msg}")
    if args.push:
        subprocess.run(["git", "-C", str(HF_CLONE), "push"], check=True)
        print("pushed to HuggingFace.")
    else:
        print("NOT pushed — review, then `git -C ~/hf/vdjdb_structure_models push` (or --push).")


if __name__ == "__main__":
    main()
