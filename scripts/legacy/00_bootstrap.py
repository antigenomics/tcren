"""Bootstrap the reformat run: diff Native2026 against the shipped dataset, build the worklist.

2026-07-19

Reads the local HF clone ``~/hf/vdjdb_structure_models`` (already cloned) and the Native2026
source, and writes into ``scripts/legacy/data_dump/`` (gitignored):

* ``worklist.tsv``            — one row per *missing* Native2026 pdbid: ``pdbid, status, vdjdb_hash``
                                (status = ``joinable`` if it maps to a VDJdb record, else ``tcren``).
* ``present_map_hashes.txt``  — hashes that already have a complementarity map (skip duplicates).

No network / GitLab access — everything is local.
"""

from __future__ import annotations

import os
import sys
import tarfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import reformat as rf  # noqa: E402

HF_CLONE = Path(os.environ.get("TCREN_HF_CLONE", os.path.expanduser("~/hf/vdjdb_structure_models")))
DATA_DUMP = Path(__file__).parent / "data_dump"


def _pdbids_in_native_archive(tgz: Path) -> set:
    ids = set()
    with tarfile.open(tgz) as t:
        for n in t.getnames():
            b = os.path.basename(n)
            if b.endswith(".pdb"):
                stem = b[:-4]
                while stem.startswith("aligned_"):
                    stem = stem[len("aligned_"):]
                ids.add(stem.lower())
    return ids


def _map_hashes(tgz: Path) -> set:
    with tarfile.open(tgz) as t:
        return {os.path.basename(n)[:-4] for n in t.getnames() if n.endswith(".svg")}


def _native2026_ids() -> set:
    ids = set()
    for f in os.listdir(rf.NATIVE2026_DIR):
        for ext in (".pdb.gz", ".cif.gz", ".pdb", ".cif"):
            if f.endswith(ext):
                ids.add(f[: -len(ext)].lower())
                break
    return ids


def main() -> None:
    DATA_DUMP.mkdir(exist_ok=True)
    data = HF_CLONE / "data"

    present = _pdbids_in_native_archive(data / "pdb_files_native.tgz")
    map_hashes = _map_hashes(data / "complementarity_maps.tgz")
    native2026 = _native2026_ids()
    missing = sorted(native2026 - present)

    print(f"present natives:   {len(present)}")
    print(f"Native2026 ids:    {len(native2026)}  (already present ∩ Native2026: {len(native2026 & present)})")
    print(f"missing to add:    {len(missing)}")

    vdjdb = rf.load_vdjdb_index()
    print(f"VDJdb structure ids indexed: {len(vdjdb)}")

    joinable = unjoinable = 0
    worklist = DATA_DUMP / "worklist.tsv"
    with open(worklist, "w") as fh:
        fh.write("pdbid\tstatus\tvdjdb_hash\n")
        for pid in missing:
            rec = vdjdb.get(pid)
            if rec is not None:
                h = (rec.get("TCR_hash") or "").strip()
                fh.write(f"{pid}\tjoinable\t{h}\n")
                joinable += 1
            else:
                fh.write(f"{pid}\ttcren\t\n")
                unjoinable += 1

    (DATA_DUMP / "present_map_hashes.txt").write_text("\n".join(sorted(map_hashes)))
    print(f"\nworklist: {joinable} joinable (VDJdb) + {unjoinable} tcren-annotated  -> {worklist}")
    print(f"present map hashes: {len(map_hashes)}  -> data_dump/present_map_hashes.txt")

    # Show one joinable VDJdb row so the join is inspectable.
    for pid in missing:
        rec = vdjdb.get(pid)
        if rec:
            print(f"\nsample joinable row  pdbid={pid}:")
            for k in ("cdr3.alpha", "v.alpha", "j.alpha", "cdr3.beta", "v.beta", "j.beta",
                      "species", "mhc.a", "mhc.b", "mhc.class", "antigen.epitope",
                      "meta.structure.id", "vdjdb.score", "TCR_hash"):
                print(f"   {k:18} {rec.get(k, '')!r}")
            break


if __name__ == "__main__":
    main()
