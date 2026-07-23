"""Validate the reformat pipeline against natives already in the dataset.

2026-07-19

For each sampled present native, runs the pipeline and compares to the shipped entry:
  * hash        — tcren-annotated vs shipped (expected to differ when arda's allele/MHC
                  nomenclature diverges from VDJdb; the joinable path avoids this).
  * num_contacts — tcren tcr_peptide[all] vs the shipped metadata value.
  * contacts     — Jaccard overlap of the CDR3–peptide/CDR3–CDR3 residue-pair sets.

Coordinates are NOT expected to match byte-for-byte: tcren's canonical frame differs from the
legacy PyMOL frame by a rigid transform (the maps/coords are self-consistent per structure).

    python 10_validate.py [--ids 1d9k,1ao7,2bnq,...]
"""

from __future__ import annotations

import argparse
import os
import sys
import tarfile
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import reformat as rf  # noqa: E402

HF_CLONE = Path(os.environ.get("TCREN_HF_CLONE", os.path.expanduser("~/hf/vdjdb_structure_models")))
DEFAULT_IDS = ["1d9k", "1ao7", "2bnq", "2gj6", "1mi5"]


def _shipped_meta(pdbid: str) -> dict | None:
    import gzip
    import csv
    with gzip.open(HF_CLONE / "vdjdb_structures_metadata.tsv.gz", "rt") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row.get("is_native") == "True" and (row.get("meta.structure.id") or "").lower() == pdbid:
                return row
    return None


def _shipped_contact_pairs(pdbid: str, hash_: str) -> set:
    name = f"{hash_}_{pdbid}_aa_contacts.tsv"
    tgz = HF_CLONE / "data" / "contacts_aa.tgz"
    with tarfile.open(tgz) as t:
        member = next((m for m in t.getmembers() if os.path.basename(m.name) == name), None)
        if not member:
            return set()
        lines = t.extractfile(member).read().decode().splitlines()[1:]
    pairs = set()
    for ln in lines:
        f = ln.split("\t")
        if len(f) >= 6:  # unordered: the shipped from/to direction differs by chain order
            pairs.add(frozenset({(f[0], f[2]), (f[5], f[4])}))
    return pairs


def _mine_contact_pairs(path: Path) -> set:
    pairs = set()
    for ln in path.read_text().splitlines()[1:]:
        f = ln.split("\t")
        if len(f) >= 6:
            pairs.add(frozenset({(f[0], f[2]), (f[5], f[4])}))
    return pairs


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ids", type=str, default=",".join(DEFAULT_IDS))
    args = ap.parse_args()

    hdr = f"{'pdbid':6} {'hash':6} {'num_contacts':16} {'contacts Jaccard':18}"
    print(hdr + "\n" + "-" * len(hdr))
    for pdbid in [x.strip().lower() for x in args.ids.split(",")]:
        meta = _shipped_meta(pdbid)
        if not meta:
            print(f"{pdbid:6} (no shipped native row)")
            continue
        ship_hash = meta["tcr_pmhc_hash"]
        ship_nc = meta.get("num_contacts", "")
        tmp = Path(tempfile.mkdtemp())
        try:
            row = rf.process(pdbid, tmp, vdjdb_rec=None, make_map=False, do_angles=False)
        except Exception as e:  # noqa: BLE001
            print(f"{pdbid:6} ERROR {type(e).__name__}: {e}")
            continue
        mine_hash = row["tcr_pmhc_hash"]
        mine_nc = row["num_contacts"]
        contacts = next(tmp.glob("*_aa_contacts.tsv"), None)
        mine_pairs = _mine_contact_pairs(contacts) if contacts else set()
        ship_pairs = _shipped_contact_pairs(pdbid, ship_hash)
        inter = len(mine_pairs & ship_pairs)
        union = len(mine_pairs | ship_pairs) or 1
        jac = f"{inter}/{union} = {inter/union:.2f}"
        hmark = "same" if mine_hash == ship_hash else "DIFF"
        ncmark = "=" if str(mine_nc) == str(float(ship_nc)) else "≠"
        print(f"{pdbid:6} {hmark:6} {str(mine_nc)+' '+ncmark+' '+ship_nc:16} {jac:18}")


if __name__ == "__main__":
    main()
