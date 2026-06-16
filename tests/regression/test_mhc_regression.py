"""MHC mapping parity: class / role / locus vs the curated dataset annotations.

Requires arda and a built MHC reference (``database/mhc/alleles.aa.fasta``); skipped
otherwise. Build the reference with ``tcren build-mhc-ref`` or ``tcren.mhc.build()``.
"""

from __future__ import annotations

import os
from pathlib import Path

import polars as pl
import pytest

pytest.importorskip("arda")

pytestmark = pytest.mark.slow  # invokes arda / mmseqs per structure

from tcren import parse_structure
from tcren.annotation import classify_chains
from tcren.mhc import map_mhc, reference

REPO = Path(__file__).resolve().parents[2]
PDB_DIR = REPO / "tests" / "assets" / "pdb"
SUMMARY = REPO / "legacy" / "data" / "summary_PDB_structures.csv"
PDB_MHC = REPO / "legacy" / "data" / "PDB_MHC_annotation.csv"

_HAVE_REF = (REPO / "database" / "mhc" / "alleles.aa.fasta").exists()
needs_ref = pytest.mark.skipif(not _HAVE_REF, reason="MHC reference not built")


def _call_class(calls) -> str:
    return "MHCII" if any(c.chain_role == "MHCb" for c in calls) else "MHCI"


@needs_ref
@pytest.mark.parametrize(
    "pdb_id,organism,exp_class,exp_locus_token",
    [
        ("1ao7", "human", "MHCI", "A*02"),
        ("5m01", "mouse", "MHCI", "H2"),
        ("4ozg", "human", "MHCII", "DQ"),
    ],
)
def test_mhc_mapping(pdb_id, organism, exp_class, exp_locus_token):
    s = parse_structure(PDB_DIR / f"{pdb_id}.pdb")
    classify_chains(s, organism=organism)
    calls = map_mhc(s)
    assert calls, "no MHC chains mapped"
    assert _call_class(calls) == exp_class
    # class I: exactly one MHCa + one B2M; class II: one MHCa + one MHCb
    roles = sorted(c.chain_role for c in calls)
    if exp_class == "MHCI":
        assert roles == ["B2M", "MHCa"]
    else:
        assert roles == ["MHCa", "MHCb"]
    mhca = next(c for c in calls if c.chain_role == "MHCa")
    assert exp_locus_token in mhca.allele


@needs_ref
@pytest.mark.skipif(not os.getenv("RUN_BENCHMARK"), reason="set RUN_BENCHMARK=1 to run")
def test_mhc_class_concordance_sweep():
    summary = pl.read_csv(SUMMARY)
    sp = {
        r["pdb.id"]: ("mouse" if r["complex.species"] == "Mouse" else "human")
        for r in summary.iter_rows(named=True)
    }
    wrong = []
    n = 0
    for r in summary.iter_rows(named=True):
        pid = r["pdb.id"]
        fp = PDB_DIR / f"{pid}.pdb"
        if not fp.exists():
            continue
        s = parse_structure(fp)
        classify_chains(s, organism=sp[pid])
        calls = map_mhc(s)
        if not calls:
            wrong.append((pid, "no calls"))
            continue
        n += 1
        if _call_class(calls) != r["mhc.class"]:
            wrong.append((pid, f"{_call_class(calls)} != {r['mhc.class']}"))
    assert len(wrong) / max(n, 1) < 0.05, f"class mismatches: {wrong[:15]}"
