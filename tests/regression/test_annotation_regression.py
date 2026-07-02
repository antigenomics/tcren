"""TCR annotation parity: arda-driven chain typing + region markup vs the oracle.

arda's IMGT region boundaries are not byte-identical to the legacy mir markup, so the
exact-match assertions use single-copy structures where the two agree on the contacting
residues; the benchmark sweep measures aggregate concordance and asserts that every
oracle TCR↔peptide contact is reproduced (region labels allowed to differ slightly).
"""

from __future__ import annotations

import os
from pathlib import Path

import polars as pl
import pytest

from tcren import ContactMap, parse_structure
from tcren.annotation import classify_chains

arda = pytest.importorskip("arda")

REPO = Path(__file__).resolve().parents[2]
PDB_DIR = REPO / "tests" / "assets" / "pdb"
CONTACT_MAPS = REPO / "legacy" / "data" / "contact_maps_PDB.csv"
SUMMARY = REPO / "legacy" / "data" / "summary_PDB_structures.csv"

# invokes arda / mmseqs per structure, and compares against the legacy mir oracle CSVs; skip
# (do not fail) when those un-fetched reference files are absent from the checkout, mirroring the
# arda importorskip above.
pytestmark = [
    pytest.mark.slow,
    pytest.mark.skipif(
        not (SUMMARY.exists() and CONTACT_MAPS.exists()),
        reason="legacy oracle data (legacy/data/summary_PDB_structures.csv, contact_maps_PDB.csv) "
        "not present in this checkout",
    ),
]

_FULL_KEYS = [
    "chain.type.from",
    "region.type.from",
    "residue.index.from",
    "residue.index.to",
    "residue.aa.from",
    "residue.aa.to",
]
_CONTACT_KEYS = ["residue.index.from", "residue.index.to", "residue.aa.from", "residue.aa.to"]


def _organism(pdb_id: str, summary: pl.DataFrame) -> str:
    row = summary.filter(pl.col("pdb.id") == pdb_id)
    if row.height == 0:
        return "human"  # not in summary; arda autodetect will retry the other organism
    return "mouse" if row["complex.species"][0] == "Mouse" else "human"


def _annotated_contact_map(pdb_id: str, summary: pl.DataFrame) -> ContactMap:
    s = parse_structure(PDB_DIR / f"{pdb_id}.pdb")
    classify_chains(s, organism=_organism(pdb_id, summary))
    return ContactMap.from_structure(s)


# Single-copy complexes where arda and mir agree on the contacting residues.
@pytest.mark.parametrize("pdb_id", ["1ao7", "1bd2", "5m01"])
def test_full_contactmap_matches_oracle(pdb_id):
    summary = pl.read_csv(SUMMARY)
    oracle = pl.read_csv(CONTACT_MAPS)
    cm = _annotated_contact_map(pdb_id, summary)
    got = set(map(tuple, cm.tcr_peptide().select(_FULL_KEYS).unique().rows()))
    want = set(
        map(tuple, oracle.filter(pl.col("pdb.id") == pdb_id).select(_FULL_KEYS).unique().rows())
    )
    assert got == want


@pytest.mark.parametrize("pdb_id", ["1ao7", "1bd2", "5m01"])
def test_cdr3_matches_summary(pdb_id):
    """arda's cdr3_aa equals the summary CDR3 minus its conserved C…F/W anchors."""
    summary = pl.read_csv(SUMMARY)
    row = summary.filter(pl.col("pdb.id") == pdb_id)
    organism = _organism(pdb_id, summary)
    s = parse_structure(PDB_DIR / f"{pdb_id}.pdb")
    classify_chains(s, organism=organism)
    for chain in s.chains:
        if chain.chain_type in ("TRA", "TRB"):
            rec = arda.annotate_sequences(
                [(chain.chain_id, chain.sequence())], seqtype="aa", organism=organism
            )[0]
            col = "cdr3a" if chain.chain_type == "TRA" else "cdr3b"
            summary_cdr3 = row[col][0]
            assert rec["cdr3_aa"] == summary_cdr3[1:-1]


def test_chain_classification_5m01():
    summary = pl.read_csv(SUMMARY)
    s = parse_structure(PDB_DIR / "5m01.pdb")
    classify_chains(s, organism="mouse")
    types = {c.chain_id: c.chain_type for c in s.chains}
    assert types == {"G": "TRA", "H": "TRB", "P": "PEPTIDE", "A": "MHC", "B": "MHC"}


@pytest.mark.skipif(not os.getenv("RUN_BENCHMARK"), reason="set RUN_BENCHMARK=1 to run")
def test_annotation_concordance_sweep():
    summary = pl.read_csv(SUMMARY)
    oracle = pl.read_csv(CONTACT_MAPS)
    pdbs = oracle["pdb.id"].unique().to_list()
    n = 0
    full_exact = 0
    missing_contacts = []  # oracle contacts not reproduced (a real failure)
    region_only_diff = 0
    for pdb_id in pdbs:
        if not (PDB_DIR / f"{pdb_id}.pdb").exists():
            continue
        n += 1
        cm = _annotated_contact_map(pdb_id, summary)
        got = cm.tcr_peptide().select(_FULL_KEYS).unique()
        want = oracle.filter(pl.col("pdb.id") == pdb_id).select(_FULL_KEYS).unique()
        g_full = set(map(tuple, got.rows()))
        o_full = set(map(tuple, want.rows()))
        g_con = set(map(tuple, got.select(_CONTACT_KEYS).unique().rows()))
        o_con = set(map(tuple, want.select(_CONTACT_KEYS).unique().rows()))
        if g_full == o_full:
            full_exact += 1
        elif g_con == o_con:
            region_only_diff += 1
        if not o_con <= g_con:
            missing_contacts.append((pdb_id, sorted(o_con - g_con)[:3]))

    print(
        f"\nconcordance over {n} structures: full-exact={full_exact} "
        f"region-only-diff={region_only_diff} not-reproduced={len(missing_contacts)}"
    )
    # Every oracle contact must be reproduced; multi-copy structures may add extras.
    assert not missing_contacts, f"{len(missing_contacts)} structures missing contacts: {missing_contacts[:10]}"
    assert full_exact / n >= 0.85
