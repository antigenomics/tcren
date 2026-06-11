"""Unit tests for arda region-coordinate projection onto structure residues."""

from __future__ import annotations

from pathlib import Path

import pytest

from tcren.annotation.arda_adapter import annotate_chain
from tcren.structure import parse_structure

arda = pytest.importorskip("arda")

REPO = Path(__file__).resolve().parents[2]
PDB_DIR = REPO / "data" / "PDB_structures"


def test_region_projection_aligns_with_arda_coords():
    s = parse_structure(PDB_DIR / "5m01.pdb")
    chain = s.chain("G")
    seq = chain.sequence()
    record = annotate_chain(chain, organism="mouse")

    assert chain.chain_type == "TRA"
    region_types = [r.region_type for r in chain.regions]
    assert region_types == ["FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "FR4"]

    # Each region's residue slice must reproduce arda's 1-based inclusive coordinates.
    arda_to_mir = {
        "fwr1": "FR1", "cdr1": "CDR1", "fwr2": "FR2", "cdr2": "CDR2",
        "fwr3": "FR3", "cdr3": "CDR3", "fwr4": "FR4",
    }
    by_type = {r.region_type: r for r in chain.regions}
    for arda_name, mir_name in arda_to_mir.items():
        start, end = record[f"{arda_name}_start"], record[f"{arda_name}_end"]
        region = by_type[mir_name]
        assert region.sequence == seq[start - 1 : end]
        assert region.start_seq_index == chain.residues[start - 1].seq_index
        assert region.end_seq_index == chain.residues[end - 1].seq_index


def test_non_tcr_chain_left_unannotated():
    s = parse_structure(PDB_DIR / "5m01.pdb")
    peptide = s.chain("P")
    record = annotate_chain(peptide, organism="mouse")
    assert peptide.chain_type is None  # arda does not call a TCR locus on the peptide
    assert peptide.regions == []
    if record is not None:
        assert record.get("locus") not in ("TRA", "TRB")
