"""Unit tests for C-gene trimming on import."""

from __future__ import annotations

from pathlib import Path

from tcren.annotation.cgene import constant_span
from tcren.structure import import_structure, parse_structure

REPO = Path(__file__).resolve().parents[2]
FULL = REPO / "tests" / "assets" / "cgene" / "1ao7_full.pdb"
VONLY = REPO / "data" / "PDB_structures" / "1ao7.pdb"


def test_constant_span_found_in_full_beta_chain():
    beta = parse_structure(FULL).chain("E")  # β chain with Cβ constant domain
    span = constant_span(beta.sequence())
    assert span is not None
    start, end = span
    assert 100 < start < len(beta.residues)  # constant is C-terminal


def test_constant_span_absent_in_variable_only():
    beta = parse_structure(FULL).chain("E")
    v_only = beta.sequence()[:110]
    assert constant_span(v_only) is None


def test_import_trims_constant_and_preserves_v_seq_index():
    full = parse_structure(FULL)
    imported = import_structure(FULL)  # trim_c_gene=True by default
    assert imported.cell_type == "ab"
    beta_full = full.chain("E")
    beta_trim = imported.chain("E")
    assert len(beta_trim.residues) < len(beta_full.residues)
    # V-domain residues keep their seq_index and sequence prefix.
    assert beta_trim.residues[0].seq_index == beta_full.residues[0].seq_index
    assert beta_full.sequence().startswith(beta_trim.sequence())


def test_keep_c_gene_retains_full_chain():
    full = parse_structure(FULL)
    kept = import_structure(FULL, keep_c_gene=True)
    assert len(kept.chain("E").residues) == len(full.chain("E").residues)


def test_trim_is_noop_on_variable_only_structure():
    full = parse_structure(VONLY)
    imported = import_structure(VONLY)
    assert [len(c.residues) for c in imported.chains] == [
        len(c.residues) for c in full.chains
    ]
