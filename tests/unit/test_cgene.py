"""Unit tests for αβ/γδ classification from the TCR constant region."""

from __future__ import annotations

from pathlib import Path

from tcren.annotation.cgene import cell_type, classify_chain_constant, classify_constants
from tcren.structure import parse_structure

REPO = Path(__file__).resolve().parents[2]
ASSETS = REPO / "tests" / "assets" / "cgene"


def test_alpha_beta_structure():
    s = parse_structure(ASSETS / "1ao7_full.pdb")
    calls = {c.chain_class: c for c in classify_constants(s)}
    assert "beta" in calls
    assert calls["beta"].gene in ("TRBC1", "TRBC2")
    assert cell_type(s) == "ab"


def test_gamma_delta_structure():
    s = parse_structure(ASSETS / "1hxm_gd.pdb")
    classes = {c.chain_class for c in classify_constants(s)}
    assert "delta" in classes and "gamma" in classes
    genes = {c.gene for c in classify_constants(s)}
    assert "TRDC" in genes and any(g.startswith("TRGC") for g in genes)
    assert cell_type(s) == "gd"


def test_variable_only_yields_no_call():
    # The β chain V domain alone (no constant) must not be mistaken for a constant.
    s = parse_structure(ASSETS / "1ao7_full.pdb")
    beta = next(c for c in s.chains if len(c.residues) > 180)
    # Truncate to the first 110 residues (variable domain only).
    v_only_seq = beta.sequence()[:110]
    assert classify_chain_constant(v_only_seq) is None


def test_unknown_when_no_constant_present():
    # A short peptide / V-only chain set yields an unknown cell type.
    s = parse_structure(ASSETS / "1ao7_full.pdb")
    # Keep only the MHC + peptide chains (drop the TCR chains with constants).
    s.chains = [c for c in s.chains if len(c.residues) < 100 or len(c.residues) > 250]
    assert cell_type(s) == "unknown"
