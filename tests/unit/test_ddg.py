"""Unit tests for the fast ΔΔG engine (S4).

Uses the same tiny hand-built contact map / potential as the scoring tests so the
energies are analytically checkable. No external oracle CSV (new method); the
checks are the ΔΔG identities and per-position structure.
"""

from __future__ import annotations

import polars as pl
import pytest

from tcren.contactmap import ContactMap
from tcren.ddg import alanine_scan, ddg, neoantigen_ddg
from tcren.potential import Potential
from tcren.scoring import score_peptides


def _toy_potential() -> Potential:
    vals = {("A", "A"): 1.0, ("A", "K"): -2.0, ("L", "A"): 0.5, ("L", "K"): 3.0,
            ("A", "G"): 0.1, ("L", "G"): 0.2}
    rows = [{"residue.aa.from": fr, "residue.aa.to": to, "value": v}
            for (fr, to), v in vals.items()]
    return Potential(name="toy", matrix=pl.DataFrame(rows), alphabet=("A", "L", "K", "G"))


def _toy_contact_map() -> ContactMap:
    # TCR 'A' contacts peptide pos 0; TCR 'L' contacts peptide pos 2.
    contacts = pl.DataFrame(
        {
            "chain.type.from": ["TRA", "TRB"],
            "chain.type.to": ["PEPTIDE", "PEPTIDE"],
            "residue.aa.from": ["A", "L"],
            "residue.aa.to": ["G", "G"],
            "region.type.from": ["CDR3", "CDR3"],
            "residue.index.from": [10, 20],
            "residue.index.to": [0, 2],
            "region.start.from": [8, 18],
            "region.start.to": [0, 0],
            "pdb.id": ["toy", "toy"],
        }
    )
    return ContactMap(pdb_id="toy", contacts=contacts, peptide_length=3)


def test_ddg_native_vs_native_is_zero():
    cm, pot = _toy_contact_map(), _toy_potential()
    assert ddg(cm, "AGK", "AGK", pot) == 0.0


def test_ddg_matches_independent_two_calls():
    cm, pot = _toy_contact_map(), _toy_potential()
    native, mutant = "AGK", "AGA"
    e_native = float(score_peptides(cm, [native], pot)["score"][0])
    e_mutant = float(score_peptides(cm, [mutant], pot)["score"][0])
    assert ddg(cm, native, mutant, pot) == pytest.approx(e_native - e_mutant)


def test_ddg_sign_and_value():
    cm, pot = _toy_contact_map(), _toy_potential()
    # native "AGK": (A,A)=1.0 + (L,K)=3.0 = 4.0
    # mutant "AGA": (A,A)=1.0 + (L,A)=0.5 = 1.5  -> ddG = 4.0 - 1.5 = 2.5 (destabilising)
    assert ddg(cm, "AGK", "AGA", pot) == pytest.approx(2.5)


def test_alanine_scan_one_row_per_position():
    cm, pot = _toy_contact_map(), _toy_potential()
    native = "AGK"
    scan = alanine_scan(cm, native, pot)
    assert scan.columns == ["pos", "wt_aa", "ddG"]
    assert scan.height == len(native)
    assert scan["pos"].to_list() == [0, 1, 2]
    # wt_aa column reproduces the native peptide.
    assert scan["wt_aa"].to_list() == list(native)
    # pos0 is already 'A' -> mutating to Ala is a no-op -> ddG 0.
    assert scan.filter(pl.col("pos") == 0)["ddG"][0] == pytest.approx(0.0)
    # pos1 'G' has no TCR contact -> no-op -> ddG 0.
    assert scan.filter(pl.col("pos") == 1)["ddG"][0] == pytest.approx(0.0)
    # pos2 'K'->'A': contributes (L,K)=3.0 native vs (L,A)=0.5 -> ddG = 2.5.
    assert scan.filter(pl.col("pos") == 2)["ddG"][0] == pytest.approx(2.5)


def test_alanine_scan_position_matches_independent_ddg():
    cm, pot = _toy_contact_map(), _toy_potential()
    native = "AGK"
    scan = alanine_scan(cm, native, pot)
    for pos in range(len(native)):
        mutant = native[:pos] + "A" + native[pos + 1:]
        expected = ddg(cm, native, mutant, pot)
        assert scan.filter(pl.col("pos") == pos)["ddG"][0] == pytest.approx(expected)


def test_neoantigen_ddg():
    cm, pot = _toy_contact_map(), _toy_potential()
    native = "AGK"
    mutants = ["AGA", "AGK"]
    df = neoantigen_ddg(cm, native, mutants, pot)
    assert df.columns == ["native", "mutant", "ddG"]
    assert df["native"].to_list() == [native, native]
    assert df["mutant"].to_list() == mutants
    assert df.filter(pl.col("mutant") == "AGA")["ddG"][0] == pytest.approx(2.5)
    assert df.filter(pl.col("mutant") == "AGK")["ddG"][0] == pytest.approx(0.0)
