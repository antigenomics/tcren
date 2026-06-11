"""Unit tests for the substitution scoring engine."""

from __future__ import annotations

import polars as pl
import pytest

from tcren.contactmap import ContactMap
from tcren.potential import Potential
from tcren.scoring import score_peptides


def _toy_potential() -> Potential:
    # 3-letter alphabet with distinct, hand-checkable values.
    rows = []
    vals = {("A", "A"): 1.0, ("A", "K"): -2.0, ("L", "A"): 0.5, ("L", "K"): 3.0,
            ("A", "G"): 0.1, ("L", "G"): 0.2}
    for (fr, to), v in vals.items():
        rows.append({"residue.aa.from": fr, "residue.aa.to": to, "value": v})
    return Potential(name="toy", matrix=pl.DataFrame(rows), alphabet=("A", "L", "K", "G"))


def _toy_contact_map() -> ContactMap:
    # Two TCR↔peptide contacts: TCR 'A' at peptide pos 0, TCR 'L' at peptide pos 2.
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


def test_hand_computed_score():
    cm = _toy_contact_map()
    pot = _toy_potential()
    # Candidate "AGK": pos0='A' (TCR 'A' -> value A,A = 1.0), pos2='K' (TCR 'L' -> L,K = 3.0)
    res = score_peptides(cm, ["AGK"], pot)
    assert res.height == 1
    assert res["score"][0] == pytest.approx(1.0 + 3.0)


def test_length_filter():
    cm = _toy_contact_map()
    pot = _toy_potential()
    res = score_peptides(cm, ["AGK", "AG", "AGKK"], pot, require_same_length=True)
    assert res["peptide"].to_list() == ["AGK"]  # only the length-3 candidate kept


def test_missing_pair_contributes_zero():
    cm = _toy_contact_map()
    pot = _toy_potential()
    # 'G' on the substituted (peptide) side at pos0 pairs with TCR 'A' -> (A,G)=0.1;
    # pos2 'A' with TCR 'L' -> (L,A)=0.5.  Total 0.6.
    res = score_peptides(cm, ["GGA"], pot)
    assert res["score"][0] == pytest.approx(0.1 + 0.5)


def test_output_schema_and_sorting():
    cm = _toy_contact_map()
    pot = _toy_potential()
    res = score_peptides(cm, ["AGK", "GGA"], pot)
    assert res.columns == ["complex.id", "peptide", "potential", "score"]
    assert res["score"].to_list() == sorted(res["score"].to_list())
    assert res["potential"].unique().to_list() == ["toy"]
