"""Unit tests for chemical contact typing (pure logic, no structure/mmseqs)."""

from __future__ import annotations

import polars as pl
import pytest

from tcren.contact_types import TYPES, _classify, classify_contacts, contact_type_counts


@pytest.mark.parametrize(
    "aa_a,aa_b,atom_a,atom_b,dist,expected",
    [
        ("K", "D", "NZ", "OD1", 3.2, "salt_bridge"),     # Lys+ ↔ Asp-
        ("R", "E", "NH1", "OE2", 3.8, "salt_bridge"),    # Arg+ ↔ Glu-
        ("K", "D", "NZ", "OD1", 4.5, "hydrogen_bond"),   # too far for salt, but polar N/O ≤ ... actually >3.5 -> other
        ("Q", "T", "NE2", "OG1", 2.9, "hydrogen_bond"),  # amide N ↔ hydroxyl O
        ("S", "N", "OG", "ND2", 3.1, "hydrogen_bond"),   # two polar
        ("F", "Y", "CZ", "CE1", 4.5, "aromatic"),        # ring–ring aromatic
        ("L", "V", "CD1", "CG1", 4.2, "hydrophobic"),    # apolar C–C
        ("G", "P", "CA", "CB", 4.0, "hydrophobic"),      # P apolar, G apolar? G not in apolar set -> other
        ("Q", "P", "OE1", "CG", 3.0, "other"),           # polar O ↔ carbon, not h-bond (C not N/O)
    ],
)
def test_classify(aa_a, aa_b, atom_a, atom_b, dist, expected):
    # correct the two annotated edge cases inline
    got = _classify(aa_a, aa_b, atom_a, atom_b, dist)
    if aa_a == "K" and dist == 4.5:
        assert got == "other"        # >3.5 Å and >4.0 Å -> neither h-bond nor salt
    elif aa_a == "G":
        assert got == "other"        # Gly not apolar side chain
    else:
        assert got == expected
    assert got in TYPES


def test_classify_contacts_adds_column():
    df = pl.DataFrame({
        "chain.id.from": ["A", "A"], "residue.index.from": [1, 2],
        "chain.id.to": ["C", "C"], "residue.index.to": [3, 3],
        "residue.aa.from": ["K", "L"], "residue.aa.to": ["D", "V"],
        "atom.from": ["NZ", "CD1"], "atom.to": ["OD1", "CG1"], "dist": [3.2, 4.2],
    })
    out = classify_contacts(df)
    assert out["contact.type"].to_list() == ["salt_bridge", "hydrophobic"]


def test_classify_contacts_empty():
    df = pl.DataFrame({c: [] for c in ["residue.aa.from", "residue.aa.to", "atom.from", "atom.to", "dist"]})
    out = classify_contacts(df)
    assert "contact.type" in out.columns and out.height == 0


class _FakeCM:
    def __init__(self, df):
        self._df = df

    def interface(self, name, tcr_regions="all"):
        return self._df


def test_contact_type_counts():
    df = pl.DataFrame({
        "chain.id.from": ["A", "A", "A"], "residue.index.from": [1, 1, 2],
        "chain.id.to": ["C", "C", "C"], "residue.index.to": [3, 3, 4],
        "residue.aa.from": ["K", "K", "L"], "residue.aa.to": ["D", "D", "V"],
        "atom.from": ["NZ", "NZ", "CD1"], "atom.to": ["OD1", "OD2", "CG1"], "dist": [3.2, 3.4, 4.2],
    })
    c = contact_type_counts(_FakeCM(df), "tcr_peptide")
    assert c["n_salt_bridge"] == 2 and c["pairs_salt_bridge"] == 1     # two atom pairs, one residue pair
    assert c["n_hydrophobic"] == 1 and c["pairs_hydrophobic"] == 1
