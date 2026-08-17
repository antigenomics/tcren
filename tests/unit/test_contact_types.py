"""Unit tests for chemical contact typing (pure logic, no structure/mmseqs)."""

from __future__ import annotations

from pathlib import Path

import polars as pl
import pytest

from tcren.contact_types import (
    TYPES_V1,
    TYPES_V2,
    _classify,
    _is_apolar,
    _types_v2,
    classify_contacts,
    contact_type_counts,
)

PDB_DIR = Path(__file__).resolve().parents[1] / "assets" / "pdb"


# --- v1, frozen: the trained recognition models were fitted on these counts -----------------------
@pytest.mark.parametrize(
    "aa_a,aa_b,atom_a,atom_b,dist,expected",
    [
        ("K", "D", "NZ", "OD1", 3.2, "salt_bridge"),     # Lys+ <-> Asp-
        ("R", "E", "NH1", "OE2", 3.8, "salt_bridge"),    # Arg+ <-> Glu-
        ("K", "D", "NZ", "OD1", 4.5, "other"),           # beyond both the salt and h-bond cutoffs
        ("Q", "T", "NE2", "OG1", 2.9, "hydrogen_bond"),  # amide N <-> hydroxyl O
        ("S", "N", "OG", "ND2", 3.1, "hydrogen_bond"),   # two polar
        ("F", "Y", "CZ", "CE1", 4.5, "aromatic"),        # ring-ring aromatic
        ("L", "V", "CD1", "CG1", 4.2, "hydrophobic"),    # apolar C-C
        ("G", "P", "CA", "CB", 4.0, "other"),            # Gly is not in v1's apolar residue set
        ("Q", "P", "OE1", "CG", 3.0, "other"),           # polar O <-> carbon: not an h-bond
    ],
)
def test_classify_v1_is_unchanged(aa_a, aa_b, atom_a, atom_b, dist, expected):
    got = _classify(aa_a, aa_b, atom_a, atom_b, dist)
    assert got == expected
    assert got in TYPES_V1


def test_v1_scheme_still_reachable_through_the_public_api():
    df = pl.DataFrame({
        "chain.id.from": ["A"], "residue.index.from": [1],
        "chain.id.to": ["C"], "residue.index.to": [3],
        "residue.aa.from": ["Y"], "residue.aa.to": ["G"],
        "atom.from": ["CE1"], "atom.to": ["CA"], "dist": [4.0],
    })
    assert classify_contacts(df, "v1")["contact.type"].to_list() == ["other"]
    assert "is_hydrophobic" not in classify_contacts(df, "v1").columns


def test_unknown_scheme_raises():
    with pytest.raises(ValueError, match="v1"):
        classify_contacts(pl.DataFrame({"dist": [1.0]}), scheme="v3")


# --- v2: atom-level apolarity ---------------------------------------------------------------------
@pytest.mark.parametrize("aa,atom,apolar", [
    ("Y", "CE1", True),    # a Tyr ring carbon IS apolar; v1's residue-level test excluded all of Tyr
    ("Y", "CZ", False),    # ...but the one bonded to OH is not
    ("R", "CG", True),     # aliphatic stem of a charged residue
    ("R", "CZ", False),    # guanidinium carbon
    ("K", "CD", True),
    ("K", "CE", False),    # bonded to NZ
    ("T", "CG2", True),
    ("T", "CB", False),    # bonded to OG1
    ("W", "CZ2", True),
    ("W", "CD1", False),   # bonded to NE1
    ("M", "SD", True),     # sulfur packs like carbon
    ("L", "CA", False),    # backbone alpha carbon is bonded to N
    ("D", "OD1", False),   # not a carbon at all
])
def test_apolarity_is_decided_per_atom(aa, atom, apolar):
    assert _is_apolar(aa, atom) is apolar


# --- v2: the classes v1 was missing ---------------------------------------------------------------
def test_mixed_polar_apolar_contacts_are_polar_not_other():
    """59% of v1's `other` rows were C-O or C-N pairs, for which it had no class."""
    assert "polar" in _types_v2("Q", "P", "OE1", "CG", 3.0)


def test_tyr_ring_against_leucine_is_hydrophobic():
    assert "hydrophobic" in _types_v2("Y", "L", "CE1", "CD1", 4.0)
    assert _classify("Y", "L", "CE1", "CD1", 4.0) == "other"      # what v1 said


def test_hbond_reaches_further_without_hydrogens():
    assert "hydrogen_bond" in _types_v2("S", "N", "OG", "OD1", 3.7)
    assert _classify("S", "N", "OG", "OD1", 3.7) == "other"       # v1's 3.5 A cutoff


def test_two_carbonyl_oxygens_are_not_a_hydrogen_bond():
    """Both are pure acceptors; v1 called any N/O pair within 3.5 A an h-bond."""
    assert "hydrogen_bond" not in _types_v2("D", "E", "OD1", "OE1", 3.0)
    assert _classify("D", "E", "OD1", "OE1", 3.0) == "hydrogen_bond"


def test_proline_backbone_nitrogen_cannot_donate():
    assert "hydrogen_bond" not in _types_v2("P", "D", "N", "OD1", 3.0)
    assert "hydrogen_bond" in _types_v2("L", "D", "N", "OD1", 3.0)


def test_cation_pi_is_detected():
    assert "cation_pi" in _types_v2("R", "W", "NH1", "CZ2", 4.5)
    assert "cation_pi" not in _types_v2("R", "L", "NH1", "CD1", 4.5)


def test_a_contact_can_carry_more_than_one_type():
    """A charged N-O pair is both a salt bridge and a hydrogen bond; one label loses that."""
    hits = _types_v2("K", "D", "NZ", "OD1", 3.2)
    assert {"salt_bridge", "hydrogen_bond"} <= hits


def test_classify_contacts_v2_emits_independent_booleans():
    df = pl.DataFrame({
        "chain.id.from": ["A", "A"], "residue.index.from": [1, 2],
        "chain.id.to": ["C", "C"], "residue.index.to": [3, 3],
        "residue.aa.from": ["K", "L"], "residue.aa.to": ["D", "V"],
        "atom.from": ["NZ", "CD1"], "atom.to": ["OD1", "CG1"], "dist": [3.2, 4.2],
    })
    out = classify_contacts(df)
    assert out["contact.type"].to_list() == ["salt_bridge", "hydrophobic"]
    assert out["is_salt_bridge"].to_list() == [True, False]
    assert out["is_hydrogen_bond"].to_list() == [True, False]      # the same row, both types
    assert all(f"is_{t}" in out.columns for t in TYPES_V2)


def test_classify_contacts_empty():
    df = pl.DataFrame({c: [] for c in
                       ["residue.aa.from", "residue.aa.to", "atom.from", "atom.to", "dist"]})
    for scheme in ("v1", "v2"):
        out = classify_contacts(df, scheme)
        assert "contact.type" in out.columns and out.height == 0


class _FakeCM:
    def __init__(self, df):
        self._df = df

    def interface(self, name, tcr_regions="all"):
        return self._df


def test_contact_type_counts_v1():
    df = pl.DataFrame({
        "chain.id.from": ["A", "A", "A"], "residue.index.from": [1, 1, 2],
        "chain.id.to": ["C", "C", "C"], "residue.index.to": [3, 3, 4],
        "residue.aa.from": ["K", "K", "L"], "residue.aa.to": ["D", "D", "V"],
        "atom.from": ["NZ", "NZ", "CD1"], "atom.to": ["OD1", "OD2", "CG1"], "dist": [3.2, 3.4, 4.2],
    })
    c = contact_type_counts(_FakeCM(df), "tcr_peptide", scheme="v1")
    assert c["n_salt_bridge"] == 2 and c["pairs_salt_bridge"] == 1   # two atom pairs, one residue pair
    assert c["n_hydrophobic"] == 1 and c["pairs_hydrophobic"] == 1


def test_contact_type_counts_v2_requires_a_structure():
    with pytest.raises(ValueError, match="structure"):
        contact_type_counts(_FakeCM(pl.DataFrame()), "tcr_peptide", scheme="v2")


# --- end to end on real interfaces ----------------------------------------------------------------
@pytest.fixture(scope="module")
def crystals():
    pytest.importorskip("arda")
    from tcren.annotation import classify_chains
    from tcren.mhc import annotate_mhc
    from tcren.structure import parse_structure

    out = []
    for pid in ("1ao7", "1bd2", "2ckb", "5m01", "6bj3"):
        s = parse_structure(PDB_DIR / f"{pid}.pdb")
        classify_chains(s, organism="human", autodetect_species=True)
        annotate_mhc(s)
        out.append(s)
    return out


@pytest.mark.slow
def test_other_no_longer_dominates(crystals):
    """The measured regression: v1 typed 72.3% of TCR:peptide contacts as `other`."""
    from tcren.contact_types import residue_pair_types

    counts = {"v1": 0, "v2": 0, "n": 0}
    for s in crystals:
        df = residue_pair_types(s, "tcr_peptide")
        counts["n"] += df.height
        counts["v2"] += int((df["contact.type"] == "other").sum())
    assert counts["n"] > 50
    assert counts["v2"] / counts["n"] < 0.25


@pytest.mark.slow
def test_every_typed_contact_is_chemically_possible(crystals):
    """No contact may be given a type its atoms cannot support."""
    from tcren.contact_types import _ANIONIC_ATOMS, _CATIONIC_ATOMS, _RING_ATOMS, _typed_atom_pairs

    for s in crystals:
        df = _typed_atom_pairs(s, "tcr_peptide", "all", 5.0)
        for row in df.iter_rows(named=True):
            aa1, aa2 = row["residue.aa.from"], row["residue.aa.to"]
            a1, a2 = row["atom.from"], row["atom.to"]
            if row["is_salt_bridge"]:
                assert ({a1, a2} & _CATIONIC_ATOMS) and ({a1, a2} & _ANIONIC_ATOMS)
                assert row["dist"] <= 4.0
            if row["is_hydrophobic"]:
                assert _is_apolar(aa1, a1) and _is_apolar(aa2, a2)
                assert row["dist"] <= 4.5
            if row["is_aromatic"]:
                assert a1 in _RING_ATOMS and a2 in _RING_ATOMS
            if row["is_hydrogen_bond"]:
                assert {a1[0], a2[0]} <= {"N", "O"}
            if row["contact.type"] == "other":
                assert row["dist"] > 4.5        # `other` now means only "distant", not "unrecognised"


@pytest.mark.slow
def test_typing_is_unchanged_by_explicit_hydrogens(crystals):
    """An AF/OpenMM-relaxed model carries H; before the fix it typed as ~100% `other`."""
    import copy

    import numpy as np

    from tcren.contact_types import residue_pair_types
    from tcren.structure.model import Atom

    import dataclasses

    s = copy.deepcopy(crystals[0])
    for chain in s.chains:
        chain.residues = [
            dataclasses.replace(res, atoms=tuple(res.atoms) + (
                Atom(name="H", element="H", coord=np.asarray(res.atoms[0].coord) + 0.4),))
            for res in chain.residues
        ]

    base = residue_pair_types(crystals[0], "tcr_peptide")
    with_h = residue_pair_types(s, "tcr_peptide")
    assert base.height == with_h.height
    assert base["contact.type"].to_list() == with_h["contact.type"].to_list()
