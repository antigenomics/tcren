"""Unit tests for atom-weighted scoring (``contact_weight="atomic"``).

``contact_weight="residue"`` (default) gives each contacting residue pair unit weight —
byte-identical to the legacy scorer. ``contact_weight="atomic"`` multiplies each pair's
potential by its ``n_atom_contacts`` heavy-atom-pair count. On a toy contact map with a
multi-atom contact the two modes must differ, and the atomic score must equal the
hand-computed weighted sum.
"""

from __future__ import annotations

import polars as pl
import pytest

from tcren.contactmap import ContactMap
from tcren.potential import Potential
from tcren.scoring import score_peptides


def _toy_potential() -> Potential:
    rows = [
        {"residue.aa.from": "A", "residue.aa.to": "G", "value": 0.1},
        {"residue.aa.from": "L", "residue.aa.to": "G", "value": 0.2},
        {"residue.aa.from": "A", "residue.aa.to": "K", "value": -2.0},
        {"residue.aa.from": "L", "residue.aa.to": "K", "value": 3.0},
    ]
    return Potential(name="toy", matrix=pl.DataFrame(rows), alphabet=("A", "L", "K", "G"))


def _toy_contact_map(n_atom_from=3, n_atom_to=5) -> ContactMap:
    # Two TCR<->peptide contacts with distinct atom-pair counts.
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
            "n_atom_contacts": [n_atom_from, n_atom_to],
            "pdb.id": ["toy", "toy"],
        }
    )
    return ContactMap(pdb_id="toy", contacts=contacts, peptide_length=3)


def test_residue_mode_matches_legacy_default():
    cm = _toy_contact_map()
    pot = _toy_potential()
    # Candidate "AGK": pos0='A' -> (A,A) absent -> 0 contribution actually? No: 'A' is the
    # substituted peptide aa at pos0, paired with TCR 'A' -> (A,A) not in potential -> the
    # pair is dropped. pos2='K' -> TCR 'L' -> (L,K)=3.0. So residue score = 3.0.
    default = score_peptides(cm, ["AGK"], pot)
    explicit = score_peptides(cm, ["AGK"], pot, contact_weight="residue")
    assert default["score"][0] == pytest.approx(explicit["score"][0])
    assert default["score"][0] == pytest.approx(3.0)


def test_atomic_weights_each_contact():
    cm = _toy_contact_map(n_atom_from=3, n_atom_to=5)
    pot = _toy_potential()
    # "AGK": only the pos2 contact (TCR 'L' x peptide 'K' = 3.0) is in the potential, with
    # n_atom_contacts=5 on the second contact -> atomic score = 3.0 * 5 = 15.0.
    atomic = score_peptides(cm, ["AGK"], pot, contact_weight="atomic")
    residue = score_peptides(cm, ["AGK"], pot, contact_weight="residue")
    assert atomic["score"][0] == pytest.approx(15.0)
    assert atomic["score"][0] != pytest.approx(residue["score"][0])


def test_atomic_both_contacts_counted():
    cm = _toy_contact_map(n_atom_from=3, n_atom_to=5)
    pot = _toy_potential()
    # "GGG": pos0 TCR 'A' x 'G' = 0.1 (w=3), pos2 TCR 'L' x 'G' = 0.2 (w=5).
    atomic = score_peptides(cm, ["GGG"], pot, contact_weight="atomic")
    assert atomic["score"][0] == pytest.approx(0.1 * 3 + 0.2 * 5)


def test_atomic_requires_count_column():
    # A contact map without n_atom_contacts must raise under atomic weighting.
    contacts = pl.DataFrame(
        {
            "chain.type.from": ["TRA"],
            "chain.type.to": ["PEPTIDE"],
            "residue.aa.from": ["A"],
            "residue.aa.to": ["G"],
            "region.type.from": ["CDR3"],
            "residue.index.from": [10],
            "residue.index.to": [0],
            "region.start.from": [8],
            "region.start.to": [0],
            "pdb.id": ["toy"],
        }
    )
    cm = ContactMap(pdb_id="toy", contacts=contacts, peptide_length=3)
    with pytest.raises(ValueError, match="n_atom_contacts"):
        score_peptides(cm, ["AGK"], _toy_potential(), contact_weight="atomic")


def test_invalid_contact_weight_rejected():
    cm = _toy_contact_map()
    with pytest.raises(ValueError, match="contact_weight"):
        score_peptides(cm, ["AGK"], _toy_potential(), contact_weight="nonsense")
