"""Unit tests for the intra-peptide energy term.

Every interface energy in the package sums over contacts between two *different* chains, so a
peptide held in its bound conformation by its own side chains costs the same as one that is not.
``intra_weight`` is that omitted term, made optional. The tests that matter are the two guards:
the term is off by default and changes nothing when it is, and the potential is symmetrised —
an intra-chain pair has no ``from``/``to`` orientation, so a directed potential such as TCRen must
not be read as though it did.
"""

from __future__ import annotations

from pathlib import Path

import polars as pl
import pytest

from tcren.contactmap import ContactMap
from tcren.potential import Potential, mj, tcren
from tcren.scoring import intra_peptide_energy, score_peptides

pytest.importorskip("arda")

from tcren.annotation import classify_chains  # noqa: E402
from tcren.structure.io import parse_structure  # noqa: E402

ASSET = Path(__file__).resolve().parents[1] / "assets" / "pdb" / "1ao7.pdb"
CANDIDATES = ("LLFGYPVYV", "LLWGYPWYV", "AAAAAAAAA")


@pytest.fixture(scope="module")
def structure():
    s = parse_structure(ASSET)
    classify_chains(s)
    return s


# --- the toy path: no structure needed, so the arithmetic is checkable by hand ------------------


def _directed_potential() -> Potential:
    """A potential that is deliberately *not* symmetric: F[A,K] = -2, F[K,A] = +4."""
    rows = [
        {"residue.aa.from": "A", "residue.aa.to": "A", "value": 1.0},
        {"residue.aa.from": "A", "residue.aa.to": "K", "value": -2.0},
        {"residue.aa.from": "K", "residue.aa.to": "A", "value": 4.0},
        {"residue.aa.from": "K", "residue.aa.to": "K", "value": 0.5},
    ]
    return Potential(name="directed", matrix=pl.DataFrame(rows), alphabet=("A", "K"))


def _toy_map(n_atom_contacts: int | None = None) -> ContactMap:
    """One intra-peptide pair: position 0 (A) touching position 4 (K)."""
    cols = {
        "chain.id.from": ["C"], "chain.id.to": ["C"],
        "residue.index.from": [0], "residue.index.to": [4],
        "residue.aa.from": ["A"], "residue.aa.to": ["K"],
        "pos.from": [0], "pos.to": [4],
    }
    if n_atom_contacts is not None:
        cols["n_atom_contacts"] = [n_atom_contacts]
    return ContactMap(pdb_id="toy", contacts=pl.DataFrame(), peptide_length=5,
                      peptide_internal=pl.DataFrame(cols))


def test_the_potential_is_symmetrised():
    """An intra-chain pair has no from/to orientation, so the two directions are averaged."""
    energy = intra_peptide_energy(_toy_map(), _directed_potential())
    assert energy == pytest.approx((-2.0 + 4.0) / 2)


def test_threading_a_candidate_reads_the_pair_positions():
    # Position 0 -> K and position 4 -> A: the same unordered pair, so the same symmetrised value.
    energy = intra_peptide_energy(_toy_map(), _directed_potential(), peptide="KAAAA")
    assert energy == pytest.approx((-2.0 + 4.0) / 2)
    # Both positions A: F[A,A] = 1.0, symmetric already.
    assert intra_peptide_energy(_toy_map(), _directed_potential(), peptide="AAAAA") == pytest.approx(1.0)


def test_atomic_weighting_multiplies_by_the_atom_pair_count():
    pot = _directed_potential()
    assert intra_peptide_energy(_toy_map(3), pot, contact_weight="atomic") == pytest.approx(3.0)


def test_atomic_weighting_without_the_count_names_the_flag():
    with pytest.raises(ValueError, match="count_atoms=True"):
        intra_peptide_energy(_toy_map(), _directed_potential(), contact_weight="atomic")


def test_a_map_without_internal_contacts_names_the_flag():
    cm = ContactMap(pdb_id="toy", contacts=pl.DataFrame())
    with pytest.raises(ValueError, match="peptide_internal=True"):
        intra_peptide_energy(cm, _directed_potential())


# --- the structure path ------------------------------------------------------------------------


def test_collecting_internal_contacts_leaves_the_contact_table_untouched(structure):
    """The interface table every score is built on must not notice that this option exists."""
    plain = ContactMap.from_structure(structure)
    with_internal = ContactMap.from_structure(structure, peptide_internal=True)
    assert plain.contacts.equals(with_internal.contacts)
    assert plain.peptide_internal is None
    assert with_internal.peptide_internal.height > 0


def test_intra_weight_zero_is_byte_identical(structure):
    cm = ContactMap.from_structure(structure, peptide_internal=True)
    pot = tcren()
    base = score_peptides(cm, CANDIDATES, pot)
    assert base.equals(score_peptides(cm, CANDIDATES, pot, intra_weight=0.0))


def test_intra_weight_adds_exactly_the_term(structure):
    cm = ContactMap.from_structure(structure, peptide_internal=True)
    pot, intra_pot, w = tcren(), mj(), 2.5
    base = score_peptides(cm, CANDIDATES, pot)
    scored = score_peptides(cm, CANDIDATES, pot, intra_weight=w, intra_potential=intra_pot)
    for row in scored.iter_rows(named=True):
        peptide = row["peptide"]
        interface = base.filter(pl.col("peptide") == peptide)["score"][0]
        term = intra_peptide_energy(cm, intra_pot, peptide=peptide)
        assert row["score"] == pytest.approx(interface + w * term)


def test_the_term_separates_candidates_the_interface_sum_cannot(structure):
    """The point of the term: two candidates can differ only in their contacts with themselves."""
    cm = ContactMap.from_structure(structure, peptide_internal=True)
    pot = mj()
    # 1ao7's peptide touches itself at positions 2 and 5 only, so these two differ there and
    # nowhere the intra term can see -- but they are not equally happy packed against themselves.
    assert intra_peptide_energy(cm, pot, peptide="LLFGYPVYV") != pytest.approx(
        intra_peptide_energy(cm, pot, peptide="LLWGYPWYV")
    )


def test_a_peptide_that_touches_nothing_of_itself_scores_zero(structure):
    """5m01's peptide makes no internal contact at all: the term is 0, not NaN or an error."""
    other = parse_structure(ASSET.with_name("5m01.pdb"))
    classify_chains(other)
    cm = ContactMap.from_structure(other, peptide_internal=True)
    assert cm.peptide_internal.height == 0
    assert intra_peptide_energy(cm, mj()) == 0.0
