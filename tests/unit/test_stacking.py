"""Unit tests for ring-stacking geometry.

The point of this module is the distinction a contact potential cannot make: two rings
face-to-face at 3.5 Å versus the same two rings edge-on. The synthetic tests build both
arrangements explicitly, so a regression that loses the angle would fail rather than merely
shift a number.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from tcren.stacking import RING_ATOMS, ring_of, ring_stacking
from tcren.structure.io import parse_structure
from tcren.structure.model import Atom, Chain, Residue, Structure

ASSET = Path(__file__).resolve().parents[1] / "assets" / "pdb" / "1ao7.pdb"


def _phe(seq, centre, normal):
    """A flat hexagon of the right atom names, centred and oriented as asked."""
    centre, normal = np.asarray(centre, float), np.asarray(normal, float)
    normal = normal / np.linalg.norm(normal)
    helper = np.array([1.0, 0.0, 0.0])
    if abs(helper @ normal) > 0.9:
        helper = np.array([0.0, 1.0, 0.0])
    u = np.cross(normal, helper); u /= np.linalg.norm(u)
    v = np.cross(normal, u)
    atoms = tuple(
        Atom(name=n, element="C",
             coord=centre + 1.4 * (np.cos(t) * u + np.sin(t) * v))
        for n, t in zip(RING_ATOMS["PHE"], np.linspace(0, 2 * np.pi, 6, endpoint=False),
                        strict=True)
    )
    return Residue(seq_index=seq, pdb_index=seq, insertion_code="", aa="F",
                   resname="PHE", atoms=atoms)


def _two_ring_chain(centre_b, normal_b):
    a = _phe(0, (0.0, 0.0, 0.0), (0.0, 0.0, 1.0))
    b = _phe(4, centre_b, normal_b)
    return Structure(pdb_id="t", chains=[Chain(chain_id="A", residues=[a, b])])


def test_parallel_stack_is_flat_and_close():
    s = _two_ring_chain((0.0, 0.0, 3.6), (0.0, 0.0, 1.0))
    row = ring_stacking(s, min_seq_sep=3).to_dicts()[0]
    assert row["centroid_distance"] == pytest.approx(3.6)
    assert row["interplanar_angle"] == pytest.approx(0.0, abs=1e-6)
    assert row["vertical"] == pytest.approx(3.6)
    assert row["lateral"] == pytest.approx(0.0, abs=1e-6)


def test_edge_to_face_is_perpendicular():
    """Same separation, same residues, entirely different interaction."""
    s = _two_ring_chain((0.0, 0.0, 3.6), (1.0, 0.0, 0.0))
    row = ring_stacking(s, min_seq_sep=3).to_dicts()[0]
    assert row["interplanar_angle"] == pytest.approx(90.0, abs=1e-6)


def test_parallel_displaced_splits_into_vertical_and_lateral():
    s = _two_ring_chain((3.0, 0.0, 4.0), (0.0, 0.0, 1.0))
    row = ring_stacking(s, min_seq_sep=3).to_dicts()[0]
    assert row["vertical"] == pytest.approx(4.0)
    assert row["lateral"] == pytest.approx(3.0)
    assert row["centroid_distance"] == pytest.approx(5.0)


def test_sequence_neighbours_can_be_excluded():
    a = _phe(0, (0.0, 0.0, 0.0), (0.0, 0.0, 1.0))
    b = _phe(1, (0.0, 0.0, 3.6), (0.0, 0.0, 1.0))
    s = Structure(pdb_id="t", chains=[Chain(chain_id="A", residues=[a, b])])
    assert ring_stacking(s, min_seq_sep=3).height == 0
    assert ring_stacking(s, min_seq_sep=1).height == 1


def test_a_residue_without_a_ring_has_none():
    res = Residue(seq_index=0, pdb_index=0, insertion_code="", aa="L", resname="LEU",
                  atoms=(Atom(name="CA", element="C", coord=np.zeros(3)),))
    assert ring_of(res, "A") is None


def test_an_incomplete_side_chain_has_no_plane():
    """Half a ring does not define one, and must not be guessed at."""
    full = _phe(0, (0.0, 0.0, 0.0), (0.0, 0.0, 1.0))
    partial = Residue(seq_index=0, pdb_index=0, insertion_code="", aa="F", resname="PHE",
                      atoms=full.atoms[:3])
    assert ring_of(partial, "A") is None


def test_runs_on_a_deposited_structure():
    s = parse_structure(ASSET)
    found = ring_stacking(s, min_seq_sep=3)
    assert found.height > 0
    assert (found["centroid_distance"] <= 7.5).all()
    assert (found["interplanar_angle"] >= 0).all()
    assert (found["interplanar_angle"] <= 90).all()
