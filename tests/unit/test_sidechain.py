"""The side-chain builder: chemistry declared, geometry measured, and the two directions invert."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from tcren.refine.sidechain import (N_CHI, TOPOLOGY, bond_angle, build_sidechain, chi_angles,
                                    dihedral, geometry, nerf, rebuild_residue)
from tcren.refine.substitute import virtual_cb

AAS = [a for a in TOPOLOGY if a not in ("A", "G")]

ASSET = Path(__file__).resolve().parents[1] / "assets" / "pdb" / "1ao7.pdb"


@pytest.fixture(scope="module")
def deposited():
    """A real structure, for the two tests that need deposited side chains rather than built ones."""
    from tcren.structure import parse_structure

    return parse_structure(ASSET, pdb_id="1ao7")


def _backbone():
    return np.array([0.0, 1.4, 0.0]), np.zeros(3), np.array([1.5, 0.0, 0.0])


def test_nerf_inverts_dihedral():
    """A sign error here is silent: every bond and angle stays right and every atom flips 180."""
    rng = np.random.default_rng(0)
    for _ in range(200):
        a, b, c, d = rng.normal(size=(4, 3)) * 3.0
        back = nerf(a, b, c, float(np.linalg.norm(d - c)), bond_angle(b, c, d),
                    dihedral(a, b, c, d))
        assert np.linalg.norm(d - back) < 1e-9


@pytest.mark.parametrize("aa", AAS)
def test_every_type_builds_its_declared_atoms(aa):
    n, ca, c = _backbone()
    built = build_sidechain(n, ca, c, virtual_cb(n, ca, c), aa, [0.5] * 4)
    assert set(built) == {t[0] for t in TOPOLOGY[aa]}
    assert all(np.isfinite(v).all() for v in built.values())


@pytest.mark.parametrize("aa", AAS)
def test_chi_round_trips_through_the_build(aa):
    """Building at a chi vector and reading it back returns the same vector."""
    n, ca, c = _backbone()
    want = [-1.1, 0.7, 2.3, -0.4][: N_CHI[aa]]
    xyz = {"N": n, "CA": ca, "C": c, "CB": virtual_cb(n, ca, c)}
    xyz |= build_sidechain(n, ca, c, xyz["CB"], aa, want + [0.0] * 4)
    got = []
    for k in range(1, N_CHI[aa] + 1):
        at, a, b, cc = next(r[:4] for r in TOPOLOGY[aa] if r[4] == k)
        got.append(dihedral(xyz[a], xyz[b], xyz[cc], xyz[at]))
    assert got == pytest.approx(want, abs=1e-9)


@pytest.mark.parametrize("aa", AAS)
def test_bond_lengths_are_chemistry(aa):
    """Every fitted bond sits in the range a covalent heavy-atom bond can take."""
    for atom, p in geometry()[aa].items():
        assert 1.2 <= p["bond"] <= 1.9, f"{aa}:{atom} bond {p['bond']}"
        assert 1.5 <= p["angle"] <= 2.5, f"{aa}:{atom} angle {p['angle']}"
        assert p["n"] > 500, f"{aa}:{atom} fitted on only {p['n']} residues"


def test_alanine_and_glycine_have_nothing_past_cb():
    n, ca, c = _backbone()
    cb = virtual_cb(n, ca, c)
    assert build_sidechain(n, ca, c, cb, "A", []) == {}
    assert build_sidechain(n, ca, c, cb, "G", []) == {}


def test_a_short_chi_vector_is_refused():
    n, ca, c = _backbone()
    with pytest.raises(ValueError, match="needs 4 chi"):
        build_sidechain(n, ca, c, virtual_cb(n, ca, c), "R", [0.1])


def test_rebuild_residue_keeps_the_backbone_and_adds_the_side_chain(deposited):
    """`substitute` stops at CB; this does not, and it moves no backbone atom doing it."""
    res = next(r for ch in deposited.chains for r in ch.residues
               if {"N", "CA", "C", "O"} <= {a.name for a in r.atoms})
    out = rebuild_residue(res, "W", [0.0, 0.0])
    before = {a.name: np.asarray(a.coord) for a in res.atoms if a.name in ("N", "CA", "C", "O")}
    after = {a.name: np.asarray(a.coord) for a in out.atoms if a.name in ("N", "CA", "C", "O")}
    assert before.keys() == after.keys()
    for k in before:
        assert np.allclose(before[k], after[k])
    assert {a.name for a in out.atoms} >= {t[0] for t in TOPOLOGY["W"]} | {"CB"}
    assert out.aa == "W" and out.resname == "TRP"


def test_chi_angles_reads_a_deposited_residue(deposited):
    """A real side chain yields as many chi as its type declares."""
    hits = 0
    for ch in deposited.chains:
        for r in ch.residues:
            if N_CHI.get(r.aa, 0) and {t[0] for t in TOPOLOGY[r.aa]} <= {a.name for a in r.atoms}:
                assert len(chi_angles(r)) == N_CHI[r.aa]
                hits += 1
    assert hits > 0, "no complete side chain in the fixture"
