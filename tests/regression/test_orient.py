"""Canonical orientation: frame geometry, chain renaming, reverse-dock flag, writer round-trip."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from tcren import parse_structure
from tcren.annotation import classify_chains

pytest.importorskip("Bio")
arda = pytest.importorskip("arda")

pytestmark = pytest.mark.slow  # arda + native DB + BioPython superposition

REPO = Path(__file__).resolve().parents[2]
# `data/Native2026` is what `setup.sh` fetches (`tcren fetch-data`), so it is the set a working
# checkout actually has. This used to read `notebooks/data/Native2022/*.pdb`, which no checkout has
# carried since the notebook data moved -- the four tests here had been failing with
# FileNotFoundError rather than skipping. The assertions are structural (chain renaming, z stacking,
# peptide N->C direction, canonical dock), so they hold for any deposit of the same PDB entry.
NATIVE_DIRS = (REPO / "data" / "Native2026", REPO / "notebooks" / "data" / "Native2022")


def _native(pid: str) -> Path:
    for d in NATIVE_DIRS:
        for name in (f"{pid}.pdb.gz", f"{pid}.pdb"):
            if (d / name).is_file():
                return d / name
    pytest.skip(f"no reference structure for {pid}; run `tcren fetch-data` to populate data/")


def _oriented(pid):
    from tcren.mhc import annotate_mhc
    from tcren.docking import canonicalize_structure

    s = parse_structure(_native(pid), pdb_id=pid)
    classify_chains(s, organism="human")
    annotate_mhc(s)
    return canonicalize_structure(s)


def _mean_z(structure, cid):
    pts = [r.ca for c in structure.chains if c.chain_id == cid for r in c.residues if r.ca is not None]
    return float(np.mean([p[2] for p in pts])) if pts else None


@pytest.mark.parametrize("pid", ["1ao7", "1fyt", "1bd2"])
def test_canonical_frame_geometry(pid):
    oriented, res = _oriented(pid)
    # chains renamed to A-E
    assert {c.chain_id for c in oriented.chains} <= {"A", "B", "C", "D", "E"}
    # z stacking: TCR (A/B) above peptide (C) above MHC floor (D)
    tcr_z = np.mean([z for z in (_mean_z(oriented, "A"), _mean_z(oriented, "B")) if z is not None])
    assert tcr_z > _mean_z(oriented, "C") > _mean_z(oriented, "D")
    # peptide N->C runs toward +y
    pep = [c for c in oriented.chains if c.chain_id == "C"][0]
    cas = [r.ca for r in pep.residues if r.ca is not None]
    assert cas[-1][1] > cas[0][1]  # C-term y > N-term y
    # 1ao7/1fyt/1bd2 are canonical docks
    assert res.reversed_dock is False
    assert res.frame == "native"


def test_oriented_pdb_round_trips(tmp_path):
    from tcren.structure.io import write_pdb

    oriented, _ = _oriented("1ao7")
    p = write_pdb(oriented, tmp_path / "1ao7_oriented.pdb")
    reparsed = parse_structure(p, pdb_id="1ao7")
    assert {c.chain_id for c in reparsed.chains} == {c.chain_id for c in oriented.chains}
    # atom coordinates preserved
    a = oriented.chains[0].residues[0].atoms[0].coord
    b = reparsed.chains[0].residues[0].atoms[0].coord
    assert np.allclose(a, b, atol=1e-3)
