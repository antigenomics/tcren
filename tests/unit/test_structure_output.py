"""Structure output: format-flag paths, PDB/mmCIF dispatch + round-trip, transform averaging."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from tcren.structure import parse_structure
from tcren.structure.io import structure_output_path, write_structure

_FIXTURE = Path(__file__).resolve().parents[1] / "assets" / "pdb" / "1ao7.pdb"


def test_output_path_format_flags():
    assert structure_output_path("/d", "x").name == "x.pdb"
    assert structure_output_path("/d", "x", mmcif=True).name == "x.cif"
    assert structure_output_path("/d", "x", compress=True).name == "x.pdb.gz"
    assert structure_output_path("/d", "x", mmcif=True, compress=True).name == "x.cif.gz"


def _atom_count(s):
    return sum(len(r.atoms) for c in s.chains for r in c.residues)


def test_pdb_and_mmcif_roundtrip(tmp_path):
    s = parse_structure(_FIXTURE, pdb_id="1ao7")
    for mmcif in (False, True):
        for compress in (False, True):
            p = write_structure(s, structure_output_path(tmp_path, "1ao7", mmcif, compress))
            s2 = parse_structure(p, pdb_id="1ao7")
            assert len(s2.chains) == len(s.chains)
            assert _atom_count(s2) == _atom_count(s)


def test_average_transform_is_orthonormal_and_exact_on_constant():
    from tcren.orient.superimpose import _average_transform

    rng = np.random.default_rng(0)
    q, _ = np.linalg.qr(rng.standard_normal((3, 3)))
    if np.linalg.det(q) < 0:
        q[:, -1] *= -1
    t = rng.standard_normal(3)
    rot, tran = _average_transform([(q, t), (q, t), (q, t)])
    assert np.allclose(rot @ rot.T, np.eye(3), atol=1e-9)        # proper rotation
    assert np.isclose(np.linalg.det(rot), 1.0, atol=1e-9)
    assert np.allclose(rot, q, atol=1e-9) and np.allclose(tran, t, atol=1e-9)
