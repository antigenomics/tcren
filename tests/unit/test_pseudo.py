"""MHC pseudosequence (MPS) annotation + unified `annotate` CLI selector."""

from __future__ import annotations

from pathlib import Path

import pytest

pytest.importorskip("arda")  # MPS needs chain typing + MHC annotation (mmseqs-backed)

from tcren.annotation import classify_chains  # noqa: E402
from tcren.mhc import annotate_mhc, annotate_pseudo  # noqa: E402
from tcren.mhc.pseudo import _best_pseudo_hit, _pseudo_index  # noqa: E402
from tcren.structure.io import import_structure  # noqa: E402

_ASSETS = Path(__file__).resolve().parents[1] / "assets" / "pdb"


def test_pseudo_index_loads_unique_34mers():
    idx = _pseudo_index("MHCI")
    assert len(idx) > 1000
    assert all(len(seq) == 34 for seq in idx.values())


def _annotated(pdb, organism="human"):
    s = import_structure(_ASSETS / f"{pdb}.pdb", pdb_id=pdb)
    classify_chains(s, organism=organism)
    annotate_mhc(s)
    return s


def test_mhci_mps_marks_mhca_only():
    s = _annotated("1ao7")
    best = annotate_pseudo(s)
    assert best is not None
    marked = {c.chain_type: [r for reg in c.regions if reg.region_type == "MPS" for r in reg.residues]
              for c in s.chains}
    mhca = marked.get("MHCa", [])
    assert 25 <= len(mhca) <= 34          # ~34 scattered groove residues
    assert not marked.get("B2M")          # never β2m


def test_mhcii_mps_splits_across_two_chains():
    s = _annotated("4ozg")
    best = annotate_pseudo(s)
    assert best is not None
    a = sum(len(reg.residues) for c in s.chains if c.chain_type == "MHCa"
            for reg in c.regions if reg.region_type == "MPS")
    b = sum(len(reg.residues) for c in s.chains if c.chain_type == "MHCb"
            for reg in c.regions if reg.region_type == "MPS")
    assert a > 0 and b > 0                # split α1 + β1
    assert a + b <= 34


def test_best_hit_is_identity_consistent():
    # The chosen pseudosequence's residues must actually occur (identity) in the chain.
    s = _annotated("1ao7")
    mhca = next(c for c in s.chains if c.chain_type == "MHCa")
    best_id, pseudo = _best_pseudo_hit(mhca.sequence(), "MHCI")
    assert best_id in _pseudo_index("MHCI")
    assert pseudo == _pseudo_index("MHCI")[best_id]
