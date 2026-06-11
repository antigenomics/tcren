"""Tests for the native-database uses: canonical alignment + potential re-derivation."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("arda")

pytestmark = pytest.mark.slow  # invokes arda / mmseqs per structure

from tcren import parse_structure
from tcren.annotation import classify_chains
from tcren.mhc import annotate_mhc
from tcren.native import NativeDatabase, align_to_native, apply_transform, derive_native_potential

REPO = Path(__file__).resolve().parents[2]
_DB = NativeDatabase()
_HAVE_REF = (REPO / "database" / "mhc" / "alleles.aa.fasta").exists()
needs_data = pytest.mark.skipif(
    not (_DB.is_present() and _HAVE_REF),
    reason="native database and/or MHC reference not built",
)


def _prep(db, pdb_id, organism="human"):
    s = parse_structure(db.cif_for(pdb_id), pdb_id=pdb_id)
    classify_chains(s, organism=organism)
    annotate_mhc(s)
    return s


@needs_data
def test_self_alignment_is_identity():
    db = NativeDatabase()
    result = align_to_native(_prep(db, "1ao7"), db, reference_id="1ao7")
    assert result.rmsd < 1e-6
    assert result.n_anchor_atoms > 100


@needs_data
def test_same_allele_aligns_tightly():
    db = NativeDatabase()
    # 1bd2 and the 1ao7 reference are both HLA-A*02 class I → groove superposes tightly.
    result = align_to_native(_prep(db, "1bd2"), db, reference_id="1ao7")
    assert result.rmsd < 2.0
    assert result.n_anchor_atoms > 50


@needs_data
def test_apply_transform_brings_into_frame():
    db = NativeDatabase()
    mobile = _prep(db, "1bd2")
    result = align_to_native(mobile, db, reference_id="1ao7")
    moved = apply_transform(mobile, result)
    # Re-deriving the transform on the moved copy must be the identity rotation + zero
    # translation (it is already optimally placed in the frame). The residual rmsd is
    # unchanged — it is the irreducible structural difference between the two grooves,
    # which a rigid transform cannot remove.
    classify_chains(moved, organism="human")
    annotate_mhc(moved)
    again = align_to_native(moved, db, reference_id="1ao7")
    assert np.allclose(again.rotation, np.eye(3), atol=1e-3)
    assert np.allclose(again.translation, 0.0, atol=1e-2)
    assert again.rmsd == pytest.approx(result.rmsd, abs=1e-6)


@needs_data
def test_derive_native_potential_shape():
    db = NativeDatabase()
    pot = derive_native_potential(
        db, include=["1ao7", "1bd2", "1mi5", "1oga", "1qrn"], use_cache=False
    )
    assert pot.matrix.height == 380  # 19 (no Cys 'from') x 20, classic
    assert pot.matrix["value"].is_finite().all()
    # 'from' axis excludes Cys (classic variant).
    assert "C" not in set(pot.matrix["residue.aa.from"])
