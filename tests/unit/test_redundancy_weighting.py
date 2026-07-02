"""Unit tests for redundancy weighting in TCRen derivation.

Covers the two new pieces:

* :func:`tcren.potential.derive.derive_tcren` ``weights`` argument — ``weights=None`` is
  byte-identical to the legacy unweighted derivation; a uniform ``weight=1.0`` map
  reproduces it numerically; non-uniform weights down-weight a structure's contributions.
* :func:`tcren.potential.redundancy.cluster_weights` — inverse-cluster-size weights that
  sum to the cluster count and assign ``1/size`` per member.
"""

from __future__ import annotations

import polars as pl
import pytest

from tcren.potential import cluster_weights, derive_tcren


def _contacts(rows: list[dict]) -> pl.DataFrame:
    return pl.DataFrame(rows)


def _two_structure_contacts() -> pl.DataFrame:
    # Structure 'x' contributes A->D pairs; 'y' contributes L->A pairs.
    return _contacts(
        [{"pdb.id": "x", "residue.aa.from": "A", "residue.aa.to": "D"}] * 4
        + [{"pdb.id": "y", "residue.aa.from": "L", "residue.aa.to": "A"}] * 2
    )


def test_weights_none_is_byte_exact():
    contacts = _two_structure_contacts()
    base = derive_tcren(contacts, variant="classic")
    again = derive_tcren(contacts, variant="classic", weights=None)
    assert again.matrix.equals(base.matrix)


def test_uniform_weight_one_matches_unweighted():
    contacts = _two_structure_contacts()
    base = derive_tcren(contacts, variant="classic")
    uniform = derive_tcren(contacts, variant="classic", weights={"x": 1.0, "y": 1.0})
    j = base.matrix.join(
        uniform.matrix, on=["residue.aa.from", "residue.aa.to"], suffix="_w"
    )
    max_abs = j.select((pl.col("value") - pl.col("value_w")).abs().max()).item()
    assert max_abs == pytest.approx(0.0, abs=1e-12)


def test_missing_id_defaults_to_weight_one():
    # 'y' absent from the weights map -> defaults to 1.0, so this equals weighting x only.
    contacts = _two_structure_contacts()
    explicit = derive_tcren(contacts, variant="classic", weights={"x": 0.5, "y": 1.0})
    defaulted = derive_tcren(contacts, variant="classic", weights={"x": 0.5})
    assert explicit.matrix.equals(defaulted.matrix)


def test_nonuniform_weights_change_potential():
    contacts = _two_structure_contacts()
    base = derive_tcren(contacts, variant="classic")
    # Down-weight 'x' heavily: its A->D contribution shrinks, so the (A,D) energy shifts.
    weighted = derive_tcren(contacts, variant="classic", weights={"x": 0.1, "y": 1.0})
    assert weighted.value("A", "D") != pytest.approx(base.value("A", "D"))


def _markup(rows: list[dict]) -> pl.DataFrame:
    return pl.DataFrame(rows)


def test_cluster_weights_inverse_size_and_sum():
    # Two tight clusters: {a1, a2} (size 2) and {b1} (size 1).
    mk = _markup([
        {"pdb.id": "a1", "cdr3a": "CASSLAPGATNEKLFF", "cdr3b": "CASSIRSSYEQYF", "peptide": "GILGFVFTL"},
        {"pdb.id": "a2", "cdr3a": "CASSLAPGATNEKLFF", "cdr3b": "CASSIRSSYEQYF", "peptide": "GILGFVFTL"},
        {"pdb.id": "b1", "cdr3a": "CAVRDGGTGNGKLTF", "cdr3b": "CASSPGQGAYEQYF", "peptide": "NLVPMVATV"},
    ])
    w = cluster_weights(mk, t=6.0)
    assert set(w) == {"a1", "a2", "b1"}
    assert w["a1"] == pytest.approx(0.5)
    assert w["a2"] == pytest.approx(0.5)
    assert w["b1"] == pytest.approx(1.0)
    # Total weight equals the number of clusters (each cluster contributes 1).
    assert sum(w.values()) == pytest.approx(2.0)


def test_cluster_weights_all_unique():
    # Three mutually distant structures -> three singletons -> every weight 1.0.
    mk = _markup([
        {"pdb.id": "a", "cdr3a": "CASSLAPGATNEKLFF", "cdr3b": "CASSIRSSYEQYF", "peptide": "GILGFVFTL"},
        {"pdb.id": "b", "cdr3a": "CAVRDGGTGNGKLTF", "cdr3b": "CASSPGQGAYEQYF", "peptide": "NLVPMVATV"},
        {"pdb.id": "c", "cdr3a": "CAASFGDNSKLIW", "cdr3b": "CASRGTGELFF", "peptide": "ELAGIGILTV"},
    ])
    w = cluster_weights(mk, t=6.0)
    assert all(v == pytest.approx(1.0) for v in w.values())
    assert sum(w.values()) == pytest.approx(3.0)


def test_cluster_weights_single_structure():
    mk = _markup([{"pdb.id": "solo", "cdr3a": "CASS", "cdr3b": "CASS", "peptide": "GIL"}])
    assert cluster_weights(mk, t=6.0) == {"solo": 1.0}
