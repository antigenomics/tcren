"""Unit tests for non-redundancy clustering (:mod:`tcren.potential.redundancy`).

A verbatim lift of the notebook-01 ``nonredundant``/``alphabeta`` helpers; the checks
here pin the three behaviours the derivation relies on: ``t=None`` is "off" (every id
kept), a small ``t`` collapses near-duplicate CDR3s, and clustering keeps exactly one
representative per cluster.
"""

from __future__ import annotations

import polars as pl

from tcren.potential import alphabeta_ids, nonredundant_ids


def _markup(rows: list[dict]) -> pl.DataFrame:
    return pl.DataFrame(rows)


def test_t_none_returns_all_ids():
    mk = _markup([
        {"pdb.id": "a", "cdr3a": "CASSL", "cdr3b": "CASSF", "peptide": "GILGFVFTL"},
        {"pdb.id": "b", "cdr3a": "CASSL", "cdr3b": "CASSF", "peptide": "GILGFVFTL"},
        {"pdb.id": "c", "cdr3a": "WWWWW", "cdr3b": "YYYYY", "peptide": "NLVPMVATV"},
    ])
    # t=None => redundancy filtering off => every id returned, in order.
    assert nonredundant_ids(mk, t=None) == ["a", "b", "c"]


def test_small_t_collapses_near_duplicates():
    # a, b are identical (distance 0); c differs by many positions in all three fields.
    mk = _markup([
        {"pdb.id": "a", "cdr3a": "CASSLAPGATNEKLFF", "cdr3b": "CASSIRSSYEQYF", "peptide": "GILGFVFTL"},
        {"pdb.id": "b", "cdr3a": "CASSLAPGATNEKLFF", "cdr3b": "CASSIRSSYEQYF", "peptide": "GILGFVFTL"},
        {"pdb.id": "c", "cdr3a": "CAVRDGGTGNGKLTF", "cdr3b": "CASSPGQGAYEQYF", "peptide": "NLVPMVATV"},
    ])
    keep = nonredundant_ids(mk, t=6.0)
    # The identical pair collapses to one representative; c stays distinct.
    assert len(keep) == 2
    assert "c" in keep
    # exactly one of the duplicate pair survives, and it is the first in row order.
    assert ("a" in keep) ^ ("b" in keep)
    assert "a" in keep


def test_one_representative_per_cluster():
    # Three tight clusters of two near-identical members each => exactly 3 kept.
    mk = _markup([
        {"pdb.id": "a1", "cdr3a": "CASSLAPGATNEKLFF", "cdr3b": "CASSIRSSYEQYF", "peptide": "GILGFVFTL"},
        {"pdb.id": "a2", "cdr3a": "CASSLAPGATNEKLFF", "cdr3b": "CASSIRSSYEQYF", "peptide": "GILGFVFTL"},
        {"pdb.id": "b1", "cdr3a": "CAVRDGGTGNGKLTF", "cdr3b": "CASSPGQGAYEQYF", "peptide": "NLVPMVATV"},
        {"pdb.id": "b2", "cdr3a": "CAVRDGGTGNGKLTF", "cdr3b": "CASSPGQGAYEQYF", "peptide": "NLVPMVATV"},
        {"pdb.id": "c1", "cdr3a": "CAASFGDNSKLIW", "cdr3b": "CASRGTGELFF", "peptide": "ELAGIGILTV"},
        {"pdb.id": "c2", "cdr3a": "CAASFGDNSKLIW", "cdr3b": "CASRGTGELFF", "peptide": "ELAGIGILTV"},
    ])
    keep = nonredundant_ids(mk, t=6.0)
    assert len(keep) == 3
    # one representative from each cluster (the first member in row order)
    assert set(keep) == {"a1", "b1", "c1"}


def test_nulls_are_treated_as_empty_strings():
    mk = _markup([
        {"pdb.id": "a", "cdr3a": None, "cdr3b": "CASSF", "peptide": "GIL"},
        {"pdb.id": "b", "cdr3a": None, "cdr3b": "CASSF", "peptide": "GIL"},
    ])
    # null cdr3a -> "" for both; identical otherwise -> collapse to one.
    assert nonredundant_ids(mk, t=6.0) == ["a"]


def test_single_row_returns_that_row():
    mk = _markup([{"pdb.id": "solo", "cdr3a": "CASS", "cdr3b": "CASS", "peptide": "GIL"}])
    assert nonredundant_ids(mk, t=6.0) == ["solo"]


def test_alphabeta_ids_keeps_only_tra_trb():
    contacts = pl.DataFrame({
        "pdb.id": ["ab", "ab", "gd", "gd", "mix", "mix"],
        "chain.type.from": ["TRA", "TRB", "TRG", "TRD", "TRA", "TRD"],
    })
    keep = alphabeta_ids(contacts)
    assert keep == ["ab"]
