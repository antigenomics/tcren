"""``epitope_weights`` — one-epitope-one-vote weighting for potential derivation."""

from __future__ import annotations

import polars as pl
import pytest

from tcren.potential import epitope_weights


def _markup(pairs: list[tuple[str, str]]) -> pl.DataFrame:
    return pl.DataFrame({"pdb.id": [p for p, _ in pairs], "peptide": [s for _, s in pairs]})


def test_singleton_peptides_get_weight_one():
    w = epitope_weights(_markup([("1abc", "GILGFVFTL"), ("2def", "NLVPMVATV")]))
    assert w == {"1abc": 1.0, "2def": 1.0}


def test_shared_peptide_splits_one_vote():
    w = epitope_weights(_markup([("a", "P"), ("b", "P"), ("c", "P"), ("d", "Q")]))
    assert w["a"] == w["b"] == w["c"] == pytest.approx(1 / 3)
    assert w["d"] == 1.0


def test_each_distinct_epitope_contributes_total_weight_one():
    m = _markup([("a", "P"), ("b", "P"), ("c", "Q"), ("d", "R"), ("e", "R"), ("f", "R")])
    w = epitope_weights(m)
    by_pep: dict[str, float] = {}
    for pid, pep in m.rows():
        by_pep[pep] = by_pep.get(pep, 0.0) + w[pid]
    assert all(v == pytest.approx(1.0) for v in by_pep.values())
    assert len(by_pep) == 3


def test_null_peptides_are_dropped_not_counted():
    m = pl.DataFrame({"pdb.id": ["a", "b", "c"], "peptide": ["P", None, "P"]})
    w = epitope_weights(m)
    assert set(w) == {"a", "c"}
    assert w["a"] == w["c"] == pytest.approx(0.5)


def test_field_is_configurable():
    m = pl.DataFrame({"pdb.id": ["a", "b"], "antigen.epitope": ["P", "P"]})
    w = epitope_weights(m, field="antigen.epitope")
    assert w == {"a": pytest.approx(0.5), "b": pytest.approx(0.5)}


def test_weighting_changes_the_derived_potential():
    """A peptide crystallized many times must not dominate the counts once weighted."""
    from tcren.potential import derive_tcren

    rows = []
    for i in range(9):                                  # nine copies of one epitope
        rows.append({"pdb.id": f"dup{i}", "residue.aa.from": "W", "residue.aa.to": "K"})
    rows.append({"pdb.id": "solo", "residue.aa.from": "W", "residue.aa.to": "D"})
    contacts = pl.DataFrame(rows)
    markup = _markup([(f"dup{i}", "PPPPPPPPP") for i in range(9)] + [("solo", "DDDDDDDDD")])

    w = epitope_weights(markup)
    assert w["dup0"] == pytest.approx(1 / 9) and w["solo"] == 1.0

    unweighted = derive_tcren(contacts, variant="classic")
    weighted = derive_tcren(contacts, variant="classic", weights=w)

    def cell(pot, a, b):
        m = pot.matrix.filter(
            (pl.col("residue.aa.from") == a) & (pl.col("residue.aa.to") == b))
        return m["value"].item()

    assert cell(unweighted, "W", "K") != pytest.approx(cell(weighted, "W", "K"), abs=1e-9)
