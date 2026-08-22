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


# --- balanced_weights: several redundancy axes at once ------------------------------


def _mk(rows: list[tuple[str, str, str, str]]) -> pl.DataFrame:
    return pl.DataFrame(
        {"pdb.id": [r[0] for r in rows], "peptide": [r[1] for r in rows],
         "cdr3a": [r[2] for r in rows], "cdr3b": [r[3] for r in rows]}
    )


def test_single_axis_is_exactly_epitope_weights():
    from tcren.potential import balanced_weights
    m = _mk([("a", "P", "X", "X"), ("b", "P", "Y", "Y"), ("c", "Q", "Z", "Z")])
    assert balanced_weights(m, axes=(("peptide",),)) == epitope_weights(m)


def test_novel_receptor_on_a_common_epitope_keeps_most_of_its_weight():
    """The case the mean exists for: 1/9 under a product rule, 0.556 under the mean."""
    from tcren.potential import balanced_weights
    rows = [(f"dup{i}", "P", f"A{i}", f"B{i}") for i in range(8)]
    rows.append(("novel", "P", "UNIQUE_A", "UNIQUE_B"))
    w = balanced_weights(_mk(rows))
    assert w["novel"] == pytest.approx((1 / 9 + 1 / 1) / 2)
    assert w["novel"] > 0.5


def test_true_resolve_duplicated_on_every_axis_gets_one_over_n():
    from tcren.potential import balanced_weights
    rows = [(f"s{i}", "P", "A", "B") for i in range(4)]
    w = balanced_weights(_mk(rows))
    assert all(v == pytest.approx(0.25) for v in w.values())


def test_unique_on_every_axis_gets_full_weight():
    from tcren.potential import balanced_weights
    w = balanced_weights(_mk([("a", "P", "A", "B"), ("b", "Q", "C", "D")]))
    assert w == {"a": pytest.approx(1.0), "b": pytest.approx(1.0)}


def test_receptor_axis_is_keyed_on_both_cdr3_loops_jointly():
    """Sharing only cdr3a is not the same receptor."""
    from tcren.potential import balanced_weights
    m = _mk([("a", "P", "A", "B1"), ("b", "Q", "A", "B2")])
    w = balanced_weights(m, axes=(("cdr3a", "cdr3b"),))
    assert w == {"a": pytest.approx(1.0), "b": pytest.approx(1.0)}
    shared = _mk([("a", "P", "A", "B"), ("b", "Q", "A", "B")])
    w2 = balanced_weights(shared, axes=(("cdr3a", "cdr3b"),))
    assert w2 == {"a": pytest.approx(0.5), "b": pytest.approx(0.5)}


def test_axis_order_does_not_matter():
    from tcren.potential import balanced_weights
    m = _mk([("a", "P", "A", "B"), ("b", "P", "C", "D"), ("c", "Q", "A", "B")])
    f = balanced_weights(m, axes=(("peptide",), ("cdr3a", "cdr3b")))
    r = balanced_weights(m, axes=(("cdr3a", "cdr3b"), ("peptide",)))
    assert f == pytest.approx(r)


def test_weight_is_bounded_by_the_per_axis_extremes():
    from tcren.potential import balanced_weights
    m = _mk([("a", "P", "A", "B"), ("b", "P", "C", "D"), ("c", "P", "A", "B")])
    w = balanced_weights(m)
    assert all(0.0 < v <= 1.0 for v in w.values())
