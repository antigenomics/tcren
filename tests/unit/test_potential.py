"""Unit tests for potential derivation and the Potential container."""

from __future__ import annotations

import numpy as np
import polars as pl
import pytest

from tcren.potential import Potential, derive_tcren, keskin, mj, tcren
from tcren.potential.derive import _AA20_CLASSIC


def _independent_classic(contacts: pl.DataFrame, pseudocount: int = 1) -> dict[tuple[str, str], float]:
    """Reference re-implementation of the classic TCRen log-odds, in plain numpy.

    Pins the orientation of total.from (sum over 'to' for fixed 'from') and total.to
    (sum over 'from' for fixed 'to'), independently of the polars implementation.
    """
    aa = list(_AA20_CLASSIC)
    idx = {a: i for i, a in enumerate(aa)}
    n = len(aa)
    counts = np.zeros((n, n))
    for row in contacts.iter_rows(named=True):
        counts[idx[row["residue.aa.from"]], idx[row["residue.aa.to"]]] += 1
    counts += pseudocount
    total_from = counts.sum(axis=1, keepdims=True)  # per 'from' row
    total_to = counts.sum(axis=0, keepdims=True)  # per 'to' column
    total = counts.sum()
    val = -np.log(counts * total / total_to / total_from)
    out = {}
    for i, fr in enumerate(aa):
        if fr == "C":
            continue
        for j, to in enumerate(aa):
            out[(fr, to)] = val[i, j]
    return out


def test_classic_formula_orientation():
    # Asymmetric synthetic contacts so a from/to transpose would be caught.
    rows = (
        [{"pdb.id": "x", "residue.aa.from": "A", "residue.aa.to": "D"}] * 5
        + [{"pdb.id": "x", "residue.aa.from": "L", "residue.aa.to": "A"}] * 3
        + [{"pdb.id": "x", "residue.aa.from": "K", "residue.aa.to": "E"}] * 2
    )
    contacts = pl.DataFrame(rows)
    pot = derive_tcren(contacts, variant="classic")
    ref = _independent_classic(contacts)

    for row in pot.matrix.iter_rows(named=True):
        key = (row["residue.aa.from"], row["residue.aa.to"])
        assert row["value"] == pytest.approx(ref[key], abs=1e-12)

    # Orientation sanity: the (A->D) cell differs from (D->A).
    assert pot.value("A", "D") != pytest.approx(pot.value("D", "A"))


def test_classic_excludes_cys_from_only():
    contacts = pl.DataFrame(
        [{"pdb.id": "x", "residue.aa.from": "A", "residue.aa.to": "D"}]
    )
    pot = derive_tcren(contacts, variant="classic")
    froms = set(pot.matrix["residue.aa.from"])
    tos = set(pot.matrix["residue.aa.to"])
    assert "C" not in froms
    assert "C" in tos


def test_as_matrix_indexing():
    pot = tcren()
    dense, index = pot.as_matrix()
    assert dense.shape == (len(pot.alphabet), len(pot.alphabet))
    fr, to = "A", "D"
    assert dense[index[fr], index[to]] == pytest.approx(pot.value(fr, to))


def test_bundled_loaders_distinct():
    p_mj, p_kes = mj(), keskin()
    assert p_mj.name == "MJ"
    assert p_kes.name == "Keskin"
    # Different potentials → different values for at least some pair.
    assert p_mj.value("A", "A") != pytest.approx(p_kes.value("A", "A"))


def test_value_missing_pair_raises():
    pot = tcren()
    with pytest.raises(KeyError):
        pot.value("C", "A")  # Cys dropped from the 'from' axis in classic
