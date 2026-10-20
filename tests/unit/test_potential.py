"""Unit tests for potential derivation and the Potential container."""

from __future__ import annotations

import numpy as np
import polars as pl
import pytest

from tcren.potential import Potential, betancourt, derive_tcren, keskin, mj, mj1996, tcren
from tcren.potential.model import AA20


def _independent_classic(contacts: pl.DataFrame, pseudocount: int = 1) -> dict[tuple[str, str], float]:
    """Reference re-implementation of the classic TCRen log-odds, in plain numpy.

    Pins the orientation of total.from (sum over 'to' for fixed 'from') and total.to
    (sum over 'from' for fixed 'to'), independently of the polars implementation.
    """
    aa = list(AA20)
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


def test_betancourt_is_thr_referenced_and_symmetric():
    """The B matrix's defining property: Thr is the reference solvent, so its row is zero.

    A transcription slip anywhere in the lower triangle would break one of these three, which is
    why the bundled file is parsed from the AAindex record rather than typed.
    """
    b = betancourt()
    assert b.name == "BT1999"
    dense, index = b.as_matrix()
    dense = np.asarray(dense, float)
    assert dense.shape == (20, 20)
    assert np.allclose(dense, dense.T)
    assert np.all(dense[index["T"], :] == 0.0)
    assert dense.min() == pytest.approx(-1.34) and dense.max() == pytest.approx(0.66)


def test_reference_states_pair_like_with_like():
    """MJ and BT are pair-contact matrices; Keskin and MJ1996 are raw contact energies.

    Comparing across the two groups compares reference states, not derivations -- the trap the
    SOURCES entry records. Guarded here so a future matrix swap cannot cross the groups silently.
    """
    pair = [np.asarray(f().as_matrix()[0], float) for f in (mj, betancourt)]
    raw = [np.asarray(f().as_matrix()[0], float) for f in (keskin, mj1996)]
    assert all(-0.5 < m.mean() < 0.5 and m.max() > 0 for m in pair)
    assert all(m.mean() < -2.0 and m.max() < 0 for m in raw)


def test_value_missing_pair_raises():
    pot = tcren()
    with pytest.raises(KeyError):
        pot.value("C", "A")  # Cys dropped from the 'from' axis in classic


def test_directed_matrix_moments_use_every_finite_cell():
    pot = Potential("directed", pl.DataFrame({
        "residue.aa.from": ["A", "A", "K", "K"],
        "residue.aa.to": ["A", "K", "A", "K"],
        "value": [1., 2., 8., float("nan")],
    }), ("A", "K"))
    assert pot.offset() == pytest.approx(np.mean([1., 2., 8.]))
    assert pot.scale() == pytest.approx(np.std([1., 2., 8.]))


# --- contact-type-conditioned derivation ----------------------------------------------------------
def test_derive_by_type_needs_the_type_column():
    from tcren.potential.derive import derive_tcren_by_type

    df = pl.DataFrame({"pdb.id": ["a"], "residue.aa.from": ["L"], "residue.aa.to": ["V"]})
    with pytest.raises(ValueError, match="contact.type"):
        derive_tcren_by_type(df)


def test_derive_by_type_reports_the_sparsity_it_creates():
    """The report is the point: splitting a fixed contact set across types is what makes it sparse."""
    from tcren.potential.derive import derive_tcren_by_type

    rng = np.random.default_rng(0)
    aas = list("LFIMVWYHAGPTSQNDERK")
    n = 600
    df = pl.DataFrame({
        "pdb.id": ["s%d" % (i % 20) for i in range(n)],
        "residue.aa.from": rng.choice(aas, n).tolist(),
        "residue.aa.to": rng.choice(aas, n).tolist(),
        "contact.type": rng.choice(["polar", "hydrogen_bond", "salt_bridge"], n,
                                   p=[0.7, 0.25, 0.05]).tolist(),
    })
    pots, report = derive_tcren_by_type(df, min_count=30)

    assert set(pots) == {"polar", "hydrogen_bond", "salt_bridge"}
    assert report["n_contacts"].sum() == n
    assert report["contact.type"].to_list()[0] == "polar"        # sorted by count, descending
    # Every type must be scoreable even where it is sparse (pseudocount fills the empty cells)...
    for p in pots.values():
        assert p.as_matrix()[0].shape == (20, 20)   # full alphabet, holes filled by pseudocount
        assert p.name.startswith("TCRen[")
    # ...and the report must say so rather than let the matrix imply otherwise.
    assert (report["frac_cells_ge_min"] < 0.1).all()


def test_derive_by_type_partitions_the_contacts():
    """Each contact goes to exactly one type, so the per-type counts sum to the pooled count."""
    from tcren.potential.derive import derive_tcren_by_type

    rng = np.random.default_rng(1)
    aas = list("LFIVAG")
    n = 200
    df = pl.DataFrame({
        "pdb.id": ["s"] * n,
        "residue.aa.from": rng.choice(aas, n).tolist(),
        "residue.aa.to": rng.choice(aas, n).tolist(),
        "contact.type": rng.choice(["polar", "hydrophobic"], n).tolist(),
    })
    _pots, report = derive_tcren_by_type(df)
    assert report["n_contacts"].sum() == df.height
