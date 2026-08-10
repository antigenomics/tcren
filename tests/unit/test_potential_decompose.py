"""Unit tests for the one-body / pair split of a contact potential.

A contact energy is not purely an interaction: burying a residue against any partner
carries a cost that depends on that residue alone. ``Potential.decompose`` separates the
two, and the reason to care is scoring, not bookkeeping -- an additive per-position model
can absorb the one-body part, so only the pair part is genuinely beyond it.
"""

from __future__ import annotations

import numpy as np
import pytest

from tcren.potential import Potential, mj, tcren

AA = "ACDEFGHIKLMNPQRSTVWY"


def test_split_is_exact():
    pot = mj()
    d = pot.decompose()
    for a in AA:
        for b in AA:
            assert d.energy(a, b) == pytest.approx(pot.value(a, b), abs=1e-12)


def test_pair_part_has_zero_marginals():
    """This is what makes the split unique, and what makes J the non-additive remainder."""
    d = mj().decompose()
    assert np.allclose(d.pair.sum(axis=0), 0.0, atol=1e-12)
    assert np.allclose(d.pair.sum(axis=1), 0.0, atol=1e-12)


def test_pair_part_is_symmetric_like_its_source():
    d = mj().decompose()
    assert np.allclose(d.pair, d.pair.T)


def test_one_body_is_the_centred_row_mean():
    pot = mj()
    dense, index = pot.as_matrix()
    d = pot.decompose()
    expected = dense.mean(axis=1) - dense.mean()
    assert np.allclose(d.one_body, expected)


def test_a_directed_potential_is_refused():
    """TCRen is TCR-to-peptide directed; splitting it this way would be meaningless."""
    with pytest.raises(ValueError, match="not symmetric"):
        tcren().decompose()


def test_a_purely_additive_potential_has_no_pair_part():
    """If e(a,b) = f(a) + f(b) exactly, every pair term must vanish."""
    f = {a: i * 0.1 for i, a in enumerate(AA)}
    rows = [{"residue.aa.from": a, "residue.aa.to": b, "value": f[a] + f[b]}
            for a in AA for b in AA]
    import polars as pl
    pot = Potential(name="additive", matrix=pl.DataFrame(rows), alphabet=tuple(AA))
    d = pot.decompose()
    assert np.allclose(d.pair, 0.0, atol=1e-12)


def test_hydrophobicity_fit_reproduces_the_mj_matrix():
    """Li, Tang and Wingreen's result, on the bundled matrix: three parameters and one
    number per residue account for most of it."""
    fit = mj().hydrophobicity_fit()
    assert fit.r2 > 0.8
    assert fit.eigenvalue_share > 0.4
    # q orders by hydrophobicity: the aliphatics and aromatics on one end, charges on the
    # other. This is the property that makes the fit interpretable rather than numerology.
    assert fit.q[fit.index["F"]] > fit.q[fit.index["K"]]
    assert fit.q[fit.index["L"]] > fit.q[fit.index["D"]]
    assert fit.q[fit.index["I"]] > fit.q[fit.index["E"]]


def test_hydrophobicity_fit_value_matches_its_own_coefficients():
    fit = mj().hydrophobicity_fit()
    qa, qb = fit.q[fit.index["A"]], fit.q[fit.index["W"]]
    assert fit.value("A", "W") == pytest.approx(
        fit.c0 + fit.c1 * (qa + qb) + fit.c2 * qa * qb
    )
    assert fit.one_body("W") == pytest.approx(fit.c1 * qb)


def test_hydrophobicity_fit_refuses_a_directed_potential():
    with pytest.raises(ValueError, match="not symmetric"):
        tcren().hydrophobicity_fit()
