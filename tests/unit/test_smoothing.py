"""Substitution-matrix pseudocounts for a sparse pair-count matrix.

The properties worth pinning are the ones a wrong inversion or a wrong blend would break: that the
recovered background matches the published BLOSUM62 one, that the prior spreads mass onto
chemically similar residues and not others, that a well-observed cell is left alone while an empty
one becomes its prior, and that the total is conserved so the log-odds downstream is on the same
scale.
"""

from __future__ import annotations

import numpy as np
import polars as pl
import pytest

from tcren.potential import (blosum_background, blosum_conditional, derive_tcren,
                             smooth_counts)
from tcren.potential.model import AA20

IDX = {a: i for i, a in enumerate(AA20)}

#: BLOSUM62 background frequencies as usually quoted, for the residues where the check is sharpest.
#: The recovered values differ by the rounding of the published scores to integers.
PUBLISHED_BACKGROUND = {"A": 0.074, "L": 0.099, "G": 0.074, "W": 0.013, "C": 0.025, "K": 0.058}


def _counts(rows: list[tuple[str, str, float]]) -> pl.DataFrame:
    return pl.DataFrame({"residue.aa.from": [r[0] for r in rows],
                         "residue.aa.to": [r[1] for r in rows],
                         "count": [r[2] for r in rows]})


# --- inverting the matrix ---------------------------------------------------------------


def test_background_is_a_distribution():
    p = blosum_background()
    assert p.shape == (20,)
    assert p.sum() == pytest.approx(1.0)
    assert (p > 0).all()


def test_background_matches_the_published_one():
    p = blosum_background()
    for aa, ref in PUBLISHED_BACKGROUND.items():
        assert abs(p[IDX[aa]] - ref) < 0.02, aa


def test_background_satisfies_its_own_marginal_condition():
    """``sum_b q_ab = p_a`` is what defines it; check the implied ``q`` really has that marginal."""
    from tcren.potential.smoothing import _scores

    p = blosum_background()
    q = np.outer(p, p) * np.exp2(_scores("BLOSUM62") / 2.0)
    q /= q.sum()
    assert np.allclose(q.sum(axis=1), p, atol=1e-9)


def test_conditional_columns_are_distributions():
    p = blosum_conditional()
    assert np.allclose(p.sum(axis=0), 1.0)
    assert (p >= 0).all()


def test_a_residue_is_its_own_most_likely_substitute():
    p = blosum_conditional()
    for aa in AA20:
        assert np.argmax(p[:, IDX[aa]]) == IDX[aa], aa


def test_conservative_substitutions_outrank_radical_ones():
    p = blosum_conditional()
    assert p[IDX["L"], IDX["I"]] > p[IDX["K"], IDX["I"]]
    assert p[IDX["V"], IDX["I"]] > p[IDX["W"], IDX["I"]]
    assert p[IDX["E"], IDX["D"]] > p[IDX["F"], IDX["D"]]


# --- the blend --------------------------------------------------------------------------


def test_beta_zero_is_the_identity_on_the_grid():
    c = _counts([("I", "F", 100.0), ("L", "F", 80.0)])
    out = smooth_counts(c, beta=0.0)
    assert out.height == 400
    nonzero = out.filter(pl.col("count") > 0).sort("count", descending=True)
    assert nonzero["count"].to_list() == [100.0, 80.0]


def test_the_total_is_conserved():
    c = _counts([("I", "F", 100.0), ("L", "F", 80.0), ("W", "D", 3.0)])
    for beta in (1.0, 20.0, 500.0):
        assert smooth_counts(c, beta=beta)["count"].sum() == pytest.approx(183.0)


def test_mass_lands_on_substitutable_neighbours_not_on_unrelated_residues():
    """Ile:Phe and Leu:Phe observations should inform Val:Phe far more than Lys:Phe."""
    out = smooth_counts(_counts([("I", "F", 100.0), ("L", "F", 80.0)]), beta=20.0)

    def cell(a, b):
        return out.filter((pl.col("residue.aa.from") == a)
                          & (pl.col("residue.aa.to") == b))["count"].item()

    assert cell("V", "F") > 5 * cell("K", "F")
    assert cell("M", "F") > cell("K", "F")
    assert cell("K", "F") > cell("W", "F")


def test_the_ordering_of_observed_cells_survives():
    """Smoothing may compress, but it must not reorder what was actually measured."""
    out = smooth_counts(_counts([("I", "F", 100.0), ("L", "F", 80.0), ("K", "E", 60.0)]),
                        beta=20.0)

    def cell(a, b):
        return out.filter((pl.col("residue.aa.from") == a)
                          & (pl.col("residue.aa.to") == b))["count"].item()

    assert cell("I", "F") > cell("L", "F") > cell("K", "E")


def test_more_smoothing_moves_more_mass_off_the_observed_cells():
    """The prior's share must grow with beta, monotonically -- that is the knob's whole meaning."""
    c = _counts([("I", "F", 100.0), ("L", "F", 80.0)])
    observed = []
    for beta in (1.0, 20.0, 200.0):
        out = smooth_counts(c, beta=beta)
        observed.append(out.filter(
            ((pl.col("residue.aa.from") == "I") | (pl.col("residue.aa.from") == "L"))
            & (pl.col("residue.aa.to") == "F"))["count"].sum())
    assert observed[0] > observed[1] > observed[2]


def test_a_residue_with_no_observations_still_gets_a_prior():
    """The whole point: valine is never seen here, and Ile/Leu should speak for it."""
    c = _counts([("I", "F", 100.0), ("L", "F", 80.0)])
    out = smooth_counts(c, beta=20.0)
    assert out.filter((pl.col("residue.aa.from") == "V")
                      & (pl.col("residue.aa.to") == "F"))["count"].item() > 0.0


def test_stronger_smoothing_compresses_the_derived_potential():
    """The blend pulls every cell toward a common prior, so the spread must shrink with beta."""
    c = _counts([("I", "F", 60.0), ("L", "F", 40.0), ("W", "D", 5.0), ("K", "E", 30.0)])
    spreads = [np.ptp(derive_tcren(
        c.with_columns(pl.lit("x").alias("pdb.id")), weight_col="count", smooth_beta=b
    ).matrix["value"].to_numpy()) for b in (0.0, 20.0, 200.0)]
    assert spreads[0] > spreads[1] > spreads[2]


def test_negative_beta_is_rejected():
    with pytest.raises(ValueError, match="non-negative"):
        smooth_counts(_counts([("A", "A", 1.0)]), beta=-1.0)


def test_smoothing_is_rejected_for_the_gapped_alphabet():
    c = _counts([("A", "A", 1.0)]).with_columns(pl.lit("x").alias("pdb.id"))
    with pytest.raises(ValueError, match="classic"):
        derive_tcren(c, variant="am", smooth_beta=10.0)
