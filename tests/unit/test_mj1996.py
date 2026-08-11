"""Unit tests for the Miyazawa--Jernigan 1996 contact energies.

The bundled ``MJ`` matrix arrived without a recorded source, and it turns out not to be raw
contact energies at all: it takes both signs. These tests pin the properties that tell the
two apart, so a future edit cannot quietly swap one for the other.
"""

from __future__ import annotations

import numpy as np
import pytest

from tcren.potential import mj, mj1996

AA = "ACDEFGHIKLMNPQRSTVWY"


def test_it_is_a_raw_contact_matrix():
    """Every published contact energy in this table is attractive."""
    dense, _ = mj1996().as_matrix()
    assert (dense < 0).all()
    assert dense.min() == pytest.approx(-7.37)
    assert dense.max() == pytest.approx(-0.12)


def test_landmark_values_match_the_published_table():
    pot = mj1996()
    assert pot.value("A", "A") == pytest.approx(-2.72)
    assert pot.value("L", "L") == pytest.approx(-7.37)


def test_it_is_symmetric_and_complete():
    dense, index = mj1996().as_matrix()
    assert dense.shape == (20, 20)
    assert set(index) == set(AA)
    assert np.allclose(dense, dense.T)
    assert not np.isnan(dense).any()


def test_it_is_distinguishable_from_the_bundled_mj_matrix():
    """The property that identifies each: one is all-attractive, the other is not."""
    raw, _ = mj1996().as_matrix()
    bundled, _ = mj().as_matrix()
    assert (raw < 0).all()
    assert (bundled > 0).any()
    assert raw.mean() < -2.0
    assert abs(bundled.mean()) < 0.5


def test_li_tang_wingreen_holds_more_strongly_on_the_raw_matrix():
    """Their result was derived on a matrix of this kind, and it shows: three parameters
    and one number per residue reproduce almost all of it."""
    fit = mj1996().hydrophobicity_fit()
    assert fit.r2 > 0.95
    assert fit.eigenvalue_share > 0.8
    assert fit.q[fit.index["F"]] > fit.q[fit.index["K"]]
