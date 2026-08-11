"""Unit tests for the Miyazawa--Jernigan effective partition energies (AAindex MIYS850101).

A pairwise matrix cannot supply a one-body term, so this scale is bundled rather than
derived. The test that matters is the last one: it checks the bundled numbers against a
hydrophobicity axis recovered from a completely different file, which is what makes a
transcription error detectable rather than invisible.
"""

from __future__ import annotations

import numpy as np
import pytest

from tcren.potential import mj, mj1996, mj_partition_energy

AA = "ACDEFGHIKLMNPQRSTVWY"


def test_covers_the_twenty_residues():
    scale = mj_partition_energy()
    assert set(scale) == set(AA)
    assert all(isinstance(v, float) for v in scale.values())


def test_landmark_values():
    scale = mj_partition_energy()
    assert scale["F"] == pytest.approx(4.37)
    assert scale["K"] == pytest.approx(1.23)
    assert scale["A"] == pytest.approx(2.36)


def test_it_orders_as_hydrophobicity():
    """Larger is more hydrophobic here -- the opposite sign convention to a contact energy."""
    scale = mj_partition_energy()
    for buried in "FMILWV":
        for exposed in "KDNEQ":
            assert scale[buried] > scale[exposed]


def test_it_agrees_with_the_axis_recovered_from_the_1996_matrix():
    """Independent check on the transcription: this scale and the leading hydrophobicity
    axis of the 1996 contact matrix come from different sources and must still agree."""
    scale = mj_partition_energy()
    order = sorted(AA)
    values = np.array([scale[a] for a in order])
    for potential in (mj1996(), mj()):
        fit = potential.hydrophobicity_fit()
        q = np.array([fit.q[fit.index[a]] for a in order])
        assert np.corrcoef(values, q)[0, 1] > 0.9
