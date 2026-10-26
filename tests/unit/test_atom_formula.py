"""Small correctness checks for the explicit atom-formula exploration."""

import numpy as np
import pytest
from sklearn.metrics import average_precision_score
from tcren.contacts.atom_formula import align30, bounded, ap_columns


def test_alignment_bounding_and_exact_tied_ap():
    assert align30([1, 3, 4, 10], [10] * 4).tolist() == [1, 3, 24, 30]
    with pytest.raises(ValueError):
        align30([4], [31])
    for kind in ["rational", "exp", "hill2"]:
        x = bounded(np.arange(100), 4, kind)
        assert x[0] == 0 and (np.diff(x) >= 0).all() and ((x >= 0) & (x <= 1)).all()
    y = np.array([1, 0, 1, 0, 1])
    s = np.array([[1, 0], [1, 1], [0, 1], [0, 0], [1, 1]])
    np.testing.assert_allclose(
        ap_columns(y, s), [average_precision_score(y, s[:, i]) for i in range(2)]
    )
    np.testing.assert_allclose(ap_columns(y, s), ap_columns(y, 3 * s + 7))


