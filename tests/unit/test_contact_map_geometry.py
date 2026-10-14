"""Independent checks of retained orientation and contact-map APIs."""
import numpy as np
import polars as pl
import pytest
from tcren.geometry import orientation_cosine, virtual_cb


def test_orientation_cosine_signs_are_the_way_round():
    """+1 is a side chain aimed at the partner. If this inverts, `g_orient` reads backwards."""
    ca = np.array([[0.0, 0.0, 0.0], [3.8, 0.0, 0.0], [6.0, 3.1, 0.0]])
    cb = virtual_cb(ca)
    direction = cb[1] - ca[1]
    toward = ca[1] + 10.0 * direction / np.linalg.norm(direction)
    away = ca[1] - 10.0 * direction / np.linalg.norm(direction)
    cos = orientation_cosine(ca, np.stack([toward, away]))
    assert cos[1, 0] == pytest.approx(1.0, abs=1e-6)
    assert cos[1, 1] == pytest.approx(-1.0, abs=1e-6)


def test_a_collinear_trace_does_not_produce_nan():
    straight = np.array([[0.0, 0, 0], [3.8, 0, 0], [7.6, 0, 0], [11.4, 0, 0]])
    cb = virtual_cb(straight)
    assert np.isfinite(cb).all()
    assert np.allclose(np.linalg.norm(cb - straight, axis=1), np.linalg.norm(cb[0] - straight[0]))


def test_map_agreement_is_one_when_the_map_matches_exactly():
    """The ceiling is a deterministic map; anything smooth reads below it, by construction."""
    from tcren.potts import map_agreement
    base = dict(**{"pdb.id": "x"})
    d = pl.DataFrame({"pdb.id": ["x"] * 4, "sigma": [1.0, 1.0, 0.0, 0.0],
                      "p_model": [0.5, 0.5, 0.0, 0.0]})
    assert abs(map_agreement(d).row(0, named=True)["js"] - 1.0) < 1e-9
    spread = pl.DataFrame({"pdb.id": ["x"] * 4, "sigma": [1.0, 1.0, 0.0, 0.0],
                           "p_model": [0.25, 0.25, 0.25, 0.25]})
    assert map_agreement(spread).row(0, named=True)["js"] < 0.7
