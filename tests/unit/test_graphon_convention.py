"""Registration and read-back are two directions of one map, and they have to agree.

They did not. ``_resample`` registered residue ``i`` at grid coordinate ``i(G-1)/(L-1)`` while the
read-back sampled ``(i+1)(G-1)/(L+1)``, so a round trip contracted a loop's footprint toward its
midpoint by ``(L-1)/(L+1)``. These tests pin the old behaviour (so the shipped priors keep
reproducing), pin the new one, and pin the size of the difference.
"""
from __future__ import annotations

import numpy as np
import pytest

from tcren.contactmap import _register
from tcren.geometry import graphon_weights, read_back, resample_graphon

GRIDS = (8, 12, 20, 32)
LENGTHS = (5, 6, 7, 9, 12, 17, 20)


def _ramp(loop_len: int, pep_len: int = 9) -> np.ndarray:
    """A loop whose Cα distance to the peptide falls linearly from 14 Å to 5 Å along the loop."""
    return np.linspace(14.0, 5.0, loop_len)[:, None] * np.ones((1, pep_len))


def _map_coordinates_reference(dmap: np.ndarray, grid: int) -> np.ndarray:
    """What ``_resample`` and ``_register`` computed before they delegated."""
    from scipy.ndimage import map_coordinates
    u, v = np.meshgrid(np.linspace(0, dmap.shape[0] - 1, grid),
                       np.linspace(0, dmap.shape[1] - 1, grid), indexing="ij")
    return map_coordinates(dmap, [u.ravel(), v.ravel()], order=1,
                           mode="nearest").reshape(grid, grid)


@pytest.mark.parametrize("grid", GRIDS)
@pytest.mark.parametrize("loop_len", LENGTHS)
def test_endpoint_reproduces_the_old_resample(grid, loop_len):
    """The delegation must not move `registered_map` or the shipped geom priors."""
    rng = np.random.default_rng(loop_len * 100 + grid)
    dmap = rng.uniform(3.0, 30.0, (loop_len, 11))
    reference = _map_coordinates_reference(dmap, grid)
    # 7e-15 Å is float round-off on a 30 Å distance, and 1e12 times finer than the 0.01 Å the
    # ensemble prior is quantised to on disk.
    assert np.allclose(resample_graphon(dmap, grid, convention="endpoint"), reference, atol=1e-12)
    assert np.allclose(_register(dmap, grid), reference, atol=1e-12)


def test_interior_registration_is_adjoint_to_the_read_back():
    """A round trip under "interior" converges to the identity; under "endpoint" it does not."""
    for loop_len in (6, 12, 17):
        dmap = _ramp(loop_len)
        coarse = np.abs(read_back(resample_graphon(dmap, 12, convention="interior"),
                                  loop_len, 9) - dmap).max()
        fine = np.abs(read_back(resample_graphon(dmap, 256, convention="interior"),
                                loop_len, 9) - dmap).max()
        assert fine < coarse, f"L={loop_len}: interior did not converge with grid"
        assert fine < 0.02, f"L={loop_len}: interior round trip off by {fine:.4f} Å at grid 256"


def test_the_endpoint_bias_does_not_shrink_with_grid():
    """The defect is a wrong coordinate map, not a resolution limit -- so refining cannot fix it."""
    for loop_len in (6, 12, 17):
        dmap = _ramp(loop_len)
        errors = [np.abs(read_back(resample_graphon(dmap, g, convention="endpoint"),
                                   loop_len, 9) - dmap).max() for g in (12, 32, 128, 256)]
        assert max(errors) - min(errors) < 1e-9, f"L={loop_len}: bias moved with grid, {errors}"
        assert errors[0] > 0.4, f"L={loop_len}: expected a bias of ~0.5-1.3 Å, got {errors[0]:.4f}"


def test_the_endpoint_bias_is_worst_for_the_shortest_loop():
    """CDR2 is a 5-7-mer and CDR3 a 12-17-mer, and this is the ordering the loss follows."""
    at_grid_32 = {}
    for loop_len in LENGTHS:
        dmap = _ramp(loop_len)
        err = read_back(resample_graphon(dmap, 32, convention="endpoint"), loop_len, 9) - dmap
        at_grid_32[loop_len] = np.abs(err).max()
        # Zero at the midpoint, largest at the termini -- a contraction, not an offset.
        assert np.abs(err[loop_len // 2, 0]) < 0.3 * at_grid_32[loop_len]
    lengths = sorted(at_grid_32)
    assert all(at_grid_32[a] > at_grid_32[b] for a, b in zip(lengths, lengths[1:])), at_grid_32
    assert at_grid_32[6] > 1.2, "a 6-mer CDR2 should be off by ~1.29 Å"
    assert at_grid_32[17] < 0.6, "a 17-mer CDR3 should be off by ~0.50 Å"


def test_the_contraction_factor_is_l_minus_one_over_l_plus_one():
    """Registering the residue index itself reads the distortion off directly."""
    for loop_len in LENGTHS:
        index = np.arange(loop_len, dtype=float)[:, None] * np.ones((1, 9))
        got = read_back(resample_graphon(index, 32, convention="endpoint"), loop_len, 9)[:, 0]
        scale = (got[-1] - got[0]) / (loop_len - 1)
        assert scale == pytest.approx((loop_len - 1) / (loop_len + 1), abs=1e-3)


@pytest.mark.parametrize("grid", GRIDS)
def test_read_back_carries_a_stack_in_one_call(grid):
    rng = np.random.default_rng(grid)
    stack = rng.uniform(3.0, 30.0, (7, grid, grid))
    one_call = read_back(stack, 13, 9)
    per_member = np.stack([read_back(s, 13, 9) for s in stack])
    assert np.array_equal(one_call, per_member)


def test_graphon_weights_are_a_partition_of_unity():
    for convention in ("interior", "endpoint"):
        for n in LENGTHS:
            w = graphon_weights(n, 32, convention=convention)
            assert w.shape == (32, n)
            assert np.allclose(w.sum(axis=1), 1.0)
            assert (w >= 0).all()


def test_graphon_weights_rejects_a_bad_convention_and_a_degenerate_grid():
    with pytest.raises(ValueError, match="convention must be"):
        graphon_weights(10, 12, convention="midpoint")
    with pytest.raises(ValueError, match="at least one residue"):
        graphon_weights(0, 12)
    with pytest.raises(ValueError, match="needs grid >= 2"):
        graphon_weights(10, 1, convention="interior")
