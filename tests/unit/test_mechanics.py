"""Unit tests for interface mechanics (stiffness tensor + rupture + coupling).

Follows the ``test_ddg.py`` style: tiny hand-built fixtures with analytically
checkable geometry. The tensor/rupture checks build :class:`InterfaceSprings`
directly (no PDB), the coupling check exercises the pure set-intersection helper,
and a single integration smoke test parses a bundled PDB (1mi5) and runs the three
public entry points end-to-end.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from tcren.mechanics import (
    InterfaceSprings,
    _coupling_counts,
    coupling_residues,
    rupture,
    stiffness_tensor,
)

REPO = Path(__file__).resolve().parents[2]


def _springs(a, b, k) -> InterfaceSprings:
    """Build an InterfaceSprings from anchors + stiffnesses (mirrors module internals)."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    k = np.asarray(k, float)
    rest = np.linalg.norm(b - a, axis=1)
    ct = np.average(a, 0, weights=k)
    cp = np.average(b, 0, weights=k)
    axis = cp - ct
    axis = axis / (np.linalg.norm(axis) + 1e-9)
    return InterfaceSprings(a, b, k, rest, axis)


# --- (1) stiffness tensor sign structure ------------------------------------

def _z_springs(n_extra=0, k_val=1.0) -> InterfaceSprings:
    """n>=3 springs purely along +z (rest length 5), spread in the xy-plane."""
    base = [[0, 0, 0], [3, 0, 0], [0, 3, 0]]
    extra = [[6, 0, 0], [0, 6, 0], [6, 6, 0]]
    a = base + extra[:n_extra]
    b = [[x, y, z + 5] for x, y, z in a]
    k = [k_val] * len(a)
    return _springs(a, b, k)


def test_stiffness_tensor_z_springs_are_tensile_dominant():
    """Springs aligned with the docking axis put all stiffness into K_tens."""
    zs = _z_springs()
    # Feed the pre-built network through the same math stiffness_tensor uses.
    u = (zs.b - zs.a) / zs.rest[:, None]
    K = (zs.k[:, None, None] * u[:, :, None] * u[:, None, :]).sum(0)
    S_tot = float(np.trace(K))
    K_tens = float(zs.axis @ K @ zs.axis)
    K_shear = S_tot - K_tens
    assert K_tens > K_shear
    # z-aligned unit springs: K is diag(0, 0, n) so shear is ~0 and tensile ~= trace.
    assert K_tens == pytest.approx(len(zs))
    assert K_shear == pytest.approx(0.0, abs=1e-6)


def test_stiffness_tensor_under_three_springs_is_nan(monkeypatch):
    """< 3 springs cannot define the tensor; every descriptor is nan, n_spring is exact."""
    import tcren.mechanics as mech

    two = _springs([[0, 0, 0], [1, 0, 0]], [[0, 0, 5], [1, 0, 5]], [1, 1])
    monkeypatch.setattr(mech, "interface_springs", lambda *a, **k: two)
    out = mech.stiffness_tensor(object())
    for key in ("S_tot", "K_tens", "K_shear", "aniso", "lam_max", "lam_min"):
        assert np.isnan(out[key])
    assert out["n_spring"] == 2.0


# --- (2) rupture monotonicity -----------------------------------------------

def _rupture_force(springs, direction="tensile", break_strain=0.5, steps=80) -> float:
    """Rupture peak force on a pre-built network (mirrors rupture()'s inner loop)."""
    from tcren.mechanics import _pull_direction

    d = _pull_direction(springs, direction)
    max_disp = 2.0 * float(springs.rest.max())
    dt = max_disp / steps
    alive = np.ones(len(springs), bool)
    best = 0.0
    for i in range(steps):
        t = (i + 1) * dt
        v = (springs.b + t * d) - springs.a
        L = np.linalg.norm(v, axis=1)
        alive &= (L - springs.rest) / springs.rest <= break_strain
        ext = np.maximum(L - springs.rest, 0.0)
        f = springs.k * ext * (v @ d) / (L + 1e-9)
        best = max(best, float((f * alive).sum()))
    return best


def test_rupture_force_monotone_in_spring_count():
    """More parallel springs resist more strongly."""
    f3 = _rupture_force(_z_springs(n_extra=0))
    f6 = _rupture_force(_z_springs(n_extra=3))
    assert f6 > f3


def test_rupture_force_monotone_in_stiffness():
    """Scaling every stiffness by 5x scales the peak force by 5x."""
    f1 = _rupture_force(_z_springs(k_val=1.0))
    f5 = _rupture_force(_z_springs(k_val=5.0))
    assert f5 == pytest.approx(5.0 * f1)
    assert f5 > f1


def test_rupture_force_monotone_in_break_strain():
    """A lower break strain ruptures sooner, so the peak force is <= a higher one."""
    zs = _z_springs()
    f_low = _rupture_force(zs, break_strain=0.2)
    f_high = _rupture_force(zs, break_strain=0.8)
    assert f_low <= f_high + 1e-9


# --- (3) coupling set-intersection correctness ------------------------------

def test_coupling_counts_set_intersections():
    """Mirror the assertion baked into mechanics._demo."""
    c = _coupling_counts(
        pep_tcr={1, 2, 3}, pep_mhc={2, 3, 4}, mhc_tcr={5}, mhc_pep={5, 6},
        tcr_pmhc={1, 2, 7}, tcr_ab={2, 8},
    )
    # couple_pep = pep_tcr & pep_mhc = {2,3} -> 2
    # couple_mhc = mhc_tcr & mhc_pep = {5}   -> 1
    # couple_tcr = tcr_pmhc & tcr_ab = {2}   -> 1
    assert (c["couple_pep"], c["couple_mhc"], c["couple_tcr"]) == (2, 1, 1)
    assert c["couple_total"] == 4
    # n_interface = |tcr_pmhc| + |pep_tcr| + |mhc_tcr| = 3 + 3 + 1 = 7
    assert c["n_interface"] == 7


def test_coupling_counts_disjoint_sets_are_zero():
    """No overlap anywhere -> all couplings zero, but n_interface still counts sizes."""
    c = _coupling_counts(
        pep_tcr={1}, pep_mhc={2}, mhc_tcr={3}, mhc_pep={4},
        tcr_pmhc={5, 6}, tcr_ab={7},
    )
    assert c["couple_pep"] == c["couple_mhc"] == c["couple_tcr"] == 0
    assert c["couple_total"] == 0
    assert c["n_interface"] == len({5, 6}) + len({1}) + len({3}) == 4


# --- (4) integration smoke test ---------------------------------------------

@pytest.mark.slow
def test_integration_1mi5_public_functions():
    """Parse + classify + annotate a bundled complex; run the three public entry points."""
    pytest.importorskip("arda")  # classify_chains is mmseqs-backed
    from tcren.annotation import classify_chains
    from tcren.mhc import annotate_mhc
    from tcren.structure import parse_structure

    pdb = REPO / "data" / "Canonical2026" / "1mi5.pdb.gz"
    if not pdb.exists():
        pytest.skip(f"bundled structure not found: {pdb}")

    s = parse_structure(pdb, pdb_id="1mi5")
    classify_chains(s, organism="human")
    annotate_mhc(s)

    st = stiffness_tensor(s)
    assert np.isfinite(st["S_tot"]) and np.isfinite(st["K_tens"])
    assert st["K_tens"] > 0.0
    assert st["n_spring"] > 3

    r = rupture(s)
    assert np.isfinite(r["rupture_force"])
    assert r["rupture_force"] > 0.0
    assert r["n_spring"] > 3

    c = coupling_residues(s)
    assert c["couple_pep"] >= 0
    assert c["couple_total"] >= 0
    assert c["n_interface"] > 0
