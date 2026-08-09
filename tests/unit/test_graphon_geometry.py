"""Graphon geometry featurisation: reachability, registered maps, binding-mode centroid.

Regression targets are the algem-monograph specs (tcren-code-sync 01/02/03). These are
``structure -> geometry`` descriptors, NOT binder scores, so the tests check the geometric contract
(shapes, invariants, the σ-split, the reachability closed form) — never binder discrimination.
"""
import math
from pathlib import Path

import numpy as np
import pytest

from tcren import binding_mode, geometry, registered_map
from tcren.annotation import classify_chains
from tcren.structure.io import import_structure

PDB = Path(__file__).resolve().parents[1] / "assets" / "pdb" / "1ao7.pdb"


@pytest.fixture(scope="module")
def complex_1ao7():
    pytest.importorskip("arda")  # classify_chains needs the arda backend; CI installs tcren without it
    s = import_structure(str(PDB))
    classify_chains(s, organism="human")
    return s


# --- 03 reachability (pure geometry, no structure) -------------------------------------------
def test_reach_max_closed_form():
    assert geometry.reach_max(13, 5.25) == pytest.approx(26.47, abs=0.1)   # spec ~26.4


def test_reachability_floor_inverts_reach_max():
    assert geometry.reachability_floor(8.1, 5.25) == 4                     # spec target
    d, R = 8.1, 5.25                                                       # floor L reaches, L-1 does not
    L = geometry.reachability_floor(d, R)
    assert geometry.reach_max(L, R) >= d > geometry.reach_max(L - 1, R)


def test_reach_max_zero_when_loop_too_short():
    assert geometry.reach_max(0, 100.0) == 0.0


def test_span_saturation_in_unit_range(complex_1ao7):
    sat = geometry.span_saturation(complex_1ao7)
    assert set(sat) == {"cdr3a", "cdr3b"}
    for v in sat.values():                                                # real loops run well below 1
        assert 0.0 < v < 1.0


# --- 01 registered map ------------------------------------------------------------------------
def test_registered_map_shape_and_grid(complex_1ao7):
    for g in (4, 8, 12):
        m = registered_map(complex_1ao7, grid=g)
        assert m.shape == (2, g, g)                                        # both loops resolve for 1ao7


def test_registered_map_contact_mode_is_binary(complex_1ao7):
    m = registered_map(complex_1ao7, grid=8, metric="contact", cutoff=5.0)
    assert set(np.unique(m)) <= {0.0, 1.0}


def test_registered_map_rejects_unknown_metric(complex_1ao7):
    with pytest.raises(ValueError, match="metric"):
        registered_map(complex_1ao7, metric="bogus")


# --- 02 binding-mode centroid -----------------------------------------------------------------
def test_binding_mode_hits_the_apex_and_sigma_involution(complex_1ao7):
    m = binding_mode(complex_1ao7, contact=8.0)
    assert m.apex_x == pytest.approx(0.47, abs=0.05)                       # the loop reads from its apex
    assert m.y_alpha < m.y_beta                                            # alpha N-terminal, beta C-terminal
    assert m.sigma_sum == pytest.approx(1.0, abs=0.1)                      # y_alpha + y_beta ~ 1
    assert m.n_contacts_alpha >= 3 and m.n_contacts_beta >= 3


def test_binding_mode_keeps_loops_separate(complex_1ao7):
    m = binding_mode(complex_1ao7, contact=8.0)
    assert not math.isclose(m.y_alpha, m.y_beta)                           # never pooled — the split is the point


def test_binding_mode_default_cutoff_resolves_a_crystal(complex_1ao7):
    # The cutoff is Cα–Cα, so a 5 Å default returned None on every real complex.
    m = binding_mode(complex_1ao7)          # no explicit cutoff — the shipped default must work
    assert m is not None and 0.3 < m.apex_x < 0.7


# --- 04 CDR3 internal coordinates -------------------------------------------------------------
def test_cdr3_internal_coords_lengths_and_physics(complex_1ao7):
    ic = geometry.cdr3_internal_coords(complex_1ao7, "cdr3b")
    n = ic.n_ca
    assert (len(ic.bonds), len(ic.angles), len(ic.torsions)) == (n - 1, n - 2, n - 3)
    assert np.allclose(ic.bonds, 3.8, atol=0.3)                        # Cα–Cα virtual bond ~3.8 Å
    assert np.all((ic.angles > np.radians(80)) & (ic.angles < np.radians(150)))  # physical Cα angle
    assert np.all((ic.torsions >= -np.pi) & (ic.torsions <= np.pi))


def test_cdr3_internal_coords_neck_matches_endpoints(complex_1ao7):
    from tcren.contactmap import _cdr3_ca
    ic = geometry.cdr3_internal_coords(complex_1ao7, "cdr3b")
    ca = _cdr3_ca(complex_1ao7, "TRB")
    assert ic.neck == pytest.approx(float(np.linalg.norm(ca[-1] - ca[0])), abs=1e-9)


def test_cdr3_internal_coords_rejects_unknown_loop(complex_1ao7):
    with pytest.raises(ValueError, match="cdr3a"):
        geometry.cdr3_internal_coords(complex_1ao7, "cdr3x")
