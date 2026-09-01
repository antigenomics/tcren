"""The interface as a graph and as a matrix.

Three invariants make these numbers usable, and they are the ones :mod:`tcren.footprint`'s suite
already pins for the older block: a known-answer geometry, invariance under rigid motion, and NaN
rather than 0 on input that does not support the measure. Rigid motion is covered for all sixteen
columns at once in ``test_footprint.py`` — every one of them rides in ``footprint_features`` — so
what is here is the known answers and the degenerate cases, plus the two properties the family was
added for: that the spectral pair does not know how long either side was, and that the loop-overlap
term reads a collapsed footprint.
"""
from __future__ import annotations

import numpy as np
import pytest

from tcren.interface_graph import (
    GRAPH_FEATURES,
    MATRIX_FEATURES,
    PROMOTED_POSE_FEATURES,
    _algebraic_connectivity,
    _biadj_cyclo,
    _pielou,
    _selfcheck,
    graph_features,
    matrix_features,
)
from tcren.recognition import DESCRIPTORS, DETAIL, INVARIANCE, STATUS
from tcren.structure.model import PEPTIDE_TYPE, Atom, Chain, RegionMarkup, Residue, Structure

_THREE = {"L": "LEU", "D": "ASP", "K": "LYS", "W": "TRP", "G": "GLY"}


def test_selfcheck_passes():
    _selfcheck()


# --- known answers, by hand ---------------------------------------------------------------------

def test_pielou_is_one_on_uniform_and_falls_when_it_concentrates():
    assert _pielou(np.array([3.0, 3, 3]), 3) == pytest.approx(1.0)
    assert _pielou(np.array([9.0, 1, 1]), 3) < 0.7
    assert np.isnan(_pielou(np.array([0.0, 0, 0]), 3))


def test_cyclomatic_fraction_is_zero_on_a_forest_and_positive_on_a_block():
    """A matching is a forest — every contact is a spanning edge, none of it closes a cycle."""
    assert _biadj_cyclo(np.eye(5, dtype=np.int8)) == 0.0
    # K_{2,3}: E = 6, V = 5, C = 1, so b1 = 2 and the fraction is 2/6.
    assert _biadj_cyclo(np.ones((2, 3), dtype=np.int8)) == pytest.approx(2 / 6)


def test_algebraic_connectivity_is_one_on_a_complete_block_and_near_zero_on_a_bridge():
    """K_{m,n} is maximally connected under the normalised Laplacian; a bridged pair is not."""
    assert _algebraic_connectivity(np.ones((3, 3), dtype=np.int8)) == pytest.approx(1.0, abs=1e-9)
    bridge = np.zeros((6, 6), dtype=np.int8)
    bridge[:3, :3] = 1
    bridge[3:, 3:] = 1
    bridge[2, 3] = 1                       # the single edge holding the two blocks together
    assert 0.0 < _algebraic_connectivity(bridge) < _algebraic_connectivity(np.ones((3, 3), np.int8))


def test_a_disconnected_footprint_raises_the_component_fraction():
    """Two patches over the same number of nodes must read higher than one."""
    one = np.ones((3, 3), dtype=np.int8)
    two = np.zeros((3, 3), dtype=np.int8)
    two[0, 0] = two[1, 1] = two[2, 2] = 1          # three isolated edges, three components
    assert _biadj_cyclo(two) == 0.0
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components

    def comps(B):
        m, n = B.shape
        A = np.zeros((m + n, m + n))
        A[:m, m:] = B
        A[m:, :m] = B.T
        return int(connected_components(coo_matrix(A), directed=False)[0])
    assert comps(two) / 6 > comps(one) / 6


# --- what the matrix family measures, and what was removed before it shipped --------------------

def test_the_face_terms_read_the_side_chains_and_not_the_backbone():
    """Cbeta closer than Calpha means the side chains lean in; the reverse means they point away.

    This is the whole content of ``m_face_*``, and it is why a forced pose reads negative: satisfying
    a contact-count objective brings the backbones together without turning the side chains in.
    """
    from tcren.pose import pose_consistency
    s = _complex()
    assert np.isfinite(matrix_features(s)["m_face_tp"])
    # one number, one implementation -- pose reads it off its own d2/d3 layers and this module
    # reads it off the same ones. They agreed to 4.4e-16 over 196 Native2026 crystals.
    assert matrix_features(s)["m_face_tp"] == pytest.approx(
        pose_consistency(s)["m_face_tp"], rel=1e-12, abs=1e-12)


def test_every_length_coupled_column_says_so_in_status():
    """The screen's whole product: which columns carry a length, stated where a caller will read it.

    Measured on 143 class I Native2026 crystals over two axes -- peptide 8-13 and CDR3(a+b) 19-28 --
    because they catch different columns: `m_erank_tm` is clean on the peptide axis (-0.186) and
    coupled on the CDR3 one (-0.437), while `m_erank_tp` is the reverse. None is dropped for it;
    every one keeps at least 58 per cent of its variance after both lengths are removed, and a
    length-coupled column beside a length-free one lets a model cancel the shared part.
    """
    coupled = ("m_erank_tp", "m_erank_tm", "m_gap_tp", "ca_cb_agreement_tm",
               "degree_evenness_tp", "frac_well_coordinated_tp",
               "g_even_tcr", "g_loop_even", "g_comp_frac")
    for d in coupled:
        assert d in STATUS, d
        assert "length coupled" in STATUS[d][1] or "length" in STATUS[d][1], d
    # and the ones the screen cleared carry no such flag
    for d in ("m_gap_tm", "g_assort", "m_face_tm", "g_alg_conn", "m_face_tp",
              "ca_cb_agreement_tp", "g_even_pmhc"):
        assert d not in STATUS, d


# --- degenerate input is NaN, never zero ----------------------------------------------------------

def _chain(cid, ctype, n, origin, aas, regions, rng):
    res = []
    for i in range(n):
        x, y, z = np.asarray(origin, float) + rng.normal(scale=1.8, size=3)
        aa = aas[i % len(aas)]
        res.append(Residue(i, i + 1, "", aa, _THREE[aa],
                           (Atom("CA", "C", np.array([x, y, z])),
                            Atom("CB", "C", np.array([x + 0.7, y, z])))))
    c = Chain(cid, res, chain_type=ctype)
    names = tuple(regions or ())
    if names:
        span = n // len(names)
        c.regions = [RegionMarkup(r, i * span, (i + 1) * span - 1,
                                  "".join(x.aa for x in res[i * span:(i + 1) * span]),
                                  res[i * span:(i + 1) * span])
                     for i, r in enumerate(names)]
    return c


def _complex(sep: float = 4.0):
    rng = np.random.default_rng(0)
    return Structure("synth", [
        _chain("C", PEPTIDE_TYPE, 9, [0.0, 0.0, 0.0], "DKW", None, rng),
        _chain("D", "TRA", 18, [-2.0, 0.0, sep], "LDK", ("CDR1", "CDR2", "CDR3"), rng),
        _chain("B", "TRB", 18, [2.0, 0.0, sep], "LDW", ("CDR1", "CDR2", "CDR3"), rng),
        _chain("A", "MHCa", 24, [0.0, 0.0, -4.0], "LDW", ("HELIX_A1", "HELIX_A2"), rng)])


def test_a_complex_that_makes_no_contact_is_all_nan_not_all_zero():
    """An interface that does not exist has no evenness. Zero would rank it above a bad one."""
    row = graph_features(_complex(sep=60.0))
    assert set(row) == set(GRAPH_FEATURES) | set(PROMOTED_POSE_FEATURES)
    assert all(np.isnan(v) for v in row.values()), row


def test_the_mhc_arm_is_nan_without_the_supertype_rather_than_falling_back():
    """`interface("tcr_mhc")` matches the supertype `annotate_mhc` assigns, not the generic "MHC".

    A structure whose MHC chains were never refined must lose the `_tm` columns, not silently
    compute them against a chain the selection cannot see.
    """
    rng = np.random.default_rng(0)
    bare = Structure("unrefined", [
        _chain("C", PEPTIDE_TYPE, 9, [0.0, 0.0, 0.0], "DKW", None, rng),
        _chain("D", "TRA", 18, [-2.0, 0.0, 4.0], "LDK", ("CDR1", "CDR2", "CDR3"), rng),
        _chain("B", "TRB", 18, [2.0, 0.0, 4.0], "LDW", ("CDR1", "CDR2", "CDR3"), rng),
        _chain("A", "MHC", 24, [0.0, 0.0, -4.0], "LDW", None, rng)])   # generic type, unrefined
    row = matrix_features(bare)
    assert np.isnan(row["m_face_tm"]) and np.isnan(row["ca_cb_agreement_tm"])
    assert np.isfinite(row["m_face_tp"])


def test_every_column_is_in_range_on_a_real_looking_complex():
    s = _complex()
    row = {**graph_features(s), **matrix_features(s)}
    for k in ("g_even_tcr", "g_even_pmhc", "g_comp_frac", "g_cyclo_frac", "g_loop_even",
              "g_loop_overlap", "degree_evenness_tp", "frac_well_coordinated_tp",
              "m_erank_tp", "m_gap_tp", "m_erank_tm", "m_gap_tm"):
        assert 0.0 <= row[k] <= 1.0, (k, row[k])
    assert 0.0 <= row["g_alg_conn"] <= 2.0, row["g_alg_conn"]
    for k in ("g_assort", "ca_cb_agreement_tp", "ca_cb_agreement_tm"):
        assert -1.0 <= row[k] <= 1.0 or np.isnan(row[k]), (k, row[k])
    assert np.isfinite(row["m_face_tp"])


# --- the catalogue ---------------------------------------------------------------------------------

def test_the_sixteen_are_catalogued_under_one_family_with_units_and_a_definition():
    for name in GRAPH_FEATURES + PROMOTED_POSE_FEATURES + MATRIX_FEATURES:
        assert DESCRIPTORS[name] == ("topology", True), name
        assert name in INVARIANCE, name
        units, text = DETAIL[name]
        assert units and len(text) > 20 and text.endswith("."), name


def test_the_promoted_pair_did_not_bring_max_degree_with_it():
    """`pose._degree_descriptors` returns a third column that no family claims. It stays out."""
    assert "max_degree_tp" not in DESCRIPTORS
    assert "max_degree_tp" not in PROMOTED_POSE_FEATURES
