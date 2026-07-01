"""Fast tests for the native interface-energy kernel (tcren._relax) + interface_energy wrapper."""

from __future__ import annotations

import numpy as np
import pytest

_relax = pytest.importorskip("tcren._relax")

# A tiny 2-class, 5-bin DOPE-like table: only table[cp=0][cq=1] is non-zero, a ramp per distance bin.
N_CLS, N_BINS, X0, DX = 2, 5, 0.0, 1.0


def _table():
    t = np.zeros((N_CLS, N_CLS, N_BINS), np.float32)
    t[0, 1] = [10, 20, 30, 40, 50]
    return t.reshape(-1)


def _ie(pep, pcl, par, qcl):
    return _relax.interface_energy(np.asarray(pep, float), np.asarray(pcl, np.int32),
                                   np.asarray(par, float), np.asarray(qcl, np.int32),
                                   _table(), N_CLS, N_BINS, X0, DX)


def test_linear_interpolation_between_bins():
    # dist 2.5, x_start=0 dx=1 -> t=2.5 -> knots[2]*0.5 + knots[3]*0.5 = 30*.5 + 40*.5 = 35
    assert _ie([[0, 0, 0]], [0], [[2.5, 0, 0]], [1]) == pytest.approx(35.0, abs=1e-5)


def test_out_of_range_contributes_zero():
    # d_max = x_start + (n_bins-1)*dx = 4.0; a pair at 10 Å is beyond range -> 0
    assert _ie([[0, 0, 0]], [0], [[10.0, 0, 0]], [1]) == 0.0


def test_negative_class_skipped():
    assert _ie([[0, 0, 0]], [-1], [[2.5, 0, 0]], [1]) == 0.0  # peptide class -1
    assert _ie([[0, 0, 0]], [0], [[2.5, 0, 0]], [-1]) == 0.0  # partner class -1


def test_short_range_capped_to_first_knot():
    # t <= 0 -> knots[0]
    assert _ie([[0, 0, 0]], [0], [[0.0, 0, 0]], [1]) == pytest.approx(10.0, abs=1e-5)


def test_sum_over_pairs():
    # two partner atoms, both class 1, at 2.5 (->35) and 1.0 (->20) -> 55
    assert _ie([[0, 0, 0]], [0], [[2.5, 0, 0], [1.0, 0, 0]], [1, 1]) == pytest.approx(55.0, abs=1e-5)


@pytest.mark.slow
def test_interface_energy_native_negative_and_separated_zero():
    pytest.importorskip("arda")
    from tcren.annotation import classify_chains
    from tcren.mhc import annotate_mhc
    from tcren.paths import reference_structure_path
    from tcren.refine.interface import interface_energy
    from tcren.structure import parse_structure
    from tcren.structure.model import PEPTIDE_TYPE, Atom, Chain, Residue, Structure

    s = parse_structure(reference_structure_path("1ao7"), pdb_id="1ao7")
    classify_chains(s, organism="human")
    annotate_mhc(s)
    assert interface_energy(s) < 0.0  # a real bound interface is favourable

    pep = next(c for c in s.chains if c.chain_type == PEPTIDE_TYPE)
    far = [Residue(r.seq_index, r.pdb_index, r.insertion_code, r.aa, r.resname,
                   tuple(Atom(a.name, a.element, a.coord + np.array([1000.0, 0, 0])) for a in r.atoms))
           for r in pep.residues]
    s2 = Structure(s.pdb_id, [Chain(pep.chain_id, far, chain_type=pep.chain_type,
                                    chain_supertype=pep.chain_supertype) if c is pep else c
                              for c in s.chains])
    assert interface_energy(s2) == 0.0  # separated -> no interaction
