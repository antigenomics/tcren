"""Potential-guided refinement: the C++ kernel (pure) + the refine_peptide wrapper (arda-gated)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

_refine = pytest.importorskip("tcren._refine")

_I = lambda *v: np.asarray(v, dtype=np.int32)  # noqa: E731


def test_kernel_relieves_clash_and_moves_away():
    # one peptide atom 1 Å from one partner atom → clash (d0=3); zero potential isolates clash+restraint.
    pep = np.array([[0.0, 0.0, 0.0]]); par = np.array([[1.0, 0.0, 0.0]])
    best, e, n_acc = _refine.refine(
        pep, _I(0), _I(0), par, _I(0), _I(0), np.zeros((1, 1)),
        cutoff=5.0, clash_d0=3.0, clash_w=10.0, restraint_w=0.1,
        n_steps=2000, trans_sigma=0.3, seed=1,
    )
    assert e < 10 * (3 - 1) ** 2          # below the start clash energy (40)
    assert np.linalg.norm(best[0] - par[0]) > 1.0   # moved away from the clashing partner
    assert n_acc > 0


def test_kernel_native_pose_barely_moves():
    # partner far away (no clash, no contact); restraint keeps the peptide at the start.
    pep = np.array([[0.0, 0.0, 0.0]]); par = np.array([[50.0, 0.0, 0.0]])
    best, e, _ = _refine.refine(
        pep, _I(0), _I(0), par, _I(0), _I(0), np.zeros((1, 1)),
        clash_w=10.0, restraint_w=1.0, n_steps=2000, trans_sigma=0.2, seed=1,
    )
    assert np.linalg.norm(best[0]) < 0.5  # stayed home
    assert e == pytest.approx(0.0, abs=1e-9)


def test_kernel_sums_potential_over_contacts():
    # peptide & partner 4 Å apart (in contact, no clash); a favourable (-2) potential lowers energy.
    pep = np.array([[0.0, 0.0, 0.0]]); par = np.array([[4.0, 0.0, 0.0]])
    e = _refine.refine(pep, _I(0), _I(0), par, _I(0), _I(0), np.array([[-2.0]]),
                       cutoff=5.0, n_steps=0, restraint_w=1.0, seed=0)[1]
    assert e == pytest.approx(-2.0)       # one contact, potential[-2], no clash/restraint at start


def test_kernel_is_deterministic():
    pep = np.array([[0.0, 0.0, 0.0]]); par = np.array([[1.0, 0.0, 0.0]])
    args = (pep, _I(0), _I(0), par, _I(0), _I(0), np.zeros((1, 1)))
    a = _refine.refine(*args, clash_w=10.0, n_steps=1000, seed=7)
    b = _refine.refine(*args, clash_w=10.0, n_steps=1000, seed=7)
    assert a[1] == b[1] and np.allclose(a[0], b[0])


# --- integration (needs arda for chain typing) ---------------------------------------------------

_FIXTURE = Path(__file__).resolve().parents[1] / "assets" / "pdb" / "1ao7.pdb"


def test_refine_peptide_native_barely_moves():
    pytest.importorskip("arda")
    from tcren.annotation import classify_chains
    from tcren.refine import refine_peptide
    from tcren.structure.io import import_structure
    from tcren.structure.model import PEPTIDE_TYPE

    s = import_structure(_FIXTURE)
    classify_chains(s, organism="human")
    ca = lambda st: np.array([r.ca for c in st.chains if c.chain_type == PEPTIDE_TYPE  # noqa: E731
                              for r in c.residues if r.ca is not None])
    before = ca(s)
    out, energy = refine_peptide(s, restraint_w=1.0, n_steps=2000, seed=1)
    rmsd = float(np.sqrt(((ca(out) - before) ** 2).sum(1).mean()))
    assert rmsd < 1.0                     # a native pose is already a local optimum
    assert isinstance(energy, float)
