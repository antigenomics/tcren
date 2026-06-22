"""Potential-guided DOPE refinement: the C++ kernel (pure) + the refine_peptide wrapper (arda-gated)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

_refine = pytest.importorskip("tcren._refine")

_I = lambda *v: np.asarray(v, dtype=np.int32)  # noqa: E731


def _table(n_cls, pairs, nbins):
    """Symmetric (n_cls, n_cls, nbins) DOPE-style table from {(i,j): [knots]}."""
    t = np.zeros((n_cls, n_cls, nbins), dtype=np.float32)
    for (i, j), v in pairs.items():
        t[i, j] = v
        t[j, i] = v
    return t


# A 1-class toy potential: knots at d = 1, 2, 3 Å (x_start=1, dx=1): repulsive, favourable, zero.
_TOY = _table(1, {(0, 0): [10.0, -2.0, 0.0]}, 3).reshape(-1)


def test_kernel_interpolates_potential():
    # atoms 2.5 Å apart → linear interp between -2 (@2 Å) and 0 (@3 Å) = -1.
    e = _refine.refine(np.array([[0.0, 0, 0]]), _I(0), np.array([[2.5, 0, 0]]), _I(0),
                       _TOY, 1, 3, 1.0, 1.0, restraint_w=0.0, n_steps=0, seed=0)[1]
    assert e == pytest.approx(-1.0)


def test_kernel_relieves_repulsion_and_moves():
    # 1 Å apart → +10 (repulsive); MC moves the peptide away to a lower-energy pose.
    best, e, n_acc = _refine.refine(np.array([[0.0, 0, 0]]), _I(0), np.array([[1.0, 0, 0]]), _I(0),
                                    _TOY, 1, 3, 1.0, 1.0, restraint_w=0.1, n_steps=2000,
                                    trans_sigma=0.3, seed=1)
    assert e < 10.0
    assert np.linalg.norm(best[0] - np.array([1.0, 0, 0])) > 1.0
    assert n_acc > 0


def test_kernel_native_barely_moves():
    # partner beyond d_max (=3 Å) → no interaction; restraint holds the peptide home.
    best, e, _ = _refine.refine(np.array([[0.0, 0, 0]]), _I(0), np.array([[50.0, 0, 0]]), _I(0),
                                _TOY, 1, 3, 1.0, 1.0, restraint_w=1.0, n_steps=2000, seed=1)
    assert np.linalg.norm(best[0]) < 0.5
    assert e == pytest.approx(0.0, abs=1e-9)


def test_kernel_skips_unmapped_class():
    # an unmapped (class -1) atom contributes no energy.
    e = _refine.refine(np.array([[0.0, 0, 0]]), _I(-1), np.array([[2.0, 0, 0]]), _I(0),
                       _TOY, 1, 3, 1.0, 1.0, restraint_w=0.0, n_steps=0, seed=0)[1]
    assert e == pytest.approx(0.0)


def test_kernel_is_deterministic():
    args = (np.array([[0.0, 0, 0]]), _I(0), np.array([[1.0, 0, 0]]), _I(0), _TOY, 1, 3, 1.0, 1.0)
    a = _refine.refine(*args, n_steps=1000, seed=7)
    b = _refine.refine(*args, n_steps=1000, seed=7)
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
    out, energy = refine_peptide(s, restraint_w=0.5, n_steps=2000, seed=1)
    rmsd = float(np.sqrt(((ca(out) - before) ** 2).sum(1).mean()))
    assert rmsd < 1.0          # the native pose is a DOPE optimum
    assert energy < 0.0        # favourable DOPE energy
