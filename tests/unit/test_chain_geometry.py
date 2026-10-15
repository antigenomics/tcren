"""``enforce_bond``: a sampled or interpolated Cα trace is not yet a chain.

Consecutive Cα atoms are a fixed virtual bond apart. Nothing in a Gaussian draw or a bilinear
read-back knows that, so every geometry that does not come off a solved structure has to be
projected before it can be called a conformation.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from tcren.geometry import BOND_CA, enforce_bond

PDB_DIR = Path(__file__).resolve().parents[1] / "assets" / "pdb"
pytestmark = pytest.mark.skipif(not PDB_DIR.exists(), reason="pdb assets missing")


def _gaps(trace: np.ndarray) -> np.ndarray:
    return np.linalg.norm(np.diff(trace, axis=-2), axis=-1)


@pytest.mark.parametrize("n", [2, 3, 8, 14, 25])
def test_projection_puts_every_gap_on_the_bond(n):
    rng = np.random.default_rng(n)
    projected = enforce_bond(rng.normal(0.0, 4.0, (n, 3)))
    assert np.allclose(_gaps(projected), BOND_CA, atol=2e-7)


def test_a_trace_that_cannot_converge_is_refused_rather_than_returned_wrong():
    with pytest.raises(ValueError, match="did not converge"):
        enforce_bond(np.random.default_rng(0).normal(0.0, 4.0, (16, 3)), max_sweeps=1)


def test_projection_is_idempotent():
    """Re-projecting a converged trace moves atoms by less than the tolerance it converged to."""
    rng = np.random.default_rng(0)
    once = enforce_bond(rng.normal(0.0, 4.0, (16, 3)))
    assert np.abs(enforce_bond(once) - once).max() < 1e-6


def test_a_stack_projects_independently():
    rng = np.random.default_rng(1)
    traces = rng.normal(0.0, 4.0, (5, 14, 3))
    assert np.allclose(enforce_bond(traces),
                       np.stack([enforce_bond(t) for t in traces]), atol=1e-9)


def test_projection_moves_atoms_as_little_as_it_can():
    """A trace already close to the bond should barely move; a wild one moves more."""
    rng = np.random.default_rng(2)
    near = np.cumsum(rng.normal(0.0, 1.0, (12, 3)), axis=0)
    near *= BOND_CA / _gaps(near).mean()
    far = rng.normal(0.0, 20.0, (12, 3))
    assert np.abs(enforce_bond(near) - near).max() < np.abs(enforce_bond(far) - far).max()


def test_a_trace_of_one_atom_has_no_bond_to_enforce():
    single = np.zeros((1, 3))
    assert np.array_equal(enforce_bond(single), single)


def test_coincident_atoms_are_refused_rather_than_silently_nan():
    trace = np.zeros((4, 3))
    with pytest.raises(ValueError, match="coincident"):
        enforce_bond(trace)


def test_bad_shape_is_refused():
    with pytest.raises(ValueError, match=r"\(\.\.\., n, 3\)"):
        enforce_bond(np.zeros((5, 2)))


def test_real_cdr_traces_are_a_virtual_bond_chain_apart_from_cis_bonds():
    """If solved loops violated the constraint, nothing downstream of it would stand.

    Measured over the 17 bundled crystals: 642 CDR Cα gaps at 3.7985 Å, s.d. 0.0668, of which 638
    (99.38 %) sit in 3.5-3.9 Å. The exceptions are *cis* peptide bonds near 2.9 Å, which is a real
    conformation and not a defect -- so this asserts the trans population and that every outlier is
    cis, rather than asserting a single value.
    """
    pytest.importorskip("arda")
    from tcren.annotation.batch import iter_annotated_set

    gaps = []
    for structure in iter_annotated_set(PDB_DIR):
        for chain in structure.chains:
            if chain.chain_type not in ("TRA", "TRB"):
                continue
            for region in chain.regions or []:
                if region.region_type not in ("CDR1", "CDR2", "CDR3"):
                    continue
                ca = [r.ca for r in region.residues if r.ca is not None]
                if len(ca) >= 2:
                    gaps.append(_gaps(np.asarray(ca, float)))
    gaps = np.concatenate(gaps)
    assert len(gaps) > 500, f"only {len(gaps)} gaps seen; the fixture set did not annotate"

    trans = np.abs(gaps - BOND_CA) < 0.35
    assert trans.mean() > 0.98, f"only {100 * trans.mean():.2f}% of gaps are trans"
    assert gaps[trans].mean() == pytest.approx(BOND_CA, abs=0.02)
    # Every non-trans gap is a cis peptide bond (Ca-Ca ~2.9 A), not an arbitrary chain break.
    assert np.all((gaps[~trans] > 2.7) & (gaps[~trans] < 3.1)), gaps[~trans].round(3)
