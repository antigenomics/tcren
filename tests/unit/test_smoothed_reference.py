"""The smoothed reference delta and its curvature (tcren.ddg.smoothed_reference)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from tcren.annotation import classify_chains
from tcren.contactmap import ContactMap
from tcren.ddg import reference_delta, smoothed_reference
from tcren.mhc import annotate_mhc
from tcren.potential import mj, tcren2
from tcren.scoring import recognition_matrix
from tcren.structure import import_structure

_ASSET = Path(__file__).resolve().parents[1] / "assets" / "cgene" / "1ao7_full.pdb"


@pytest.fixture(scope="module")
def cm():
    pytest.importorskip("arda")  # classify_chains needs the arda backend; CI installs tcren without it
    s = import_structure(_ASSET)
    classify_chains(s)
    annotate_mhc(s)
    return ContactMap.from_structure(s, cutoff=5.0)


@pytest.fixture(scope="module")
def pot():
    return tcren2()


@pytest.fixture(scope="module")
def peptide():
    s = import_structure(_ASSET)
    classify_chains(s)
    return next(c.sequence() for c in s.chains if str(c.chain_type).upper().startswith("PEP"))


def test_the_third_interface_drops_out_of_each_direction(cm, pot, peptide):
    """A substitution on one chain leaves whole interfaces untouched, and they cancel exactly.

    Varying the peptide cannot change the TCR:MHC energy -- no peptide residue is in that
    interface, and neither partner moves. Varying the TCR cannot change peptide:MHC for the same
    reason. This holds on the virtual path (the contact map is frozen) and on the structural one
    (a truncated side chain moves no atom of the other two chains), so the two-term Hamiltonian
    each direction reduces to is exact rather than an approximation.
    """
    assert reference_delta(cm, peptide, mj(), interface="tcr_mhc") == pytest.approx(0.0, abs=1e-12)
    assert reference_delta(cm, peptide, pot, interface="tcr_peptide") != pytest.approx(0.0, abs=1e-6)

    pm = cm.interface("peptide_mhc")
    sides = set(pm["chain.type.from"].to_list()) | set(pm["chain.type.to"].to_list())
    assert not (sides & {"TRA", "TRB"}), sides       # no TCR residue to substitute


def test_the_cold_limit_is_the_arithmetic_mean_field(cm, pot):
    """As beta -> 0 the free energy of the background collapses to its mean.

    That limit IS the reference state a combinatorial peptide library realises: every position
    other than the one being read is held at an equimolar mixture, so the completing term is the
    twenty-residue mean rather than any single residue.
    """
    got = smoothed_reference(cm, pot, side="peptide", beta=1e-6)
    want = _mean_field(cm, pot, "peptide")
    assert got["dPhi"] == pytest.approx(want, abs=1e-3)


def test_the_hot_limit_is_the_distance_from_the_best_residue(cm, pot):
    """As beta grows the tilted weight becomes a point mass on the best residue at each position.

    The limit carries the background's own entropy: the free energy is
    ``min_a phi(a) - log p(a*) / beta``, not ``min_a phi(a)``, so at finite beta the delta sits
    ``n * log 20 / beta`` below the minimum for an equimolar background. The test asserts that
    identity rather than the bare limit, which is what tells the two apart; the tolerance leaves
    room for the runner-up residue's ``exp(-beta * gap)``, which is not identically zero at any
    finite beta.
    """
    beta = 500.0
    got = smoothed_reference(cm, pot, side="peptide", beta=beta)
    n = got["n_positions"]
    want = _min_field(cm, pot, "peptide") - n * np.log(20.0) / beta
    assert got["dPhi"] == pytest.approx(want, abs=1e-4)
    assert got["varPhi"] == pytest.approx(0.0, abs=1e-6)      # the tilted weight is a point mass


def test_the_two_chains_partition_the_receptor_direction(cm, pot):
    """TRA and TRB are disjoint position sets, so their deltas add to the pooled one exactly.

    This is why both are emitted: a linear model handed the two parts can form any contrast
    between them -- including the TRB - TRA difference the benchmark's contact-subset search found
    -- while one handed only their sum cannot.
    """
    both = smoothed_reference(cm, pot, side="tcr")
    a = smoothed_reference(cm, pot, side="tcr", chain="TRA")
    b = smoothed_reference(cm, pot, side="tcr", chain="TRB")
    assert a["dPhi"] + b["dPhi"] == pytest.approx(both["dPhi"], abs=1e-9)
    assert a["n_positions"] + b["n_positions"] == both["n_positions"]


def test_the_curvature_is_a_sum_of_variances(cm, pot):
    for side in ("peptide", "tcr"):
        assert smoothed_reference(cm, pot, side=side)["varPhi"] >= 0.0


def test_the_interface_weights_are_the_native_scales(cm, pot):
    """Each surviving interface enters divided by its own Native2026 spread, not raw.

    Passing both weights as 1 must change the answer, or the coefficients are not doing anything.
    """
    from tcren.pipeline import _phi_scale

    assert _phi_scale("tcr_peptide", pot) == pytest.approx(1.6389607864045381)
    assert _phi_scale("peptide_mhc", mj()) == pytest.approx(4.301317384887072)

    weighted = smoothed_reference(cm, pot, side="peptide")
    raw = smoothed_reference(cm, pot, side="peptide",
                             weights={"tcr_peptide": 1.0, "peptide_mhc": 1.0})
    assert weighted["dPhi"] != pytest.approx(raw["dPhi"], rel=1e-3)


def test_the_side_and_beta_are_validated(cm, pot):
    with pytest.raises(ValueError, match="side must be"):
        smoothed_reference(cm, pot, side="mhc")
    with pytest.raises(ValueError, match="beta must be positive"):
        smoothed_reference(cm, pot, beta=0.0)


def _fields(cm, pot, side):
    """The same two-interface local field the function builds, assembled independently here."""
    from tcren.ddg import SMOOTH_INTERFACES
    from tcren.pipeline import _phi_scale

    out: dict = {}
    aa: tuple = ()
    for iface, which in SMOOTH_INTERFACES[side]:
        p = pot if iface == "tcr_peptide" else mj()
        c = 1.0 / _phi_scale(iface, p)
        rm = recognition_matrix(cm, p, interface=iface, side=which)
        aa = aa or rm.aa
        for i, key in enumerate(rm.positions):
            v = c * np.asarray(rm.energy, float)[i]
            out[key] = v if key not in out else out[key] + v
    keys = list(out)
    return np.vstack([out[k] for k in keys]), np.array([aa.index(k[3]) for k in keys])


def _mean_field(cm, pot, side):
    phi, native = _fields(cm, pot, side)
    return float(np.nansum(phi[np.arange(len(native)), native] - np.nanmean(phi, axis=1)))


def _min_field(cm, pot, side):
    phi, native = _fields(cm, pot, side)
    return float(np.nansum(phi[np.arange(len(native)), native] - np.nanmin(phi, axis=1)))
