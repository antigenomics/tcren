"""The amino-acid background taken in 3D: three outcomes, and the kernel matches its reference."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from tcren.energetics import background as bg

ASSET = Path(__file__).resolve().parents[1] / "assets" / "pdb" / "1ao7.pdb"
pytestmark = pytest.mark.slow


@pytest.fixture(scope="module")
def structure():
    pytest.importorskip("arda")
    from tcren.annotation import classify_chains
    from tcren.structure import parse_structure

    s = parse_structure(ASSET, pdb_id="1ao7")
    classify_chains(s, organism="human", autodetect_species=True)
    return s


@pytest.fixture(scope="module")
def scanned(structure):
    return bg.scan_structure(structure, "tcr_peptide", "to", step=60.0, clash_tol=0.8)


def test_the_native_kernel_matches_the_numpy_reference(structure):
    """The only way a hand-written C++ inner loop is trustworthy."""
    nat, _p, _q = bg.scan_structure(structure, "tcr_peptide", "to")
    ref, _p, _q = bg.scan_structure(structure, "tcr_peptide", "to", engine="numpy")
    assert np.abs(nat.p - ref.p).max() < 1e-9
    assert (nat.admissible == ref.admissible).all()
    assert (nat.n_clash == ref.n_clash).all()
    assert (nat.n_rotamer == ref.n_rotamer).all()


def test_the_scan_covers_the_peptide(scanned):
    r, pos, par = scanned
    assert [x[1].aa for x in pos] == list("LLFGYPVYV"), "1ao7 carries the Tax nonamer"
    assert r.p.shape == (len(r.pair_i), 20)
    assert r.admissible.shape == (len(pos), 20)


def test_probabilities_are_probabilities(scanned):
    r, _pos, _par = scanned
    assert np.isfinite(r.p).all()
    assert (r.p >= 0.0).all() and (r.p <= 1.0).all()


def test_alanine_and_glycine_are_always_admissible(scanned):
    """Neither reaches past Cbeta, so neither can collide on a backbone the crystal already holds."""
    r, _pos, _par = scanned
    for aa in ("A", "G"):
        assert r.admissible[:, bg.AA.index(aa)].all()


def test_a_deposited_residue_is_admissible_at_its_own_position(scanned):
    """The invariant the grid is calibrated by: the crystal proves the residue fits."""
    r, pos, _par = scanned
    for i, (_c, res) in enumerate(pos):
        assert r.admissible[i, bg.AA.index(res.aa)], f"{res.aa} at position {i} rejected itself"


def test_all_three_outcomes_occur(scanned):
    """A scan that never collides, or never contacts, is not measuring what it claims to."""
    r, _pos, _par = scanned
    adm = r.admissible[r.pair_i]
    assert (~adm).any(), "nothing collided"
    assert (adm & (r.p <= 0)).any(), "nothing was a clean non-contact"
    assert (adm & (r.p > 0)).any(), "nothing contacted"


def test_a_collision_is_every_rotamer(scanned):
    """Inadmissible means every rotamer clashed, not merely most of them."""
    r, _pos, _par = scanned
    bad = ~r.admissible
    assert (r.n_clash[bad] == r.n_rotamer[bad]).all()
    assert (r.n_clash[r.admissible] < r.n_rotamer[r.admissible]).all()


def test_an_inadmissible_candidate_leaves_the_free_energy(scanned):
    """The point of the exercise: an impossible substitution is absent, not merely unfavourable."""
    from tcren.potential import tcren2

    r, pos, par = scanned
    out = bg.reference_free_energy(r, pos, par, tcren2(), side="to")
    assert np.isnan(out["phi"][~r.admissible]).all()
    assert (out["n_admissible"] == r.admissible.sum(1)).all()
    assert np.isfinite(out["F"]).all()
    # F is a free energy over the admissible set, so it sits inside that set's range.
    for i in range(len(pos)):
        live = out["phi"][i][r.admissible[i]]
        assert live.min() - 1e-9 <= out["F"][i] <= live.max() + 1e-9


def test_the_potential_direction_is_not_symmetric(scanned):
    """TCRen2 is directional; reading it the wrong way round transposes it silently."""
    from tcren.potential import tcren2

    r, pos, par = scanned
    a = bg.reference_free_energy(r, pos, par, tcren2(), side="to")["dPhi"]
    b = bg.reference_free_energy(r, pos, par, tcren2(), side="from")["dPhi"]
    assert not np.allclose(np.nan_to_num(a), np.nan_to_num(b))


def test_an_unknown_interface_or_side_is_refused(structure):
    with pytest.raises(ValueError, match="unknown interface"):
        bg.scan_structure(structure, "nonsense", "to")
    with pytest.raises(ValueError, match="side must be"):
        bg.scan_structure(structure, "tcr_peptide", "sideways")
    with pytest.raises(ValueError, match="engine must be"):
        bg.scan_structure(structure, "tcr_peptide", "to", engine="cuda")


def test_native_scan_validates_packed_inputs():
    """Synthetic empty scan tests the binding without structures or scientific reruns."""
    from tcren import _relax

    args = dict(bb=np.empty((0, 4, 3)), n_atom=np.zeros(20, np.int32),
                n_chi=np.zeros(20, np.int32), chi_vals=np.array([0.0]),
                env_xyz=np.empty((0, 3)), env_cls=np.empty(0, np.int32),
                env_rad=np.empty(0), env_res=np.empty(0, np.int32),
                env_ptr=np.array([0], np.int32), env_idx=np.empty(0, np.int32),
                par_xyz=np.empty((0, 3)), par_ptr=np.array([0], np.int32),
                pair_i=np.empty(0, np.int32), pair_j=np.empty(0, np.int32),
                dope_table=np.zeros((1, 1, 2), np.float32), n_cls=1, n_bins=2,
                x_start=0.0, dx=1.0, cutoff=5.0, clash_tol=0.8,
                temperature=1.0, n_threads=1)
    for key in ('t_a', 't_b', 't_c', 't_chi', 't_cls'):
        args[key] = np.zeros((20, 1), np.int32)
    for key in ('t_bond', 't_ang', 't_off', 't_rad'):
        args[key] = np.zeros((20, 1))
    assert _relax.substitution_scan(**args)['p'].shape == (0, 20)
    for change in ({'temperature': 0.0}, {'dx': float('nan')},
                   {'bb': np.empty((0, 3))}, {'n_chi': np.full(20, 5, np.int32)},
                   {'env_ptr': np.array([1], np.int32)},
                   {'pair_i': np.array([0], np.int32), 'pair_j': np.array([0], np.int32)},
                   {'t_a': np.zeros((20, 2), np.int32)[:, ::2]}):
        with pytest.raises(ValueError):
            _relax.substitution_scan(**(args | change))
    # A single synthetic backbone exercises active topology and the released-GIL path.
    args.update(bb=np.array([[[0., 1., 0.], [0., 0., 0.], [1., 0., 0.],
                              [0., 0., 1.]]]), env_ptr=np.array([0, 0], np.int32),
                par_xyz=np.array([[0., 0., 2.]]), par_ptr=np.array([0, 1], np.int32),
                pair_i=np.array([0], np.int32), pair_j=np.array([0], np.int32),
                n_atom=np.ones(20, np.int32), n_chi=np.ones(20, np.int32),
                t_b=np.ones((20, 1), np.int32), t_c=np.full((20, 1), 3, np.int32),
                t_chi=np.ones((20, 1), np.int32), t_bond=np.ones((20, 1)),
                t_ang=np.ones((20, 1)))
    native = _relax.substitution_scan(**args)
    reference = bg.scan_numpy(**{k: v for k, v in args.items() if k != 'n_threads'})
    for key in native:
        np.testing.assert_allclose(native[key], reference[key])
    for change in ({'t_a': np.full((20, 1), 4, np.int32)},
                   {'t_chi': np.full((20, 1), 2, np.int32)},
                   {'t_cls': np.ones((20, 1), np.int32)},
                   {'pair_j': np.array([1], np.int32)}):
        with pytest.raises(ValueError):
            _relax.substitution_scan(**(args | change))
