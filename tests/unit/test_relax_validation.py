"""Native relax boundaries reject malformed arrays before unchecked access."""
import numpy as np
import pytest
from tcren import _relax


def _args(kind):
    xyz = np.array([[0., 0., 0.], [1., 0., 0.], [0., 1., 0.]])
    ints = lambda x: np.array(x, dtype=np.int32)
    common = dict(dope_table=np.ones((1, 1, 3), np.float32), n_cls=1,
                  n_bins=3, x_start=0., dx=1.)
    if kind == 'interface_energy':
        return common | dict(pep_xyz=xyz[:1], pep_class=ints([0]),
                             par_xyz=xyz[1:], par_class=ints([0, 0]))
    common |= dict(xyz=xyz, atom_class=ints([0, 0, 0]), n_steps=1, temperature=1.)
    if kind == 'repack':
        return common | dict(res_lo=ints([0]), res_hi=ints([2]), res_chi_ptr=ints([0, 1]),
                             chi_a=ints([0]), chi_b=ints([1]), chi_mov_ptr=ints([0, 1]),
                             moving=ints([1]), env_atom=ints([2]), env_ptr=ints([0, 1]))
    return common | dict(pep_lo=0, pep_hi=2, pep_res=ints([0, 0]), tor_a=ints([0]),
                         tor_b=ints([1]), tor_mov_ptr=ints([0, 1]), tor_mov=ints([1]),
                         par_atom=ints([2]), anchor_atom=ints([0]), burn_in=0)


@pytest.mark.parametrize('kind', ['interface_energy', 'repack', 'relax_interface'])
def test_relax_rejects_invalid_spacing(kind):
    with pytest.raises(ValueError):
        getattr(_relax, kind)(**(_args(kind) | {'dx': -1.}))


@pytest.mark.parametrize('kind', ['interface_energy', 'repack', 'relax_interface'])
def test_relax_rejects_shape_class_and_table(kind):
    args = _args(kind)
    coord, cls = ('pep_xyz', 'pep_class') if kind == 'interface_energy' else ('xyz', 'atom_class')
    for change in ({coord: np.zeros((1, 2))}, {cls: np.array([], np.int32)},
                   {cls: np.full(args[cls].shape, 1, np.int32)},
                   {'dope_table': np.zeros(1, np.float32)},
                   {'dope_table': np.zeros(6, np.float32)[::2]},
                   {coord: np.full(args[coord].shape, np.nan)}):
        with pytest.raises(ValueError):
            getattr(_relax, kind)(**(args | change))


@pytest.mark.parametrize('kind', ['repack', 'relax_interface'])
def test_relax_rejects_indices_and_csr(kind):
    args = _args(kind)
    ptr = 'chi_mov_ptr' if kind == 'repack' else 'tor_mov_ptr'
    axis = 'chi_a' if kind == 'repack' else 'tor_a'
    changes = [{ptr: np.array([0, 2], np.int32)}, {axis: np.array([3], np.int32)},
               {'temperature': 0.}, {'n_steps': 0}]
    changes += ([{'res_hi': np.array([4], np.int32)}, {'n_steps': 2**31 - 1,
                  'res_chi_ptr': np.array([0, 2], np.int32),
                  'chi_a': np.array([0, 0], np.int32), 'chi_b': np.array([1, 1], np.int32),
                  'chi_mov_ptr': np.array([0, 0, 1], np.int32)}]
                if kind == 'repack' else [{'anchor_atom': np.array([2], np.int32)},
                                          {'pep_hi': 4}, {'pep_res': np.array([0], np.int32)}])
    for change in changes:
        with pytest.raises(ValueError):
            getattr(_relax, kind)(**(args | change))


@pytest.mark.parametrize('kind', ['interface_energy', 'repack', 'relax_interface'])
def test_relax_valid_minimal_input(kind):
    out = getattr(_relax, kind)(**_args(kind))
    if kind == 'interface_energy':
        assert out == pytest.approx(2.)
    else:
        assert np.isfinite(out['energy']).all()
        assert np.isfinite(out['xyz']).all()
