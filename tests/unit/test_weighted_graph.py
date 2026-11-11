"""Hand-counted checks for energy-weighted recognition graphs."""
import numpy as np
import pytest

from tcren.topology.weighted import _weighted_features


def test_weighted_square_strength_cut_and_spectrum():
    edges = {(0, 2): 1., (0, 3): 2., (1, 2): 3., (1, 3): 4.}
    row = _weighted_features(edges, {0, 1}, {0: 'a', 1: 'b'}, {2: 'p', 3: 'm'})
    assert row['w_total'] == 10
    assert row['s_max'] == 7
    assert row['s_mean'] == 5
    assert row['clustering'] == row['transitivity'] == 0
    assert row['beta1'] == 1
    assert row['mincut_tcr_pmhc'] == 4
    assert row['wmincut_tcr_pmhc'] == 10
    assert row['w_cyclo'] == pytest.approx(.1)
    a = np.array([[0, 0, 1, 2], [0, 0, 3, 4], [1, 3, 0, 0], [2, 4, 0, 0]])
    assert row['graph_energy'] == pytest.approx(abs(np.linalg.eigvalsh(a)).sum())
    assert 0 <= row['alg_conn_w_bounded'] <= 1
    assert row['mincut_frac'] == 1
    assert row['partcoef_tcr_w'] == pytest.approx(((1-(1/3)**2-(2/3)**2)+(1-(3/7)**2-(4/7)**2))/2)


def test_triangle_and_missing_weights():
    row = _weighted_features({(0, 1): 1., (0, 2): 1., (1, 2): 1.}, {0}, {0: 'a'}, {1: 'p', 2: 'm'})
    assert row['clustering'] == row['transitivity'] == 1
    row = _weighted_features({(0, 1): None}, {0}, {0: 'a'}, {1: 'p'})
    assert row['missing_weight'] == 1
    assert np.isnan(row['w_total'])
    assert row['edges'] == 1


def test_zero_weights_and_empty_graph():
    row = _weighted_features({(0, 1): 0.}, {0}, {0: 'a'}, {1: 'p'})
    assert row['w_total'] == 0
    assert np.isnan(row['w_cyclo'])
    assert np.isnan(row['wmincut_frac'])
    assert not _weighted_features({}, set(), {}, {})['available']


def test_betweenness_tied_paths_and_capacity_boundary():
    from tcren.topology.weighted import _betweenness, _cut
    adjacency = np.array([[0., 0., 1., 1.], [0., 0., 1., 1.],
                          [1., 1., 0., 0.], [1., 1., 0., 0.]])
    np.testing.assert_allclose(_betweenness(adjacency), np.full(4, 1/6))
    with pytest.raises(ValueError, match='capacity'):
        _cut(np.array([[0, 3_000_000_000], [3_000_000_000, 0]], dtype=np.int64), np.array([True, False]))


def test_directed_potentials_ignore_physical_chain_order(monkeypatch):
    from types import SimpleNamespace
    import polars as pl
    import tcren.exports
    import tcren.potential
    from tcren.topology.weighted import weighted_graph_features
    from tcren.structure.model import Chain, RegionMarkup, Residue, Structure
    def chain(key, aa, kind):
        residue = Residue(0, 1, '', aa, 'ALA', ())
        c = Chain(key, [residue], chain_type=kind)
        if kind == 'TRA':
            c.regions = [RegionMarkup('CDR3', 0, 0, aa, [residue])]
        return c
    structure = Structure('direction', [chain('Z', 'A', 'TRA'), chain('A', 'G', 'PEPTIDE'), chain('B', 'R', 'MHCa')])
    parts = {'tcr_peptide': {(('A', 0), ('Z', 0))}, 'tcr_mhc': {(('B', 0), ('Z', 0))},
             'peptide_mhc': {(('A', 0), ('B', 0))}}
    monkeypatch.setattr(tcren.exports, '_recognition_graphs', lambda *a, **k:
        ([('tp_pm_tm', set().union(*parts.values()))], parts, set(), set(), 0))
    def potential(rows):
        return SimpleNamespace(matrix=pl.DataFrame(rows, schema=['residue.aa.from', 'residue.aa.to', 'value'], orient='row'))
    monkeypatch.setattr(tcren.potential, 'tcren2', lambda: potential([('A','G',-1.),('G','A',100.),('R','G',-2.),('G','R',100.)]))
    monkeypatch.setattr(tcren.potential, 'mj', lambda: potential([('A','R',-3.),('R','A',100.)]))
    row = weighted_graph_features(structure, None).row(0, named=True)
    assert row['w_total'] == 6
    assert row['missing_weight'] == 0
    # Explicit asymmetric PM weights must use peptide G -> MHC R, not its transpose.
    custom = weighted_graph_features(structure, None,
        peptide_potential=potential([('A','G',-4.),('G','A',100.)]),
        mhc_potential=potential([('A','R',-5.),('R','A',100.)]),
        pmhc_potential=potential([('G','R',-6.),('R','G',-100.)])).row(0,named=True)
    assert custom['w_total'] == 15
    assert custom['missing_weight'] == 0
    boltzmann = weighted_graph_features(structure, None,
        peptide_potential=potential([('A','G',1.),('G','A',100.)]),
        mhc_potential=potential([('A','R',0.),('R','A',100.)]),
        pmhc_potential=potential([('G','R',-1.),('R','G',-100.)]),
        weighting='boltzmann', beta=1.).row(0,named=True)
    assert boltzmann['w_total'] == pytest.approx(np.exp(-1)+1+np.exp(1))
    assert boltzmann['w_max'] == pytest.approx(np.exp(1))
    assert boltzmann['missing_weight'] == 0
    for scheme in ['positive','boltzmann']:
        missing = weighted_graph_features(structure,None,
            peptide_potential=potential([('A','G',float('inf'))]),
            mhc_potential=potential([('A','R',0.)]),
            pmhc_potential=potential([('G','R',-1.)]),
            weighting=scheme).row(0,named=True)
        assert missing['missing_weight']==1
        assert np.isnan(missing['w_total'])
    with pytest.raises(ValueError,match='overflow/underflow'):
        weighted_graph_features(structure,None,weighting='boltzmann',beta=1000)


def test_per_contacting_residue_normalization_and_missing_weights():
    row = _weighted_features({(0,2):1.,(0,3):2.,(1,2):3.,(1,3):4.},
                             {0,1},{0:'a',1:'b'},{2:'p',3:'m'})
    assert row['edges_per_vertex'] == 1
    assert row['edges_per_vertex_bounded'] == .5
    assert row['w_per_vertex'] == 2.5
    assert row['w_per_vertex_bounded'] == pytest.approx(10/14)
    row = _weighted_features({(0,1):None},{0},{0:'a'},{1:'p'})
    assert row['edges_per_vertex'] == .5
    assert np.isnan(row['w_per_vertex'])


@pytest.mark.parametrize("beta", [0.,-1.,float("nan"),float("inf")])
def test_boltzmann_rejects_invalid_beta(beta):
    from tcren.topology.weighted import weighted_graph_features
    with pytest.raises(ValueError,match="finite beta"):
        weighted_graph_features(None,None,weighting="boltzmann",beta=beta)
