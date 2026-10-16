"""Spatial chemistry requires co-located fields, not products of marginal means."""
import numpy as np
import polars as pl
import pytest


def test_spatial_covariance_oracle_missingness_and_export():
    from tcren.topology.pmhc import spatial_features

    cells = pl.DataFrame({'pdb.id': ['a']*6, 'iy': [0]*6, 'ix': list(range(6)),
                         'peptide': ['AL']*6, 'mhc_class': ['I']*6,
                         'source': ['peptide']*4+['mhc_helix_a1']*2,
                         'h': [0., 1., 2., 3., 1., 1.],
                         'phobic': [-1., 0., 1., 2., 0., 0.]})
    r = spatial_features(cells).row(0, named=True)
    assert r['ridge_height_A'] == 1
    assert r['ridge_excess_hydropathy_cov'] == .875
    assert r['height_hydropathy_cov'] == 1.25
    assert r['ridge_excess_hydropathy_mean'] == 1.25
    assert r['n_above_ridge_cells'] == 2
    assert r['above_ridge_hydropathy'] == 1.5
    assert r['ridge_excess_hydropathy_mean'] == r['ridge_excess_mean_A']*r['phobic_mean']+r['ridge_excess_hydropathy_cov']
    changed = cells.with_columns(pl.Series('phobic', [0., 2., -1., 1., 0., 0.]))
    q = spatial_features(changed).row(0, named=True)
    assert q['phobic_mean'] == r['phobic_mean']
    assert q['ridge_excess_hydropathy_cov'] != r['ridge_excess_hydropathy_cov']
    no_ridge = spatial_features(cells.filter(pl.col('source') == 'peptide')).row(0, named=True)
    assert no_ridge['ridge_excess_hydropathy_cov'] is None
    buried = cells.with_columns(pl.when(pl.col('source') != 'peptide').then(10.).otherwise(pl.col('h')).alias('h'))
    q = spatial_features(buried).row(0, named=True)
    assert q['ridge_excess_hydropathy_cov'] == 0 and q['above_ridge_hydropathy'] is None
    for bad in (pl.concat([cells,cells]), cells.with_columns(pl.lit('III').alias('mhc_class'))):
        with pytest.raises(ValueError):
            spatial_features(bad)
    from test_exports import typed_structure
    from tcren.exports import export_structure
    s = typed_structure()
    a = export_structure(s, 'pmhc-cells', None, None)
    for c in s.by_type('TRA', 'TRB'):
        for res in c.residues:
            for atom in res.atoms:
                atom.coord[:] += np.array([100., -100., -200.])
    b = export_structure(s, 'pmhc-cells', None, None)
    assert a.equals(b)
    assert spatial_features(a).height == 1
    both = pl.concat([cells, cells.with_columns(pl.lit('b').alias('pdb.id'), pl.lit('II').alias('mhc_class'))])
    assert spatial_features(both)['ridge_excess_hydropathy_cov'].to_list() == [.875, .875]
    pp = a.filter(pl.col('source') == 'peptide')
    hh = pp['h'].to_numpy(); kk = pp['phobic'].to_numpy()
    oracle = np.mean((hh-hh.mean())*(kk-kk.mean()))
    assert spatial_features(a)['height_hydropathy_cov'][0] == pytest.approx(oracle, abs=1e-12)
