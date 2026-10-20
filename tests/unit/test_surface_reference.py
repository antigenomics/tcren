"""Surface normalization preserves direction and cancels table offsets exactly."""
from types import SimpleNamespace

import polars as pl
import pytest

from tcren.contacts.surface_reference import surface_reference_score


def test_directed_surface_partition_offset_missingness_and_empty_sum():
    atoms = pl.DataFrame({'pdb.id':['x']*4,'chain.id':['B','B','T','T'],
        'residue.index':[0,1,0,1],'atom.name':['CA']*4,'residue.aa':['D','K','D','K'],
        'role':['binder']*2+['target']*2,'region.type':['CDR3']*2+['PEPTIDE']*2,
        'chain.type':['TRA']*2+['PEPTIDE']*2,'isolated_sasa_A2':[100.,50.,25.,75.]})
    edges = pl.DataFrame({'pdb.id':['x']*2,'binder_atom':[0,1],'target_atom':[3,2],
        'role.from':['binder']*2,'role.to':['target']*2,'region.type.from':['CDR3']*2,
        'chain.type.to':['PEPTIDE']*2,'residue.aa.from':['D','K'],'residue.aa.to':['K','D'],
        'binder_area_A2':[10.,20.]})
    matrix = pl.DataFrame({'residue.aa.from':['D','D','K','K'],
        'residue.aa.to':['D','K','D','K'],'value':[0.,-2.,1.,0.]})
    def score(m=matrix, e=edges):
        return surface_reference_score(atoms,e,SimpleNamespace(matrix=m))
    result = score()
    assert result['surface_reference_factor_mean'] == pytest.approx(1.,abs=1e-14)
    assert result['surface_raw_sum'] == 0.
    assert result['surface_context_sum'] == pytest.approx(30*result['surface_logZ'])
    shifted = score(matrix.with_columns(pl.col('value')+7.3))
    assert shifted['surface_context_sum'] == pytest.approx(result['surface_context_sum'],abs=1e-12)
    transposed = matrix.rename({'residue.aa.from':'residue.aa.to','residue.aa.to':'residue.aa.from'})
    assert score(transposed)['surface_context_sum'] != pytest.approx(result['surface_context_sum'])
    missing = score(matrix.filter(~((pl.col('residue.aa.from')=='D')&(pl.col('residue.aa.to')=='K'))))
    assert missing['surface_missing_potential_area_A2'] == 10.
    assert missing['surface_supported_area_A2'] == 20.
    empty = score(e=edges.head(0))
    assert empty['surface_context_sum'] == empty['surface_raw_sum'] == 0.
    assert empty['surface_context_mean'] is None and empty['surface_status'] == 'no_supported_area'
