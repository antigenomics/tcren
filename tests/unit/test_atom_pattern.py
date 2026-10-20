"""The reference measure must match the per-edge query measure."""
import math

import polars as pl
import pytest

from tcren.contacts.atom_patterns import derive_atom_patterns, validate_resource
from tcren.contacts.atom_counts import TYPES, score_atom_counts


def test_edge_mass_equivalence_and_frozen_scorer():
    rows = [dict(edge=e, aa_from='Y', atom_from=i, aa_to='P', atom_to=j)
            for e, i, j in [('dense','OH','O'),('dense','CZ','CB'),('single','OH','O')]]
    data = pl.DataFrame(rows)
    r = derive_atom_patterns(data, data.head(0))
    cell = {tuple(c[k] for k in TYPES): c for c in r['cells']}
    assert cell['Y','OH','P','O']['mass'] == 1.5
    assert cell['Y','CZ','P','CB']['mass'] == .5
    assert cell['Y','OH','P','O']['value'] == pytest.approx(-math.log(4/3))
    # Duplicating the entire reference changes support, not its distribution.
    duplicate = data.with_columns((pl.col('edge')+'-copy').alias('edge'))
    twice = derive_atom_patterns(pl.concat([data,duplicate]),data.head(0))
    assert [c['value'] for c in twice['cells']] == pytest.approx([c['value'] for c in r['cells']])
    with pytest.raises(ValueError, match='duplicate'):
        derive_atom_patterns(pl.concat([data,data.head(1)]),data.head(0))
    bad = dict(r, cells=[dict(c) for c in r['cells']]); bad['cells'][0]['mass'] += .1
    with pytest.raises(ValueError):
        validate_resource(bad)
    contacts = pl.DataFrame({'chain.id.from':['A'], 'residue.index.from':[1],
        'chain.id.to':['C'], 'residue.index.to':[2], 'region.type.from':['CDR3'], 'original_contribution':[-.2]})
    obs = pl.DataFrame({'chain.id.from':['A']*2, 'residue.index.from':[1]*2,
        'chain.id.to':['C']*2, 'residue.index.to':[2]*2, 'role.from':['tcr']*2,
        'chem.from':['Y:OH','Y:XX'], 'chem.to':['P:O']*2})
    detail, score = score_atom_counts(contacts,obs,r,-.2)
    assert score['atom_pattern'] == pytest.approx(-.2-math.log(4/3)/2)
    assert score['atom_pattern_unsupported_atoms'] == 1
    assert 'atom_count_increment' not in detail.columns
    # Chemical equivalents retain both physical observations and one edge unit.
    eq = {'D:OD':['D:OD1','D:OD2']}
    d = pl.DataFrame(dict(edge=['x','x'],aa_from=['D']*2,atom_from=['OD1','OD2'],aa_to=['R']*2,atom_to=['NE']*2))
    merged = derive_atom_patterns(d,d.head(0),eq)
    assert len(merged['cells']) == 1
    assert merged['cells'][0]['count'] == 2 and merged['cells'][0]['mass'] == 1.


def test_freeze_with_string_paths_and_equivalence_provenance(tmp_path):
    import json
    from tcren.contacts.atom_patterns import _freeze
    n = pl.DataFrame({'pdb.id':['x'], 'mhc_class':['I'], 'interface':['tcr_peptide'],
        'role.from':['tcr'], 'chain.id.from':['A'], 'chain.id.to':['C'],
        'residue.index.from':[1], 'residue.index.to':[2],
        'chem.from':['D:OD1'],'chem.to':['R:NE'],
        'region.type.from':['CDR3'],'region.type.to':['PEPTIDE']})
    n.write_parquet(tmp_path/'n.parquet')
    e = pl.DataFrame(dict(pdb=['y'],instance=['i'],model=['1'],ab_chain=['H'],
        ab_auth=['1'],ab_ins=[''],target_chain=['T'],target_auth=['2'],target_ins=[''],
        target_type=['PEPTIDE'],ab_aa=['D'],ab_atom=['OD2'],target_aa=['R'],target_atom=['NE']))
    e.write_parquet(tmp_path/'e.parquet')
    (tmp_path/'eq.json').write_text(json.dumps({'chemical_equivalence':{'D:OD':['D:OD1','D:OD2']}}))
    _freeze(*[str(tmp_path/x) for x in ['n.parquet','e.parquet','out.json','eq.json']])
    result = json.loads((tmp_path/'out.json').read_text())
    assert result['reference_counts'] == dict(native_atom_pairs=1,native_edges=1,external_atom_pairs=1,external_edges=1)
    assert result['cells'][0]['mass'] == 2 and result['sources']['equivalence_resource_sha256']
    with pytest.raises(FileExistsError):
        _freeze(*[str(tmp_path/x) for x in ['n.parquet','e.parquet','out.json','eq.json']])
