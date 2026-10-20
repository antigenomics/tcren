"""Hand-counted chemical core and its physical-edge application."""
import copy
import math

import numpy as np
import polars as pl
import pytest

from tcren.contacts.chemical_core import derive_chemical_core, score_chemical_core, validate_resource


def reference(last=5):
    rows=[dict(aa_from='A',atom_from=i,aa_to='A',atom_to=j,count=n)
          for i,j,n in [('CA','CA',5),('CA','C',3),('C','CA',7),('C','C',last)] if n]
    n=pl.DataFrame(rows)
    return derive_chemical_core(n,n.head(0))


def physical(pairs):
    edge=pl.DataFrame({'chain.id.from':['A'],'residue.index.from':[1],
        'chain.id.to':['C'],'residue.index.to':[2],'region.type.from':['CDR3'],
        'original_contribution':[-.2]})
    obs=pl.DataFrame({'chain.id.from':['A']*len(pairs),'residue.index.from':[1]*len(pairs),
        'chain.id.to':['C']*len(pairs),'residue.index.to':[2]*len(pairs),
        'role.from':['tcr']*len(pairs),'role.to':['peptide']*len(pairs),
        'chem.from':[a for a,b in pairs],'chem.to':[b for a,b in pairs]})
    return edge,obs


def test_hand_counts_partition_gauge_and_original_conservation():
    r=reference();validate_resource(r)
    cells=r['cells'];assert len(cells)==3
    assert sorted(x['count'] for x in cells)==[5,5,10]
    assert sorted(x['value'] for x in cells)==pytest.approx(sorted([-math.log(25/24)]*2+[-math.log(25/26)]))
    assert r['normalizers'][0]['logZ']==pytest.approx(0.,abs=1e-14)
    edge,obs=physical([('A:CA','A:CA'),('A:CA','A:C'),('A:CB','A:C')])
    detail,s=score_chemical_core(edge,obs,r,-.2)
    expected=(-math.log(25/24)-math.log(25/26))/3
    assert s['chemical_core']==pytest.approx(-.2+expected)
    assert detail['original_contribution'].item()==-.2
    assert s['chemical_core_unsupported_atoms']==1
    assert s['chemical_core_unseen_atoms']==0
    _,s2=score_chemical_core(edge.with_columns(pl.lit(.7).alias('original_contribution')),obs,r,.7)
    assert s2['chemical_core_increment']==s['chemical_core_increment']
    # Constant gauge on every completed core cell, not just observed cells.
    q_b=np.array([x['probability'] for x in r['binder']]);q_t=np.array([x['probability'] for x in r['target']])
    w=np.empty((2,2))
    for c in cells:
        w[c['type_a'],c['type_b']]=w[c['type_b'],c['type_a']]=math.exp(-c['value'])
    u=-np.log(w);z=q_b@w@q_t
    assert np.allclose(u+math.log(z),(u+.8)+math.log(q_b@(w*math.exp(-.8))@q_t),atol=1e-14)


def test_unseen_core_is_not_unsupported_zero_and_physical_multiplicity():
    r=reference(0);edge,obs=physical([('A:C','A:C')])
    _,s=score_chemical_core(edge,obs,r,-.2)
    assert s['chemical_core_unseen_atoms']==1 and s['chemical_core_unsupported_atoms']==0
    assert s['chemical_core_increment']==pytest.approx(math.log(1+7/15*3/15))
    with pytest.raises(ValueError,match='duplicate'):
        score_chemical_core(edge,pl.concat([obs,obs]),r,-.2)
    with pytest.raises(ValueError,match='missing'):
        score_chemical_core(edge,obs.head(0),r,-.2)
    with pytest.raises(ValueError,match='class I only'):
        score_chemical_core(edge,obs,r,-.2,mhc_class='II')
    with pytest.raises(ValueError,match='original contribution'):
        score_chemical_core(edge,obs,r,5.)
    _,empty=score_chemical_core(edge.head(0),obs.head(0),r,0.)
    assert empty['chemical_core']==empty['chemical_core_pairs']==0


def test_validation_and_symmetric_pooling_of_roles():
    r=reference();bad=copy.deepcopy(r);bad['normalizers'][0]['logZ']+=.1
    with pytest.raises(ValueError,match='normalizer'):
        validate_resource(bad)
    bad=copy.deepcopy(r);bad['cells'][0]['value']=float('nan')
    with pytest.raises(ValueError,match='finite'):
        validate_resource(bad)
    bad=copy.deepcopy(r);bad['mapping'].append(bad['mapping'][0])
    with pytest.raises(ValueError,match='mapping'):
        validate_resource(bad)
    bad=copy.deepcopy(r);bad['mapping'][0]['type']=1-bad['mapping'][0]['type']
    with pytest.raises(ValueError,match='mapping'):
        validate_resource(bad)
    # A normalized but impossible endpoint/core incidence: shift one binder end.
    bad=copy.deepcopy(r)
    bad['binder'][0].update(count=9,probability=.45)
    bad['binder'][1].update(count=11,probability=.55)
    for cell,q,w in zip(bad['cells'],[.27,.51,.22],[25/27,50/51,25/22],strict=True):
        cell.update(reference_probability=q,value=-math.log(w))
    with pytest.raises(ValueError,match='incidence'):
        validate_resource(bad)
    n=pl.DataFrame([dict(aa_from='D',atom_from='OD1',aa_to='A',atom_to='CA',count=3)])
    with pytest.raises(ValueError,match='positive integer'):
        derive_chemical_core(n.with_columns(pl.lit(.5).alias('count')),n.head(0))
    rr=derive_chemical_core(n,n.head(0))
    edge,obs=physical([('D:OD1','A:CA'),('A:CA','D:OD1')])
    _,s=score_chemical_core(edge,obs,rr,-.2)
    assert s['chemical_core_unsupported_atoms']==1  # Reverse roles have no conditional support.


def test_all_unknown_projection_and_equivalent_physical_atoms():
    r=reference();edge,obs=physical([('A:CB','A:CB')])
    _,s=score_chemical_core(edge,obs,r,-.2)
    assert s['chemical_core']==-.2 and s['chemical_core_unsupported_atoms']==1
    n=pl.DataFrame([dict(aa_from='D',atom_from=a,aa_to='A',atom_to=b,count=2)
                    for a in ['OD1','OD2'] for b in ['CA','C']])
    r=derive_chemical_core(n,n.head(0))
    edge,obs=physical([('D:OD1','A:CA'),('D:OD2','A:CA')])
    _,s=score_chemical_core(edge,obs,r,-.2)
    assert s['chemical_core_pairs']==2 and s['chemical_core_unsupported_atoms']==0


def test_fit_cli_owns_pool_and_freezes_identical_counts(tmp_path):
    import json
    from typer.testing import CliRunner
    from tcren.contacts.chemical_core import fit_command
    from tcren.runtime import resolve_threads
    import typer
    native=tmp_path/'native.parquet';external=tmp_path/'external.parquet'
    rows=dict(aa_from='A',atom_from='CA',aa_to='A',atom_to='C',count=1)
    pl.DataFrame({'pdb.id':['native'],'mhc_class':['I'],'interface':['tcr_peptide'],'role.from':['tcr'],
                  'chem.from':['A:CA'],'chem.to':['A:C']}).write_parquet(native)
    pl.DataFrame({'pdb':['external'],'instance':['external:1'],'target_type':['PEPTIDE'],'ab_aa':['A'],'ab_atom':['C'],
                  'target_aa':['A'],'target_atom':['CA']}).write_parquet(external)
    app=typer.Typer();app.command()(fit_command)
    outputs=[]
    for threads in [1,0]:
        dest=tmp_path/f'resource-{threads}.json'
        args=['--native',str(native),'--external',str(external),'-o',str(dest),'-t',str(threads)]
        result=CliRunner().invoke(app,args)
        assert result.exit_code==0,result.exception
        r=json.loads(dest.read_text());validate_resource(r)
        assert r['polars_threads']==resolve_threads(threads)
        assert r['reference_counts']==dict(native_atom_pairs=1,external_atom_pairs=1)
        assert len(r['sources']['native_sha256'])==64
        outputs.append(r)
        result=CliRunner().invoke(app,args)
        assert result.exit_code!=0 and isinstance(result.exception,FileExistsError)
    for key in ['cells','binder','target','normalizers','mapping','types']:
        assert outputs[0][key]==outputs[1][key]
