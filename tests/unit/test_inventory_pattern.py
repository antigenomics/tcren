"""Inventory corrections retain physical zeros and normalize nonempty states."""
import itertools

import numpy as np
import polars as pl
import pytest

from tcren.contacts.inventory_patterns import (
    derive_inventory_patterns, pattern_increment, score_inventory_patterns,
    validate_resource, prepare_resource,
)


def reference():
    return pl.DataFrame(dict(aa_from=['A','A'],atom_from=['N','CB'],aa_to=['A','A'],
        atom_to=['CB','CB'],count=[8,2],opportunity=[10,10]))


def test_nonempty_normalization_and_physical_permutations():
    p=np.array([.8,.2,.4]);q=np.array([.3,.3,.3])
    total=0.
    for state in itertools.product([0,1],repeat=3):
        if not any(state):
            continue
        z=np.array(state)
        qz=np.prod(q**z*(1-q)**(1-z))/(1-np.prod(1-q))
        total += qz*np.exp(-pattern_increment(z,p,q))
        assert pattern_increment(z,p,q) == pytest.approx(pattern_increment(z[::-1],p[::-1],q[::-1]))
    assert total == pytest.approx(1.)
    assert pattern_increment([1],[.2],[.8]) == pytest.approx(0.)
    assert pattern_increment([1,0],[.8,.2],[.5,.5]) != pytest.approx(pattern_increment([1,0],[.8,.8],[.5,.5]))


def test_inventory_resource_and_whole_edge_abstention():
    data=reference(); empty=data.head(0)
    resource=derive_inventory_patterns(data,empty)
    assert resource['cells'][0]['q'] == pytest.approx(.5)
    validate_resource(resource)
    contacts=pl.DataFrame({'chain.id.from':['A'],'residue.index.from':[1],
        'chain.id.to':['P'],'residue.index.to':[2],'region.type.from':['CDR3'],
        'original_contribution':[-2.]})
    obs=pl.DataFrame({'role.from':['tcr'],'chain.id.from':['A'],'residue.index.from':[1],
        'chem.from':['A:N'],'chain.id.to':['P'],'residue.index.to':[2],'chem.to':['A:CB']})
    inventory=pl.DataFrame({'chain.id':['A','A','P'],'residue.index':[1,1,2],
        'residue.aa':['A','A','A'],'atom.name':['N','CB','CB'],'contact_eligible':[True]*3})
    _, summary=score_inventory_patterns(contacts,obs,inventory,resource,-2.,prepared=prepare_resource(resource))
    expected=pattern_increment([1,0],[8.5/11,2.5/11],[.5,.5])
    assert summary['inventory_pattern'] == pytest.approx(-2+expected)
    assert summary['inventory_pattern_pairs']==2
    assert sum(summary['inventory_pattern_'+part+'_increment'] for part in ['contact','noncontact','normalizer']) == pytest.approx(expected)
    missing=inventory.with_columns(pl.when(pl.col('atom.name')=='CB').then(pl.lit('O')).otherwise(pl.col('atom.name')).alias('atom.name'))
    # Preserve contacted target; only a noncontact binder cell becomes unsupported.
    missing=missing.with_columns(pl.when(pl.col('chain.id')=='P').then(pl.lit('CB')).otherwise(pl.col('atom.name')).alias('atom.name'))
    _, summary=score_inventory_patterns(contacts,obs,missing,resource,-2.)
    assert summary['inventory_pattern']==-2.
    assert summary['inventory_pattern_unsupported_edges']==1
    bad=derive_inventory_patterns(data,empty)
    bad['cells'][0]['p']=.999
    with pytest.raises(ValueError,match='probability'):
        validate_resource(bad)


def test_equivalence_pools_counts_not_physical_opportunities():
    data=pl.DataFrame(dict(aa_from=['R','R'],atom_from=['NH1','NH2'],aa_to=['A','A'],
        atom_to=['O','O'],count=[3,1],opportunity=[5,5]))
    resource=derive_inventory_patterns(data,data.head(0),{'R:NH':['R:NH1','R:NH2']})
    cell=resource['cells'][0]
    assert cell['count']==4
    assert cell['opportunity']==10
    assert cell['p']==pytest.approx(4.5/11)
