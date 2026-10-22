"""Directed atom residuals preserve the original residue scorer."""
import math

import polars as pl
import pytest

from tcren.contacts.atom_counts import derive_atom_counts, score_atom_counts


def test_directed_conditional_odds_and_edge_average():
    # One residue pair: endpoint frequencies 1/2 each; joint diagonal3/8.
    rows = [dict(aa_from='Y', atom_from=i, aa_to='P', atom_to=j, count=n)
            for i,j,n in [('OH','O',3),('OH','CB',1),('CZ','O',1),('CZ','CB',3)]]
    native = pl.DataFrame(rows)
    resource = derive_atom_counts(native, native.head(0))
    lookup = {(r['aa_from'],r['atom_from'],r['aa_to'],r['atom_to']):r['value'] for r in resource['cells']}
    assert lookup['Y','OH','P','O'] == pytest.approx(-math.log(1.5))
    assert ('P','O','Y','OH') not in lookup
    contacts = pl.DataFrame({'chain.id.from':['A'], 'residue.index.from':[1],
        'chain.id.to':['C'], 'residue.index.to':[2], 'region.type.from':['CDR3'], 'original_contribution':[-.2]})
    obs = pl.DataFrame({'chain.id.from':['A','A','A'], 'residue.index.from':[1,1,1],
        'chain.id.to':['C','C','C'], 'residue.index.to':[2,2,2], 'role.from':['tcr']*3,
        'role.to':['peptide']*3, 'chem.from':['Y:OH','Y:CZ','Y:XX'], 'chem.to':['P:O','P:O','P:O']})
    edges, summary = score_atom_counts(contacts,obs,resource, -.2)
    expected = (-math.log(1.5)+math.log(2))/3
    assert edges['atom_count_increment'].item() == pytest.approx(expected)
    assert summary['atom_count'] == pytest.approx(-.2+expected)
    assert summary['atom_count_unsupported_atoms'] == 1
    with pytest.raises(ValueError,match='duplicate'):
        score_atom_counts(contacts,pl.concat([obs,obs.head(1)]),resource,-.2)
    with pytest.raises(ValueError,match='missing'):
        score_atom_counts(contacts,obs.head(0),resource,-.2)
    with pytest.raises(ValueError,match='positive integer'):
        derive_atom_counts(native.with_columns(pl.lit(-1).alias('count')),native.head(0))
    _, empty = score_atom_counts(contacts.head(0),obs.head(0),resource,0.)
    assert empty['atom_count'] == empty['atom_count_pairs'] == 0
    with pytest.raises(ValueError,match='class I only'):
        score_atom_counts(contacts.head(0),obs.head(0),resource,0.,mhc_class='II')


def test_equivalent_chemistry_preserves_physical_atom_observations():
    data = pl.DataFrame(dict(aa_from=['D','D'],atom_from=['OD1','OD2'],
        aa_to=['R','R'],atom_to=['NH1','NH2'],count=[2,2]))
    resource = derive_atom_counts(data,data.head(0),{'D:OD':['D:OD1','D:OD2'],'R:NH':['R:NH1','R:NH2']})
    assert len(resource['cells']) == 1
    assert resource['cells'][0]['count'] == 4
    assert resource['cells'][0]['value'] == 0
