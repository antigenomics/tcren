"""Historical RNG order, anchor positions, observed labels and identity-preserving shuffle."""
import numpy as np
import polars as pl
import pytest
from tcren.historical_candidates import AA, generate_candidates


def test_anchored_random_keeps_exact_rng_and_shuffle_identities(tmp_path):
    path=tmp_path/'markup.csv'
    pl.DataFrame({'pdb.id':['a','b'],'peptide':['AAAA','AAAA'],'cdr3a':['SS','SS'],'cdr3b':['TT','TT']}).write_csv(path)
    spec={'markup':str(path),'ids':['a','b']}
    rows,_=generate_candidates('anchored-random',spec)
    assert rows.height==2002
    rng=np.random.default_rng(0);rng.shuffle(np.arange(400))
    first=rng.integers(0,20,size=(1000,4));first[:,1]=AA.index('A');first[:,-1]=AA.index('A')
    assert rows['peptide'][1:1001].to_list()==[''.join(AA[k] for k in row) for row in first]
    shuffled,_=generate_candidates('receptor-shuffle',spec)
    assert shuffled.select('complex.id','candidate.id','y').rows()==[('a','a',1),('a','b',0),('b','a',0),('b','b',1)]


def test_observed_union_keeps_cognate_once_and_labels_exactly(tmp_path):
    observations=tmp_path/'observations.csv';cognates=tmp_path/'cognates.csv'
    pl.DataFrame({'peptide':['AAAA','CCCC','bad']}).write_csv(observations)
    pl.DataFrame({'id':['a'],'cognate':['AAAA']}).write_csv(cognates)
    rows,_=generate_candidates('observed-union',dict(observations=str(observations),cognates=str(cognates),
        peptide_column='peptide',id_column='id',cognate_column='cognate',ids=['a']))
    assert rows.select('peptide','y').rows()==[('AAAA',1),('CCCC',0)]
    with pytest.raises(ValueError,match='Unknown'):
        generate_candidates('new-control',{})


@pytest.mark.parametrize('suffix', ['.csv', '.tsv.gz', '.parquet'])
def test_observed_union_preserves_numeric_looking_ids_and_protocol_order(tmp_path, suffix):
    import gzip

    observations=tmp_path/'observations.csv'
    pl.DataFrame({'peptide':['AAAA','DDDD','bad']}).write_csv(observations)
    cognates=tmp_path/('cognates'+suffix)
    table=pl.DataFrame({'id':['26990','002'],'cognate':['AAAA','CCCC']})
    if suffix == '.parquet':
        table.write_parquet(cognates)
    elif suffix == '.tsv.gz':
        cognates.write_bytes(gzip.compress(table.write_csv(separator='\t').encode()))
    else:
        table.write_csv(cognates)
    spec=dict(observations=str(observations),cognates=str(cognates),peptide_column='peptide',
              id_column='id',cognate_column='cognate',ids=['002','26990'])
    rows,_=generate_candidates('observed-union',spec)
    assert rows.rows()==[
        ('002','0','AAAA',0),('002','1','CCCC',1),('002','2','DDDD',0),
        ('26990','3','AAAA',1),('26990','4','DDDD',0)]
    spec['ids']=[26990]
    rows,_=generate_candidates('observed-union',spec)
    assert rows.rows()==[(26990,'0','AAAA',1),(26990,'1','DDDD',0)]
