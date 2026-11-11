"""Saturation must count atoms once and never reward incomplete inventories."""
import math

import polars as pl
import pytest

from tcren.contacts.saturation import saturation_counts, fit_reference, score_saturation


def test_atom_coverage_and_probability_bounds():
    inventory = pl.DataFrame({
        'chain.id': ['A'] * 11, 'residue.index': [1]*7+[2,3,3,4],
        'residue.aa': ['R']*7+['A','R','R','G'],
        'atom.name': ['CB','CG','CD','NE','CZ','NH1','NH2','CB','CB','NH1','CA'],
    })
    contacts = pl.DataFrame({
        'chain.id.from': ['A']*9, 'residue.index.from': [1]*8+[2],
        'atom.from': ['CB','CG','CD','NE','CZ','NH1','NH2','NH2','CB'],
        'interface': ['peptide']*9,
    })
    rows = saturation_counts(inventory, contacts)
    r = rows.filter((pl.col('interface')=='union') & (pl.col('residue.index')==1)).row(0,named=True)
    assert (r['n'],r['k'],r['complete']) == (7,7,True)
    assert not rows.filter(pl.col('residue.index')==3)['complete'].any()
    assert rows.filter(pl.col('residue.index')==4)['n'].to_list()==[0,0,0]
    # Two Arg observations, only one fully engaged: Jeffreys p_full = 1/2.
    native = pl.DataFrame({'residue.aa':['R','R','A'], 'interface':['union']*3,
        'n':[7,7,1], 'k':[7,1,1], 'complete':[True]*3})
    ref = fit_reference(native)
    assert ref['cells']['union:R']['p_full']==.5
    query=native.with_columns(pl.lit('x').alias('pdb.id'))
    residues, scores=score_saturation(query,ref)
    assert residues['full_bits'][0]==1
    expected=(1-math.log2(.75))/(2-math.log2(.75))
    assert scores['full'][0]==pytest.approx(expected)
    assert 0 <= scores['binomial'][0] <= 1
    assert scores['full_sum'][0]==pytest.approx((1-math.log2(.75))/(2-math.log2(.75)))
    # Conditional Binomial(1,p) given K>0 is deterministic, contributes no bits.
    assert residues['binomial_bits'][2]==pytest.approx(0,abs=1e-12)
    assert residues['binomial_max_bits'][2]==pytest.approx(0,abs=1e-12)
    zero=query.with_columns(pl.lit(0).alias('k'))
    _, z=score_saturation(zero,ref)
    assert z['full'][0]==z['binomial'][0]==0
    with pytest.raises(ValueError):
        score_saturation(query.with_columns(pl.lit(99).alias('k')),ref)


def test_reference_zeros_and_unsupported_amino_acids():
    native=pl.DataFrame({'residue.aa':['R','R'],'interface':['union']*2,
        'n':[7,7],'k':[7,0],'complete':[True,True]})
    ref=fit_reference(native)
    assert ref['cells']['union:R']['N']==1
    assert ref['atom_probability']['union']['p']==.5
    query=native.with_columns(pl.lit('x').alias('pdb.id'),pl.lit('W').alias('residue.aa'))
    _, scores=score_saturation(query,ref)
    assert scores['unsupported'][0]==1
    assert scores['full'][0]==0


def test_bounded_mixture():
    from tcren.contacts.saturation import bounded_mixture
    assert bounded_mixture([0,1],[0,1],[0,1],2,4).tolist()==[0,1]
    assert bounded_mixture([.5],[.8],[.2],1,2)[0]==pytest.approx(.425)
    with pytest.raises(ValueError):
        bounded_mixture([.5],[.8],[.2],-1,2)


def test_structural_extraction_ignores_alternate_conformer():
    from dataclasses import replace
    from tcren.contacts.saturation import residue_saturation
    from tcren.structure.model import Atom, Chain, RegionMarkup, Residue, Structure
    import numpy as np
    def residue(aa, atoms):
        return Residue(0,1,'',aa,'ALA' if aa=='A' else 'GLY',tuple(atoms))
    def atom(name,x,alt=''):
        return Atom(name,'C',np.array([x,0.,0.]),alt)
    a=residue('A',[atom('CA',20),atom('CB',20,'A'),atom('CB',0,'B')])
    g=residue('G',[atom('CA',25)])
    peptide=residue('A',[atom('CA',1),atom('CB',1)])
    mhc=residue('A',[atom('CA',40),atom('CB',40)])
    chains=[]
    for name,kind,r in [('A','TRA',a),('B','TRB',g),('C','PEPTIDE',peptide),('D','MHCa',mhc)]:
        c=Chain(name,[r],chain_type=kind,chain_supertype='MHCI' if kind=='MHCa' else kind)
        c.regions=[RegionMarkup('CDR3' if kind in ['TRA','TRB'] else kind,0,0,r.aa,[r])]
        chains.append(c)
    s=Structure('fixture',chains)
    out=residue_saturation(s)
    assert out.filter(pl.col('chain.id')=='A')['k'].to_list()==[0,0,0]
    # Moving the selected conformer into range must create complete Ala contact.
    chains[0].residues=[replace(a,atoms=(atom('CA',20),atom('CB',0,'A'),atom('CB',20,'B')))]
    out=residue_saturation(s)
    assert out.filter((pl.col('chain.id')=='A') & (pl.col('interface')=='union'))['k'][0]==1


def test_mock_balance_and_upper_tail():
    from tcren.contacts.saturation import balance_mock_residues
    rows=pl.DataFrame({'pdb.id':['a','b','c'], 'cdr3_pair':['X|Y','X|Y','Z|W'],
        'chain.type':['TRA']*3,'epitope':['E1','E1','E2'],
        'residue.aa':['R']*3,'interface':['union']*3,'n':[7]*3,'k':[1,1,7], 'complete':[True]*3})
    weighted=balance_mock_residues(rows)
    assert weighted['weight'].to_list()==pytest.approx([.75,.75,1.5])
    ref=fit_reference(weighted)
    assert ref['cells']['union:R']['p_full']==.5
    query=rows.with_columns(pl.Series('k',[1,4,7]))
    scored, summary=score_saturation(query,ref)
    assert scored['tail_bits'][0]==pytest.approx(0)
    assert 0 < scored['tail_bits'][1] < scored['tail_bits'][2]
    assert summary['tail_sum'].is_between(0,1).all()
    with pytest.raises(ValueError):
        fit_reference(rows.with_columns(pl.lit(-1.).alias('weight')))


def test_tail_opportunity_normalization_counts_noncontacts():
    ref=fit_reference(pl.DataFrame({'residue.aa':['R','R'],'interface':['union']*2,
        'n':[7,7],'k':[7,0],'complete':[True,True]}))
    query=pl.DataFrame({'pdb.id':['one','two','two','ala','zero'],
        'residue.aa':['R','R','R','A','R'],'interface':['union']*5,
        'n':[7,7,7,1,7],'k':[7,7,0,1,0],'complete':[True]*5})
    detail,summary=score_saturation(query,ref)
    result=dict(summary.select('pdb.id','tail').iter_rows())
    assert result==pytest.approx({'one':1.,'two':.5,'ala':0.,'zero':0.})
    assert detail.filter(pl.col('pdb.id')=='ala')['tail_bits'][0]==0
    assert detail.filter(pl.col('pdb.id')=='zero')['tail_max_bits'][0]>0
