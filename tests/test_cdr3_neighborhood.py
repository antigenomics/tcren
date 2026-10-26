"""All-atom witnesses, separate interfaces and fixed prepared corrections."""
import numpy as np
import polars as pl
import pytest

from tcren.topology.neighborhood import neighborhood_counts, neighborhood_scores, RULES


def atoms():
    return pl.DataFrame({
        'chain.id.from':['A','A','B'], 'residue.index.from':[95,95,96],
        'chain.id.to':['C','C','D'], 'residue.index.to':[6,6,64],
        'chain.type.from':['TRA','TRA','TRB'], 'region.type.from':['CDR3']*3,
        'residue.aa.from':['V','V','R'], 'residue.aa.to':['A','A','K'],
        'atom.from':['O','CG1','NH1'], 'atom.to':['O','CB','NZ'],
        'dist':[3.0,3.5,3.8], 'interface':['peptide','peptide','mhc']})


def test_all_atom_support_and_separate_mhc_conflict():
    counts=neighborhood_counts(atoms())
    assert counts['peptide_edges']==1
    assert counts['peptide_supported']==1  # The closest O–O witness alone misses this.
    assert counts['mhc_conflicted']==1
    assert counts['mhc_unanchored_conflicted']==1
    assert counts['pmhc_chains']==2
    assert counts['chemical_pmhc_chains']==1


def test_peptide_anchor_protects_only_local_mhc_rule():
    rows=atoms().to_dicts()
    rows.append(dict(rows[1],**{'chain.id.from':'B','residue.index.from':96,
                               'chain.type.from':'TRB'}))
    counts=neighborhood_counts(pl.DataFrame(rows))
    assert counts['mhc_conflicted']==1
    assert counts['mhc_unanchored_conflicted']==0


def test_zero_interfaces_are_finite_and_prepared_is_reproduced():
    counts=neighborhood_counts(atoms().head(0))
    models=pl.DataFrame([{'pdb.id':'empty','iptm':.8,'topology':.4,'confidence_weight':.55,**counts}])
    result=neighborhood_scores(models)
    assert result['prepared'][0]==pytest.approx(.62)
    assert result['gate_no_cdr3_peptide'][0]==0
    assert result['pep_conflict'][0]==pytest.approx(.62)
    assert np.isfinite(result.select(['prepared']+list(RULES)).to_numpy()).all()


def test_invalid_atom_distance_fails_explicitly():
    with pytest.raises(ValueError,match='distance'):
        neighborhood_counts(atoms().with_columns(pl.lit(float('nan')).alias('dist')))


def test_fixed_penalties_preserve_confidence_component():
    counts=neighborhood_counts(atoms())
    row={'pdb.id':'one','iptm':.8,'topology':.4,'confidence_weight':.55,**counts}
    result=neighborhood_scores(pl.DataFrame([row]))
    assert result['mhc_conflict'][0]==pytest.approx(.44)
    assert result['pep_support'][0]==pytest.approx(.62)
    assert result['pep_packing'][0]<=result['prepared'][0]


def test_cli_preserves_models_and_records_fixed_rules(tmp_path):
    from typer.testing import CliRunner
    from tcren.cli import app
    from tcren.provenance import stamp
    features=tmp_path/'features.tsv'
    models=tmp_path/'models.tsv'
    output=tmp_path/'scores.tsv'
    pl.DataFrame([{'pdb.id':'one',**neighborhood_counts(atoms())}]).write_csv(features,separator='\t')
    stamp(features,command='export',extra={'quantity':'cdr3-neighborhood'})
    pl.DataFrame({'pdb.id':['one'],'iptm':[.8],'topology':[.4],'confidence_weight':[.55]}).write_csv(models,separator='\t')
    result=CliRunner().invoke(app,['neighborhood-rules','--features',str(features),'--models',str(models),'-o',str(output)])
    assert result.exit_code==0,result.output
    assert pl.read_csv(output,separator='\t')['prepared'][0]==pytest.approx(.62)


def test_atom_order_does_not_change_rule_counts():
    frame=atoms()
    assert neighborhood_counts(frame.reverse())==neighborhood_counts(frame)
