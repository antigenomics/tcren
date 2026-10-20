"""Mechanism accounting must explain score changes without changing membership."""
import numpy as np
import polars as pl


def test_pair_rank_accounting_and_exact_signed_decomposition():
    from tcren.contacts.mechanism_summary import rank_changes, signed_terms
    # Original energies rank the first positive correctly and the second wrongly.
    base = np.array([-2., 2., 0.])
    changed = np.array([1., -1., 0.])
    result = rank_changes(np.array([1, 1, 0]), base, changed)
    assert result['n_pairs'] == 2
    assert result['gained_pairs'] == result['lost_pairs'] == 1
    assert result['delta_auc'] == 0
    assert rank_changes(np.array([1]),np.array([0.]),np.array([1.]))['delta_auc'] is None
    assert rank_changes(np.array([1, 0]), np.array([0., 0.]), np.array([-1., 0.]))['delta_auc'] == .5
    edges = pl.DataFrame({'original_contribution': [-2., 3., None], 'unnaturalness': [.8, .9, 1.],
        'geometry_status': ['supported_cdr_geometry']*3})
    terms = signed_terms(edges)
    assert terms['weighted_favourable_magnitude'].sum() == 1.6
    assert terms['weighted_positive_magnitude'].sum() == 2.7
    np.testing.assert_allclose(2*(terms['weighted_favourable_magnitude']-terms['weighted_positive_magnitude']).sum(), -2.2)


def test_end_to_end_audit_retains_zero_contact_models(tmp_path):
    from test_contact_observations import complex_fixture
    from tcren.contacts.observations import contact_observations
    from tcren.contactmap import ContactMap
    from tcren.contacts.plausibility import fit_resource, score_contact_penalties, comparison_command
    from tcren.contacts.mechanism import geometry_diagnostics
    from tcren.contacts.mechanism_summary import audit_tables
    from tcren.potential import tcren2
    structure = complex_fixture()
    obs = contact_observations(structure, interface='tcr_peptide')
    ref = pl.concat([obs, obs.with_columns((pl.col('dist')*.95).alias('dist'))])
    resource = fit_resource(ref,obs)
    contacts = ContactMap.from_structure(structure).interface('tcr_peptide')
    edge,summary = score_contact_penalties(contacts,obs,tcren2(),resource)
    _,zero = score_contact_penalties(contacts.clear(),obs.clear(),tcren2(),resource)
    run,geo = tmp_path/'run',tmp_path/'geometry'
    (run/'scoring').mkdir(parents=True)
    geo.mkdir()
    labels = pl.DataFrame(dict(id=['a','b','c'],dataset=['test']*3,task=['ranking']*3,
        cohort=['case']*3,mhc_class=['I']*3,template_coverage=['unknown']*3,y=[1,0,1],iptm=[.9,.2,.5]))
    labels.write_csv(run/'labels.tsv',separator='\t')
    ledger=pl.DataFrame(dict(id=['a','b','c'],sha256=['raw-a','raw-b','raw-c'],status=['ok']*3))
    ledger.write_csv(run/'scoring/accounting.tsv',separator='\t')
    ledger.write_csv(geo/'accounting.tsv',separator='\t')
    pl.DataFrame([dict(id=i,mhc_class='I',**(zero if i=='c' else summary)) for i in ['a','b','c']]).write_csv(run/'scoring/scores.tsv',separator='\t')
    import gzip
    with gzip.open(run/'scoring/contacts.tsv.gz','wt') as stream:
        stream.write(pl.concat([edge.with_columns(pl.lit(i).alias('pdb.id')) for i in ['a','b']]).write_csv(separator='\t'))
    with gzip.open(geo/'geometry.tsv.gz','wt') as stream:
        stream.write(pl.concat([geometry_diagnostics(obs,resource).with_columns(pl.lit(i).alias('pdb.id')) for i in ['a','b']]).write_csv(separator='\t'))
    comparison_command(run/'labels.tsv',run/'scoring/scores.tsv',run/'comparison')
    audit_tables(run,geo,tmp_path/'audit')
    result = pl.read_csv(tmp_path/'audit/models.tsv',separator='\t')
    assert result.height==3
    assert result.filter(pl.col('id')=='c')['n_contacts'].item()==0
    assert result.filter(pl.col('id')=='c')['geometry_virtual_equivalent'].item()==0

    labels.with_columns(pl.lit(1).alias('y')).write_csv(run/'labels.tsv',separator='\t')
    comparison_command(run/'labels.tsv',run/'scoring/scores.tsv',run/'single-comparison')
    import shutil
    for filename in ['metrics.tsv','membership.tsv']:
        shutil.copy2(run/'single-comparison'/filename,run/'comparison'/filename)
    audit_tables(run,geo,tmp_path/'single-audit')
    contrasts=pl.read_csv(tmp_path/'single-audit/motif_contrasts.tsv',separator='\t')
    assert contrasts['positive_minus_negative_per_model_n_edges'].null_count()==contrasts.height
    ledger.head(2).write_csv(run/'scoring/accounting.tsv',separator='\t')
    import pytest
    with pytest.raises(ValueError,match='complete nonnull'):
        audit_tables(run,geo,tmp_path/'invalid-ledger')
