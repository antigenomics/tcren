"""Frozen opportunity ranks preserve signed entries, symmetry and missing support."""
import numpy as np
import polars as pl

from tcren.contacts.weighted import STRATA, VERSION, calibrated_weights, multiplicity_edges


def test_inclusive_rank_native_precedence_and_explicit_neutral():
    edges = pl.DataFrame(dict(n_distance=[2,2,2,3,4,2],n_omega=[0]*6,n_theta=[0]*6,n_phi=[0]*6,
        minimum_tail=[.4,.3,.9,.2,.1,None]))
    records = [dict(source=source,n_distance=n,n_omega=0,n_theta=0,n_phi=0,
        minimum_tails=tails,n_edges=len(tails),n_structures=1) for source,n,tails in
        [('native_tcr_peptide',2,[.4,.4,.8]),('external_antibody_antigen',2,[.1]),('external_antibody_antigen',3,[.1,.5])]]
    result = calibrated_weights(edges,dict(version=VERSION,strata=STRATA,records=records))
    np.testing.assert_allclose(result['contact_weight'],[.75,.25,1.,2/3,1.,1.])
    assert result['calibration_source'].to_list() == ['native_tcr_peptide']*3+['external_antibody_antigen','unseen_multiplicity','unavailable_geometry']


def test_angular_multiplicity_and_reversal_are_invariant():
    from test_contact_observations import complex_fixture
    from tcren.contacts.observations import contact_observations
    from tcren.contacts.plausibility import fit_resource
    obs = contact_observations(complex_fixture(),interface='tcr_peptide')
    ref = pl.concat([obs,obs.with_columns((pl.col('dist')*.9).alias('dist'))])
    resource = fit_resource(ref,obs)
    original = multiplicity_edges(obs,resource)
    repeated = multiplicity_edges(pl.concat([obs,obs]),resource)
    np.testing.assert_array_equal(original.select(STRATA[1:]),repeated.select(STRATA[1:]))
    assert original['n_distance'].sum()*2 == repeated['n_distance'].sum()
    pairs = [(c,c[:-5]+'.to') for c in obs.columns if c.endswith('.from') and c[:-5]+'.to' in obs.columns]
    swapped = obs.with_columns(*[pl.col(b).alias(a) for a,b in pairs],*[pl.col(a).alias(b) for a,b in pairs],
        pl.col('theta.to_rad').alias('theta.from_rad'),pl.col('theta.from_rad').alias('theta.to_rad'),
        pl.col('phi.to_rad').alias('phi.from_rad'),pl.col('phi.from_rad').alias('phi.to_rad'))
    reversed_edges = multiplicity_edges(swapped,resource)
    np.testing.assert_array_equal(original.select(STRATA),reversed_edges.select(STRATA))
    np.testing.assert_allclose(original['minimum_tail'],reversed_edges['minimum_tail'])


def test_weighted_signed_identity_and_zero_contacts():
    from test_contact_observations import complex_fixture
    from tcren.contacts.observations import contact_observations
    from tcren.contacts.plausibility import fit_resource, score_contact_penalties
    from tcren.contacts.weighted import fit_calibration, apply_weighted
    from tcren.contactmap import ContactMap
    from tcren.potential import tcren2
    obs = contact_observations(complex_fixture(),interface='tcr_peptide')
    ref = pl.concat([obs,obs.with_columns((pl.col('dist')*.9).alias('dist'))])
    resource = fit_resource(ref,obs)
    # External schema is deliberately different; do not assume author/native numbering.
    external=ref.rename({'chain.id.from':'source_chain.from','residue.index.from':'source_resid.from',
        'chain.id.to':'source_chain.to','residue.index.to':'source_resid.to'})
    calibration=fit_calibration(external,obs,resource)
    for observations in [obs,obs.clear()]:
        contacts=ContactMap.from_structure(complex_fixture()).interface('tcr_peptide')
        if observations.is_empty():
            contacts=contacts.clear()
        table,summary=score_contact_penalties(contacts,observations,tcren2(),resource)
        weights=calibrated_weights(multiplicity_edges(observations,resource),calibration)
        result,summary=apply_weighted(table,weights,summary)
        np.testing.assert_allclose(result['weighted_contribution'],result['original_contribution']*result['contact_weight'])
        assert summary['weighted']==result['weighted_contribution'].sum()
        assert summary['uncalibrated_cdr_contacts']==0


def test_missing_observation_edge_is_explicitly_neutral():
    from tcren.contacts.weighted import apply_weighted, KEYS
    table=pl.DataFrame({'pdb.id':['a']*3,'chain.id.from':['A']*3,'residue.index.from':[1,2,3],
        'chain.id.to':['P']*3,'residue.index.to':[1]*3,
        'geometry_status':['unavailable_cdr_geometry','unassessed_framework','unavailable_cdr_geometry'],
        'original_contribution':[-2.,5.,None]})
    weights=table.select(KEYS).clear().with_columns(pl.lit(1.).alias('contact_weight'),pl.lit('native_tcr_peptide').alias('calibration_source'))
    result,summary=apply_weighted(table,weights,{})
    assert result['calibration_source'].to_list()==['unavailable_geometry','unassessed_framework','unavailable_geometry']
    assert result['weighted_contribution'].to_list()==[-2.,5.,None]
    assert summary['weighted']==3. and summary['uncalibrated_cdr_contacts']==2
