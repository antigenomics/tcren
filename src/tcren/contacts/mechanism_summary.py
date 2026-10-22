"""CLI-owned descriptive audit of accepted contact decompositions and new geometry.

No new scorer, fitted parameter, threshold, significance test or cohort pooling.
The original metric ledger is read only to verify exact rank-change accounting.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import polars as pl
import typer

from .mechanism import EDGE
from .plausibility import ARMS
from ..associations import GROUPS, _unique
from ..runtime import single_threaded


def signed_terms(edges):
    e, u = pl.col('original_contribution'), pl.col('unnaturalness')
    return edges.with_columns(
        (u*(-e).clip(lower_bound=0)).alias('weighted_favourable_magnitude'),
        (u*e.clip(lower_bound=0)).alias('weighted_positive_magnitude'),
        (pl.col('geometry_status')=='supported_cdr_geometry').alias('assessed'),
        e.is_not_null().alias('known_potential'))


def rank_changes(y, original_energy, new_energy):
    """Exact concordance accounting in 32-positive blocks; ties carry half credit."""
    positive, negative = np.flatnonzero(y==1), np.flatnonzero(y==0)
    n = len(positive)*len(negative)
    result = dict(n_pairs=n, gained_pairs=0, lost_pairs=0, gain_credit=0., loss_credit=0.,
                  old_credit=0., new_credit=0.)
    best = {}
    for offset in range(0, len(positive), 32):
        ii = positive[offset:offset+32]
        old = original_energy[negative][None, :]-original_energy[ii, None]
        new = new_energy[negative][None, :]-new_energy[ii, None]
        a = (old>0).astype(float)+.5*(old==0)
        b = (new>0).astype(float)+.5*(new==0)
        change = b-a
        result['old_credit'] += float(a.sum())
        result['new_credit'] += float(b.sum())
        for name, mask in [('gain',change>0),('loss',change<0)]:
            result['gained_pairs' if name=='gain' else 'lost_pairs'] += int(mask.sum())
            result[name+'_credit'] += float(np.abs(change[mask]).sum())
            if mask.any():
                # Explain the strongest new separation among pairs whose ranking changed.
                weight = np.where(mask, np.abs(new), -1.)
                i,j = np.unravel_index(np.argmax(weight), weight.shape)
                candidate = (float(weight[i,j]), int(ii[i]), int(negative[j]))
                if name not in best or candidate[0]>best[name][0]:
                    best[name] = candidate
    result['delta_auc'] = (result['new_credit']-result['old_credit'])/n if n else None
    result['examples'] = {name: [a,b] for name,(_,a,b) in best.items()}
    return result


def _read(path):
    return pl.read_csv(path, separator='\t', infer_schema_length=None)


def _mean_profiles(models, variables):
    return models.unpivot(on=variables, index=[*GROUPS,'potential','y','id'],
        variable_name='quantity', value_name='value').group_by(*GROUPS,'potential','y','quantity').agg(
            pl.len().alias('n_models'), pl.col('value').count().alias('n_available'),
            pl.col('value').sum().alias('total'), pl.col('value').mean().alias('mean'),
            pl.col('value').median().alias('median')).sort(*GROUPS,'potential','quantity','y')


@single_threaded
def audit_tables(run: Path, geometry: Path, out: Path):
    """Audit every accepted row and physical edge; never recalculate structures/scores."""
    from ..exports import _sha256
    labels, scores = _read(run/'labels.tsv'), _read(run/'scoring/scores.tsv')
    edges = _read(run/'scoring/contacts.tsv.gz')
    diagnosed = _read(geometry/'geometry.tsv.gz')
    metrics = _read(run/'comparison/metrics.tsv')
    _unique(labels,['id'],'labels')
    _unique(scores,['id','potential'],'scores')
    _unique(edges,[*EDGE,'potential'],'contacts')
    _unique(diagnosed,EDGE,'geometry')
    if set(labels['id']) != set(scores['id']):
        raise ValueError('accepted model population differs')
    if not set(labels['y']) <= {0,1} or labels['y'].null_count():
        raise ValueError('finite binary labels required')
    for _,part in scores.group_by('potential'):
        if set(part['id']) != set(labels['id']):
            raise ValueError('each potential must retain the complete accepted population')
    if not scores.select(pl.all_horizontal(pl.col(['original',*ARMS]).is_finite()).all()).item() or not labels['iptm'].is_finite().all():
        raise ValueError('accepted energies and ipTM must be finite')

    ledger=_read(geometry/'accounting.tsv')
    accepted_sources=_read(run/'scoring/accounting.tsv')
    _unique(ledger,['id'],'geometry accounting')
    _unique(accepted_sources,['id'],'accepted source accounting')
    if set(accepted_sources['id']) != set(labels['id']) or accepted_sources['sha256'].null_count() or ledger['sha256'].null_count():
        raise ValueError('complete nonnull raw-source identities required in both ledgers')
    if set(ledger['id']) != set(labels['id']) or set(ledger['status']) != {'ok'}:
        raise ValueError('geometry must account successfully for every accepted model, including zero contacts')
    raw_check=ledger.select('id','sha256').join(accepted_sources.select('id','sha256'),on='id',validate='1:1',suffix='_accepted')
    if raw_check.filter(pl.col('sha256')!=pl.col('sha256_accepted')).height:
        raise ValueError('raw geometry sources differ from accepted coordinates')
    if not _read(run/'comparison/membership.tsv')['used'].all():
        raise ValueError('this audit requires the completed all-row matched comparison')
    joined = edges.join(diagnosed,on=EDGE,how='left',validate='m:1')
    assessed = joined.filter(pl.col('geometry_status')=='supported_cdr_geometry')
    if assessed['audit_unnaturalness'].null_count() or not np.allclose(
        assessed['unnaturalness'], assessed['audit_unnaturalness'],rtol=0,atol=2e-15):
        raise ValueError('new diagnostics do not exactly reproduce frozen edge weights')
    source_keys = assessed.select(EDGE).unique()
    if diagnosed.join(source_keys,on=EDGE,how='anti').height:
        raise ValueError('diagnostics contain edges outside accepted CDR graph')
    joined = signed_terms(joined).rename({'pdb.id':'id'}).join(labels,on='id',validate='m:1')
    joined = joined.with_columns(
        (pl.col('assessed') & pl.col('chemical_supported').fill_null(False)).alias('chemical_real'),
        (pl.col('assessed') & ~pl.col('chemical_supported').fill_null(False)).alias('chemical_proximity'),
        pl.when(pl.col('assessed')).then(1-pl.col('unnaturalness')).otherwise(None).alias('geometry_real_equivalent'),
        pl.when(pl.col('assessed')).then(pl.col('unnaturalness')).otherwise(None).alias('geometry_virtual_equivalent'),
        (pl.col('assessed') & pl.col('chemical_supported') & ~pl.col('closest_chemical_supported')).fill_null(False).alias('closest_misses_chemical'),
        pl.when(pl.col('known_potential') & pl.col('assessed')).then(pl.col('unnaturalness')).otherwise(0.).alias('scored_geometry_virtual_equivalent'),
        (pl.col('geometry_status')=='unassessed_framework').alias('unassessed_framework'),
        (pl.col('original_contribution')<0).alias('favourable_entries'),
        (pl.col('original_contribution')>0).alias('positive_entries'))
    sums = ['chemical_real','chemical_proximity','geometry_real_equivalent','geometry_virtual_equivalent',
        'closest_misses_chemical','scored_geometry_virtual_equivalent','weighted_favourable_magnitude','weighted_positive_magnitude',
        'unassessed_framework','assessed','known_potential','favourable_entries','positive_entries']
    summaries = joined.group_by('id','potential').agg(*[pl.col(c).sum().alias(c) for c in sums],
        *[pl.col(c).sum().alias(c) for c in ['n_atom_pairs','n_distinct_marginals']],
        pl.col('unnaturalness').filter(pl.col('assessed')).mean().alias('mean_u'),
        pl.col('unnaturalness').filter(pl.col('assessed')).median().alias('median_u'),
        *[(pl.col('winner_family')==f).sum().alias('winner_'+f+'_edges') for f in ['distance','omega','theta','phi']],
        *[pl.col('original_contribution').filter(pl.col('is_'+t).fill_null(False)).sum().alias('overlapping_energy_'+t)
          for t in ['hydrogen_bond','salt_bridge','hydrophobic','aromatic','polar','vdw','other']])
    models = scores.join(labels,on='id',validate='m:1',suffix='_label').join(summaries,on=['id','potential'],how='left',validate='1:1')
    fill = [*sums,'n_atom_pairs','n_distinct_marginals',*[f'winner_{f}_edges' for f in ['distance','omega','theta','phi']],
            *[c for c in summaries.columns if c.startswith('overlapping_energy_')]]
    models = models.with_columns(pl.col(fill).fill_null(0))
    checks = []
    for arm, expected in [('fixed',pl.col('original')+pl.col('scored_geometry_virtual_equivalent')),
        ('matrix',pl.col('original')+pl.col('matrix_max_abs')*pl.col('scored_geometry_virtual_equivalent')),
        ('flip_favourable',pl.col('original')+2*pl.col('weighted_favourable_magnitude')),
        ('flip',pl.col('original')+2*pl.col('weighted_favourable_magnitude')-2*pl.col('weighted_positive_magnitude'))]:
        error = models.select((pl.col(arm)-expected).abs().max()).item()
        if error>1e-10:
            raise ValueError('signed sum identity failed: '+arm)
        checks.append(dict(identity=arm,max_abs_error=error))
    out.mkdir(parents=True,exist_ok=False)
    models.write_csv(out/'models.tsv',separator='\t')
    variables = ['n_contacts',*fill,'mean_u','median_u','original',*ARMS,'iptm']
    _mean_profiles(models,variables).write_csv(out/'label_profiles.tsv',separator='\t')
    classification = joined.group_by(*GROUPS,'potential','y','union_type','winner_family','geometry_status').agg(
        pl.len().alias('n_edges'),pl.col('id').n_unique().alias('n_models_with_edges'),
        pl.col('unnaturalness').sum().alias('sum_u'),pl.col('original_contribution').sum().alias('sum_original'),
        pl.col('fixed_increment').sum().alias('sum_fixed_increment'),
        pl.col('flip_increment').sum().alias('sum_flip_increment'),
        pl.col('flip_favourable_increment').sum().alias('sum_favourable_flip_increment'))
    classification.write_csv(out/'classification.tsv',separator='\t')
    joined.filter(pl.col('assessed')).group_by(*GROUPS,'potential','y','winner_family','winner_parameter_source').agg(
        pl.len().alias('n_edges'),pl.col('id').n_unique().alias('n_models_with_edges'),
        pl.col('unnaturalness').mean().alias('mean_u'),pl.col('unnaturalness').median().alias('median_u'),
        pl.col('min_distance_A').median().alias('median_closest_distance_A'),
        pl.col('max_distance_A').median().alias('median_largest_distance_A'),
        pl.col('winner_dist').median().alias('median_winner_atom_distance_A'),
        pl.col('n_atom_pairs').median().alias('median_atom_pairs'),
        pl.col('n_distinct_marginals').median().alias('median_distinct_marginals')).write_csv(out/'geometry_profiles.tsv',separator='\t')

    cells = joined.group_by(*GROUPS,'potential','y','residue.aa.from','residue.aa.to').agg(
        pl.len().alias('n_edges'),pl.col('assessed').sum().alias('n_assessed_edges'),
        pl.col('known_potential').sum().alias('n_known_potential_edges'),
        pl.col('chemical_real').sum().alias('n_chemical_retained'),
        pl.col('geometry_real_equivalent').sum().alias('geometry_retained_equivalent'),
        pl.col('original_contribution').first().alias('unchanged_matrix_value'),
        pl.col('unnaturalness').filter(pl.col('assessed')).mean().alias('mean_u'))
    cells=cells.with_columns(
        pl.when(pl.col('n_assessed_edges')>0).then(pl.col('n_chemical_retained')/pl.col('n_assessed_edges')).otherwise(None).alias('chemical_retention_fraction'),
        pl.when(pl.col('n_assessed_edges')>0).then(pl.col('geometry_retained_equivalent')/pl.col('n_assessed_edges')).otherwise(None).alias('geometry_mean_weight'))
    cells.write_csv(out/'matrix_selection.tsv',separator='\t')
    motifs = joined.group_by(*GROUPS,'potential','y','chain.type.from','region.type.from','pos.to','residue.aa.from','residue.aa.to').agg(
        pl.len().alias('n_edges'),pl.col('id').n_unique().alias('n_models_with_motif'),
        *[pl.col(c).sum().alias(c) for c in ['original_contribution','fixed_increment','flip_increment','flip_favourable_increment',
            'weighted_positive_magnitude','weighted_favourable_magnitude','chemical_real','chemical_proximity','geometry_virtual_equivalent']])
    denominators = models.group_by(*GROUPS,'potential','y').len(name='n_models')
    motifs=motifs.join(denominators,on=[*GROUPS,'potential','y'],validate='m:1')
    motif_values=['original_contribution','fixed_increment','flip_increment','flip_favourable_increment','weighted_positive_magnitude','weighted_favourable_magnitude','chemical_real','chemical_proximity','geometry_virtual_equivalent']
    motifs=motifs.with_columns(*[(pl.col(c)/pl.col('n_models')).alias('per_model_'+c) for c in ['n_edges',*motif_values]])
    motifs.write_csv(out/'motifs.tsv',separator='\t')
    motif_keys=[*GROUPS,'potential','chain.type.from','region.type.from','pos.to','residue.aa.from','residue.aa.to']
    means=['per_model_'+c for c in ['n_edges',*motif_values]]
    positive=motifs.filter(pl.col('y')==1).select(*motif_keys,*[pl.col(c).alias('positive_'+c) for c in means])
    negative=motifs.filter(pl.col('y')==0).select(*motif_keys,*[pl.col(c).alias('negative_'+c) for c in means])
    contrasts=positive.join(negative,on=motif_keys,how='full',coalesce=True,validate='1:1').fill_null(0)
    class_support=denominators.group_by(*GROUPS,'potential').agg((pl.col('y')==1).any().alias('has_positive_class'),(pl.col('y')==0).any().alias('has_negative_class'))
    contrasts=contrasts.join(class_support,on=[*GROUPS,'potential'],validate='m:1')
    contrasts=contrasts.with_columns(*[pl.when(pl.col('has_'+side+'_class')).then(pl.col(side+'_'+c)).otherwise(None).alias(side+'_'+c) for side in ['positive','negative'] for c in means])
    contrasts=contrasts.with_columns(*[(pl.col('positive_'+c)-pl.col('negative_'+c)).alias('positive_minus_negative_'+c) for c in means])
    contrasts.write_csv(out/'motif_contrasts.tsv',separator='\t')

    rank_rows, examples = [], []
    for key, part in models.sort('id').group_by(*GROUPS,'potential',maintain_order=True):
        context = dict(zip([*GROUPS,'potential'],key))
        yy = part['y'].to_numpy(); e0=part['original'].to_numpy()
        for reference,arm in [('original',a) for a in [*ARMS,'iptm']]+[('flip_favourable','flip')]:
            old=part[reference].to_numpy()
            new = -part[arm].to_numpy() if arm=='iptm' else part[arm].to_numpy()
            result = rank_changes(yy,old,new)
            for kind,(a,b) in result.pop('examples').items():
                examples.append(dict(**context,reference_arm=reference,arm=arm,change=kind,positive_id=part['id'][a],negative_id=part['id'][b],
                    original_positive=e0[a],original_negative=e0[b],reference_positive_energy=old[a],reference_negative_energy=old[b],
                    new_positive_energy=new[a],new_negative_energy=new[b]))
            ledger=metrics.filter(pl.all_horizontal([pl.col(c)==v for c,v in context.items()]) & (pl.col('metric')=='roc_auc'))
            target=ledger.filter(pl.col('arm')==arm)['value'].item()
            baseline=ledger.filter(pl.col('arm')==reference)['value'].item()
            expected=None if target is None or baseline is None else target-baseline
            actual=result['delta_auc']
            matches = (actual is None and expected is None) or (actual is not None and expected is not None and np.isclose(actual,expected,rtol=0,atol=2e-15))
            if not matches:
                raise ValueError('rank accounting differs from accepted AUC delta')
            rank_rows.append(dict(**context,reference_arm=reference,arm=arm,**result,recorded_delta_auc=expected))
    pl.DataFrame(rank_rows).write_csv(out/'rank_changes.tsv',separator='\t')
    pl.DataFrame(examples).write_csv(out/'case_pairs.tsv',separator='\t')
    joined.write_parquet(out/'contacts.parquet',compression='zstd')
    (out/'VERIFICATION.json').write_text(json.dumps(dict(complete=True,rows=labels.height,
        model_potential_rows=models.height,physical_edges=edges.select(EDGE).unique().height,
        assessed_physical_edges=source_keys.height,diagnosed_edges=diagnosed.height,
        exact_frozen_weights=True,identities=checks,rank_accounting_rows=len(rank_rows),
        unchanged_scores_sha256=_sha256(run/'scoring/scores.tsv'),
        unchanged_labels_sha256=_sha256(run/'labels.tsv'),geometry_sha256=_sha256(geometry/'geometry.tsv.gz'),
        no_fitting=True,no_new_accuracy_metrics=True),indent=2)+'\n')


def mechanism_command(
    run: Path = typer.Option(...,'--run',exists=True,file_okay=False),
    geometry: Path = typer.Option(...,'--geometry',exists=True,file_okay=False),
    out: Path = typer.Option(...,'-o','--out'),
):
    """Explain contact classifications and exact rank changes on accepted rows."""
    audit_tables(run,geometry,out)
