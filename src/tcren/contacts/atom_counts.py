"""Frozen directed residue-conditional atom-contact enrichment, no assay fitting.

Raw native and antibody counts are pooled with their actual observation counts.
Endpoint frequencies describe observed contact opportunities, not exposed atoms.
One mean residual per physical residue edge prevents linear atom multiplicity.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import polars as pl
import typer

from ..runtime import single_threaded

VERSION = 'directed-atom-count-1'
TYPES = ['aa_from', 'atom_from', 'aa_to', 'atom_to']
EDGE = ['chain.id.from', 'residue.index.from', 'chain.id.to', 'residue.index.to']


def derive_atom_counts(native, external, equivalences=None):
    """Freeze D=-ln[P(i,j|a,b)/(qB(i|a)qT(j|b))] on observed directed cells.

    Inputs contain TYPES and positive integer count. Unseen cells have no
    estimated value; application uses an explicit neutral missing-support policy.
    Chemical equivalence changes category labels, never physical observation IDs.
    """
    equivalences = equivalences or {}
    aliases = {raw: key for key, values in equivalences.items() for raw in values}
    frames = []
    for label, data in [('native', native), ('external', external)]:
        if not set([*TYPES, 'count']) <= set(data.columns):
            raise ValueError('directed atom types and count required')
        if data.select(pl.any_horizontal(pl.col([*TYPES, 'count']).is_null()).any()).item():
            raise ValueError('nonnull directed atom counts required')
        for side in ['from','to']:
            if data.filter(~pl.col('aa_'+side).is_in(list('ACDEFGHIKLMNPQRSTVWY')) |
                           (pl.col('atom_'+side).str.len_chars()==0) |
                           pl.col('atom_'+side).str.contains(':')).height:
                raise ValueError('standard amino acids and nonempty atom names required')
        values = data['count'].to_numpy()
        if not np.isfinite(values).all() or (values <= 0).any() or (values != np.floor(values)).any():
            raise ValueError('positive integer counts required')
        data = data.with_columns(*[(pl.col('aa_'+s)+':'+pl.col('atom_'+s)).replace(aliases).str.split_exact(':',1)
            .struct.field('field_1').alias('atom_'+s) for s in ['from','to']])
        frames.append(data.group_by(TYPES).agg(pl.col('count').sum().cast(pl.Int64).alias(label+'_count')))
    counts = frames[0].join(frames[1], on=TYPES, how='full', coalesce=True, validate='1:1').with_columns(
        pl.col('native_count','external_count').fill_null(0)).with_columns(
            (pl.col('native_count')+pl.col('external_count')).alias('count'))
    if counts.is_empty():
        raise ValueError('nonempty reference required')
    counts = counts.with_columns(
        pl.col('count').sum().over('aa_from','aa_to').alias('residue_pair_count'),
        (pl.col('count').sum().over('aa_from','atom_from')/pl.col('count').sum().over('aa_from')).alias('q_binder'),
        (pl.col('count').sum().over('aa_to','atom_to')/pl.col('count').sum().over('aa_to')).alias('q_target'))
    counts = counts.with_columns((-(pl.col('count')/pl.col('residue_pair_count')/
        (pl.col('q_binder')*pl.col('q_target'))).log()).alias('value')).sort(TYPES)
    return dict(version=VERSION, cutoff_A=5., direction='binder_to_target', mhc_class='I',
        formula='-ln[(n_xy/N_ab)/(qB(i|a)*qT(j|b))]',
        pooling='native+external actual integer atom-contact counts; correlated observations, no assay weighting',
        aggregation='mean over all physical atom pairs per known CDR residue edge; one unit coefficient',
        missing='unseen atom cells: zero correction, retained denominator; not established compatibility',
        interpretation='residue-conditional atom-contact enrichment, not isolated chemistry, joint likelihood or free energy',
        chemical_equivalence=equivalences, reference_counts=dict(native_atom_pairs=int(native['count'].sum()),
            external_atom_pairs=int(external['count'].sum())), cells=counts.to_dicts())


def score_atom_counts(contacts, observations, resource, original, *, mhc_class="I"):
    """Add one mean atom residual to each supported original CDR contribution."""
    pattern = resource.get('version') == 'within-residue-atom-pattern-1'
    if (not pattern and resource.get('version') != VERSION) or resource.get('cutoff_A') != 5.:
        raise ValueError('unsupported directed atom-count resource')
    if mhc_class != resource.get('mhc_class'):
        raise ValueError('atom-count reference supports class I only')
    if contacts.select(EDGE).is_duplicated().any():
        raise ValueError('duplicate physical residue contacts')
    if not np.isclose(float(contacts['original_contribution'].sum()), original, rtol=0, atol=1e-12):
        raise ValueError('original contribution mismatch')
    swap = pl.col('role.from') != 'tcr'
    atoms = observations.select(*[
        pl.when(swap).then(pl.col(k+'.'+other)).otherwise(pl.col(k+'.'+side)).alias(k+'.'+side)
        for side,other in [('from','to'),('to','from')] for k in ['chain.id','residue.index','chem']])
    if atoms.select(*EDGE,'chem.from','chem.to').is_duplicated().any():
        raise ValueError('duplicate physical atom observations')
    aliases = {raw: key for key, values in resource['chemical_equivalence'].items() for raw in values}
    atoms = atoms.with_columns(*[
        pl.col('chem.'+side).replace(aliases).str.split_exact(':',1).struct.field(field).alias(dest+'_'+side)
        for side in ['from','to'] for field,dest in [('field_0','aa'),('field_1','atom')]])
    matrix = (pl.DataFrame(resource['cells']).select(*TYPES,'value') if resource['cells'] else
              pl.DataFrame(schema={**dict.fromkeys(TYPES,pl.String),'value':pl.Float64}))
    if matrix.select(TYPES).is_duplicated().any() or not matrix['value'].is_finite().all():
        raise ValueError('unique finite atom-count cells required')
    atoms = atoms.join(matrix,on=TYPES,how='left',validate='m:1')
    reduced = atoms.group_by(EDGE).agg(pl.len().alias('atom_count_pairs'),
        pl.col('value').null_count().alias('atom_count_unsupported_pairs'),
        pl.col('value').fill_null(0).mean().alias('atom_count_increment'))
    result = contacts.join(reduced,on=EDGE,how='left',validate='1:1')
    cdr = pl.col('region.type.from').is_in(['CDR1','CDR2','CDR3']).fill_null(False)
    known = pl.col('original_contribution').is_finite().fill_null(False)
    if result.filter(cdr & known & pl.col('atom_count_increment').is_null()).height:
        raise ValueError('missing atom observations for scored CDR edge')
    result = result.with_columns(pl.when(cdr & known).then(pl.col('atom_count_increment')).otherwise(0.).alias('atom_count_increment'),
        pl.col('atom_count_pairs','atom_count_unsupported_pairs').fill_null(0))
    result = result.with_columns((pl.col('original_contribution')+pl.col('atom_count_increment')).alias('atom_count_contribution'))
    increment = float(result['atom_count_increment'].sum())
    summary = dict(atom_count=original+increment, atom_count_increment=increment,
        atom_count_pairs=int(result['atom_count_pairs'].sum()),
        atom_count_unsupported_atoms=int(result['atom_count_unsupported_pairs'].sum()))
    if pattern:
        rename = {c:c.replace('atom_count','atom_pattern',1) for c in result.columns if c.startswith('atom_count')}
        result = result.rename(rename)
        summary = {c.replace('atom_count','atom_pattern',1):v for c,v in summary.items()}
    return result, summary


@single_threaded
def fit_command(
    native: Path = typer.Option(..., '--native', exists=True, dir_okay=False),
    external: Path = typer.Option(..., '--external', exists=True, dir_okay=False),
    out: Path = typer.Option(..., '-o', '--out'),
    equivalence_resource: Path | None = typer.Option(None, '--equivalence-resource', exists=True, dir_okay=False),
):
    """Freeze directed atom-count residuals from existing class-I/antibody tables."""
    from ..exports import _sha256
    from datetime import datetime, timezone
    if out.exists():
        raise FileExistsError(out)
    n = pl.read_parquet(native).filter((pl.col('mhc_class')=='I') & (pl.col('interface')=='tcr_peptide'))
    swap = pl.col('role.from') != 'tcr'
    n = n.select(*[pl.when(swap).then(pl.col('chem.'+other)).otherwise(pl.col('chem.'+side)).alias('chem_'+side)
        for side,other in [('from','to'),('to','from')]])
    n = n.with_columns(*[pl.col('chem_'+side).str.split_exact(':',1).struct.field(field).alias(dest+'_'+side)
        for side in ['from','to'] for field,dest in [('field_0','aa'),('field_1','atom')]])
    e = pl.read_parquet(external).filter(pl.col('target_type')=='PEPTIDE').select(
        pl.col('ab_aa').alias('aa_from'),pl.col('ab_atom').alias('atom_from'),
        pl.col('target_aa').alias('aa_to'),pl.col('target_atom').alias('atom_to'))
    eq = json.loads(equivalence_resource.read_text())['chemical_equivalence'] if equivalence_resource else {}
    resource = derive_atom_counts(n.group_by(TYPES).len(name='count'),e.group_by(TYPES).len(name='count'),eq)
    resource['sources'] = dict(native_sha256=_sha256(native),external_sha256=_sha256(external),
        equivalence_resource_sha256=_sha256(equivalence_resource) if equivalence_resource else None,
        producer_sha256=_sha256(Path(__file__)))
    resource['frozen_utc'] = datetime.now(timezone.utc).isoformat()
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(resource,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(resource_sha256=_sha256(out),cells=len(resource['cells']),**resource['reference_counts'])),flush=True)
