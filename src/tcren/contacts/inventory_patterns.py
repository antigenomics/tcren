"""Frozen inventory-aware Bernoulli whole-pattern ratios on contacting edges.

Cell independence is a geometric composite-model assumption. Probabilities are
estimated on contacting residue edges, then both distributions are conditioned
on a nonempty pattern. They are neither binding free energies nor assay fits.
"""
from __future__ import annotations

import math
from pathlib import Path

import typer

import numpy as np
import polars as pl

from .atom_counts import EDGE, TYPES

VERSION = 'inventory-whole-pattern-1'


def _aliases(equivalences):
    aliases = {}
    for canonical, names in equivalences.items():
        for raw in [canonical, *names]:
            if raw.split(':')[0] != canonical.split(':')[0] or (raw in aliases and aliases[raw] != canonical):
                raise ValueError('ambiguous or residue-changing chemical equivalence')
            aliases[raw] = canonical
    return aliases


def _canonical(data, aliases):
    return data.with_columns(*[
        (pl.col('aa_'+s)+':'+pl.col('atom_'+s)).replace(aliases).str.split_exact(':',1)
        .struct.field('field_1').alias('atom_'+s) for s in ['from','to']])


def derive_inventory_patterns(native_counts, external_counts, equivalences=None):
    """Pool integer physical contacts/opportunities, including zero-contact cells."""
    equivalences = equivalences or {}
    aliases = _aliases(equivalences)
    frames = []
    for source, data in [('native',native_counts),('external',external_counts)]:
        data = data.select(*TYPES,'count','opportunity')
        if data.null_count().sum_horizontal().item():
            raise ValueError('nonnull inventory counts required')
        for side in ['from','to']:
            if data.filter(~pl.col('aa_'+side).is_in(list('ACDEFGHIKLMNPQRSTVWY')) |
                           ~pl.col('atom_'+side).str.contains(r'^[A-Z][A-Z0-9]*$')).height:
                raise ValueError('standard residues and atom names required')
        for key in ['count','opportunity']:
            x = data[key].to_numpy()
            if not np.isfinite(x).all() or (x < 0).any() or (x != np.floor(x)).any():
                raise ValueError('nonnegative integer physical counts required')
        if data.filter((pl.col('opportunity') <= 0) | (pl.col('count') > pl.col('opportunity'))).height:
            raise ValueError('positive opportunities and count <= opportunity required')
        frames.append(_canonical(data,aliases).group_by(TYPES).agg(*[
            pl.col(k).sum().cast(pl.Int64).alias(source+'_'+k) for k in ['count','opportunity']]))
    c = frames[0].join(frames[1],on=TYPES,how='full',coalesce=True,validate='1:1').fill_null(0)
    if c.is_empty():
        raise ValueError('nonempty reference required')
    c = c.with_columns(*[(pl.col('native_'+k)+pl.col('external_'+k)).alias(k) for k in ['count','opportunity']])
    c = c.with_columns(((pl.col('count')+.5)/(pl.col('opportunity')+1.)).alias('p'),
        ((pl.col('count').sum().over('aa_from','aa_to')+.5)/
         (pl.col('opportunity').sum().over('aa_from','aa_to')+1.)).alias('q')).sort(TYPES)
    resource = dict(version=VERSION,cutoff_A=5.,mhc_class='I',chemical_equivalence=equivalences,
        formula='-ln[P(z|inventory,K>=1)/Q(z|inventory,K>=1)]',
        estimator='Jeffreys Beta(0.5,0.5) posterior means; independent physical Bernoulli cells',
        pooling='native+external actual integer opportunities and contacts; no source or assay weights',
        support='contacting residue edges only; unsupported whole edge abstains',
        aggregation='one whole-pattern correction per original CDR edge, coefficient one',
        reference_counts={source+'_'+key:int(c[source+'_'+key].sum()) for source in ['native','external'] for key in ['count','opportunity']},cells=c.to_dicts())
    validate_resource(resource)
    return resource


def validate_resource(resource):
    """Reconstruct probabilities from source counts; reject malformed resources."""
    try:
        if resource['version'] != VERSION or resource['cutoff_A'] != 5. or resource['mhc_class'] != 'I':
            raise ValueError('unsupported inventory-pattern resource')
        _aliases(resource['chemical_equivalence'])
        c = pl.DataFrame(resource['cells'])
        if c.is_empty() or c.select(TYPES).is_duplicated().any() or c.null_count().sum_horizontal().item():
            raise ValueError('unique nonnull inventory cells required')
        for side in ['from','to']:
            if c.filter(~pl.col('aa_'+side).is_in(list('ACDEFGHIKLMNPQRSTVWY')) |
                        ~pl.col('atom_'+side).str.contains(r'^[A-Z][A-Z0-9]*$')).height:
                raise ValueError('invalid inventory atom labels')
        for key in ['count','opportunity']:
            for prefix in ['native_','external_','']:
                x = c[prefix+key]
                if x.dtype not in (pl.Int64,pl.Int32) or (x < 0).any():
                    raise ValueError('integer inventory counts required')
            if not (c[key] == c['native_'+key]+c['external_'+key]).all():
                raise ValueError('pooled inventory counts mismatch')
        for source in ['native','external']:
            for key in ['count','opportunity']:
                if resource['reference_counts'][source+'_'+key] != int(c[source+'_'+key].sum()):
                    raise ValueError('reference count conservation failed')
        for prefix in ['native_','external_','']:
            if (c[prefix+'count'] > c[prefix+'opportunity']).any():
                raise ValueError('contacts exceed opportunities')
        if (c['opportunity'] <= 0).any():
            raise ValueError('positive pooled opportunities required')
        expected = c.with_columns(((pl.col('count').sum().over('aa_from','aa_to')+.5)/
            (pl.col('opportunity').sum().over('aa_from','aa_to')+1.)).alias('expected_q'))
        if not np.allclose(c['p'],(c['count']+.5)/(c['opportunity']+1.),rtol=0,atol=1e-12) or not np.allclose(c['q'],expected['expected_q'],rtol=0,atol=1e-12):
            raise ValueError('inventory probability mismatch')
    except (KeyError,TypeError,pl.exceptions.PolarsError) as exc:
        raise ValueError('malformed inventory-pattern resource') from exc


def pattern_increment(z, p, q):
    """Exact nonempty-conditioned ratio; each entry is one physical atom pair."""
    z,p,q = np.asarray(z),np.asarray(p,dtype=float),np.asarray(q,dtype=float)
    if z.ndim != 1 or p.shape != z.shape or q.shape != z.shape or not len(z) or not np.isin(z,[0,1]).all() or not z.any() or not (np.isfinite(p).all() and np.isfinite(q).all()) or not ((p>0)&(p<1)&(q>0)&(q<1)).all():
        raise ValueError('nonempty binary pattern and probabilities strictly between zero and one required')
    return float(-np.sum(z*np.log(p/q)+(1-z)*np.log((1-p)/(1-q)))+
                 math.log(-math.expm1(float(np.log1p(-p).sum())))-
                 math.log(-math.expm1(float(np.log1p(-q).sum()))))


def prepare_resource(resource):
    """Validate once and compile immutable reference inputs for a score batch."""
    validate_resource(resource)
    return pl.DataFrame(resource['cells']).select(*TYPES,'p','q'), _aliases(resource['chemical_equivalence'])


def score_inventory_patterns(contacts, observations, inventory, resource, original, *, mhc_class='I', prepared=None):
    """Add whole-pattern corrections using every observed eligible endpoint atom."""
    matrix, aliases = prepare_resource(resource) if prepared is None else prepared
    if mhc_class != resource['mhc_class']:
        raise ValueError('inventory reference supports class I only')
    if contacts.select(EDGE).is_duplicated().any() or not np.isclose(contacts['original_contribution'].sum(),original,rtol=0,atol=1e-12):
        raise ValueError('unique original edges and conserved contributions required')
    inv = inventory.filter(pl.col('contact_eligible')).select('chain.id','residue.index','residue.aa','atom.name')
    if inv.null_count().sum_horizontal().item() or inv.group_by('chain.id','residue.index').agg(pl.col('residue.aa').n_unique().alias('n')).filter(pl.col('n')!=1).height:
        raise ValueError('nonnull inventory with one amino acid per residue required')
    if inv.select('chain.id','residue.index','atom.name').is_duplicated().any():
        raise ValueError('unique physical inventory atoms required')
    swap = pl.col('role.from') != 'tcr'
    obs = observations.select(*[
        pl.when(swap).then(pl.col(k+'.'+other)).otherwise(pl.col(k+'.'+side)).alias(k+'.'+side)
        for side,other in [('from','to'),('to','from')] for k in ['chain.id','residue.index','chem']])
    if obs.is_duplicated().any():
        raise ValueError('duplicate physical atom contacts')
    obs = obs.with_columns(*[pl.col('chem.'+s).str.split_exact(':',1).struct.field('field_1').alias('raw_'+s) for s in ['from','to']])
    edges = contacts.filter(pl.col('region.type.from').is_in(['CDR1','CDR2','CDR3']) & pl.col('original_contribution').is_finite()).select(EDGE)
    cells = edges
    for side in ['from','to']:
        cells = cells.join(inv.rename({'chain.id':'chain.id.'+side,'residue.index':'residue.index.'+side,
            'residue.aa':'aa_'+side,'atom.name':'raw_'+side}),on=['chain.id.'+side,'residue.index.'+side],how='inner')
    keys = [*EDGE,'raw_from','raw_to']
    selected_obs = obs.join(edges,on=EDGE,how='inner').select(keys)
    if selected_obs.join(cells.select(keys),on=keys,how='anti').height:
        raise ValueError('contact atom absent from eligible inventory')
    observed_cells = obs.join(cells,on=keys,how='inner')
    if observed_cells.filter((pl.col('chem.from') != pl.col('aa_from')+':'+pl.col('raw_from')) |
                             (pl.col('chem.to') != pl.col('aa_to')+':'+pl.col('raw_to'))).height:
        raise ValueError('contact chemistry disagrees with physical inventory')
    cells = cells.join(selected_obs.with_columns(pl.lit(True).alias('z')),on=keys,how='left',validate='1:1').with_columns(pl.col('z').fill_null(False))
    cells = cells.with_columns(*[pl.col('raw_'+s).alias('atom_'+s) for s in ['from','to']])
    cells = _canonical(cells,aliases).join(matrix,on=TYPES,how='left',validate='m:1')
    cells = cells.with_columns((-(pl.col('z').cast(pl.Float64)*(pl.col('p')/pl.col('q')).log())).alias('contact_term'),
        (-((~pl.col('z')).cast(pl.Float64)*((1-pl.col('p'))/(1-pl.col('q'))).log())).alias('noncontact_term'))
    reduced = cells.group_by(EDGE).agg(pl.len().alias('inventory_pattern_pairs'),pl.col('z').sum().alias('K'),
        pl.col('p').null_count().alias('missing'),
        pl.col('contact_term').sum().alias('contact'),pl.col('noncontact_term').sum().alias('noncontact'),
        (1-pl.col('p')).log().sum().alias('lp0'),(1-pl.col('q')).log().sum().alias('lq0'))
    if edges.join(reduced.select(EDGE),on=EDGE,how='anti').height or reduced.filter(pl.col('K')==0).height:
        raise ValueError('missing or empty physical pattern on scored edge')
    supported = reduced['missing'].to_numpy()==0
    normalizer = np.zeros(reduced.height)
    normalizer[supported] = np.log(-np.expm1(reduced['lp0'].to_numpy()[supported])) - np.log(-np.expm1(reduced['lq0'].to_numpy()[supported]))
    reduced = reduced.with_columns(pl.Series('inventory_pattern_contact_increment',np.where(supported,reduced['contact'].to_numpy(),0.)),
        pl.Series('inventory_pattern_noncontact_increment',np.where(supported,reduced['noncontact'].to_numpy(),0.)),
        pl.Series('inventory_pattern_normalizer_increment',normalizer))
    correction = reduced['inventory_pattern_contact_increment'].to_numpy()+reduced['inventory_pattern_noncontact_increment'].to_numpy()+normalizer
    reduced = reduced.with_columns((pl.col('missing')>0).cast(pl.Int64).alias('inventory_pattern_unsupported_edges'),
        pl.Series('inventory_pattern_increment',correction))
    result = contacts.join(reduced.select(*EDGE,'inventory_pattern_pairs','inventory_pattern_unsupported_edges','inventory_pattern_increment','inventory_pattern_contact_increment','inventory_pattern_noncontact_increment','inventory_pattern_normalizer_increment',pl.col('K').alias('inventory_pattern_contact_pairs'),pl.col('missing').alias('inventory_pattern_unsupported_pairs')),on=EDGE,how='left',validate='1:1').with_columns(pl.col('inventory_pattern_pairs','inventory_pattern_contact_pairs','inventory_pattern_unsupported_pairs','inventory_pattern_unsupported_edges').fill_null(0),pl.col('inventory_pattern_increment','inventory_pattern_contact_increment','inventory_pattern_noncontact_increment','inventory_pattern_normalizer_increment').fill_null(0.))
    result = result.with_columns((pl.col('original_contribution')+pl.col('inventory_pattern_increment')).alias('inventory_pattern_contribution'))
    increment = float(result['inventory_pattern_increment'].sum())
    return result,dict(inventory_pattern=original+increment,inventory_pattern_increment=increment,
        **{key:float(result[key].sum()) for key in ['inventory_pattern_contact_increment','inventory_pattern_noncontact_increment','inventory_pattern_normalizer_increment']},
        inventory_pattern_pairs=int(result['inventory_pattern_pairs'].sum()),
        inventory_pattern_contact_pairs=int(result['inventory_pattern_contact_pairs'].sum()),
        inventory_pattern_unsupported_pairs=int(result['inventory_pattern_unsupported_pairs'].sum()),
        inventory_pattern_unsupported_edges=int(result['inventory_pattern_unsupported_edges'].sum()))


def _freeze(native, external, out, equivalence_resource=None):
    """Freeze reference-only P/Q before any assay evaluation."""
    from datetime import datetime, timezone
    import json
    from pathlib import Path
    from ..exports import _sha256

    native, external, out = Path(native), Path(external), Path(out)
    if out.exists():
        raise FileExistsError(out)
    eqpath = Path(equivalence_resource) if equivalence_resource else None
    sources = dict(native_sha256=_sha256(native),external_sha256=_sha256(external),
        equivalence_resource_sha256=_sha256(eqpath) if eqpath else None)
    eq = json.loads(eqpath.read_text())['chemical_equivalence'] if eqpath else {}
    resource = derive_inventory_patterns(pl.read_parquet(native),pl.read_parquet(external),eq)
    if sources != dict(native_sha256=_sha256(native),external_sha256=_sha256(external),
        equivalence_resource_sha256=_sha256(eqpath) if eqpath else None):
        raise ValueError('reference input changed during freeze')
    resource.update(sources=dict(**sources,producer_sha256=_sha256(Path(__file__))),
        frozen_utc=datetime.now(timezone.utc).isoformat(),polars_threads=pl.thread_pool_size())
    out.parent.mkdir(parents=True,exist_ok=True)
    with out.open('x') as stream:
        stream.write(json.dumps(resource,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(resource_sha256=_sha256(out),cells=len(resource['cells']))),flush=True)


def fit_command(native: Path=typer.Option(...,'--native',exists=True,dir_okay=False),
                external: Path=typer.Option(...,'--external',exists=True,dir_okay=False),
                out: Path=typer.Option(...,'-o','--out'),
                equivalence_resource: Path|None=typer.Option(None,'--equivalence-resource',exists=True,dir_okay=False),
                threads: int=typer.Option(0,'-t','--threads')):
    """Freeze inventory-aware native/antibody P/Q, with no assay inputs."""
    import json
    import os
    import subprocess
    import sys
    from ..runtime import resolve_threads

    budget = resolve_threads(threads)
    if out.exists():
        raise FileExistsError(out)
    env = os.environ.copy()
    for key in ['OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS','RAYON_NUM_THREADS']:
        env[key] = '1'
    env.update(POLARS_MAX_THREADS=str(budget),POLARS_ASYNC_THREAD_COUNT='1',POLARS_MAX_BLOCKING_THREAD_COUNT='4')
    args = dict(native=str(native.resolve()),external=str(external.resolve()),out=str(out.resolve()),
        equivalence_resource=str(equivalence_resource.resolve()) if equivalence_resource else None)
    code = ('import polars as pl; pl.thread_pool_size(); import json,sys; '
            'from tcren.contacts.inventory_patterns import _freeze; _freeze(**json.loads(sys.argv[1]))')
    subprocess.run([sys.executable,'-c',code,json.dumps(args)],env=env,check=True)
