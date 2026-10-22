"""Edge-balanced, directed within-residue atom marks from fixed references."""
from __future__ import annotations

import json
import math
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
import polars as pl
import typer

from .atom_counts import TYPES

VERSION = 'within-residue-atom-pattern-1'


def derive_atom_patterns(native, external, equivalences=None):
    """Every physical residue edge contributes one unit, split over its atom pairs."""
    equivalences = equivalences or {}
    aliases = {raw:key for key, values in equivalences.items() for raw in values}
    if any(raw.split(':')[0] != key.split(':')[0] for raw,key in aliases.items()):
        raise ValueError('chemical equivalence must preserve residue identity')
    frames, reference = [], {}
    for name, data in [('native',native),('external',external)]:
        data = data.select('edge',*TYPES)
        if data.null_count().sum_horizontal().item():
            raise ValueError('nonnull physical edges and atom labels required')
        for s in ['from','to']:
            if data.filter(~pl.col('aa_'+s).is_in(list('ACDEFGHIKLMNPQRSTVWY')) |
                           ~pl.col('atom_'+s).str.contains(r'^[A-Z][A-Z0-9]*$')).height:
                raise ValueError('standard amino acids and atom names required')
        if data.select('edge',*TYPES).is_duplicated().any():
            raise ValueError('duplicate physical atom observations')
        if data.select('edge','aa_from','aa_to').unique().group_by('edge').len().filter(pl.col('len') != 1).height:
            raise ValueError('one residue identity pair per physical edge required')
        reference[name+'_atom_pairs'] = data.height
        reference[name+'_edges'] = data['edge'].n_unique()
        data = data.with_columns((1./pl.len().over('edge')).alias('mass'),*[
            (pl.col('aa_'+s)+':'+pl.col('atom_'+s)).replace(aliases).str.split_exact(':',1)
            .struct.field('field_1').alias('atom_'+s) for s in ['from','to']])
        frames.append(data.group_by(TYPES).agg(pl.len().cast(pl.Int64).alias(name+'_count'),
                                               pl.col('mass').sum().alias(name+'_mass')))
    counts = frames[0].join(frames[1],on=TYPES,how='full',coalesce=True,validate='1:1').with_columns(
        pl.col('native_count','external_count').fill_null(0),pl.col('native_mass','external_mass').fill_null(0.))
    if counts.is_empty():
        raise ValueError('nonempty reference required')
    counts = counts.with_columns((pl.col('native_count')+pl.col('external_count')).alias('count'),
                                 (pl.col('native_mass')+pl.col('external_mass')).alias('mass'))
    counts = counts.with_columns(pl.col('mass').sum().over('aa_from','aa_to').alias('residue_pair_mass'),
        (pl.col('mass').sum().over('aa_from','atom_from')/pl.col('mass').sum().over('aa_from')).alias('q_binder'),
        (pl.col('mass').sum().over('aa_to','atom_to')/pl.col('mass').sum().over('aa_to')).alias('q_target'))
    counts = counts.with_columns((-(pl.col('mass')/pl.col('residue_pair_mass')/
        (pl.col('q_binder')*pl.col('q_target'))).log()).alias('value')).sort(TYPES)
    resource = dict(version=VERSION,cutoff_A=5.,mhc_class='I',direction='binder_to_target',
        formula='-ln[(sum_edges n_abij/K_edge / N_edges_ab)/(qB(i|a)*qT(j|b))]',
        pooling='native+external actual physical residue edges, each total mass one; no assay/source weights',
        aggregation='mean over physical atom pairs per known CDR residue edge; unit coefficient',
        missing='unseen atom cells: neutral zero abstention, denominator retained',
        interpretation='edge-balanced categorical atom-mark log ratio, not occupancy, joint likelihood or free energy',
        chemical_equivalence=equivalences,reference_counts=reference,cells=counts.to_dicts())
    validate_resource(resource)
    return resource


def validate_resource(resource):
    """Reconstruct every finite log ratio and source mass/count conservation."""
    if resource.get('version') != VERSION or resource.get('cutoff_A') != 5. or resource.get('mhc_class') != 'I':
        raise ValueError('unsupported atom-pattern resource')
    try:
        c = pl.DataFrame(resource['cells'])
        if c.is_empty() or c.select(TYPES).is_duplicated().any() or c.null_count().sum_horizontal().item():
            raise ValueError('nonempty unique nonnull pattern cells required')
        for s in ['from','to']:
            if c.filter(~pl.col('aa_'+s).is_in(list('ACDEFGHIKLMNPQRSTVWY')) |
                        ~pl.col('atom_'+s).str.contains(r'^[A-Z][A-Z0-9]*$')).height:
                raise ValueError('invalid atom-pattern labels')
        for source in ['native','external']:
            ns, ms = c[source+'_count'],c[source+'_mass']
            if ns.dtype not in (pl.Int64,pl.Int32) or (ns < 0).any() or not ms.is_finite().all() or (ms < 0).any() or (ms > ns+1e-12).any() or ((ms==0)!=(ns==0)).any():
                raise ValueError('invalid physical counts or fractional masses')
            ref = resource['reference_counts']
            if type(ref[source+'_atom_pairs']) is not int or type(ref[source+'_edges']) is not int or ref[source+'_edges'] < 0:
                raise ValueError('integer reference counts required')
            if ns.sum() != ref[source+'_atom_pairs'] or not math.isclose(ms.sum(),ref[source+'_edges'],rel_tol=0,abs_tol=1e-8):
                raise ValueError('source count/mass conservation failed')
        if not ((c['count']==c['native_count']+c['external_count']) & (c['count']>0)).all() or not np.allclose(c['mass'],c['native_mass']+c['external_mass'],rtol=0,atol=1e-12):
            raise ValueError('pooled pattern count/mass mismatch')
        expected = c.with_columns(pl.col('mass').sum().over('aa_from','aa_to').alias('M'),
            (pl.col('mass').sum().over('aa_from','atom_from')/pl.col('mass').sum().over('aa_from')).alias('B'),
            (pl.col('mass').sum().over('aa_to','atom_to')/pl.col('mass').sum().over('aa_to')).alias('T'))
        for saved, want in [('residue_pair_mass','M'),('q_binder','B'),('q_target','T')]:
            if not np.allclose(expected[saved],expected[want],rtol=0,atol=1e-12):
                raise ValueError('conditional mass/probability mismatch')
        values = -(expected['mass']/expected['M']/(expected['B']*expected['T'])).log()
        if not c['value'].is_finite().all() or not np.allclose(c['value'],values,rtol=0,atol=1e-12):
            raise ValueError('pattern log-ratio mismatch')
    except (KeyError,TypeError,pl.exceptions.PolarsError) as exc:
        raise ValueError('malformed atom-pattern resource') from exc


def _freeze(native,external,out,equivalence_resource=None):
    from datetime import datetime,timezone
    from ..exports import _sha256
    native,external,out = Path(native),Path(external),Path(out)
    if out.exists():
        raise FileExistsError(out)
    equivalence_resource = Path(equivalence_resource) if equivalence_resource else None
    source = dict(native_sha256=_sha256(native),external_sha256=_sha256(external))
    n = pl.scan_parquet(native).filter((pl.col('mhc_class')=='I')&(pl.col('interface')=='tcr_peptide'))
    swap = pl.col('role.from') != 'tcr'
    n = n.select('pdb.id',*[
        pl.when(swap).then(pl.col(key+'.'+other)).otherwise(pl.col(key+'.'+side)).alias(key+'.'+side)
        for side,other in [('from','to'),('to','from')] for key in ['chain.id','residue.index','chem','region.type']]).filter(
            pl.col('region.type.from').is_in(['CDR1','CDR2','CDR3'])).collect()
    membership = dict(native_pdb_entries=n['pdb.id'].n_unique())
    n = n.select(pl.struct('pdb.id','chain.id.from','residue.index.from','chain.id.to','residue.index.to').alias('edge'),*[
        pl.col('chem.'+s).str.split_exact(':',1).struct.field(f).alias(dest+'_'+s)
        for s in ['from','to'] for f,dest in [('field_0','aa'),('field_1','atom')]])
    e = pl.read_parquet(external).filter(pl.col('target_type')=='PEPTIDE')
    membership.update(external_pdb_entries=e['pdb'].n_unique(),external_interfaces=e['instance'].n_unique())
    e = e.select(pl.struct('pdb','instance','model','ab_chain','ab_auth','ab_ins','target_chain','target_auth','target_ins').alias('edge'),
        pl.col('ab_aa').alias('aa_from'),pl.col('ab_atom').alias('atom_from'),pl.col('target_aa').alias('aa_to'),pl.col('target_atom').alias('atom_to'))
    eq = json.loads(Path(equivalence_resource).read_text())['chemical_equivalence'] if equivalence_resource else {}
    r = derive_atom_patterns(n,e,eq)
    assert source == dict(native_sha256=_sha256(native),external_sha256=_sha256(external))
    r.update(sources=dict(**source,producer_sha256=_sha256(Path(__file__)),
        equivalence_resource_sha256=_sha256(equivalence_resource) if equivalence_resource else None),
        reference_membership=membership,frozen_utc=datetime.now(timezone.utc).isoformat(),polars_threads=pl.thread_pool_size())
    out.parent.mkdir(parents=True,exist_ok=True)
    with out.open('x') as stream:
        stream.write(json.dumps(r,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(resource_sha256=_sha256(out),cells=len(r['cells']),**r['reference_counts'])),flush=True)


def fit_command(native: Path=typer.Option(...,'--native',exists=True,dir_okay=False),
                external: Path=typer.Option(...,'--external',exists=True,dir_okay=False),
                out: Path=typer.Option(...,'-o','--out'),
                equivalence_resource: Path|None=typer.Option(None,'--equivalence-resource',exists=True,dir_okay=False),
                threads: int=typer.Option(0,'-t','--threads')):
    """Freeze one edge-balanced atom-pattern reference with no assay input."""
    from ..runtime import resolve_threads
    budget = resolve_threads(threads)
    if out.exists():
        raise FileExistsError(out)
    code = ('import polars as pl; pl.thread_pool_size(); import json,sys; '
            'from tcren.contacts.atom_patterns import _freeze; _freeze(**json.loads(sys.argv[1]))')
    env = os.environ.copy()
    for key in ['OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS','RAYON_NUM_THREADS']:
        env[key] = '1'
    env.update(POLARS_MAX_THREADS=str(budget),POLARS_ASYNC_THREAD_COUNT='1',POLARS_MAX_BLOCKING_THREAD_COUNT='4')
    args = dict(native=str(native.resolve()),external=str(external.resolve()),out=str(out.resolve()),
                equivalence_resource=str(equivalence_resource.resolve()) if equivalence_resource else None)
    subprocess.run([sys.executable,'-c',code,json.dumps(args)],env=env,check=True)
