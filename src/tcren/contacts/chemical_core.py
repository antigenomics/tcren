"""Small symmetric contact-composition core with directed residue context.

The explicit map and empirical reference are frozen inputs. Neither the core,
conditional partition nor physical-edge aggregation is a binding free energy.
"""
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

from .. import contact_types as chemistry
from .atom_counts import EDGE, TYPES, VERSION as ATOM_VERSION, score_atom_counts

VERSION = 'chemical-core-count-1'
_AA = set('ACDEFGHIKLMNPQRSTVWY')


def _chemical_type(aa, atom):
    element = chemistry._elem(atom)
    if element not in {'C', 'N', 'O', 'S'}:
        raise ValueError('standard heavy-atom chemistry required')
    if element == 'S':
        if (aa, atom) not in {('M', 'SD'), ('C', 'SG')}:
            raise ValueError('unknown sulfur atom')
        return 'S|' + ('methionine_thioether' if aa == 'M' else 'cysteine_state_unknown')
    if aa == 'H' and atom in {'ND1', 'NE2'}:
        return 'N|histidine_ring_protonation_unknown'
    donor, acceptor = chemistry._donor_acceptor(aa, atom)
    cation, anion = chemistry._charged(aa, atom)
    flags = [atom in {'N', 'CA', 'C', 'O', 'OXT'}, chemistry._is_apolar(aa, atom),
             donor, acceptor, cation, anion,
             aa in chemistry._AROMATIC_RES and atom in chemistry._RING_ATOMS]
    state = element + '|' + '|'.join(str(int(x)) for x in flags)
    if element == 'C':
        state += '|carbonyl=' + str(int(atom == 'C' or (aa, atom) in
                                      {('D', 'CG'), ('E', 'CD'), ('N', 'CG'), ('Q', 'CD')}))
    return state


def _counts(data):
    if not set([*TYPES, 'count']) <= set(data.columns):
        raise ValueError('directed atom types and count required')
    data = data.select(*TYPES, 'count')
    if data.null_count().sum_horizontal().item():
        raise ValueError('nonnull directed atom counts required')
    for side in ['from', 'to']:
        if data.filter(~pl.col('aa_'+side).is_in(sorted(_AA)) |
                       ~pl.col('atom_'+side).str.contains(r'^[A-Z][A-Z0-9]*$')).height:
            raise ValueError('standard amino acids and atom names required')
    if data['count'].dtype == pl.Boolean:
        raise ValueError('positive integer counts required')
    values = data['count'].to_numpy()
    if (not np.isfinite(values).all() or (values <= 0).any() or
            (values != np.floor(values)).any() or sum(data['count'].to_list()) > 2**53):
        raise ValueError('positive integer counts within exact float64 count range required')
    return data.with_columns(pl.col('count').cast(pl.Int64))


def derive_chemical_core(native, external):
    """Freeze unordered enrichment and role-specific conditional compositions."""
    native, external = _counts(native), _counts(external)
    total = int(native['count'].sum()+external['count'].sum())
    if not total or total > 2**53:
        raise ValueError('nonempty reference within exact float64 count range required')
    raw = sorted({f'{a}:{i}' for data in [native, external] for side in ['from', 'to']
                  for a, i in data.select('aa_'+side, 'atom_'+side).unique().iter_rows()})
    signatures = {key: _chemical_type(*key.split(':')) for key in raw}
    names = sorted(set(signatures.values()))
    if len(names) > 19:
        raise ValueError('reference exceeds the versioned 19-type partition')
    index = {name: i for i, name in enumerate(names)}
    mapping = [dict(chem=key, type=index[signatures[key]]) for key in raw]
    types = pl.DataFrame(mapping)
    frames = []
    for label, data in [('native', native), ('external', external)]:
        for side in ['from', 'to']:
            data = data.with_columns((pl.col('aa_'+side)+':'+pl.col('atom_'+side)).alias('chem'))
            data = data.join(types, on='chem', how='left', validate='m:1').rename({'type': 'type_'+side}).drop('chem')
        frames.append(data.with_columns(pl.lit(label).alias('source')))
    data = pl.concat(frames)
    endpoints = {}
    for side, name in [('from', 'binder'), ('to', 'target')]:
        e = data.group_by(pl.col('aa_'+side).alias('aa'), pl.col('type_'+side).alias('type')).agg(pl.col('count').sum())
        e = e.with_columns((pl.col('count')/pl.col('count').sum().over('aa')).alias('probability')).sort('aa', 'type')
        endpoints[name] = e.to_dicts()
    k = len(names)
    global_b = np.bincount(data['type_from'].to_numpy(), weights=data['count'].to_numpy(), minlength=k)/total
    global_t = np.bincount(data['type_to'].to_numpy(), weights=data['count'].to_numpy(), minlength=k)/total
    unordered = data.with_columns(pl.min_horizontal('type_from', 'type_to').alias('type_a'),
                                 pl.max_horizontal('type_from', 'type_to').alias('type_b'))
    counts = {(a,b,s): n for a,b,s,n in unordered.group_by('type_a', 'type_b', 'source').agg(pl.col('count').sum()).iter_rows()}
    cells = []
    w = np.ones((k, k))
    for a in range(k):
        for b in range(a, k):
            n, e = int(counts.get((a,b,'native'), 0)), int(counts.get((a,b,'external'), 0))
            q = float(global_b[a]*global_t[b] + (global_b[b]*global_t[a] if a != b else 0))
            value = -math.log((n+e)/(total*q)) if n+e else 0.
            w[a,b] = w[b,a] = math.exp(-value)
            cells.append(dict(type_a=a,type_b=b,native_count=n,external_count=e,count=n+e,
                              reference_probability=q,value=value))
    qb, aa_b = _endpoint_array(endpoints['binder'], k)
    qt, aa_t = _endpoint_array(endpoints['target'], k)
    z = np.log(qb @ w @ qt.T)
    normalizers = [dict(aa_from=a,aa_to=b,logZ=float(z[i,j]))
                   for i,a in enumerate(aa_b) for j,b in enumerate(aa_t)]
    resource = dict(version=VERSION,cutoff_A=5.,mhc_class='I',mapping=mapping,types=names,
        cells=cells,**endpoints,normalizers=normalizers,
        reference_counts=dict(native_atom_pairs=int(native['count'].sum()),external_atom_pairs=int(external['count'].sum())),
        estimator='unordered actual-count enrichment; global directed endpoint reference; unseen W=1',
        conditional_reference='directed contact-endpoint composition conditional on residue identity; not exposure',
        missing='unseen supported core: logZ; unmapped/zero conditional support: zero abstention, denominator retained',
        aggregation='one mean per known original CDR residue edge, unit coefficient',
        interpretation='dimensionless compositional correction; not a joint structure likelihood or free energy')
    validate_resource(resource)
    return resource


def _endpoint_array(rows, k):
    aa = sorted({r['aa'] for r in rows})
    array = np.zeros((len(aa), k))
    index = {a:i for i,a in enumerate(aa)}
    seen = set()
    for r in rows:
        key = (r['aa'], r['type'])
        if (key in seen or r['aa'] not in _AA or not isinstance(r['type'], int) or
                isinstance(r['type'],bool) or not 0 <= r['type'] < k or not math.isfinite(r['probability']) or
                r['probability'] <= 0 or not isinstance(r['count'], int) or isinstance(r['count'],bool) or r['count'] <= 0):
            raise ValueError('unique finite positive endpoint probabilities/counts required')
        seen.add(key)
        array[index[r['aa']],r['type']] = r['probability']
    if not aa or not np.allclose(array.sum(axis=1), 1., rtol=0, atol=1e-12):
        raise ValueError('normalized endpoint probabilities required')
    for a in aa:
        part = [r for r in rows if r['aa']==a]
        total = sum(r['count'] for r in part)
        if any(not math.isclose(r['probability'],r['count']/total,rel_tol=0,abs_tol=1e-12) for r in part):
            raise ValueError('endpoint count/probability mismatch')
    return array, aa


def validate_resource(resource):
    """Reject malformed finite models rather than classify them as missing support."""
    if resource.get('version') != VERSION or resource.get('cutoff_A') != 5. or resource.get('mhc_class') != 'I':
        raise ValueError('unsupported chemical-core resource')
    try:
        names, mapping = resource['types'], resource['mapping']
        k = len(names)
        if not 0 < k <= 19 or any(not isinstance(n,str) or not n for n in names) or len(set(names)) != k:
            raise ValueError('unique chemical types required')
        if not mapping or len({x['chem'] for x in mapping}) != len(mapping):
            raise ValueError('unique nonempty atom mapping required')
        for row in mapping:
            aa,atom = row['chem'].split(':')
            if (aa not in _AA or not atom or not isinstance(row['type'], int) or
                    isinstance(row['type'],bool) or not 0 <= row['type'] < k or
                    names[row['type']] != _chemical_type(aa,atom)):
                raise ValueError('invalid atom mapping')
        allowed = {(r['chem'].split(':')[0],r['type']) for r in mapping}
        if any((r['aa'],r['type']) not in allowed for side in ['binder','target'] for r in resource[side]):
            raise ValueError('endpoint type absent from atom mapping')
        qb, aa_b = _endpoint_array(resource['binder'], k)
        qt, aa_t = _endpoint_array(resource['target'], k)
        ref = resource['reference_counts']
        totals = [ref['native_atom_pairs'], ref['external_atom_pairs']]
        if any(not isinstance(n,int) or isinstance(n,bool) or n < 0 for n in totals) or not 0 < sum(totals) <= 2**53:
            raise ValueError('invalid reference counts')
        total = sum(totals)
        for side in ['binder','target']:
            if sum(r['count'] for r in resource[side]) != total:
                raise ValueError('endpoint/reference count mismatch')
        b = np.zeros(k);t = np.zeros(k)
        for r in resource['binder']:
            b[r['type']] += r['count']/total
        for r in resource['target']:
            t[r['type']] += r['count']/total
        cells = resource['cells']; expected = {(a,b) for a in range(k) for b in range(a,k)}
        if len(cells) != len(expected) or {(x['type_a'],x['type_b']) for x in cells} != expected:
            raise ValueError('complete unique unordered core required')
        w = np.ones((k,k));n_native = n_external = 0
        incidence = [0]*k
        for c in cells:
            a,d = c['type_a'],c['type_b']
            if type(a) is not int or type(d) is not int:
                raise ValueError('integer core indices required')
            ns = [c['native_count'],c['external_count'],c['count']]
            if any(not isinstance(n,int) or isinstance(n,bool) or n < 0 for n in ns) or ns[0]+ns[1] != ns[2]:
                raise ValueError('invalid core counts')
            value,q = c['value'],c['reference_probability']
            if not math.isfinite(value) or not math.isfinite(q) or q < 0:
                raise ValueError('finite core values and reference probabilities required')
            want_q = b[a]*t[d]+(b[d]*t[a] if a != d else 0.)
            if not math.isclose(q,want_q,rel_tol=0,abs_tol=1e-12) or (c['count'] and q <= 0):
                raise ValueError('core reference probability mismatch')
            want = -math.log(c['count']/(total*q)) if c['count'] else 0.
            if not math.isclose(value,want,rel_tol=0,abs_tol=1e-12):
                raise ValueError('core enrichment mismatch')
            w[a,d] = w[d,a] = math.exp(-value)
            n_native += c['native_count'];n_external += c['external_count']
            incidence[a] += c['count'];incidence[d] += c['count']
        if [n_native,n_external] != totals:
            raise ValueError('core/reference count mismatch')
        endpoints = [sum(r['count'] for side in ['binder','target'] for r in resource[side] if r['type']==i) for i in range(k)]
        if incidence != endpoints:
            raise ValueError('endpoint/core incidence mismatch')
        z = np.log(qb @ w @ qt.T)
        normal = resource['normalizers'];expected = {(a,b) for a in aa_b for b in aa_t}
        if len(normal) != len(expected) or {(x['aa_from'],x['aa_to']) for x in normal} != expected:
            raise ValueError('complete unique conditional normalizers required')
        lookup = {(a,b):float(z[i,j]) for i,a in enumerate(aa_b) for j,b in enumerate(aa_t)}
        if any(not math.isfinite(x['logZ']) or not math.isclose(x['logZ'],lookup[x['aa_from'],x['aa_to']],rel_tol=0,abs_tol=1e-12) for x in normal):
            raise ValueError('conditional normalizer mismatch')
    except (KeyError, TypeError, OverflowError, IndexError) as exc:
        raise ValueError('malformed chemical-core resource') from exc


def score_chemical_core(contacts, observations, resource, original, *, mhc_class='I'):
    """Project only query chemical pairs, reusing the existing physical-edge scorer."""
    validate_resource(resource)
    if mhc_class != 'I':
        raise ValueError('chemical-core reference supports class I only')
    swap = pl.col('role.from') != 'tcr'
    atoms = observations.select(*[pl.when(swap).then(pl.col(key+'.'+other)).otherwise(pl.col(key+'.'+side)).alias(key+'.'+side)
        for side,other in [('from','to'),('to','from')] for key in ['chain.id','residue.index','chem']])
    atoms = atoms.with_columns(*[pl.col('chem.'+s).str.split_exact(':',1).struct.field(f).alias(dest+'_'+s)
        for s in ['from','to'] for f,dest in [('field_0','aa'),('field_1','atom')]])
    query = atoms.select(*TYPES,'chem.from','chem.to').unique()
    for side,endpoint in [('from','binder'),('to','target')]:
        query = query.join(pl.DataFrame(resource['mapping']).rename({'chem':'chem.'+side,'type':'type_'+side}),on='chem.'+side,how='left',validate='m:1')
        q = pl.DataFrame(resource[endpoint]).select(pl.col('aa').alias('aa_'+side),pl.col('type').alias('type_'+side))
        query = query.join(q.with_columns(pl.lit(True).alias('supported_'+side)),on=['aa_'+side,'type_'+side],how='left',validate='m:1')
    query = query.with_columns(pl.min_horizontal('type_from','type_to').alias('type_a'),pl.max_horizontal('type_from','type_to').alias('type_b'))
    query = query.join(pl.DataFrame(resource['cells']).select('type_a','type_b','count','value'),on=['type_a','type_b'],how='left',validate='m:1')
    query = query.join(pl.DataFrame(resource['normalizers']),on=['aa_from','aa_to'],how='left',validate='m:1')
    known = pl.col('supported_from').fill_null(False)&pl.col('supported_to').fill_null(False)&pl.col('logZ').is_finite().fill_null(False)
    query = query.filter(known).with_columns((pl.col('value')+pl.col('logZ')).alias('value'),(pl.col('count')==0).alias('unseen'))
    # Execution-owned projection, not another estimated table or persistent cache.
    projection = dict(version=ATOM_VERSION,cutoff_A=5.,mhc_class='I',chemical_equivalence={},cells=query.select(*TYPES,'value').to_dicts())
    base = contacts.select(*EDGE,'region.type.from','original_contribution')
    detail,summary = score_atom_counts(base,observations,projection,original,mhc_class=mhc_class)
    extra = atoms.join(query.select(*TYPES,'unseen'),on=TYPES,how='left',validate='m:1').group_by(EDGE).agg(pl.col('unseen').fill_null(False).sum().alias('chemical_core_unseen_pairs'))
    detail = detail.rename({c:c.replace('atom_count','chemical_core',1) for c in detail.columns if c.startswith('atom_count')})
    detail = contacts.join(detail.select(*EDGE,*[c for c in detail.columns if c.startswith('chemical_core')]),on=EDGE,how='left',validate='1:1')
    detail = detail.join(extra,on=EDGE,how='left',validate='1:1').with_columns(pl.col('chemical_core_unseen_pairs').fill_null(0))
    summary = {c.replace('atom_count','chemical_core',1):v for c,v in summary.items()}
    summary['chemical_core_unseen_atoms'] = int(detail['chemical_core_unseen_pairs'].sum())
    return detail,summary


def _freeze(native, external, out):
    from datetime import datetime, timezone
    from ..exports import _sha256
    native, external, out = Path(native),Path(external),Path(out)
    if out.exists():
        raise FileExistsError(out)
    source_hashes = dict(native_sha256=_sha256(native),external_sha256=_sha256(external))
    n = pl.scan_parquet(native).filter((pl.col('mhc_class')=='I')&(pl.col('interface')=='tcr_peptide'))
    swap = pl.col('role.from') != 'tcr'
    n = n.select(pl.col('pdb.id').alias('pdb'),pl.col('pdb.id').alias('instance'),
        *[pl.when(swap).then(pl.col('chem.'+other)).otherwise(pl.col('chem.'+side)).alias('chem_'+side)
          for side,other in [('from','to'),('to','from')]]).collect()
    membership = dict(native=dict(pdb_entries=n['pdb'].n_unique(),interfaces=n['instance'].n_unique(),physical_atom_pairs=n.height))
    n = n.with_columns(*[pl.col('chem_'+s).str.split_exact(':',1).struct.field(f).alias(dest+'_'+s)
        for s in ['from','to'] for f,dest in [('field_0','aa'),('field_1','atom')]]).group_by(TYPES).len(name='count')
    e = pl.scan_parquet(external).filter(pl.col('target_type')=='PEPTIDE').select('pdb','instance',
        pl.col('ab_aa').alias('aa_from'),pl.col('ab_atom').alias('atom_from'),pl.col('target_aa').alias('aa_to'),pl.col('target_atom').alias('atom_to')).collect()
    membership['external'] = dict(pdb_entries=e['pdb'].n_unique(),interfaces=e['instance'].n_unique(),physical_atom_pairs=e.height)
    e = e.group_by(TYPES).len(name='count')
    if source_hashes != dict(native_sha256=_sha256(native),external_sha256=_sha256(external)):
        raise ValueError('reference source changed during loading')
    resource = derive_chemical_core(n,e)
    resource['reference_membership'] = membership
    resource['sources'] = dict(**source_hashes,producer_sha256=_sha256(Path(__file__)),
        chemical_source_sha256=_sha256(Path(chemistry.__file__)))
    resource['frozen_utc'] = datetime.now(timezone.utc).isoformat()
    resource['polars_threads'] = pl.thread_pool_size()
    out.parent.mkdir(parents=True,exist_ok=True)
    with out.open('x') as stream:
        stream.write(json.dumps(resource,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(resource_sha256=_sha256(out),types=len(resource['types']),cells=len(resource['cells']),**resource['reference_counts'])),flush=True)


def fit_command(native: Path=typer.Option(...,'--native',exists=True,dir_okay=False),
                external: Path=typer.Option(...,'--external',exists=True,dir_okay=False),
                out: Path=typer.Option(...,'-o','--out'),
                threads: int=typer.Option(0,'-t','--threads')):
    """Freeze the chemical core using one allocation-sized Polars compute pool."""
    from ..runtime import resolve_threads
    budget = resolve_threads(threads)
    if out.exists():
        raise FileExistsError(out)
    # Polars' pool cannot resize. Initialize it before tcren in this dedicated
    # table-only child; BLAS/OpenMP stay serial and no outer workers are started.
    code = ('import polars as pl; pl.thread_pool_size(); '
            'import json,sys; from tcren.contacts.chemical_core import _freeze; '
            '_freeze(**json.loads(sys.argv[1]))')
    env = os.environ.copy()
    for key in ['OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS',
                'NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS','RAYON_NUM_THREADS']:
        env[key] = '1'
    env['POLARS_MAX_THREADS'] = str(budget)
    env['POLARS_ASYNC_THREAD_COUNT'] = '1';env['POLARS_MAX_BLOCKING_THREAD_COUNT'] = '4'
    subprocess.run([sys.executable,'-c',code,json.dumps(dict(native=str(native.resolve()),external=str(external.resolve()),out=str(out.resolve())))],env=env,check=True)
