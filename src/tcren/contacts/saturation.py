"""Native-frequency CDR3 atom saturation; experimental, not binding probabilities."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import polars as pl
import typer
from scipy.stats import binom

from ..refine.sidechain import TOPOLOGY
from ..runtime import parallel_map, resolve_threads, single_threaded

EXPECTED = {aa: sorted(({'CB'} | {a[0] for a in atoms}) if aa != 'G' else set())
            for aa, atoms in TOPOLOGY.items()}
KEY = ['chain.id', 'residue.index']
INTERFACES = ('peptide', 'mhc', 'union')


def bounded_mixture(iptm, prepared, surprise, w1: float, w2: float):
    """Nonnegative convex mixture; no fitting or implicit missing-data imputation."""
    values=np.asarray([iptm,prepared,surprise],dtype=float)
    weights=np.asarray([1.,w1,w2])
    if not np.isfinite(values).all() or np.any((values<0)|(values>1)):
        raise ValueError('finite component scores in [0,1] required')
    if not np.isfinite(weights).all() or np.any(weights<0):
        raise ValueError('finite nonnegative weights required')
    return np.average(values,axis=0,weights=weights)


def saturation_counts(inventory: pl.DataFrame, contacts: pl.DataFrame) -> pl.DataFrame:
    """Complete canonical side chains only; duplicate atom witnesses count once."""
    inv = inventory.group_by(KEY+['residue.aa']).agg(pl.col('atom.name').unique().alias('observed'))
    expected = pl.DataFrame({'residue.aa': list(EXPECTED), 'expected': list(EXPECTED.values())})
    inv = inv.join(expected, on='residue.aa', how='left', validate='m:1').with_columns(
        pl.col('expected').list.len().fill_null(0).alias('n'),
        (pl.col('expected').list.set_difference(pl.col('observed')).list.len()==0)
        .fill_null(False).alias('complete'))
    parts = []
    for interface in INTERFACES:
        hits = contacts if interface=='union' else contacts.filter(pl.col('interface')==interface)
        hits = hits.select(pl.col('chain.id.from').alias('chain.id'),
                           pl.col('residue.index.from').alias('residue.index'),
                           pl.col('atom.from').alias('atom.name')).group_by(KEY).agg(
                               pl.col('atom.name').unique().alias('hit'))
        part = inv.join(hits,on=KEY,how='left',validate='1:1').with_columns(
            pl.col('expected').list.set_intersection(pl.col('hit')).list.len().fill_null(0).alias('k'),
            pl.lit(interface).alias('interface'))
        parts.append(part.drop('expected','observed','hit'))
    return pl.concat(parts).sort('interface',*KEY)


def fit_reference(rows: pl.DataFrame) -> dict:
    """Jeffreys full-contact frequencies and pooled atom-contact probabilities."""
    usable = rows.filter(pl.col('complete') & (pl.col('n')>0))
    if usable.is_empty():
        raise ValueError('no complete native side chains')
    if 'weight' not in usable.columns:
        usable=usable.with_columns(pl.lit(1.).alias('weight'))
    weights=usable['weight'].to_numpy()
    if not np.isfinite(weights).all() or np.any(weights<=0):
        raise ValueError('positive finite reference weights required')
    cells = usable.filter(pl.col('k')>0).group_by('interface','residue.aa').agg(
        pl.len().alias('N'), (pl.col('k')==pl.col('n')).sum().alias('F'),
        pl.col('weight').sum().alias('weighted_N'),
        (pl.col('weight')*(pl.col('k')==pl.col('n'))).sum().alias('weighted_F'))
    pooled = usable.group_by('interface').agg(pl.col('k').sum(),pl.col('n').sum(),
        (pl.col('weight')*pl.col('k')).sum().alias('weighted_k'),
        (pl.col('weight')*pl.col('n')).sum().alias('weighted_n'))
    return dict(schema='tcren-saturation-1',cutoff_A=5.,
        cells={r['interface']+':'+r['residue.aa']: dict(N=r['N'],F=r['F'],weighted_N=r['weighted_N'],weighted_F=r['weighted_F'],
               p_full=(r['weighted_F']+.5)/(r['weighted_N']+1))
               for r in cells.to_dicts()},
        atom_probability={r['interface']:dict(contacted=r['k'],available=r['n'],weighted_contacted=r['weighted_k'],
                          weighted_available=r['weighted_n'],p=(r['weighted_k']+.5)/(r['weighted_n']+1))
                          for r in pooled.to_dicts()})


def score_saturation(rows: pl.DataFrame, reference: dict) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Information sums / maximum attainable sums, independently in each interface.

    Binomial K|K>0 uses a pooled native atom probability including noncontacting
    native residues in estimation. Full-only score is a weighted saturation fraction.
    Empty/unsupported information contributes zero with explicit support counts.
    """
    if reference.get('schema')!='tcren-saturation-1':
        raise ValueError('unknown saturation reference')
    n, k = rows['n'].to_numpy().astype(int), rows['k'].to_numpy().astype(int)
    if np.any(n<0) or np.any(k<0) or np.any(k>n):
        raise ValueError('atom counts require 0 <= k <= n')
    full = np.zeros(len(n)); maximum = np.zeros(len(n))
    tail = np.zeros(len(n)); tmax = np.zeros(len(n))
    bits = np.zeros(len(n)); bmax = np.zeros(len(n)); supported = np.zeros(len(n),dtype=bool)
    eligible = rows['complete'].to_numpy() & (n>0) & (k>0)
    for interface in INTERFACES:
        mask = eligible & (rows['interface'].to_numpy()==interface)
        opportunity = rows['complete'].to_numpy() & (n>1) & (rows['interface'].to_numpy()==interface)
        if not (mask.any() or opportunity.any()):
            continue
        if interface not in reference['atom_probability']:
            raise ValueError(f'missing native atom reference for {interface}')
        p = reference['atom_probability'][interface]['p']
        if not 0<p<1:
            raise ValueError('reference probabilities must lie strictly between zero and one')
        no = n[opportunity]
        tmax[opportunity] = np.maximum(0,(-no*np.log(p)+np.log(-np.expm1(no*np.log1p(-p))))/np.log(2))
        nn, kk = n[mask], k[mask]
        lognorm = np.log(-np.expm1(nn*np.log1p(-p)))
        tail[mask] = np.maximum(0,(-binom.logsf(kk-1,nn,p)+lognorm)/np.log(2))
        bits[mask] = np.maximum(0,(-binom.logpmf(kk,nn,p)+lognorm)/np.log(2))
        bmax[mask] = np.maximum(0,(-np.minimum(binom.logpmf(1,nn,p),binom.logpmf(nn,nn,p))+lognorm)/np.log(2))
        # n=1 given K>0 has exactly one outcome, independent of floating arithmetic.
        tail[mask & (n==1)] = bits[mask & (n==1)] = bmax[mask & (n==1)] = 0
        for aa in EXPECTED:
            cell = reference['cells'].get(interface+':'+aa)
            amask = mask & (rows['residue.aa'].to_numpy()==aa)
            if cell is None:
                continue
            if not 0<cell['p_full']<1:
                raise ValueError('invalid full-contact reference probability')
            supported[amask] = True
            maximum[amask] = -np.log2(cell['p_full'])
            full[amask & (k==n)] = -np.log2(cell['p_full'])
    scored = rows.with_columns(pl.Series('eligible_contact',eligible),pl.Series('full_supported',supported),
        pl.Series('full_bits',full),pl.Series('full_max_bits',maximum),
        pl.Series('tail_bits',tail),pl.Series('tail_max_bits',tmax),pl.Series('binomial_bits',bits),pl.Series('binomial_max_bits',bmax))
    summary = scored.group_by('pdb.id','interface').agg(
        pl.col('full_bits','full_max_bits','binomial_bits','binomial_max_bits','tail_bits','tail_max_bits').sum(),
        pl.col('eligible_contact').sum().alias('contacting_residues'),
        (pl.col('eligible_contact') & (pl.col('k')==pl.col('n'))).sum().alias('fully_contacting'),
        (pl.col('eligible_contact') & ~pl.col('full_supported')).sum().alias('unsupported'),
        ((pl.col('n')>0) & ~pl.col('complete')).sum().alias('incomplete'))
    for kind in ('full','binomial','tail'):
        summary = summary.with_columns((pl.col(kind+'_bits')/(1+pl.col(kind+'_bits'))).alias(kind+'_sum'))
        summary = summary.with_columns(pl.when(pl.col(kind+'_max_bits')>0)
            .then(pl.col(kind+'_bits')/pl.col(kind+'_max_bits')).otherwise(0.).alias(kind))
        if not summary[kind].is_between(-1e-12,1+1e-12).all():
            raise RuntimeError('surprisal normalization violated [0,1]')
        summary = summary.with_columns(pl.col(kind).clip(0,1))
    return scored, summary.sort('pdb.id','interface')


def residue_saturation(structure) -> pl.DataFrame:
    """Fresh inventory and atom contacts; all CDR3 residues remain in the output."""
    from dataclasses import replace
    from .observations import atom_inventory, _residue_atoms
    from .geometry import BACKBONE_ATOMS
    from ..contactmap import ContactMap
    from ..exports import residue_identity
    inv = atom_inventory(structure).filter(pl.col('chain.type').is_in(['TRA','TRB']) & (pl.col('region.type')=='CDR3'))
    # Reuse the inventory's conformer policy on both endpoints, retaining whole MHC.
    selected={c.chain_id:[] for c in structure.chains}
    for chain,residue,_,_,_,_,_,kept in _residue_atoms(structure):
        selected[chain.chain_id].append(replace(residue,atoms=tuple(a for _,a in kept)))
    contact_structure=replace(structure,chains=[replace(c,residues=selected[c.chain_id]) for c in structure.chains])
    cm = ContactMap.from_structure(contact_structure,cutoff=5.,atom_pairs=True)
    hits = pl.concat([cm.interface(name).filter(pl.col('region.type.from')=='CDR3')
        .with_columns(pl.lit(label).alias('interface'))
        for name,label in [('tcr_peptide','peptide'),('tcr_mhc','mhc')]])
    result = saturation_counts(inv,hits)
    identity = residue_identity(structure).filter(pl.col('chain.type').is_in(['TRA','TRB']) & (pl.col('region.type')=='CDR3'))
    identity = identity.with_columns((pl.col('residue.index')-pl.col('residue.index').min().over('chain.id')+2).alias('junction_position'),
        (pl.len().over('chain.id')+2).alias('junction_length'))
    result = result.join(identity.drop('residue.aa'),on=KEY,validate='m:1')
    seqs = {c.chain_type: ''.join(r.aa for reg in c.regions if reg.region_type=='CDR3' for r in reg.residues)
            for c in structure.chains if c.chain_type in ('TRA','TRB')}
    # Anchor-excluded CDR3 pair is used only to audit exact reference/query overlap.
    result = result.with_columns(pl.lit(seqs.get('TRA','')+'|'+seqs.get('TRB','')).alias('cdr3_pair'))
    target_identity=residue_identity(structure).select(pl.col('chain.id').alias('chain.id.to'),
        pl.col('residue.index').alias('residue.index.to'),pl.col('residue.index.pdb').alias('target_position'))
    targets = hits.filter(~pl.col('atom.from').is_in(BACKBONE_ATOMS)).join(target_identity,on=['chain.id.to','residue.index.to'],validate='m:1').select(
        pl.col('chain.id.from').alias('chain.id'),pl.col('residue.index.from').alias('residue.index'),
        'interface',pl.concat_str('chain.type.to','residue.aa.to',pl.col('target_position').cast(pl.String),separator=':').alias('target'))
    parts = [targets,targets.with_columns(pl.lit('union').alias('interface'))]
    target = pl.concat(parts).group_by(KEY+['interface']).agg(pl.col('target').unique().sort().str.join(';').alias('targets'))
    return result.join(target,on=KEY+['interface'],how='left',validate='1:1').with_columns(pl.col('targets').fill_null(''))


def _extract_chunk(structures):
    return pl.concat([residue_saturation(s) for s in structures])


@single_threaded
def export_command(
    structures: list[str] = typer.Option(...,'-s','--structures'),
    output: Path = typer.Option(...,'-o','--output'),
    threads: int = typer.Option(0,'-t','--threads'),
    native: bool = typer.Option(False,'--native',help='Retain nonwithdrawn class-I alpha-beta native complexes; record exclusions'),
):
    """Export fresh complete CDR3 atom coverage, using one batched annotation."""
    from ..cli import _typed_batch
    from ..provenance import stamp, annotation_indexes
    threads=resolve_threads(threads)
    typed=_typed_batch(structures,'human',threads=threads)
    excluded=[]
    if native:
        from ..potential import withdrawn_ids
        withdrawn=set(withdrawn_ids())
        retained=[]
        for s in typed:
            types={c.chain_type for c in s.chains}
            classes={c.chain_supertype for c in s.chains if c.chain_type in {'MHCa','MHCb'}}
            peptides=s.by_type('PEPTIDE')
            reasons=[]
            if s.pdb_id.lower() in withdrawn:
                reasons.append('withdrawn')
            if not {'TRA','TRB','PEPTIDE','MHCa'}<=types or types & {'TRG','TRD','IGH','IGK','IGL'}:
                reasons.append('not alpha-beta TCR-pMHC')
            if classes!={'MHCI'}:
                reasons.append('not unambiguous class I')
            if len(peptides)!=1 or set(peptides[0].sequence())-set(EXPECTED):
                reasons.append('not one standard-amino-acid peptide')
            if reasons:
                excluded.append(dict(id=s.pdb_id,reasons=reasons))
            else:
                retained.append(s)
        typed=retained
    if not typed or len({s.pdb_id for s in typed})!=len(typed):
        raise ValueError('nonempty structures with unique IDs required')
    count=min(threads,len(typed))
    chunks=[typed[i*len(typed)//count:(i+1)*len(typed)//count] for i in range(count)]
    table=pl.concat(parallel_map(_extract_chunk,chunks,threads)).sort('pdb.id','interface',*KEY)
    table.write_parquet(output)
    stamp(output,command='saturation-export',extra={'models':len(typed),'threads':threads,'cutoff_A':5.,'region':'CDR3',
        'excluded_native_models':excluded,'annotation_indexes':annotation_indexes(typed)})
    typer.echo(f'Exported {len(typed)} models; {table.height} residue/interface rows')


def score_command(
    native: Path = typer.Option(...,'--native'),
    features: Path = typer.Option(...,'--features'),
    output: Path = typer.Option(...,'-o','--output'),
):
    """Freeze empirical native probabilities, excluding exact query receptor pairs."""
    reference_rows, rows = pl.read_parquet(native),pl.read_parquet(features)
    pairs=rows.select('cdr3_pair').unique()
    excluded=reference_rows.join(pairs,on='cdr3_pair',how='semi').select('pdb.id').unique()
    reference_rows=reference_rows.join(excluded,on='pdb.id',how='anti')
    reference=fit_reference(reference_rows)
    reference.update(excluded_native_models=excluded['pdb.id'].to_list(),native_models=reference_rows['pdb.id'].n_unique())
    scored,summary=score_saturation(rows,reference)
    output.mkdir(parents=True,exist_ok=False)
    (output/'reference.json').write_text(json.dumps(reference,indent=2)+'\n')
    scored.write_parquet(output/'residues.parquet')
    summary.write_csv(output/'scores.tsv',separator='\t')
    typer.echo(f'Native models {reference["native_models"]}; query models {rows["pdb.id"].n_unique()}')



def balance_mock_residues(rows: pl.DataFrame) -> pl.DataFrame:
    """Equal epitope mass per chain; equal distinct CDR3 mass within epitope.

    Each model/chain has one weight, irrespective of length or contact count.
    Duplicate CDR3s split their mass. Normalize model/chain mean weight to one
    so the fixed half-count smoothing has the same scale in both references.
    """
    keys=['pdb.id','chain.type']
    models=rows.select(*keys,'epitope','cdr3_pair').unique()
    if models.select(keys).n_unique()!=models.height:
        raise ValueError('one epitope and receptor per model/chain required')
    models=models.with_columns(pl.when(pl.col('chain.type')=='TRA')
        .then(pl.col('cdr3_pair').str.split('|').list.get(0))
        .otherwise(pl.col('cdr3_pair').str.split('|').list.get(1)).alias('cdr3'))
    models=models.with_columns((1./(pl.len().over('epitope','chain.type','cdr3') *
        pl.col('cdr3').n_unique().over('epitope','chain.type'))).alias('weight'))
    models=models.with_columns(pl.col('weight')/pl.col('weight').mean())
    return rows.join(models.select(*keys,'weight'),on=keys,validate='m:1',maintain_order='left')


def mock_command(
    reference_rows: Path = typer.Option(...,'--reference'),
    metadata: Path = typer.Option(...,'--metadata'),
    features: Path = typer.Option(...,'--features'),
    queries: Path = typer.Option(...,'--queries'),
    output: Path = typer.Option(...,'-o','--output'),
):
    """Fit raw and balanced mock references after epitope/exact-chain exclusions."""
    refrows=pl.read_parquet(reference_rows)
    rows=pl.read_parquet(features)
    meta=pl.read_csv(metadata,separator='\t')
    querymeta=pl.read_parquet(queries)
    if meta['pdb.id'].n_unique()!=meta.height or set(refrows['pdb.id'])!=set(meta['pdb.id']):
        raise ValueError('reference metadata must match all unique models exactly')
    if not (meta['y']==0).all():
        raise ValueError('mock reference requires exclusively labelled mispairings')
    refrows=refrows.join(meta.select('pdb.id','epitope'),on='pdb.id',validate='m:1')
    pairs=rows['cdr3_pair'].unique().to_list()
    a={p.split('|')[0] for p in pairs}; b={p.split('|')[1] for p in pairs}
    excluded=refrows.filter(pl.col('epitope').is_in(querymeta['epitope'].unique().to_list()) |
        pl.col('cdr3_pair').str.split('|').list.get(0).is_in(a) |
        pl.col('cdr3_pair').str.split('|').list.get(1).is_in(b)).select('pdb.id').unique()
    refrows=refrows.join(excluded,on='pdb.id',how='anti')
    output.mkdir(parents=True,exist_ok=False)
    refrows.write_parquet(output/'reference_residues.parquet')
    for name,table in [('raw',refrows),('balanced',balance_mock_residues(refrows))]:
        ref=fit_reference(table)
        ref.update(reference_kind='mock',balance=name,models=table['pdb.id'].n_unique(),
                   excluded=excluded['pdb.id'].sort().to_list())
        (output/(name+'.json')).write_text(json.dumps(ref,indent=2)+'\n')
        detail,summary=score_saturation(rows,ref)
        detail.write_parquet(output/(name+'_residues.parquet'))
        summary.write_csv(output/(name+'_scores.tsv'),separator='\t')
    typer.echo(f'Mock reference: {refrows["pdb.id"].n_unique()} retained; {excluded.height} excluded')
