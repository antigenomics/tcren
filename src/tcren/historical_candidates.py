"""Explicit historical candidate protocols; never invoked by ordinary scoring."""
import json
from pathlib import Path

import numpy as np
import polars as pl
import typer

from .associations import read_table as _read_table
from .potential.redundancy import alphabeta_ids

AA = 'LFIMVWYCHAGPTSQNDERK'


def read_table(path, *, schema_overrides=None):
    if schema_overrides is not None:
        if Path(path).suffix == '.parquet':
            return _read_table(path).cast(schema_overrides)
        return pl.read_csv(path, separator='\t' if str(path).endswith(('.tsv', '.tsv.gz')) else ',',
                           infer_schema_length=None, schema_overrides=schema_overrides)
    return pl.read_csv(path, separator='\t', infer_schema_length=None) if str(path).endswith('.tsv.gz') else _read_table(path)


def markup_table(residues):
    records = []
    for (pid,), part in residues.group_by('pdb.id', maintain_order=True):
        row = {'pdb.id': pid}
        for typ, region, column in [('TRA', 'CDR3', 'cdr3a'), ('TRB', 'CDR3', 'cdr3b'), ('PEPTIDE', None, 'peptide')]:
            selected = part.filter(pl.col('chain.type') == typ)
            if region:
                selected = selected.filter(pl.col('region.type') == region)
            row[column] = ''.join(selected.sort('residue.index')['residue.aa']) if selected.height and selected['chain.id'].n_unique() == 1 else None
        records.append(row)
    return pl.DataFrame(records)


def native_accounting(markup, contacts, ids):
    info = {r['pdb.id']: r for r in markup.to_dicts()}
    has_contacts, alpha_beta = set(contacts['pdb.id']), set(alphabeta_ids(contacts))
    ledger = []
    for pid in ids:
        row = info[pid]
        pep = row['peptide']
        reason = '' if pep and len(pep) >= 4 and set(pep) <= set(AA) and pid in has_contacts else 'missing/nonstandard peptide or no TCR contacts'
        if pid not in alpha_beta:
            reason = 'outside alpha-beta contact scope'
        ledger.append({'id': pid, 'cognate_eligible': not reason, 'reason': reason,
                       'cdr3a_observed': bool(row['cdr3a']), 'cdr3b_observed': bool(row['cdr3b'])})
    return pl.DataFrame(ledger)


def normalize_mhc(value):
    if not value:
        return None
    first = value.split(':')[0]
    return 'H-2' if first.startswith(('H2-', 'H-2')) else first[:6] if first.startswith('HLA-D') else first


def generate_candidates(protocol, spec):
    """Execute established ordering, RNG consumption and fixed candidate counts."""
    rows, accounting = [], []

    def append(pid, peptide, label):
        rows.append({'complex.id': pid, 'candidate.id': str(len(rows)), 'peptide': peptide, 'y': label})

    if protocol in {'anchored-random', 'matched-observed', 'receptor-shuffle'}:
        markup = read_table(Path(spec['markup']))
        info = {r['pdb.id']: r for r in markup.to_dicts()}
        ids = spec['ids']
    if protocol == 'anchored-random':
        rng = np.random.default_rng(0)
        rng.shuffle(np.arange(400))  # Historical preceding notebook RNG consumption.
        for pid in ids:
            cog = info[pid]['peptide']
            indices = np.array([AA.index(c) for c in cog])
            dec = rng.integers(0, 20, size=(1000, len(cog)))
            dec[:, 1], dec[:, -1] = indices[1], indices[-1]
            append(pid, cog, 1)
            for row in dec:
                append(pid, ''.join(AA[k] for k in row), 0)
    elif protocol == 'matched-observed':
        observed = read_table(Path(spec['observations']))
        for column, value in spec['equals'].items():
            observed = observed.filter(pl.col(column) == value)
        pep_column = spec['peptide_column']
        observed = observed.filter(pl.col(pep_column).str.contains('^[LFIMVWYCHAGPTSQNDERK]+$'))
        pools = {}
        for gene, pep in observed.select(spec['mhc_column'], pep_column).iter_rows():
            pools.setdefault((normalize_mhc(gene), len(pep)), set()).add(pep)
        rng = np.random.default_rng(1)
        for pid in ids:
            cog = info[pid]['peptide']
            pool = sorted(pools.get((normalize_mhc(spec['mhc'].get(pid)), len(cog)), set()) - {cog})
            accounting.append({'id': pid, 'pool_size': len(pool), 'eligible': len(pool) > 100,
                'reason': '' if len(pool) > 100 else 'at_most_100_MHC_length_matched_candidates'})
            if len(pool) <= 100:
                continue
            if len(pool) > 1000:
                pool = [pool[k] for k in rng.choice(len(pool), 1000, replace=False)]
            append(pid, cog, 1)
            for pep in pool:
                append(pid, pep, 0)
    elif protocol == 'ranked-assay':
        observed = read_table(Path(spec['observations']))
        rng = np.random.default_rng(0)
        for receptor, pid in spec['templates'].items():
            selected = observed.filter(pl.col(spec['receptor_column']) == receptor).sort(spec['rank_column'], descending=True)
            real = [v for v in selected[spec['peptide_column']] if set(v) <= set(AA)][:50]
            if len(real) != 50 or len({len(v) for v in real}) != 1:
                raise ValueError('Ranked assay requires 50 same-length standard peptides per template')
            cog = real[0]
            dec = rng.integers(0, 20, size=(1000, len(cog)))
            dec[:, 1], dec[:, -1] = AA.index(cog[1]), AA.index(cog[-1])
            for pep in real:
                append(pid, pep, 1)
            for row in dec:
                append(pid, ''.join(AA[k] for k in row), 0)
    elif protocol == 'receptor-shuffle':
        groups = {}
        for pid in ids:
            row = info[pid]
            if all(row[k] and set(row[k]) <= set(AA) for k in ('cdr3a', 'cdr3b', 'peptide')):
                groups.setdefault(tuple(len(row[k]) for k in ('cdr3a', 'cdr3b', 'peptide')), []).append(pid)
        for members in groups.values():
            if len(members) < 2:
                continue
            for pid in members:
                for donor in members:
                    rows.append({'complex.id': pid, 'candidate.id': donor, 'cdr3a': info[donor]['cdr3a'],
                                 'cdr3b': info[donor]['cdr3b'], 'y': int(pid == donor)})
    elif protocol == 'observed-matched':
        observed = read_table(Path(spec['observations'])).select(peptide=spec['peptide_column'], HLA=spec['mhc_column']).unique(maintain_order=True)
        templates = {(r['mhc'], r['length']): r['id'] for r in spec['templates']}
        for pep, hla in observed.iter_rows():
            key = (hla.replace(':', '').lower(), len(pep))
            valid = key in templates and set(pep) <= set(AA)
            accounting.append({'peptide': pep, 'HLA': hla, 'eligible': valid,
                'reason': '' if valid else 'no_matching_HLA_length_template_or_nonstandard_peptide'})
            if valid:
                append(templates[key], pep, int(pep in spec['positive_peptides']))
    elif protocol == 'observed-union':
        candidates = read_table(Path(spec['observations']))[spec['peptide_column']].to_list()
        cognates = read_table(Path(spec['cognates']), schema_overrides={spec['id_column']: pl.String})
        cognate = {str(pid): peptide for pid, peptide in cognates.select(spec['id_column'], spec['cognate_column']).iter_rows()}
        for pid in spec['ids']:
            cog = cognate[str(pid)]
            for pep in sorted(set(candidates) | {cog}):
                if set(pep) <= set(AA):
                    append(pid, pep, int(pep == cog))
    else:
        raise ValueError('Unknown historical protocol')
    if not rows:
        raise ValueError('No candidates satisfy the declared historical protocol')
    return pl.DataFrame(rows), pl.DataFrame(accounting)


def command(
    protocol: str = typer.Option(..., help='anchored-random, matched-observed, ranked-assay, receptor-shuffle, observed-matched or observed-union'),
    config: Path = typer.Option(..., exists=True, dir_okay=False),
    output: Path = typer.Option(..., '--output', '-o'),
):
    """Regenerate explicitly requested historical candidates; no scoring or fitting."""
    from . import __version__
    from .exports import _sha256
    if output.exists():
        raise ValueError('Output exists; choose a new directory')
    spec = json.loads(config.read_text())
    if not isinstance(spec, dict):
        raise ValueError("Config must be a JSON object")
    try:
        rows, accounting = generate_candidates(protocol, spec)
    except KeyError as exc:
        raise ValueError(f"Missing protocol input or column: {exc}") from exc
    output.mkdir(parents=True)
    rows.write_csv(output/'candidates.tsv', separator='\t')
    if accounting.height:
        accounting.write_csv(output/'accounting.tsv', separator='\t')
    sources = [config, *[Path(spec[k]) for k in ('markup', 'observations', 'cognates') if k in spec]]
    (output/'SOURCES.json').write_text(json.dumps({'version': __version__, 'protocol': protocol,
        'parameters': spec, 'inputs': {str(p): _sha256(p) for p in sources}}, indent=2)+'\n')


def markup_command(
    residues: Path = typer.Option(..., exists=True, dir_okay=False),
    contacts: Path = typer.Option(..., exists=True, dir_okay=False),
    membership: Path = typer.Option(..., exists=True, dir_okay=False, help='Frozen membership CSV with pdb.id'),
    output: Path = typer.Option(..., '--output', '-o'),
):
    """Reduce fresh annotated residues to ordered CDR3/peptide strings and eligibility."""
    from .exports import _sha256
    if output.exists():
        raise ValueError('Output exists; choose a new directory')
    markup = markup_table(read_table(residues))
    ids = sorted(read_table(membership)['pdb.id'].unique().to_list())
    ledger = native_accounting(markup, read_table(contacts), ids)
    output.mkdir(parents=True)
    markup.write_csv(output/'markup.csv')
    ledger.write_csv(output/'accounting.tsv', separator='\t')
    (output/'SOURCES.json').write_text(json.dumps({'inputs': {str(p): _sha256(p) for p in [residues, contacts, membership]}}, indent=2)+'\n')
