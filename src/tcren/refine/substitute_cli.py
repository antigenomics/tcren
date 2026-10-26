"""Batch backbone-preserving substitutions from an explicit template/variant table."""
from pathlib import Path
import json
import re

import polars as pl
import typer

from .substitute import substitute_peptide


def command(
    structures: str = typer.Option(..., '-s', '--structures'),
    variants: Path = typer.Option(..., exists=True, help='TSV columns template, id, peptide'),
    output: Path = typer.Option(..., '-o', '--output'),
    organism: str = typer.Option('human'),
    threads: int = typer.Option(0, '-t', '--threads'),
):
    """Write equal-length peptide variants without refinement or per-variant searches.

    Preserves backbone/C-beta and existing substitution semantics; non-alanine
    substitutions are side-chain stubs, not reconstructed full side chains.
    """
    from .. import __version__
    from ..annotation import classify_chains
    from ..annotation.batch import annotate_batch
    from ..exports import _sha256
    from ..runtime import resolve_threads
    from ..structure.io import iter_structures, parse_structure, resolve_sources, structure_paths, write_structure

    if output.exists():
        raise ValueError('Output exists; choose a new directory')
    rows = pl.read_csv(variants, separator='\t', schema_overrides={c: pl.String for c in ('template','id','peptide')})
    if not rows.height or set(('template','id','peptide'))-set(rows.columns):
        raise ValueError('Nonempty template, id, peptide table required')
    if any(rows[c].null_count() for c in ('template','id','peptide')) or rows['id'].n_unique()!=rows.height:
        raise ValueError('Present template/peptide and unique variant IDs required')
    if any(not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', key) for key in rows['id']):
        raise ValueError('Variant IDs must be safe filenames containing letters, digits, _, . or -')
    if any(not re.fullmatch(r'[ACDEFGHIKLMNPQRSTVWY]+', peptide) for peptide in rows['peptide']):
        raise ValueError('Peptides must use uppercase standard amino acids')
    sources = resolve_sources(structures)
    raw = [p for source in sources for p in structure_paths(source)]
    templates = [s for source in sources for _,s in iter_structures(source, importer=parse_structure)]
    by_id = {s.pdb_id:s for s in templates}
    if len(by_id)!=len(templates) or set(rows['template'])-set(by_id):
        raise ValueError('Template identifiers must be unique and contain every requested template')
    threads = resolve_threads(threads)
    annotations = annotate_batch(templates, organisms=(organism,), threads=threads)
    for structure, annotation in zip(templates, annotations):
        classify_chains(structure, autodetect_species=False, precomputed_records=annotation)
    # Validate every substitution before writing any partial output.
    prepared = []
    for row in rows.iter_rows(named=True):
        variant = substitute_peptide(by_id[row['template']], row['peptide'])
        variant.pdb_id = row['id']
        prepared.append((row, variant))
    output.mkdir(parents=True)
    accounting = []
    for row, variant in prepared:
        path = output/(row['id']+'.cif')
        write_structure(variant, path)
        accounting.append({**row, 'status':'ok', 'file':path.name, 'sha256':_sha256(path)})
    pl.DataFrame(accounting).write_csv(output/'accounting.tsv', separator='\t')
    (output/'SOURCES.json').write_text(json.dumps(dict(version=__version__, threads=threads,
        inputs={str(p):_sha256(p) for p in [variants,*raw]},
        n_templates=len(templates), n_variants=len(prepared), refinement=False),indent=2)+'\n')
