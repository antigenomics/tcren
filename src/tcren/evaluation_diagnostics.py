"""Explicit descriptive reducers for observed quantities; never fit recognition scores."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import t
import typer

from .associations import read_table


def regression_bands(rows, x, y):
    """OLS line and 95% Student-t mean-response/prediction bands on 60 grid points."""
    a, b = rows[x].cast(pl.Float64).to_numpy(), rows[y].cast(pl.Float64).to_numpy()
    keep = np.isfinite(a) & np.isfinite(b)
    a, b = a[keep], b[keep]
    if len(a) < 3 or np.ptp(a) == 0:
        raise ValueError("Regression bands require three finite pairs and variable x")
    slope, intercept = np.polyfit(a, b, 1)
    residual = b - (intercept + slope*a)
    scale = np.sqrt(residual @ residual / (len(a)-2))
    grid = np.linspace(a.min(), a.max(), 60)
    leverage = 1/len(a) + (grid-a.mean())**2 / ((a-a.mean())**2).sum()
    half = t.ppf(.975, len(a)-2) * scale
    fit = intercept + slope*grid
    return pl.DataFrame(dict(x=grid, fit=fit, lo=fit-half*np.sqrt(leverage),
                            hi=fit+half*np.sqrt(leverage), plo=fit-half*np.sqrt(1+leverage),
                            phi=fit+half*np.sqrt(1+leverage), n_used=[len(a)]*60))


def wild_type_affinity(rows):
    """ATLAS supplied free-energy WT medians; preserve the leading-number source policy."""
    required = ['true_PDB', 'DeltaG_kcal_per_mol', 'TCR_mut', 'MHC_mut', 'PEP_mut']
    if set(required)-set(rows.columns):
        raise ValueError(f"Required source columns: {required}")
    mutation = [~pl.col(c).cast(pl.String).fill_null('').str.strip_chars().str.to_lowercase()
                .is_in(['wt', '', 'nan', r'\n', 'none', '-']) for c in required[2:]]
    return (rows.with_columns(
        pl.col('true_PDB').cast(pl.String).str.strip_chars().str.to_lowercase().alias('structure_name'),
        pl.col('DeltaG_kcal_per_mol').cast(pl.String).str.extract(r'^\s*(-?[0-9.]+)', 1)
        .cast(pl.Float64, strict=False).alias('atlas_dG'), pl.any_horizontal(mutation).alias('mutated'))
        .filter(~pl.col('mutated') & pl.col('atlas_dG').is_finite() & (pl.col('structure_name') != ''))
        .group_by('structure_name').agg(pl.col('atlas_dG').median()).sort('structure_name'))


def backbone_deformation(rows, groups):
    """Summarise matched intrinsic windows, wrapping ca_torsion differences to [-pi, pi]."""
    rows = rows.with_columns((pl.col('value_model')-pl.col('value_ref')).alias('delta'),
                            (pl.col('value_model').is_finite() & pl.col('value_ref').is_finite())
                            .fill_null(False).alias('used'))
    delta = rows['delta'].cast(pl.Float64).to_numpy().copy()
    torsion = rows['metric'].to_numpy() == 'ca_torsion'
    delta[torsion] = np.arctan2(np.sin(delta[torsion]), np.cos(delta[torsion]))
    rows = rows.with_columns(pl.Series('delta', delta).fill_nan(None))
    summary = rows.group_by(groups, maintain_order=True).agg(
        pl.col('used').sum().alias('n_windows'),
        pl.col('value_ref').is_finite().fill_null(False).sum().alias('n_ref'),
        pl.col('value_model').is_finite().fill_null(False).sum().alias('n_model'),
        (~pl.col('used')).sum().alias('n_missing'),
        pl.col('delta').filter(pl.col('used')).pow(2).mean().sqrt().alias('rms_difference'),
        pl.col('delta').filter(pl.col('used')).abs().max().alias('max_difference'))
    return {'windows': rows, 'summary': summary}


def grouped_summary(rows, columns, groups, statistic='mean'):
    """Explicit group aggregates or within-group centering with observed-value means."""
    if not groups or not columns or statistic not in {'mean', 'median', 'min', 'max', 'center', 'standardize'}:
        raise ValueError('Explicit groups, columns and mean/median/min/max/center/standardize required')
    expressions = [pl.col(c).cast(pl.Float64).fill_nan(None) for c in columns]
    if statistic == 'standardize':
        return rows.with_columns([pl.when(expr.is_null().any().over(groups))
            .then(float('nan')).otherwise((expr-expr.mean().over(groups))/expr.std(ddof=1).over(groups))
            .alias(name) for name, expr in zip(columns, expressions)])
    if statistic == 'center':
        return rows.with_columns([(expr-expr.mean().over(groups)).alias(name)
                                  for name, expr in zip(columns, expressions)])
    return rows.group_by(groups, maintain_order=True).agg(
        [getattr(expr, statistic)().alias(name) for name, expr in zip(columns, expressions)],
        n_rows=pl.len())


def rank_selection(rows, groups, x, y):
    """Report first-tied maximum-confidence rank against maximum observed outcome."""
    result = []
    for group, part in rows.group_by(groups, maintain_order=True):
        score, value = part[x].to_numpy(), part[y].to_numpy()
        if not np.isfinite(score).all() or not np.isfinite(value).all():
            raise ValueError('Rank selection requires finite confidence and outcome')
        selected, best = int(np.argmax(score)), int(np.argmax(value))
        result.append({**dict(zip(groups, group)), 'n_models': part.height,
            'selected_rank': part['rank'][selected], 'selected_value': value[selected],
            'best_rank': part['rank'][best], 'best_value': value[best],
            'loss_to_best': value[best]-value[selected]})
    return pl.DataFrame(result)


def command(
    results: Path = typer.Option(..., exists=True, dir_okay=False),
    operation: str = typer.Option(..., help='regression-bands | wild-type-affinity | backbone-deformation | grouped-summary | rank-selection | backbone-extrema | distribution-summary | context-correlations'),
    output: Path = typer.Option(..., '--output', '-o'),
    groups: str = typer.Option('', help='Comma-separated summary groups'),
    x: str = typer.Option('x'), y: str = typer.Option('y'),
    columns: str = typer.Option(''), statistic: str = typer.Option('mean'),
    contexts: str = typer.Option(''), keys: str = typer.Option('id'),
    reference: str = typer.Option('iptm'),
):
    """Export explicitly selected descriptive summaries from supplied observations."""
    from . import __version__
    from .exports import _sha256
    if output.exists():
        raise ValueError('Output exists; choose a new directory')
    rows = read_table(results)
    if operation == 'regression-bands':
        tables = {'bands': regression_bands(rows, x, y)}
    elif operation == 'wild-type-affinity':
        tables = {'affinity': wild_type_affinity(rows)}
    elif operation == 'grouped-summary':
        tables = {'summary': grouped_summary(rows, columns.split(',') if columns else [],
                                             groups.split(',') if groups else [], statistic)}
    elif operation == 'backbone-extrema':
        from .diagnostics import backbone_extrema
        tables = backbone_extrema(rows, keys.split(','))
    elif operation == 'distribution-summary':
        from .diagnostics import grouped_summary as distribution_summary
        if not groups or not columns:
            raise ValueError('Explicit groups and columns required')
        tables = {'summary': distribution_summary(rows, columns.split(','), groups.split(','))}
    elif operation == 'context-correlations':
        from .diagnostics import context_correlations
        if not groups or not columns or not contexts:
            raise ValueError('Explicit groups, columns and contexts required')
        tables = {'associations': context_correlations(rows, columns.split(','), contexts.split(','),
            groups.split(','), keys.split(','), reference)}
    elif operation == 'rank-selection':
        if not groups:
            raise ValueError('Explicit rank selection groups required')
        tables = {'selection': rank_selection(rows, groups.split(','), x, y)}
    elif operation == 'backbone-deformation':
        if not groups:
            raise ValueError('Explicit backbone summary groups required')
        tables = backbone_deformation(rows, groups.split(','))
    else:
        raise ValueError('Unknown diagnostic operation')
    output.mkdir(parents=True)
    for name, table in tables.items():
        table.write_csv(output/f'{name}.tsv', separator='\t')
    (output/'SOURCES.json').write_text(json.dumps(dict(version=__version__, operation=operation,
        inputs={str(results): _sha256(results)}, groups=groups, x=x, y=y, columns=columns, statistic=statistic,
        contexts=contexts, keys=keys, reference=reference), indent=2)+'\n')
