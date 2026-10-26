"""Exact established R ranking summaries from an explicit, portable table contract."""
from importlib.resources import files
import json
from pathlib import Path
import re
import shutil
import subprocess

import typer


def command(
    config: Path = typer.Option(..., exists=True, dir_okay=False),
    output: Path = typer.Option(..., '--output', '-o'),
):
    """Summarize supplied panels; requires optional Rscript, data.table, jsonlite and Rcpp.

    JSON panels declare name, rows, arms {column: +/-1}, display columns and
    response_name. Optional response_column supplies continuous outcomes; matrices
    maps group names to raw positional matrices (position plus amino-acid columns).
    Rows use id, epitope, y, iptm and screen_score; optional native_structure_exists
    and flag_* columns retain context and reporting counts. No rows are filtered
    implicitly: callers supply the scientifically reviewed population explicitly.
    """
    from . import __version__
    from .exports import _sha256

    if output.exists():
        raise ValueError('Output exists; choose a new directory')
    runtime = shutil.which('Rscript')
    if not runtime:
        raise RuntimeError('ranking-summary requires optional Rscript with data.table, jsonlite and Rcpp')
    check = subprocess.run([runtime, '-e', "stopifnot(requireNamespace('data.table',quietly=TRUE),requireNamespace('jsonlite',quietly=TRUE),requireNamespace('Rcpp',quietly=TRUE))"], capture_output=True, text=True)
    if check.returncode:
        raise RuntimeError('Install optional R packages data.table, jsonlite and Rcpp for ranking-summary')
    spec = json.loads(config.read_text())
    if not isinstance(spec, dict) or not isinstance(spec.get('panels'), list) or not spec['panels']:
        raise ValueError('config must contain nonempty panels')
    sources = [config]
    names = []
    for panel in spec['panels']:
        required = {'name', 'rows', 'arms', 'display', 'response_name'}
        if not isinstance(panel, dict) or required - panel.keys():
            raise ValueError('Each panel requires name, rows, arms, display and response_name')
        name = panel['name']
        if not isinstance(name, str) or not re.fullmatch(r'[A-Za-z0-9_-]+', name) or name in names:
            raise ValueError('Panel names must be unique safe file stems')
        names.append(name)
        if not panel.get('response_name') or not panel.get('arms') or not panel.get('display'):
            raise ValueError('Each panel requires response_name, arms and display')
        if not isinstance(panel['arms'], dict) or any(type(v) is not int or v not in (-1, 1) for v in panel['arms'].values()):
            raise ValueError('arms must map score columns to integer +1/-1')
        if not isinstance(panel['display'], list) or set(panel['display']) - panel['arms'].keys():
            raise ValueError('display must list columns declared in arms')
        if not isinstance(panel.get('matrices', {}), dict):
            raise ValueError('matrices must map group names to input paths')
        for option in ('required_true', 'required_false'):
            if option in panel and (not isinstance(panel[option], list) or any(not isinstance(x, str) for x in panel[option])):
                raise ValueError(f'{option} must list boolean column names')
        paths = [panel['rows'], *panel.get('matrices', {}).values()]
        for value in paths:
            path = Path(value)
            if not path.is_absolute() or not path.is_file():
                raise ValueError('Panel input paths must be existing absolute files')
            sources.append(path)
    data = files('tcren.data')
    script = Path(str(data.joinpath('ranking_summary.R')))
    subprocess.run([runtime, str(script), str(config.resolve()), str(output.resolve()), str(script.parent)], check=True)
    (output / 'SOURCES.json').write_text(json.dumps({'version': __version__,
        'inputs': {str(p): _sha256(p) for p in sources},
        'resources': {name: _sha256(Path(str(data.joinpath(name)))) for name in
                      ['ranking_summary.R', 'ranking_metrics.R', 'tsv_exact.R']},
        'test_null': 'Pearson correlation equals zero', 'parameters': spec}, indent=2)+'\n')
