"""Explicit ranking protocol contract and exact tied groups."""
import json
import shutil
import subprocess

import polars as pl
import pytest
from typer.testing import CliRunner

from tcren.cli import app


def test_ranking_summary_raw_matrix_and_tied_scores(tmp_path):
    runtime = shutil.which('Rscript')
    if not runtime or subprocess.run([runtime, '-e', "stopifnot(all(vapply(c('Rcpp','data.table','jsonlite'),requireNamespace,logical(1),quietly=TRUE)))"], capture_output=True).returncode:
        pytest.skip('Optional R protocol dependencies unavailable')
    rows = tmp_path/'rows.tsv'
    pl.DataFrame({'id': ['AA', 'AC', 'CA', 'CC'], 'epitope': ['group']*4,
                  'y': [1, 0, 1, 0], 'iptm': [1., 1., 0., 0.],
                  'screen_score': [1., 1., 0., 0.]}).write_csv(rows, separator='\t')
    matrix = tmp_path/'matrix.tsv'
    pl.DataFrame({'position': [1, 2], 'A': [0, 0], 'C': [1, 2]}).write_csv(matrix, separator='\t')
    config = tmp_path/'config.json'
    config.write_text(json.dumps({'panels': [{'name': 'panel', 'rows': str(rows),
        'arms': {'screen_score': 1, 'iptm': 1}, 'display': ['screen_score', 'iptm'],
        'matrices': {'group': str(matrix)}, 'response_name': 'matrix_response'}]}))
    out = tmp_path/'out'
    result = CliRunner().invoke(app, ['ranking-summary', '--config', str(config), '-o', str(out)])
    assert result.exit_code == 0, result.exception
    table = pl.read_csv(out/'metrics.tsv', separator='\t')
    assert table.filter(pl.col('endpoint') == 'call')['value'].to_list() == [.5]*8
    graded = pl.read_csv(out/'panel_rows.tsv', separator='\t')['graded'].to_list()
    assert graded == pytest.approx([0., 1.584962500721156, 1., 2.])
