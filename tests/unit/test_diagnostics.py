"""Preserve tied extrema and missing-label descriptive groups."""
import polars as pl
import pytest
from tcren.diagnostics import backbone_extrema, grouped_summary, context_correlations


def test_extrema_keep_ties_and_separate_regions():
    rows = pl.DataFrame({'id': ['a']*4, 'chain.type': ['PEPTIDE']*3+['TRA'],
        'region.type': ['PEPTIDE']*3+['CDR3'], 'metric': ['ca_angle']*4,
        'units': ['degrees']*4, 'value': [10., 20., 20., 40.]})
    out = backbone_extrema(rows)
    assert out['extreme_windows'].height == 4
    assert out['geometry_wide']['peptide__ca_angle__maximum'][0] == 20.
    assert out['geometry_wide']['CDR3alpha__ca_angle__minimum'][0] == 40.


def test_summary_and_correlations_retain_unknown_label_group():
    rows = pl.DataFrame(dict(id=list('abcdef'), y=[0.,0.,0.,None,None,None],
        feature=[1.,2.,3.,3.,2.,1.], context=[1.,2.,3.,1.,2.,3.], iptm=[.1,.2,.3,.1,.2,None]))
    summary = grouped_summary(rows, ['feature'], ['y'])
    assert summary['n_models'].to_list() == [3,3]
    assert summary['q25'].to_list() == [1.5,1.5]
    correlations = context_correlations(rows, ['feature'], ['context'], ['y'])
    assert correlations['n_matched'].to_list() == [3,2]
    assert correlations['pearson_R'][0] == pytest.approx(1.)
    assert correlations['pearson_R'][1] is None
    assert correlations['y'].null_count() == 1


def test_distribution_cli_and_provenance(tmp_path):
    import json
    from typer.testing import CliRunner
    from tcren.cli import app
    rows = pl.DataFrame({'id': ['a','b'], 'label': [None,None], 'value': [1.,3.]})
    path = tmp_path/'rows.tsv'
    rows.write_csv(path, separator='\t')
    args = ['evaluate-diagnostics', '--operation', 'distribution-summary', '--results', str(path),
            '--groups', 'label', '--columns', 'value', '-o', str(tmp_path/'output')]
    result = CliRunner().invoke(app,args)
    assert result.exit_code == 0, result.output
    summary = pl.read_csv(tmp_path/'output/summary.tsv', separator='\t')
    assert summary['median'].to_list() == [2.]
    provenance = json.loads((tmp_path/'output/SOURCES.json').read_text())
    assert provenance['operation'] == 'distribution-summary'
    assert CliRunner().invoke(app,args).exit_code != 0
