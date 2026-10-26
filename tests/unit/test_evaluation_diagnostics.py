"""Exact descriptive conventions, separate from recognition-score fitting."""
import numpy as np
import polars as pl
from scipy.stats import t
from tcren.evaluation_diagnostics import regression_bands, wild_type_affinity, backbone_deformation


def test_regression_band_matches_ols_definition():
    x = np.array([1.,2.,4.,5.]); y = np.array([2.,1.,4.,7.])
    out = regression_bands(pl.DataFrame({'x':x,'y':y}), 'x','y')
    design = np.column_stack([np.ones(4), x])
    beta = np.linalg.lstsq(design,y,rcond=None)[0]
    grid = np.column_stack([np.ones(60),out['x'].to_numpy()])
    residual = y-design@beta
    variance = residual@residual/2
    half = t.ppf(.975,2)*np.sqrt(variance*np.einsum('ij,jk,ik->i',grid,np.linalg.inv(design.T@design),grid))
    np.testing.assert_allclose(out['fit'],grid@beta)
    np.testing.assert_allclose(out['lo'],grid@beta-half)
    assert (out['plo'] <= out['lo']).all()


def test_wild_type_affinity_median_and_leading_token_policy():
    rows=pl.DataFrame({'true_PDB':[' ABC ','abc','abc','def'],
        'DeltaG_kcal_per_mol':['-6.0','-8.0','-100','>5'],
        'TCR_mut':['WT','-', 'A1B',None], 'MHC_mut':['']*4,'PEP_mut':['']*4})
    assert wild_type_affinity(rows).to_dicts()==[{'structure_name':'abc','atlas_dG':-7.}]


def test_circular_deformation_keeps_missing_windows():
    rows=pl.DataFrame({'group':['a']*3,'metric':['ca_torsion']*3,
                      'value_ref':[3.,0.,None], 'value_model':[-3.,1.,2.]})
    out=backbone_deformation(rows,['group','metric'])
    row=out['summary'].row(0,named=True)
    assert (row['n_windows'],row['n_missing'],row['n_ref'],row['n_model'])==(2,1,2,3)
    np.testing.assert_allclose(row['rms_difference'],np.sqrt(((2*np.pi-6)**2+1)/2))


def test_grouped_summary_cli(tmp_path):
    from typer.testing import CliRunner
    from tcren.cli import app
    source=tmp_path/'rows.tsv'
    pl.DataFrame({'group':['a','a','b'],'value':[1.,3.,5.]}).write_csv(source,separator='\t')
    out=tmp_path/'out'
    result=CliRunner().invoke(app,['evaluate-diagnostics','--results',str(source),
        '--operation','grouped-summary','--groups','group','--columns','value','-o',str(out)])
    assert result.exit_code==0,result.output
    assert pl.read_csv(out/'summary.tsv',separator='\t')['value'].to_list()==[2.,5.]


def test_rank_selection_retains_first_tie():
    from tcren.evaluation_diagnostics import rank_selection
    rows=pl.DataFrame({'group':['a']*3,'rank':[0,1,2],'x':[.8,.8,.4],'y':[.3,.9,.6]})
    row=rank_selection(rows,['group'],'x','y').row(0,named=True)
    assert row['selected_rank']==0 and row['best_rank']==1
    np.testing.assert_allclose(row['loss_to_best'],.6)
