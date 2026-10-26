"""Historical potential recipes must be expressible exactly through the CLI."""
import numpy as np
import polars as pl
import pytest
from typer.testing import CliRunner

from tcren.cli import app
from tcren.potential import Potential, derive_tcren


def test_derive_cli_matches_weighted_cysteine_retaining_api(tmp_path):
    contacts = pl.DataFrame({
        "pdb.id": ["a", "a", "b", "b"],
        "residue.aa.from": ["C", "A", "A", "D"],
        "residue.aa.to": ["A", "D", "D", "A"],
        "vote": [0.5, 0.5, 1.0, 1.0],
    })
    source, out = tmp_path / "contacts.csv", tmp_path / "potential.csv"
    contacts.write_csv(source)
    result = CliRunner().invoke(app, [
        "derive-potential", "-i", str(source), "-o", str(out),
        "--weight-column", "vote", "--keep-cys",
    ])
    assert result.exit_code == 0, result.output
    expected = derive_tcren(contacts, weight_col="vote", drop_cys=False)
    actual = Potential.from_csv(out)
    keys = ["residue.aa.from", "residue.aa.to"]
    np.testing.assert_allclose(actual.matrix.sort(keys)["value"],
                               expected.matrix.sort(keys)["value"], atol=1e-14)
    assert actual.matrix.filter(pl.col("residue.aa.from") == "C").height == 20


@pytest.mark.parametrize("bad", [None, -1.0, float("nan"), float("inf")])
def test_derive_cli_rejects_invalid_contact_weights(tmp_path, bad):
    source, out = tmp_path / "contacts.csv", tmp_path / "potential.csv"
    pl.DataFrame({"pdb.id": ["a"], "residue.aa.from": ["A"],
                  "residue.aa.to": ["D"], "vote": [bad]},
                 schema_overrides={"vote": pl.Float64}).write_csv(source)
    result = CliRunner().invoke(app, [
        "derive-potential", "-i", str(source), "-o", str(out),
        "--weight-column", "vote",
    ])
    assert result.exit_code != 0
    assert "finite nonnegative" in result.output
    assert not out.exists()


def test_selected_derivation_requires_a_positive_weight(tmp_path):
    source, summary, out = (tmp_path / n for n in ["contacts.csv", "summary.csv", "out.csv"])
    pl.DataFrame({"pdb.id": ["a", "b"], "residue.aa.from": ["A", "D"],
                  "residue.aa.to": ["D", "A"], "vote": [0.0, 1.0]}).write_csv(source)
    pl.DataFrame({"pdb.id": ["a", "b"], "nonred": [True, False]}).write_csv(summary)
    result = CliRunner().invoke(app, [
        "derive-potential", "-i", str(source), "-o", str(out),
        "--summary", str(summary), "--nonred", "--weight-column", "vote",
    ])
    assert result.exit_code != 0
    assert "positive weight" in result.output
    assert not out.exists()


@pytest.mark.parametrize("interface,keep_cys", [("tcr_mhc",False),("peptide_mhc",True)])
def test_mhc_derivation_uses_directed_contacts_and_balanced_weights(monkeypatch, tmp_path, interface, keep_cys):
    import tcren.paper
    contacts = pl.DataFrame({
        'pdb.id':['a','a','b'], 'residue.aa.from':['A','D','A'],
        'residue.aa.to':['D','A','D'], 'chain.type.from':['TRA','TRB','TRA']})
    markup = pl.DataFrame({'pdb.id':['a','b'], 'cdr3a':['AAA','AAA'],
                           'cdr3b':['BBB','BBB'], 'peptide':['PEP','PEP']})
    captured = {}
    def annotate(path, **kwargs):
        captured.update(kwargs)
        return contacts, markup
    monkeypatch.setattr(tcren.paper, 'annotate_structure_set', annotate)
    out, used = tmp_path/'TCRen_m.csv', tmp_path/'contacts.csv'
    result = CliRunner().invoke(app, ['derive-potential','--structure-dir',str(tmp_path),
        '--interface',interface,'--balance','both','-t','2',
        '--keep-cys' if keep_cys else '--drop-cys',
        '-o',str(out),'--contacts-out',str(used)])
    assert result.exit_code == 0, result.output
    assert captured == {'interface':interface,'threads':2,'on_error':'raise'}
    table = pl.read_csv(used)
    assert table['structure_weight'].to_list() == [.5,.5,.5]
    actual = Potential.from_csv(out)
    expected = derive_tcren(contacts, weights={'a':.5,'b':.5}, drop_cys=not keep_cys)
    keys = ['residue.aa.from','residue.aa.to']
    np.testing.assert_allclose(actual.matrix.sort(keys)['value'],expected.matrix.sort(keys)['value'])
    assert actual.matrix.filter(pl.col('residue.aa.from')=='C').height == (20 if keep_cys else 0)
    values = {(a,b):v for a,b,v in actual.matrix.iter_rows()}
    assert values['A','D'] != values['D','A']
