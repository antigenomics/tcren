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
