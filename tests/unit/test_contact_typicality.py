from math import log

import numpy as np
import polars as pl
import pytest

from tcren.contact_typicality import contact_typicality
from tcren.native_reference import load_native_reference


def test_typicality_keeps_zero_and_unknown_contact_rows():
    probability = load_native_reference()["contact_profiles"]["atom/tcr_peptide"]["probability"]
    cell = max(probability, key=probability.get)
    observations = pl.DataFrame({"id": ["known", "known", "unknown"],
                                 "cell": [cell, cell, "unsupported"]})
    result = contact_typicality(observations, ["unknown", "empty", "known"])
    assert result["id"].to_list() == ["unknown", "empty", "known"]
    assert result["contact_count"].to_list() == [1, 0, 2]
    assert result["typicality_reason"].to_list() == ["unknown_cell", "no_contacts", None]
    assert result["contact_action"][2] == pytest.approx(-2 * log(probability[cell]))
    assert result["contact_typicality"][2] < 0
    assert result["contact_typicality_mcse"][2] > 0
    assert result.equals(contact_typicality(observations.reverse(), result["id"].to_list()))


def test_typicality_requires_declared_structures():
    observations = pl.DataFrame({"id": ["unknown"], "cell": ["unsupported"]})
    with pytest.raises(ValueError, match="structure IDs"):
        contact_typicality(observations, ["other"])


def test_cli_typicality_matches_the_frozen_contact_definition(tmp_path, monkeypatch):
    from typer.testing import CliRunner

    from tcren import cli
    from tcren.structure.model import Atom, Chain, Residue, Structure

    def residue(i, x, name, element="C"):
        coord = np.asarray([x, 0.0, 0.0])
        atom = Atom(name=name, element=element, coord=coord)
        deuterium = Atom(name="D", element="D", coord=coord)
        return Residue(i, i, "", "A", "ALA", (atom, deuterium))

    structure = Structure("pose", [
        Chain("A", [residue(0, 0, "CA"), residue(1, 0, "CG")], "TRA"),
        Chain("B", [residue(0, 3, "CA"), residue(1, 20, "CB")], "PEPTIDE"),
    ])
    output = tmp_path / "typicality.csv"
    monkeypatch.setattr(cli, "_typed_batch", lambda *args, **kwargs: [structure])
    result = CliRunner().invoke(cli.app, [
        "contacts", "-s", str(tmp_path), "-o", str(output), "--typicality", "--atom-pairs", "all",
    ])
    assert result.exit_code == 0, result.output
    row = pl.read_csv(output).row(0, named=True)
    assert row["id"] == "pose"
    assert row["contact_count"] == 1
    assert row["typicality_reason"] is None
