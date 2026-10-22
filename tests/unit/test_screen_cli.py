"""Screen CLI wiring retains the interface count needed to discard failed poses."""
from types import SimpleNamespace

import polars as pl
import pytest
from typer.testing import CliRunner

from tcren.screen import CohortKind


@pytest.mark.parametrize("task", list(CohortKind))
def test_screen_cli_discards_zero_contact_models_and_reports_flags(task, monkeypatch, tmp_path):
    import tcren.cli as cli
    import tcren.contactmap as contactmap
    import tcren.energetics.mutation as mutation
    import tcren.exports as exports
    import tcren.recognition as recognition
    import tcren.screen as screen

    structures = [SimpleNamespace(pdb_id=f"m{i}", chains=[
        SimpleNamespace(chain_type="PEPTIDE", sequence=lambda: "AAA")]) for i in range(2)]
    monkeypatch.setattr(cli, "_typed_batch", lambda *a, **kw: structures)
    peps, tcrs = (["AAA"] * 2, [("a", "b"), ("c", "d")])
    if task is CohortKind.EPITOPE_RANKING:
        peps, tcrs = ["AAA", "BBB"], [("a", "b")] * 2
    monkeypatch.setattr(screen, "identities", lambda s: (peps, tcrs))
    monkeypatch.setattr(recognition, "recognition_table", lambda *a, **kw: [
        {"complex.id": s.pdb_id, "sc_shape": 0.6, "g_cyclo_frac": 0.2} for s in structures])
    monkeypatch.setattr(contactmap.ContactMap, "from_structure", lambda s, **kw:
        SimpleNamespace(tcr_peptide=lambda: pl.DataFrame({"edge": list(range(
            0 if s.pdb_id == "m0" else 20))})))
    monkeypatch.setattr(exports, "receptor_reference", lambda *a, **kw:
        {"dFT_tcr_pep_cdr": 1.0})
    monkeypatch.setattr(mutation, "reference_delta", lambda *a, **kw: 1.0)
    monkeypatch.setattr(cli, "_screen_confidence", lambda *a, **kw:
        pl.DataFrame({"complex.id": ["m0", "m1"], "iptm": [0.4, 0.5]}))
    out = tmp_path / "screen.tsv"
    result = CliRunner().invoke(cli.app, ["screen", "-s", "models", "-o", str(out),
                                         "--min-cohort", "1"])
    assert result.exit_code == 0, result.output
    table = pl.read_csv(out, separator="\t")
    assert table["no_interface"].to_list() == [True, False]
    assert table["screen_score"][0] is None
    assert table["screen_score"][1] is not None
    assert table["cohort_n_no_interface"].to_list() == [1, 1]
    assert "discarded with no TCR:peptide interface" in result.output
