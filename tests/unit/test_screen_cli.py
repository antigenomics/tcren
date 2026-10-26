"""Screen CLI wiring retains the interface count needed to discard failed poses."""
from pathlib import Path
from types import SimpleNamespace

import polars as pl
import pytest
from typer.testing import CliRunner
from click import unstyle

from tcren.screen import CohortKind


@pytest.mark.parametrize("mode", ["metadata", "status"])
@pytest.mark.parametrize("gate", [None, True, False])
@pytest.mark.parametrize("task", list(CohortKind))
def test_screen_cli_discards_zero_contact_models_and_reports_flags(task, gate, mode, monkeypatch, tmp_path):
    import tcren.cli as cli
    import tcren.contactmap as contactmap
    import tcren.energetics.mutation as mutation
    import tcren.exports as exports
    import tcren.recognition as recognition
    import tcren.screen as screen
    import tcren.topology.graph as graph
    import tcren.contacts.motif as motif
    monkeypatch.setattr(cli, "_screen_residue_chunk", lambda structs: pl.DataFrame({"id": [s.pdb_id for s in structs]}))
    monkeypatch.setattr(motif, "bundled_enrichment", lambda residues: residues.with_columns(enrichment=pl.lit(.6)))
    monkeypatch.setattr(graph, "recognition_graph_features", lambda cm: {"g_recognition_cyclo_frac": .3})

    structures = [SimpleNamespace(pdb_id=f"m{i}", chains=[
        SimpleNamespace(chain_type="PEPTIDE", sequence=lambda: "AAA")]) for i in range(2)]
    monkeypatch.setattr(cli, "_typed_batch", lambda *a, **kw: structures)
    peps, tcrs = (["AAA"] * 2, [("a", "b"), ("c", "d")])
    if task is CohortKind.EPITOPE_RANKING:
        peps, tcrs = ["AAA", "BBB"], [("a", "b")] * 2
    monkeypatch.setattr(screen, "identities", lambda s: (peps, tcrs))
    def descriptors_without_reannotation(*args, **kwargs):
        assert kwargs["annotate"] is False
        return [{"complex.id": s.pdb_id, "sc_shape": 0.6, "g_cyclo_frac": 0.2}
                for s in structures]

    monkeypatch.setattr(recognition, "recognition_table", descriptors_without_reannotation)
    monkeypatch.setattr(contactmap.ContactMap, "from_structure", lambda s, **kw:
        SimpleNamespace(tcr_peptide=lambda: pl.DataFrame({"edge": list(range(
            0 if s.pdb_id == "m0" else 20))})))
    monkeypatch.setattr(exports, "receptor_reference", lambda *a, **kw:
        {"dFT_tcr_pep_cdr": 1.0})
    monkeypatch.setattr(mutation, "reference_delta", lambda *a, **kw: 1.0)
    monkeypatch.setattr(cli, "_screen_confidence", lambda *a, **kw:
        pl.DataFrame({"complex.id": ["m0", "m1"], "iptm": [0.4, 0.5], "plddt": [93., 93.]}))
    out = tmp_path / "screen.tsv"
    args = ["screen", "-s", "models", "-o", str(out), "--min-cohort", "1", "-t", "1"]
    if gate is not None:
        native = tmp_path / "native.tsv"
        native.write_text("epitope\n" + ("AAA" if gate else "CCC") + "\n")
        args += (["--template-metadata", str(native)] if mode == "metadata"
                 else ["--template-status", "present" if gate else "absent"])
    result = CliRunner().invoke(cli.app, args)
    if gate is not None and task is CohortKind.EPITOPE_RANKING:
        assert result.exit_code != 0
        assert "only to TCR ranking" in result.output
        assert not out.exists()
        return
    assert result.exit_code == 0, result.output
    table = pl.read_csv(out, separator="\t")
    assert table["no_interface"].to_list() == [True, False]
    assert table["screen_score"][0] is None
    assert table["screen_score"][1] is not None
    assert table["cohort_n_no_interface"].to_list() == [1, 1]
    assert "discarded with no TCR:peptide interface" in result.output

    if gate is not None:
        import json
        provenance = json.loads(Path(str(out) + ".provenance.json").read_text())
        assert table["native_template"].to_list() == [gate] * 2
        assert table["confidence_weight"].to_list() == pytest.approx([.666] * 2)
        assert provenance["template_gate"]["native_template"] is gate
        assert provenance["effective_weight"] == pytest.approx(.666)
        assert table["screen_score"][1] == pytest.approx(.666*.5+.246*.3+.088*.6)


def test_invalid_template_status_before_loading_structures(tmp_path):
    import tcren.cli as cli
    result = CliRunner().invoke(cli.app, ["screen", "-s", "missing", "--template-status", "unknown"])
    assert result.exit_code != 0
    assert "must be present or absent" in unstyle(result.output)


def test_screen_confidence_preserves_plddt(tmp_path):
    from tcren.cli import _screen_confidence
    table = tmp_path / "confidence.tsv"
    table.write_text("complex.id\tiptm\tplddt\na\t0.8\t93\n")
    assert _screen_confidence(table, []).row(0) == ("a", .8, 93.)


def test_screen_confidence_merges_inputs_and_rejects_conflicts(monkeypatch):
    from tcren.cli import _screen_confidence
    import tcren.metadata as metadata
    tables = {
        "one": pl.DataFrame({"id": ["a"], "iptm": [.8], "plddt": [93.]}),
        "two": pl.DataFrame({"id": ["b"], "iptm": [.7], "plddt": [92.]}),
        "conflict": pl.DataFrame({"id": ["a"], "iptm": [.9], "plddt": [93.]}),
    }
    monkeypatch.setattr(metadata, "read_metadata", tables.get)
    assert _screen_confidence(None, ["one", "two", "one"]).height == 2
    with pytest.raises(Exception, match="conflicting confidence metadata"):
        _screen_confidence(None, ["one", "conflict"])


def test_single_structure_nearest_context_is_explicit_and_checks_pmhc(monkeypatch, tmp_path):
    import tcren.cli as cli
    import tcren.screen as screen

    structs = [SimpleNamespace(pdb_id="one", chains=[SimpleNamespace(chain_type="MHCa", sequence=lambda: "MHC")])]
    monkeypatch.setattr(cli, "_typed_batch", lambda *a, **kw: structs)
    monkeypatch.setattr(screen, "identities", lambda s: (["AAA"] * len(s), [("a", "b")] * len(s)))
    monkeypatch.setattr(cli, "_screen_descriptors", lambda *a, **kw: pl.DataFrame({
        "complex.id": ["one"], "g_recognition_cyclo_frac": [.4], "enrichment": [.6],
        "dFT_tcr_pep_cdr": [0.], "sc_shape": [.5], "n_contacts_tp": [10]}))
    monkeypatch.setattr(cli, "_screen_confidence", lambda *a, **kw: pl.DataFrame({"complex.id": ["one"], "iptm": [.8], "plddt": [93.]}))
    output = tmp_path / "single.tsv"
    args = ["screen", "-s", "one.pdb", "--task", "tcr", "-o", str(output)]
    default = CliRunner().invoke(cli.app, args)
    assert default.exit_code != 0 and not output.exists()
    explicit = CliRunner().invoke(cli.app, args + ["--confidence-context", "nearest-plddt", "--min-cohort", "1"])
    assert explicit.exit_code == 0, explicit.output
    result = pl.read_csv(output, separator="\t")
    assert result["screen_score"][0] == pytest.approx(.666*.8 + .246*.4 + .088*.6)
    assert result["proxy_context_n"].to_list() == [1]
    assert result["confidence_context"].to_list() == ["nearest-plddt"]
    structs.append(SimpleNamespace(pdb_id="two", chains=[SimpleNamespace(chain_type="MHCa", sequence=lambda: "OTHER")]))
    mixed = CliRunner().invoke(cli.app, args + ["--confidence-context", "nearest-mean-plddt", "--min-cohort", "1"])
    assert mixed.exit_code != 0 and "MHC" in unstyle(mixed.output) and "sequences" in unstyle(mixed.output)
