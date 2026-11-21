from pathlib import Path

import pytest
from typer.testing import CliRunner

from tcren.cli import app


@pytest.mark.parametrize("task, kind, label", [("tcr", "binary", "2"), ("epitope", "continuous", "nan")])
def test_invalid_labels_fail_before_output(tmp_path: Path, task, kind, label):
    (tmp_path / "model.pdb").touch()
    metadata = tmp_path / "models.tsv"
    metadata.write_text(
        f"structure_id\tstructure_path\tiptm\tplddt\tlabel\nm1\tmodel.pdb\t0.7\t80\t{label}\n"
    )
    output = tmp_path / "dashboard"
    result = CliRunner().invoke(app, [
        "dashboard", "import", str(metadata), "--task", task, "--label-kind", kind,
        "--output", str(output),
    ])
    assert result.exit_code == 2
    assert "invalid label" in result.output.lower()
    assert not output.exists()


@pytest.mark.parametrize("task, kind, labels", [
    ("tcr", "binary", [1, 0]), ("epitope", "continuous", [1.25, -2.5]),
])
def test_import_uses_screen_and_exports_the_declared_ids(tmp_path, monkeypatch, task, kind, labels):
    import json
    from types import SimpleNamespace

    import polars as pl
    from polars.testing import assert_frame_equal

    from tcren import cli, screen
    from tcren.complementarity import import_command
    from tcren.complementarity.dashboard import write_dashboard

    from test_complementarity_views import example

    ids = ["001", "<model & 2>"]
    metadata = tmp_path / "models.tsv"
    for folder in ("first", "second"):
        (tmp_path / folder).mkdir()
        (tmp_path / folder / "same.pdb").touch()
    pl.DataFrame({
        "structure_id": ids, "structure_path": ["first/same.pdb", "second/same.pdb"],
        "iptm": [0.7, 0.8], "plddt": [82, 94], "label": labels,
    }).write_csv(metadata, separator="\t")

    def typed_batch(sources, *args, **kwargs):
        return [SimpleNamespace(pdb_id=path.stem) for path in sorted(Path(sources[0]).glob("*.pdb"))]

    monkeypatch.setattr(cli, "_typed_batch", typed_batch)
    monkeypatch.setattr(screen, "identities", lambda models: (
        (["AAA", "AAA"], [("a", "b"), ("c", "d")]) if task == "tcr"
        else (["AAA", "AAB"], [("a", "b"), ("a", "b")])
    ))

    def descriptors(models, *args, **kwargs):
        return pl.DataFrame({
            "complex.id": [s.pdb_id for s in models],
            "n_contacts_tp": [0, 25],
            "dFT_tcr_pep_cdr": [0.0, 2.0],
            "dPhi_tcr_pep": [0.0, 2.0],
            "g_recognition_cyclo_frac": [0.0, 0.5],
            "enrichment": [0.0, 0.6],
            "sc_shape": [0.0, 0.6],
        })

    monkeypatch.setattr(cli, "_screen_descriptors", descriptors)

    def export_synthetic(**kwargs):
        assert kwargs["task"] == task and kwargs["portable"]
        scene_inputs = pl.read_csv(kwargs["metadata"], separator="\t", infer_schema=False)
        assert scene_inputs.columns == ["structure_id", "structure_path"]
        assert scene_inputs["structure_id"].to_list() == ids
        scenes = [dict(example(task), id=identity) for identity in ids]
        write_dashboard(scenes, kwargs["output"])

    monkeypatch.setattr(import_command, "complementarity_command", export_synthetic)
    output = tmp_path / "dashboard"
    result = CliRunner().invoke(app, [
        "dashboard", "import", str(metadata), "--task", task, "--label-kind", kind,
        "--output", str(output), "-t", "1",
    ])
    assert result.exit_code == 0, (result.output, result.exception)
    actual = pl.read_csv(output / "scores.tsv", separator="\t", schema_overrides={"structure_id": pl.String})
    assert actual["structure_id"].to_list() == ids
    assert actual["label"].to_list() == labels
    assert actual["screen_score"][0] == 0
    assert actual["score_unavailable"][1]

    direct_inputs = tmp_path / "direct"
    direct_inputs.mkdir()
    for i in range(2):
        (direct_inputs / f"model{i:08d}.pdb").touch()
    confidence = tmp_path / "confidence.tsv"
    pl.DataFrame({"complex.id": ["model00000000", "model00000001"],
                  "iptm": [0.7, 0.8], "plddt": [82., 94.]}).write_csv(confidence, separator="\t")
    direct = tmp_path / "direct.tsv"
    result = CliRunner().invoke(app, ["screen", "-s", str(direct_inputs), "--task", task,
                                    "--confidence", str(confidence), "-o", str(direct), "-t", "1"])
    assert result.exit_code == 0, (result.output, result.exception)
    expected = pl.read_csv(direct, separator="\t")
    assert_frame_equal(actual.select("screen_score", "score_unavailable"),
                       expected.select("screen_score", "score_unavailable"))
    for view in ("maps", "views"):
        assert (output / view / "index.html").is_file()
        index = json.loads((output / view / "index.json").read_text())
        assert [row["id"] for row in index] == ids
    for directory in (output, output / "maps", output / "views"):
        assert (directory / "dashboard.css").read_bytes() == (output / "dashboard.css").read_bytes()
        assert (directory / "index.html").read_text().count('href="dashboard.css"') == 1
    home = (output / "index.html").read_text()
    assert home.index('name="color-scheme"') < home.index('href="dashboard.css"')
    assert 'data-theme-only' in home
    assert (output / "view-controls.js").read_bytes() == (output / "maps" / "view-controls.js").read_bytes()
    assert "&lt;model &amp; 2&gt;" in home
    assert "views/?group=reference&amp;id=%3Cmodel+%26+2%3E" in home
