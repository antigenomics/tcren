"""Per-template frozen matrices preserve directed scores and one annotation batch."""

import polars as pl
import pytest
from typer.testing import CliRunner
from click import unstyle

import tcren.cli as cli
from tcren.potential import Potential
from tcren.runtime import available_cpus
from test_cli_score_receptors import _paired


def inputs(tmp_path, command):
    candidates = tmp_path / "candidates.tsv"
    if command == "score":
        candidates.write_text("complex.id\tpeptide\nfirst\tD\nsecond\tD\n")
    else:
        candidates.write_text(
            "complex.id\tcandidate.id\tcdr3a\tcdr3b\nfirst\tnative\tK\tA\nsecond\tnative\tK\tA\n"
        )
    directory = tmp_path / "maps"
    directory.mkdir()
    for name, value in [("low", 1.0), ("high", 3.0)]:
        Potential(
            name,
            pl.DataFrame(
                {
                    "residue.aa.from": ["K", "A"],
                    "residue.aa.to": ["D", "D"],
                    "value": [value, value],
                }
            ),
            ("A", "D", "K"),
        ).to_csv(directory / (name + ".csv"))
    mapping = directory / "potentials.tsv"
    mapping.write_text("complex.id\tpotential\nFIRST\tlow.csv\nsecond\thigh.csv\n")
    return candidates, mapping


@pytest.mark.parametrize("command", ["score", "score-receptors"])
def test_template_matrices_serial_parallel_and_one_annotation_batch(tmp_path, monkeypatch, command):
    if available_cpus() < 2:
        pytest.skip("requires two allocated CPUs")
    candidates, mapping = inputs(tmp_path, command)
    calls = []

    def typed(*args, **kwargs):
        calls.append(kwargs)
        return [_paired("first"), _paired("second")]

    monkeypatch.setattr(cli, "_typed_batch", typed)
    for threads in (1, 2):
        result = CliRunner().invoke(
            cli.app,
            [
                command,
                "-s",
                "unused.pdb",
                "-c",
                str(candidates),
                "--potential-map",
                str(mapping),
                "-t",
                str(threads),
                "-o",
                str(tmp_path / f"{threads}.csv"),
            ],
        )
        assert result.exit_code == 0, (result.output, result.exception)
    assert [x["threads"] for x in calls] == [1, 2]
    assert (tmp_path / "1.csv").read_bytes() == (tmp_path / "2.csv").read_bytes()
    scores = pl.read_csv(tmp_path / "1.csv")
    assert scores["score"].to_list() == [2.0, 6.0]


@pytest.mark.parametrize("command", ["score", "score-receptors"])
@pytest.mark.parametrize("defect", ["duplicate", "missing", "extra", "null", "file", "conflict"])
def test_bad_template_map_cannot_write_scores(tmp_path, monkeypatch, command, defect):
    candidates, mapping = inputs(tmp_path, command)
    monkeypatch.setattr(cli, "_typed_batch", lambda *a, **k: [_paired("first"), _paired("second")])
    text = mapping.read_text()
    if defect == "duplicate":
        text += "first\tlow.csv\n"
    if defect == "missing":
        text = text.replace("second\thigh.csv\n", "")
    if defect == "extra":
        text += "third\tlow.csv\n"
    if defect == "null":
        text = text.replace("low.csv", "")
    if defect == "file":
        text = text.replace("low.csv", "absent.csv")
    mapping.write_text(text)
    output = tmp_path / "scores.csv"
    result = CliRunner().invoke(
        cli.app,
        [
            command,
            "-s",
            "unused.pdb",
            "-c",
            str(candidates),
            "--potential-map",
            str(mapping),
            "-t",
            "1",
            "-o",
            str(output),
            *(["-p", "mj"] if defect == "conflict" else []),
        ],
    )
    assert result.exit_code != 0
    assert "potential" in unstyle(result.output) + str(result.exception)
    assert not output.exists()


@pytest.mark.parametrize("command", ["score", "score-receptors"])
def test_bundled_map_matches_single_potential(tmp_path, monkeypatch, command):
    candidates, mapping = inputs(tmp_path, command)
    mapping.write_text("complex.id\tpotential\nfirst\tmj\nsecond\tmj\n")
    monkeypatch.setattr(cli, "_typed_batch", lambda *a, **k: [_paired("first"), _paired("second")])
    for name, options in [("map", ["--potential-map", str(mapping)]), ("single", ["-p", "mj"])]:
        result = CliRunner().invoke(
            cli.app,
            [
                command,
                "-s",
                "unused.pdb",
                "-c",
                str(candidates),
                "-t",
                "1",
                "-o",
                str(tmp_path / (name + ".csv")),
                *options,
            ],
        )
        assert result.exit_code == 0, (result.output, result.exception)
    assert (tmp_path / "map.csv").read_bytes() == (tmp_path / "single.csv").read_bytes()
