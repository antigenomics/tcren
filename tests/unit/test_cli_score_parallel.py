"""Candidate CLI workers retain each template's ranking, options and failures."""
from dataclasses import replace

import polars as pl
import pytest
from typer.testing import CliRunner

import tcren.cli as cli
from tcren.runtime import available_cpus
from test_contact_observations import complex_fixture
from tcren.structure.model import RegionMarkup


def _complex():
    structure = complex_fixture()
    peptide = structure.chain("P")
    peptide.regions = [RegionMarkup("PEPTIDE", 0, 0, "D", peptide.residues)]
    return structure


@pytest.mark.parametrize("options", [[], ["--delta", "--reference", "virtual"],
    ["--delta", "--reference", "structural"],
    ["--drop-untyped", "--position-weights", "central", "--intra-weight", ".5"],
    ["--soft", "--delta"]])
@pytest.mark.parametrize("manifest", [False, True])
def test_candidate_score_serial_parallel_exact(tmp_path, monkeypatch, options, manifest):
    if available_cpus() < 2:
        pytest.skip("requires two allocated CPUs")
    structures = [replace(_complex(), pdb_id=f"synthetic-{i}") for i in (3, 1, 2)]
    calls = []
    def typed(*args, **kwargs):
        calls.append(kwargs)
        return structures
    monkeypatch.setattr(cli, "_typed_batch", typed)
    candidates = tmp_path/"candidates.txt"
    if manifest:
        candidates.write_text("complex.id\tpeptide\nSYNTHETIC-3\tD\nsynthetic-3\tA\nsynthetic-3\tDD\nsynthetic-1\tV\nsynthetic-1\tD\n")
    else:
        candidates.write_text("D\nA\nV\nA\nDD\n")
    for threads in (1, 2):
        result = CliRunner().invoke(cli.app, ["score", "-s", "unused.pdb", "-c", str(candidates),
            "-o", str(tmp_path/f"{threads}.csv"), "-t", str(threads), *options])
        assert result.exit_code == 0, (result.output, result.exception)
    assert calls == [dict(mhc="--soft" in options, threads=count) for count in (1, 2)]
    assert (tmp_path/"1.csv").read_bytes() == (tmp_path/"2.csv").read_bytes()
    result = pl.read_csv(tmp_path/"2.csv")
    assert result["complex.id"].unique(maintain_order=True).to_list() == [s.pdb_id for s in structures[:2 if manifest else 3]]
    assert ("dphi" in result.columns) == ("--delta" in options and "--soft" not in options)


def test_candidate_parallel_empty_and_pool_failure(tmp_path, monkeypatch):
    if available_cpus() < 2:
        pytest.skip("requires two allocated CPUs")
    structures = [replace(_complex(), pdb_id=f"synthetic-{i}") for i in range(2)]
    monkeypatch.setattr(cli, "_typed_batch", lambda *args, **kwargs: structures)
    candidates, output = tmp_path/"candidates.txt", tmp_path/"scores.csv"
    candidates.write_text("")
    args = ["score", "-s", "unused.pdb", "-c", str(candidates), "-o", str(output), "-t", "2"]
    result = CliRunner().invoke(cli.app, args)
    assert result.exit_code == 0, result.output
    assert output.read_text() == "\n"
    output.unlink()
    candidates.write_text("D\n")
    def fail(*args):
        raise RuntimeError("deliberate pool failure")
    monkeypatch.setattr(cli, "parallel_map", fail)
    result = CliRunner().invoke(cli.app, args)
    assert result.exit_code != 0
    assert str(result.exception) == "deliberate pool failure"
    assert not output.exists()
