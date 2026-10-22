"""Paired CDR3 threading preserves directed sums and candidate identities."""
from dataclasses import replace

import polars as pl
import pytest
from typer.testing import CliRunner
from click import unstyle

import tcren.cli as cli
from tcren.potential import Potential
from tcren.runtime import available_cpus
from tcren.structure.model import RegionMarkup
from test_cli_score_parallel import _complex


def _paired(name):
    structure = replace(_complex(), pdb_id=name)
    beta = structure.chain("B")
    residue = replace(structure.chain("T").residues[0], aa="A")
    beta.residues = [residue]
    beta.regions = [RegionMarkup("CDR3", 0, 0, "A", [residue])]
    return structure


def test_receptor_candidates_directed_scores_ids_and_spawn(tmp_path, monkeypatch):
    if available_cpus() < 2:
        pytest.skip("requires two allocated CPUs")
    structures = [_paired("first"), _paired("second")]
    monkeypatch.setattr(cli, "_typed_batch", lambda *a, **k: structures)
    # Opposite entries deliberately differ: substitution must use the 'from' axis.
    potential = Potential("directed", pl.DataFrame({"residue.aa.from": ["A", "K", "D", "D"],
                                       "residue.aa.to": ["D", "D", "A", "K"],
                                       "value": [2., 5., 20., 50.]}), ("A", "D", "K"))
    monkeypatch.setattr(cli, "_load_potential", lambda *a: potential)
    candidates = tmp_path / "candidates.tsv"
    candidates.write_text("complex.id\tcandidate.id\tcdr3a\tcdr3b\n"
                          "FIRST\tdonor1\tK\tA\nfirst\tdonor2\tK\tA\n"
                          "first\tdonor3\tA\tK\nsecond\tdonor1\tA\tA\n")
    for threads in (1, 2):
        result = CliRunner().invoke(cli.app, ["score-receptors", "-s", "unused.pdb",
            "-c", str(candidates), "-o", str(tmp_path / f"{threads}.csv"), "-t", str(threads)])
        assert result.exit_code == 0, (result.output, result.exception)
    assert (tmp_path / "1.csv").read_bytes() == (tmp_path / "2.csv").read_bytes()
    output = pl.read_csv(tmp_path / "1.csv")
    assert output["candidate.id"].to_list() == ["donor1", "donor2", "donor3", "donor1"]
    assert output["score.tra"].to_list() == [5., 5., 2., 2.]
    assert output["score.trb"].to_list() == [2., 2., 5., 2.]
    assert output["score"].to_list() == [7., 7., 7., 4.]
    explicit = tmp_path / "explicit-cdr3.csv"
    result = CliRunner().invoke(cli.app, ["score-receptors", "-s", "unused.pdb", "-c",
        str(candidates), "-o", str(explicit), "-t", "1", "--receptor-region", "cdr3"])
    assert result.exit_code == 0, (result.output, result.exception)
    assert explicit.read_bytes() == (tmp_path / "1.csv").read_bytes()


@pytest.mark.parametrize("rows,error", [
    ("first\td1\tKK\tA\n", "length 1"),
    ("missing\td1\tK\tA\n", "templates not found"),
    ("first\td1\tK\tA\nfirst\td1\tA\tK\n", "unique"),
    ("first\td1\tX\tA\n", "standard uppercase"),
])
def test_receptor_candidate_errors(tmp_path, monkeypatch, rows, error):
    monkeypatch.setattr(cli, "_typed_batch", lambda *a, **k: [_paired("first")])
    candidates, output = tmp_path / "candidates.tsv", tmp_path / "scores.csv"
    candidates.write_text("complex.id\tcandidate.id\tcdr3a\tcdr3b\n" + rows)
    result = CliRunner().invoke(cli.app, ["score-receptors", "-s", "unused.pdb", "-c",
        str(candidates), "-o", str(output), "-t", "1"])
    assert result.exit_code != 0
    assert error in result.output + str(result.exception)
    assert not output.exists()


def test_receptor_requires_each_physical_cdr3_and_rejects_peptide_options(tmp_path, monkeypatch):
    structure = _paired("first")
    structure.chain("B").regions = []
    monkeypatch.setattr(cli, "_typed_batch", lambda *a, **k: [structure])
    candidates = tmp_path / "candidates.tsv"
    candidates.write_text("complex.id\tcandidate.id\tcdr3a\tcdr3b\nfirst\td1\tK\tA\n")
    args = ["score-receptors", "-s", "unused.pdb", "-c", str(candidates), "-t", "1"]
    result = CliRunner().invoke(cli.app, args)
    assert result.exit_code != 0
    assert "exactly one observed TRB CDR3" in str(result.exception)
    result = CliRunner().invoke(cli.app, [*args, "--soft"])
    assert result.exit_code != 0 and "No such option" in result.output


@pytest.mark.parametrize("peptide_case", ["missing", "ambiguous", "nonstandard"])
def test_receptor_requires_one_standard_peptide(tmp_path, monkeypatch, peptide_case):
    structure = _paired("first")
    peptide = structure.chain("P")
    if peptide_case == "missing":
        structure.chains.remove(peptide)
    elif peptide_case == "ambiguous":
        structure.chains.append(replace(peptide, chain_id="Q"))
    else:
        peptide.residues = [replace(peptide.residues[0], aa="X")]
    monkeypatch.setattr(cli, "_typed_batch", lambda *a, **k: [structure])
    candidates, output = tmp_path / "candidates.tsv", tmp_path / "scores.csv"
    candidates.write_text("complex.id\tcandidate.id\tcdr3a\tcdr3b\nfirst\td1\tK\tA\n")
    result = CliRunner().invoke(cli.app, ["score-receptors", "-s", "unused.pdb", "-c",
        str(candidates), "-o", str(output), "-t", "1"])
    assert result.exit_code != 0
    assert "exactly one observed standard-amino-acid peptide" in str(result.exception)
    assert not output.exists()


@pytest.mark.parametrize("cutoff", ["0", "-1", "nan", "inf"])
def test_receptor_cutoff_rejected_before_annotation(tmp_path, monkeypatch, cutoff):
    def annotation_must_not_run(*args, **kwargs):
        pytest.fail("invalid cutoff reached annotation")
    monkeypatch.setattr(cli, "_typed_batch", annotation_must_not_run)
    result = CliRunner().invoke(cli.app, ["score-receptors", "-s", "unused.pdb", "-c",
        str(tmp_path / "unused.tsv"), "--cutoff", cutoff])
    assert result.exit_code != 0
    assert "cutoff must be finite and positive" in unstyle(result.output)


def _junction(name):
    structure = _paired(name)
    for kind, anchor in (("TRA", "F"), ("TRB", "W")):
        chain = structure.by_type(kind)[0]
        cdr = replace(chain.residues[0], seq_index=2)
        left = replace(cdr, seq_index=1, aa="C")
        right = replace(cdr, seq_index=3, aa=anchor)
        chain.residues = [left, cdr, right]
        chain.regions = [RegionMarkup("FR3", 1, 1, "C", [left]),
                         RegionMarkup("CDR3", 2, 2, cdr.aa, [cdr]),
                         RegionMarkup("FR4", 3, 3, anchor, [right])]
    return structure


def test_junction_observed_w_anchor_directed_sums_and_spawn(tmp_path, monkeypatch):
    if available_cpus() < 2:
        pytest.skip("requires two allocated CPUs")
    structures = [_junction(name) for name in ("first", "second")]
    monkeypatch.setattr(cli, "_typed_batch", lambda *a, **k: structures)
    potential = Potential("directed", pl.DataFrame({
        "residue.aa.from": ["C", "K", "A", "F", "W", "D"],
        "residue.aa.to": ["D", "D", "D", "D", "D", "W"],
        "value": [3., 5., 2., 7., 11., 110.]}), ("A", "C", "D", "F", "K", "W"))
    monkeypatch.setattr(cli, "_load_potential", lambda *a: potential)
    candidates = tmp_path / "junctions.tsv"
    candidates.write_text("complex.id\tcandidate.id\tjunctiona\tjunctionb\n"
                          "first\tnative\tCKF\tCAW\nfirst\tduplicate\tCKF\tCAW\n"
                          "second\tnative\tCKW\tCAF\n")
    for threads in (1, 2):
        result = CliRunner().invoke(cli.app, ["score-receptors", "-s", "unused.pdb", "-c",
            str(candidates), "--receptor-region", "junction", "-t", str(threads),
            "-o", str(tmp_path / f"junction-{threads}.csv")])
        assert result.exit_code == 0, (result.output, result.exception)
    assert (tmp_path / "junction-1.csv").read_bytes() == (tmp_path / "junction-2.csv").read_bytes()
    output = pl.read_csv(tmp_path / "junction-1.csv")
    assert output["junctionb"].to_list() == ["CAW", "CAW", "CAF"]
    assert output["score.tra"].to_list() == [15., 15., 19.]
    assert output["score.trb"].to_list() == [16., 16., 12.]
    assert output["score"].to_list() == [31., 31., 31.]
    assert structures[0].by_type("TRB")[0].regions[-1].residues[0].aa == "W"


@pytest.mark.parametrize("defect", ["missing", "ambiguous", "gap", "wrong-anchor"])
def test_junction_requires_actual_adjacent_anchors(tmp_path, monkeypatch, defect):
    structure = _junction("first")
    alpha = structure.by_type("TRA")[0]
    if defect == "missing":
        alpha.regions = alpha.regions[1:]
    elif defect == "ambiguous":
        alpha.regions.append(alpha.regions[0])
    elif defect == "gap":
        alpha.regions[0].end_seq_index = 0
    else:
        alpha.regions[0].residues[0] = replace(alpha.regions[0].residues[0], aa="A")
    monkeypatch.setattr(cli, "_typed_batch", lambda *a, **k: [structure])
    candidates, output = tmp_path / "junctions.tsv", tmp_path / "output.csv"
    candidates.write_text("complex.id\tcandidate.id\tjunctiona\tjunctionb\nfirst\tnative\tCKF\tCAW\n")
    result = CliRunner().invoke(cli.app, ["score-receptors", "-s", "unused.pdb", "-c",
        str(candidates), "--receptor-region", "junction", "-t", "1", "-o", str(output)])
    assert result.exit_code != 0 and "anchors" in str(result.exception)
    assert not output.exists()


@pytest.mark.parametrize("mode,columns,alpha,error", [
    ("junction", "cdr3a\tcdr3b", "K", "junctiona,junctionb"),
    ("cdr3", "junctiona\tjunctionb", "CKF", "cdr3a,cdr3b"),
    ("junction", "junctiona\tjunctionb", "AKF", "start with C"),
    ("junction", "junctiona\tjunctionb", "CKA", "end with F or W"),
    ("junction", "junctiona\tjunctionb", "CKKF", "length 3"),
])
def test_junction_mode_specific_input(tmp_path, monkeypatch, mode, columns, alpha, error):
    monkeypatch.setattr(cli, "_typed_batch", lambda *a, **k: [_junction("first")])
    candidates = tmp_path / "junctions.tsv"
    candidates.write_text(f"complex.id\tcandidate.id\t{columns}\nfirst\tnative\t{alpha}\tCAW\n")
    result = CliRunner().invoke(cli.app, ["score-receptors", "-s", "unused.pdb", "-c",
        str(candidates), "--receptor-region", mode, "-t", "1"])
    assert result.exit_code != 0 and error in result.output + str(result.exception)
