"""Installed MD exports preserve the existing API quantities and explicit input contract."""
import hashlib
import json

import polars as pl
import pytest
from typer.testing import CliRunner

from tcren.cli import app
from tcren.validation.md import inspect_md_trajectory


def inputs(tmp_path):
    text = "".join(
        f"ATOM  {i:5d}  CA  ALA {chain}{j:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           C\n"
        for i, (chain, j, x, y, z) in enumerate([
            ("A", 1, 0, 0, 0), ("A", 2, 1, 0, 0), ("A", 3, 0, 1, 0),
            ("B", 1, 0, 0, 2), ("B", 2, 1, 0, 2), ("B", 3, 0, 1, 2),
            ("C", 1, 0, 0, 4), ("C", 2, 1, 0, 4), ("C", 3, 0, 1, 4),
        ], 1)
    )
    paths = {name: tmp_path / name for name in ["reference.pdb", "trajectory.pdb", "potential.xvg", "components.json"]}
    paths["reference.pdb"].write_text(text + "END\n")
    paths["trajectory.pdb"].write_text("".join(
        f"TITLE     t= {time}\nMODEL        {i}\n{text}ENDMDL\n"
        for i, time in enumerate([0, 10], 1)
    ))
    paths["potential.xvg"].write_text("0 -10\n10 -20\n")
    paths["components.json"].write_text(json.dumps({"TCR": ["A"], "MHC": ["B"], "peptide": ["C"]}))
    args = ["md", "inspect", "--trajectory", str(paths["trajectory.pdb"]),
            "--reference", str(paths["reference.pdb"]), "--potential", str(paths["potential.xvg"]),
            "--components", str(paths["components.json"]), "--mhc-class", "I"]
    return paths, args


def test_md_cli_exports_exact_api_and_source_identity(tmp_path):
    paths, args = inputs(tmp_path)
    output = tmp_path / "export"
    result = CliRunner().invoke(app, [*args, "--pbc-whole", "-o", str(output)])
    assert result.exit_code == 0, result.output
    expected = inspect_md_trajectory(paths["trajectory.pdb"], paths["reference.pdb"],
                                    {"TCR": ["A"], "MHC": ["B"], "peptide": ["C"]},
                                    paths["potential.xvg"], contact_cutoff=5)
    assert pl.read_csv(output / "rmsd.tsv", separator="\t").to_dicts() == expected["rmsd"]
    assert json.loads((output / "minimum-potential-frame.json").read_text()) == expected["minimum_potential"]
    assert (output / "minimum-potential-frame.pdb").read_text() == expected["selected_frame_pdb"]
    record = json.loads((output / "provenance.json").read_text())
    assert record["state"] == "complete" and record["mhc_class"] == "I"
    assert record["sources"]["trajectory"]["sha256"] == hashlib.sha256(paths["trajectory.pdb"].read_bytes()).hexdigest()
    assert record["cpu_budget"] == 1
    assert record["annotation_provenance"] == {}
    identities = pl.read_csv(output / "backbone_atoms.tsv", separator="\t")
    assert identities.height == 9 and identities["component"].unique().len() == 3
    assert (output / "docking.tsv").read_text().startswith("frame_index\ttime_ps")
    assert CliRunner().invoke(app, [*args, "--pbc-whole", "-o", str(output)]).exit_code != 0


def test_md_cli_requires_pbc_and_records_failure(tmp_path):
    paths, args = inputs(tmp_path)
    output = tmp_path / "export"
    assert CliRunner().invoke(app, [*args, "-o", str(output)]).exit_code != 0
    assert not output.exists()
    paths["potential.xvg"].write_text("0 -10\n9 -20\n")
    assert CliRunner().invoke(app, [*args, "--pbc-whole", "-o", str(output)]).exit_code != 0
    record = json.loads((output / "provenance.json").read_text())
    assert record["state"] == "failed" and "timestamp" in record["error"]


def test_md_quality_exports_existing_atom_pairs(tmp_path):
    paths, _ = inputs(tmp_path)
    output = tmp_path / "quality"
    result = CliRunner().invoke(app, ["md", "quality", "-s", str(paths["reference.pdb"]),
                                     "--mhc-class", "II", "--cutoff", "2.1", "-o", str(output)])
    assert result.exit_code == 0, result.output
    from tcren.validation.structures import close_atom_pairs
    assert pl.read_csv(output / "close_atom_pairs.tsv", separator="\t").to_dicts() == close_atom_pairs(paths["reference.pdb"], cutoff=2.1)
    assert json.loads((output / "provenance.json").read_text())["mhc_class"] == "II"


@pytest.mark.parametrize("mapping", [{}, {"all": ["A"]}, {"TCR": "A"},
                                    {"TCR": ["A", "A"]}, {"TCR": ["A"], "MHC": ["A"]}])
def test_md_components_reject_ambiguous_mapping(tmp_path, mapping):
    paths, _ = inputs(tmp_path)
    with pytest.raises(ValueError, match="[Cc]omponent"):
        inspect_md_trajectory(paths["trajectory.pdb"], paths["reference.pdb"], mapping,
                              paths["potential.xvg"])


def test_md_atom_correspondence_includes_author_identity(tmp_path):
    paths, _ = inputs(tmp_path)
    paths["trajectory.pdb"].write_text(paths["trajectory.pdb"].read_text().replace("ALA A   1", "ALA A 101"))
    with pytest.raises(ValueError, match="atom identities"):
        inspect_md_trajectory(paths["trajectory.pdb"], paths["reference.pdb"],
                              {"TCR": ["A"]}, paths["potential.xvg"])


def test_md_budget_validated_before_parsing(tmp_path, monkeypatch):
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "1")
    with pytest.raises(ValueError, match="exceeds"):
        inspect_md_trajectory("missing", "missing", {}, "missing", threads=2)


def test_md_cli_compares_same_coordinate_reruns(tmp_path):
    paths, _ = inputs(tmp_path)
    for name, text in {"ref-energy.xvg": "0 -100\n", "energy.xvg": "0 -99\n",
                       "ref-forces.xvg": "0 3 4 0\n", "forces.xvg": "0 3 4 1\n"}.items():
        (tmp_path / name).write_text(text)
    args = ["md", "compare-forces", "--reference-energy", str(tmp_path / "ref-energy.xvg"),
            "--reference-forces", str(tmp_path / "ref-forces.xvg"), "--energy", str(tmp_path / "energy.xvg"),
            "--forces", str(tmp_path / "forces.xvg"), "--atoms", "1", "--mhc-class", "I",
            "--topology", str(paths["components.json"]), "--coordinates", str(paths["reference.pdb"])]
    output = tmp_path / "forces"
    assert CliRunner().invoke(app, [*args, "-o", str(output)]).exit_code != 0
    result = CliRunner().invoke(app, [*args, "--same-system", "-o", str(output)])
    assert result.exit_code == 0, result.output
    record = json.loads((output / "comparison.json").read_text())
    assert record["force_error_relative_rms"] == pytest.approx(0.2)
    assert record["potential_difference_kJ_mol"] == 1
