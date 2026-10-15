"""Exact exports preserve existing quantities, physical identity and failed inputs."""
import json
import tarfile
from pathlib import Path

import polars as pl
import pytest
from polars.testing import assert_frame_equal
from typer.testing import CliRunner

from tcren import exports
from tcren.cli import app
from tcren.contactmap import ContactMap
from tcren.potential import mj, tcren as classic
from tcren.structure import parse_structure
from tcren.structure.model import RegionMarkup

PDB = Path(__file__).resolve().parents[1] / "assets/pdb/1ao7.pdb"


def typed_structure():
    structure = parse_structure(PDB)
    for chain in structure.chains:
        chain.chain_type = {"A": "MHCa", "B": "B2M", "C": "PEPTIDE",
                            "D": "TRA", "E": "TRB"}[chain.chain_id]
        chain.chain_supertype = "MHCI" if chain.chain_type in {"MHCa", "B2M"} else chain.chain_type
        if chain.chain_type == "PEPTIDE":
            chain.regions = [RegionMarkup("PEPTIDE", 0, len(chain.residues)-1, chain.sequence(), chain.residues)]
        if chain.chain_type in {"TRA", "TRB"}:
            residues = chain.residues
            chain.regions = [RegionMarkup("CDR3", 0, len(residues)-1, chain.sequence(), residues)]
    return structure


def test_exact_export_quantities_and_author_identity():
    from tcren.energetics.scoring import position_profile
    from tcren.pipeline import _interface_energy
    from tcren.viz.pymol import residue_importance

    structure = typed_structure()
    pot, other = classic(), mj()
    cm = ContactMap.from_structure(structure)
    before = cm.contacts.clone()
    row = exports.receptor_reference(cm, pot, other)
    for interface, suffix, potential in (("tcr_peptide", "tcr_pep", pot), ("tcr_mhc", "tcr_mhc", other)):
        contacts = cm.interface(interface)
        expected = _interface_energy(contacts, potential) - _interface_energy(
            contacts.with_columns(pl.lit("A").alias("residue.aa.from")), potential)
        assert row[f"dFT_{suffix}_all"] == expected
        assert row[f"dFT_{suffix}_cdr3"] == expected
        assert row[f"dFT_{suffix}_cdr12"] == 0
    for name in exports.RECEPTOR_REGIONS:
        assert row[f"dFT_total_{name}"] == row[f"dFT_tcr_pep_{name}"] + row[f"dFT_tcr_mhc_{name}"]
    assert_frame_equal(cm.contacts, before)
    for quantity, expected in (("position", position_profile(cm, pot, structure=structure)),
                               ("residue", residue_importance(structure, potential=pot))):
        actual = exports.export_structure(structure, quantity, pot, other)
        assert_frame_equal(actual.select(expected.columns).sort(expected.columns),
                           expected.sort(expected.columns))
        assert actual["residue.index.pdb"].null_count() == 0
        assert "insertion_code" in actual.columns
    identity = exports.residue_identity(structure)
    assert identity.height == sum(len(c.residues) for c in structure.chains)
    placement = exports.export_structure(structure, "placement", pot, other)
    assert placement.height == 1
    assert placement["tcr_low"][0] <= placement["tcr_high"][0]


@pytest.mark.parametrize("archived", [False, True])
def test_batch_accounts_for_parse_failure_and_zero_contact_rows(tmp_path, monkeypatch, archived):
    import tcren.annotation
    import tcren.annotation.batch
    import tcren.mhc
    import tcren.structure

    good = tmp_path / "good_full_stem.pdb"
    good.write_bytes(PDB.read_bytes())
    bad = tmp_path / "missing.pdb"
    calls = []
    original_parse = tcren.structure.parse_structure

    def parse(path, pdb_id):
        if pdb_id == "missing":
            raise ValueError("invalid fixture")
        original_parse(path, pdb_id=pdb_id)
        structure = typed_structure()
        structure.pdb_id = pdb_id
        return structure

    monkeypatch.setattr(tcren.structure, "parse_structure", parse)
    monkeypatch.setattr(tcren.annotation.batch, "annotate_batch",
                        lambda structures, **kw: calls.append(len(structures)) or [{} for s in structures])
    monkeypatch.setattr(tcren.annotation, "classify_chains", lambda *a, **kw: None)
    monkeypatch.setattr(tcren.mhc, "annotate_mhc_batch", lambda *a, **kw: None)
    monkeypatch.setattr(exports, "annotation_lineage", lambda organisms: {"organisms": list(organisms)})
    out = tmp_path / "reference.tsv"
    sources = [good, bad]
    if archived:
        bad.write_text("broken structure")
        archive = tmp_path / "inputs.tar.gz"
        with tarfile.open(archive, "w:gz") as stream:
            for path in sources:
                stream.add(path, arcname=path.name)
        sources = [archive]
    ledger = exports.export_batch(sources, out, quantity="contacts", peptide_potential=classic(),
                                   mhc_potential=mj(), cutoff=0.00001)
    assert ledger["status"].to_list() == ["ok", "error"]
    assert ledger["n_rows"].to_list() == [0, 0]
    assert calls == [1]
    assert ledger["pdb.id"].to_list() == ["good_full_stem", "missing"]
    if archived:
        assert ledger["archive_sha256"].null_count() == 0
        assert ledger["source"].str.contains("inputs.tar.gz#").all()
    meta = json.loads(Path(str(out) + ".provenance.json").read_text())
    assert (meta["n_input"], meta["n_ok"], meta["n_error"]) == (2, 1, 1)
    assert meta["potentials"]["peptide"]["sha256"]
    assert Path(str(out) + ".residues.tsv").exists()


def test_all_failed_cli_retains_accounting(tmp_path):
    out = tmp_path / "failed.tsv"
    annotation_path = Path(str(out) + ".annotation.json")
    annotation_path.write_text('{"old": "stale"}')
    result = CliRunner().invoke(app, ["export", "-s", str(tmp_path / "absent.pdb"),
                                    "-o", str(out), "--potential", "karnaukhov2022"])
    assert result.exit_code == 1, result.output
    ledger = pl.read_csv(str(out) + ".accounting.tsv", separator="\t")
    assert ledger["status"].to_list() == ["error"]
    assert "FileNotFoundError" in ledger["error"][0]
    assert json.loads(annotation_path.read_text()) == {}


def test_export_rejects_invalid_cutoff():
    with pytest.raises(ValueError, match="cutoff"):
        exports.export_structure(typed_structure(), "contacts", classic(), mj(), cutoff=float("nan"))
