"""Frozen native-reference shape, integrity and inspection."""

from hashlib import sha256
import json
from pathlib import Path

import pytest
from typer.testing import CliRunner

from tcren.cli import app
from tcren.native_reference import contact_probability, load_native_reference


def test_native_reference_contains_recorded_counts_and_scales():
    ref = load_native_reference()
    assert ref["version"] == "native387-v1"
    assert ref["native_structures"] == 387
    assert ref["size"] == pytest.approx({"n": 386, "centre": 174.0, "spread": 59.304})
    assert len(ref["energy_by_peptide"]) == 239
    assert sum(group["n"] for group in ref["energy_by_peptide"].values()) == 387
    assert sum(group["spread"] is None for group in ref["energy_by_peptide"].values()) == 168
    assert len(ref["geometry"]) == 18
    assert all(term["n"] == 387 for term in ref["geometry"].values())
    profiles = ref["contact_profiles"]
    assert {key: len(value["probability"]) for key, value in profiles.items()} == {
        "atom/complex": 3843, "atom/tcr_peptide": 1285,
        "residue/complex": 1178, "residue/tcr_peptide": 393,
    }
    q = profiles["atom/tcr_peptide"]["probability"]
    assert sum(q.values()) == pytest.approx(1.0)
    assert min(q.values()) > 0
    assert contact_probability(ref, "atom", "tcr_peptide", "CA|CA") == q["CA|CA"]
    assert contact_probability(ref, "atom", "tcr_peptide", "unseen") is None


def test_native_reference_rejects_a_changed_resource(tmp_path):
    path = tmp_path / "reference.json"
    raw = b'{"format":"wrong"}'
    path.write_bytes(raw)
    with pytest.raises(ValueError, match="SHA-256"):
        load_native_reference(path)
    with pytest.raises(ValueError, match="format"):
        load_native_reference(path, expected_sha256=sha256(raw).hexdigest())
    with pytest.raises(FileNotFoundError):
        load_native_reference(Path(tmp_path, "missing.json"))


def test_cli_reference_show_agrees_with_installed_resource():
    result = CliRunner().invoke(app, ["reference", "--show"])
    assert result.exit_code == 0, result.output
    assert json.loads(result.stdout) == load_native_reference()
