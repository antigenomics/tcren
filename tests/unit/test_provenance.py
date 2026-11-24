"""OBSOLETE.md must list exactly what is marked in the source, or it is worse than no list."""
from __future__ import annotations

import re
from pathlib import Path

from tcren._provenance import marked, not_in_tcren2

OBSOLETE = Path(__file__).resolve().parents[2] / "docs/OBSOLETE.md"


def test_obsolete_md_matches_the_markings():
    listed = set(re.findall(r"^- `([\w.]+)`", OBSOLETE.read_text(), flags=re.M))
    assert listed == {name for name, _ in marked()}, (
        "OBSOLETE.md is stale -- regenerate with `python -m tcren._provenance`"
    )


def test_marker_records_the_reason_and_keeps_the_docstring():
    @not_in_tcren2("measured and not adopted")
    def f():
        """Do a thing."""

    assert f.__not_in_tcren2__ == "measured and not adopted"
    assert "measured and not adopted" in f.__doc__
    assert "Do a thing." in f.__doc__


def test_additive_contact_count_preserves_historical_stamps(tmp_path):
    import json

    import pytest

    from tcren.provenance import StaleTableError, check, registry_digest, sidecar_path

    table = tmp_path / "historic.csv"
    old = "222b4733759b81c4b26dd8024ae3517732a4875b4a83809168b60eadde96285b"
    assert registry_digest() != old
    metadata = {"tcren": "4.9.0", "registry_digest": old, "columns": ["n_contacts_tp"]}
    sidecar = sidecar_path(table)
    sidecar.write_text(json.dumps(metadata))
    assert check(table) == metadata
    metadata["columns"].append("n_atom_contacts_tp")
    sidecar.write_text(json.dumps(metadata))
    with pytest.raises(StaleTableError, match="different descriptor catalogue"):
        check(table)
