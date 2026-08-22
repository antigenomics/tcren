"""OBSOLETE.md must list exactly what is marked in the source, or it is worse than no list."""
from __future__ import annotations

import re
from pathlib import Path

from tcren._provenance import marked, not_in_tcren2

OBSOLETE = Path(__file__).resolve().parents[2] / "OBSOLETE.md"


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
