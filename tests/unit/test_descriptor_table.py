"""The catalogue documents itself, and the docs table cannot drift from it."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

from tcren.recognition import DESCRIPTORS, DETAIL, INVARIANCE

REPO = Path(__file__).resolve().parents[2]


def test_every_descriptor_carries_units_and_a_definition():
    """A descriptor with no DETAIL entry would reach a feature table undocumented."""
    assert set(DETAIL) == set(DESCRIPTORS)
    for name, (units, definition) in DETAIL.items():
        assert units, f"{name} has no units"
        assert len(definition) > 20, f"{name}'s definition is too short to be one"
        assert definition.endswith("."), f"{name}'s definition is not a sentence"


def test_units_come_from_a_closed_vocabulary():
    """Units are what a transform has to respect, so they are a fixed set, not free text."""
    allowed = {
        "A", "A^2", "A^3", "deg", "rad", "cosine", "kT", "kT/site", "N", "N/m", "J",
        "count", "fraction", "signed fraction", "ratio", "log-odds", "log-odds^2",
        "class I / II",
    }
    seen = {u for u, _ in DETAIL.values()}
    assert seen <= allowed, f"unknown units: {sorted(seen - allowed)}"


def test_every_descriptor_is_classified():
    assert set(INVARIANCE) == set(DESCRIPTORS)


def test_the_generated_docs_table_is_current():
    """`python scripts/gen_descriptor_table.py` brings it back into line."""
    result = subprocess.run(
        [sys.executable, str(REPO / "scripts" / "gen_descriptor_table.py"), "--check"],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr


def test_the_generated_family_graph_is_current():
    """`python scripts/gen_family_graph.py` brings it back into line."""
    result = subprocess.run(
        [sys.executable, str(REPO / "scripts" / "gen_family_graph.py"), "--check"],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr
