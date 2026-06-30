"""CLI regression for ``tcren rank`` without ``-c/--candidates`` (slow — arda/mmseqs).

Guards the no-candidates default path: when no candidate peptides are passed, ``rank``
must default to each structure's *native* peptide. Regression for the v2.0.0 bug where
the default-peptide branch referenced ``c.sequence`` (the bound method) instead of calling
``c.sequence()``, crashing with ``TypeError: object of type 'method' has no len()``.
"""

from __future__ import annotations

from pathlib import Path

import polars as pl
import pytest

pytestmark = pytest.mark.slow  # invokes arda / mmseqs to classify chains

ASSET = Path(__file__).resolve().parents[1] / "assets" / "pdb" / "1ao7.pdb"


def test_rank_defaults_to_native_peptide(tmp_path):
    from typer.testing import CliRunner

    from tcren.cli import app

    out = tmp_path / "rank.csv"
    result = CliRunner().invoke(app, ["rank", "-s", str(ASSET), "-o", str(out)])
    assert result.exit_code == 0, result.output

    df = pl.read_csv(out)
    # No candidates passed -> exactly the native peptide is ranked.
    assert df.height == 1
    assert df["complex.id"][0] == "1ao7"
    assert df["peptide"][0] == "LLFGYPVYV"  # 1ao7 native (HTLV-1 Tax)
    assert 0.0 <= df["rank_pct"][0] <= 1.0
