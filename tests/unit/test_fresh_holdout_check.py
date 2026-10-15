"""The full reproduction gate must reject a reduced hold-out before fitting."""
import subprocess
import sys
from pathlib import Path

import polars as pl


def test_fresh_holdout_check_rejects_missing_manifest_rows(tmp_path):
    features = tmp_path / "features.parquet"
    pl.DataFrame({"id": ["not-a-manifest-id"]}).write_parquet(features)
    script = Path(__file__).resolve().parents[2] / "scripts/check_fresh_holdout.py"
    p = subprocess.run([sys.executable, str(script), "--features", str(features),
                        "--out", str(tmp_path / "report.json")],
                       capture_output=True, text=True, timeout=30)
    assert p.returncode != 0
    assert "missing hold-out structures" in p.stderr
