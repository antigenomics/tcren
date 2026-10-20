"""The guard exports real native-pool limits without claiming measured concurrency."""
import json
from pathlib import Path
import subprocess
import sys


def test_bounded_guard_thread_configuration(tmp_path):
    report = tmp_path / "guard.json"
    guard = Path(__file__).resolve().parents[2] / "scripts/run_bounded.py"
    child = "import os; assert os.environ['POLARS_MAX_BLOCKING_THREAD_COUNT']=='4'; assert os.environ['OMP_NUM_THREADS']=='1'"
    subprocess.run([sys.executable, str(guard), str(report), "10", "1", "--",
                    sys.executable, "-c", child], check=True, capture_output=True, timeout=15)
    record = json.loads(report.read_text())
    assert record["exit_code"] == 0 and record["guard_stop"] is None
    assert record["configured_native_thread_limit"] == 1
    assert "compute_threads" not in record
    assert record["environment_limits"]["POLARS_MAX_BLOCKING_THREAD_COUNT"] == "4"
