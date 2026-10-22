"""The guard exports real native-pool limits without claiming measured concurrency."""
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time


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


def test_bounded_guard_reaps_sampled_separate_session(tmp_path):
    report = tmp_path / "guard.json"
    guard = Path(__file__).resolve().parents[2] / "scripts/run_bounded.py"
    # Keep the parent alive long enough for the sampled guard to observe ownership.
    child = ("import subprocess,sys,time; "
             "subprocess.Popen([sys.executable,'-c','import time; time.sleep(30)'], "
             "start_new_session=True); time.sleep(1)")
    result = subprocess.run([sys.executable, str(guard), str(report), "10", "1", "--",
                             sys.executable, "-c", child], capture_output=True, timeout=25)
    assert result.returncode == 124, result.stderr
    record = json.loads(report.read_text())
    assert record["guard_stop"] == "surviving_descendants"
    assert record["surviving_sampled_pids"] == []
    assert record["cleanup_errors"] == []


def test_sigterm_cleans_up_children_and_writes_cancellation_report(tmp_path):
    report, ready = tmp_path / "guard.json", tmp_path / "ready.json"
    guard = Path(__file__).resolve().parents[2] / "scripts/run_bounded.py"
    child = (
        "import json,os,subprocess,sys,time; from pathlib import Path; "
        "p=subprocess.Popen([sys.executable,'-c','import time; time.sleep(60)'],"
        "start_new_session=True); "
        f"Path({str(ready)!r}).write_text(json.dumps([os.getpid(),p.pid])); time.sleep(60)"
    )
    process = subprocess.Popen([sys.executable, str(guard), str(report), "30", "1", "--",
                                sys.executable, "-c", child],
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    pids = []
    try:
        deadline = time.monotonic() + 10
        while not ready.exists() and process.poll() is None and time.monotonic() < deadline:
            time.sleep(.02)
        assert ready.exists(), process.communicate(timeout=5)
        pids = json.loads(ready.read_text())
        time.sleep(.3)  # let the sampled guard observe the separately sessioned descendant
        process.send_signal(signal.SIGTERM)
        time.sleep(.01)
        if process.poll() is None:
            process.send_signal(signal.SIGTERM)
        _, error = process.communicate(timeout=12)
        assert process.returncode == 124, error
        record = json.loads(report.read_text())
        assert record["guard_stop"] == "cancelled"
        assert record["cancellation_signal"] == signal.SIGTERM
        assert record["cleanup_errors"] == [] and record["surviving_sampled_pids"] == []
        live = subprocess.run(["/bin/ps", "-o", "stat=", "-p", ",".join(map(str, pids))],
                              capture_output=True, text=True, timeout=5)
        assert all(state.strip().startswith("Z") for state in live.stdout.splitlines())
    finally:
        if process.poll() is None:
            process.kill()
            process.communicate(timeout=5)
        for pid in pids:
            try:
                os.kill(pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
