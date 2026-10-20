"""Run an authorized contact job with a sampled process-tree RSS/time guard.

Usage: python scripts/run_bounded.py REPORT.json SECONDS GIB -- COMMAND ...
The guard is sampled, not an OS allocation limit. All outputs must be new paths.
"""
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time

report = Path(sys.argv[1])
if report.exists():
    raise FileExistsError(report)
seconds, cap = float(sys.argv[2]), float(sys.argv[3]) * 1024**3
command = sys.argv[5:]
assert sys.argv[4] == "--" and command and seconds > 0 and cap > 0
env = os.environ.copy()
for name in ("POLARS_MAX_THREADS", "POLARS_ASYNC_THREAD_COUNT", "OPENBLAS_NUM_THREADS",
             "OMP_NUM_THREADS", "OMP_THREAD_LIMIT", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS",
             "VECLIB_MAXIMUM_THREADS", "RAYON_NUM_THREADS"):
    env[name] = "1"
env["POLARS_MAX_BLOCKING_THREADS"] = "4"
start, peak, reason = time.monotonic(), 0, None
child = subprocess.Popen(command, env=env, start_new_session=True)

def group_alive():
    try:
        os.killpg(child.pid, 0)
        return True
    except ProcessLookupError:
        return False

try:
    while True:
        status = child.poll()
        raw = subprocess.run(["/bin/ps", "-axo", "pid=,ppid=,pgid=,rss="], capture_output=True, text=True, check=True)
        processes = [tuple(map(int, row.split())) for row in raw.stdout.splitlines()]
        members = {pid for pid, _, pgid, _ in processes if pgid == child.pid}
        if status is not None:
            if members:
                reason = "surviving_descendants"
            break
        selected = members | {os.getpid()}
        while True:
            expanded = selected | {pid for pid, parent, _, _ in processes if parent in selected}
            if expanded == selected:
                break
            selected = expanded
        rss = sum(kib * 1024 for pid, _, _, kib in processes if pid in selected)
        peak = max(peak, rss)
        if rss > cap or time.monotonic() - start > seconds:
            reason = "rss" if rss > cap else "timeout"
            break
        time.sleep(.1)
finally:
    if group_alive():
        try:
            os.killpg(child.pid, signal.SIGTERM)
        except ProcessLookupError:
            pass
        deadline = time.monotonic() + 5
        while group_alive() and time.monotonic() < deadline:
            child.poll()  # Reap the root even while surviving children ignore SIGTERM.
            time.sleep(.05)
        if group_alive():
            try:
                os.killpg(child.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
    child.wait(timeout=5)
    result = dict(command=command, exit_code=child.returncode, guard_stop=reason,
                  peak_sampled_tree_rss_bytes=peak, rss_cap_bytes=cap,
                  timeout_s=seconds, elapsed_s=time.monotonic()-start, poll_s=.1,
                  compute_threads=1, environment_limits={k: env[k] for k in env if k in
                      ("POLARS_MAX_THREADS", "POLARS_MAX_BLOCKING_THREADS", "OMP_NUM_THREADS")})
    report.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result), flush=True)
raise SystemExit(124 if reason else child.returncode)
