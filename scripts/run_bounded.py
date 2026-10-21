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
env["POLARS_MAX_BLOCKING_THREAD_COUNT"] = "4"
start, peak, reason = time.monotonic(), 0, None
child = subprocess.Popen(command, env=env, start_new_session=True)

# Retain sampled descendants across reparenting and separate sessions. Start times
# prevent a recycled PID from acquiring ownership during a later sample.
owned = {}
cleanup_errors = []


def sample_owned():
    raw = subprocess.run(["/bin/ps", "-axo", "pid=,ppid=,pgid=,rss=,stat=,lstart="],
                         capture_output=True, text=True, check=True)
    processes = {}
    for row in raw.stdout.splitlines():
        pid, parent, group, rss, state, started = row.split(maxsplit=5)
        processes[int(pid)] = (int(parent), int(group), int(rss), state, started)
    selected = {pid for pid, info in processes.items()
                if owned.get(pid) == info[4] or info[1] == child.pid}
    while True:
        expanded = selected | {pid for pid, info in processes.items() if info[0] in selected}
        if expanded == selected:
            break
        selected = expanded
    for pid in selected:
        owned[pid] = processes[pid][4]
    live = {pid for pid in selected if not processes[pid][3].startswith("Z")}
    return live, sum(processes[pid][2] * 1024 for pid in live)


try:
    while True:
        status = child.poll()
        members, rss = sample_owned()
        peak = max(peak, rss)
        if status is not None:
            if members:
                reason = "surviving_descendants"
            break
        if rss > cap or time.monotonic() - start > seconds:
            reason = "rss" if rss > cap else "timeout"
            break
        time.sleep(.1)
finally:
    survivors = set()
    try:
        for sig in (signal.SIGTERM, signal.SIGKILL):
            deadline = time.monotonic() + 5
            sent = set()
            while True:
                child.poll()
                survivors, _ = sample_owned()
                if not survivors:
                    break
                for pid in survivors - sent:
                    try:
                        os.kill(pid, sig)
                    except ProcessLookupError:
                        pass
                    except OSError as error:
                        cleanup_errors.append(f"pid {pid}: {error}")
                    sent.add(pid)
                if time.monotonic() >= deadline:
                    break
                time.sleep(.05)
            if not survivors:
                break
        child.wait(timeout=5)
    except (OSError, subprocess.SubprocessError) as error:
        cleanup_errors.append(str(error))
    if cleanup_errors or survivors:
        reason = reason or "cleanup_failed"
    result = dict(command=command, exit_code=child.returncode, guard_stop=reason,
                  peak_sampled_tree_rss_bytes=peak, rss_cap_bytes=cap,
                  timeout_s=seconds, elapsed_s=time.monotonic()-start, poll_s=.1,
                  cleanup_errors=cleanup_errors, surviving_sampled_pids=sorted(survivors),
                  configured_native_thread_limit=1, environment_limits={k: env[k] for k in env if k in
                      ("POLARS_MAX_THREADS", "POLARS_MAX_BLOCKING_THREAD_COUNT", "OMP_NUM_THREADS")})
    report.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result), flush=True)
raise SystemExit(124 if reason else child.returncode)
