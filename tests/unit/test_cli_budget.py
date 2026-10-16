"""CLI failures and annotation/native workers must obey one CPU budget."""
from pathlib import Path

import polars as pl
from typer.testing import CliRunner

from tcren.cli import app

FIXTURE = Path(__file__).resolve().parents[1] / 'assets/pdb/1ao7.pdb'


def test_surface_grid_is_checked_before_annotation(tmp_path, monkeypatch):
    import tcren.annotation.batch as batch
    monkeypatch.setattr(batch, 'annotate_batch', lambda *a, **k: (_ for _ in ()).throw(AssertionError('searched')))
    result = CliRunner().invoke(app, ['surface', '-s', str(FIXTURE), '--grid', '0x0', '-o', str(tmp_path/'bad.csv')])
    assert result.exit_code == 2, result.output
    assert 'positive' in result.output


def test_surface_batch_keeps_parse_failures_and_budget(tmp_path, monkeypatch):
    import tcren.annotation as annotation
    import tcren.annotation.batch as batch
    import tcren.mhc as mhc
    from tcren.topology import surface
    import numpy as np

    calls = []
    monkeypatch.setattr(batch, 'annotate_batch', lambda items, **kw:
                        calls.append(('tcr', len(items), kw['threads'])) or [{} for _ in items])
    monkeypatch.setattr(annotation, 'classify_chains', lambda *a, **kw: None)
    monkeypatch.setattr(mhc, 'annotate_mhc_batch', lambda items, **kw: calls.append(('mhc', len(items), kw['threads'])))
    def map_one(s, **kw):
        return surface.SurfaceMap(s.pdb_id, (2, 2), surface.DEFAULT_EXTENT,
                                  {'h': np.ones((2, 2)), 'phobic': np.zeros((2, 2)), 'charge': np.zeros((2, 2))},
                                  np.ones((2, 2), dtype=int))
    monkeypatch.setattr(surface, 'surface_map', map_one)
    for stem in ('one', 'two'):
        (tmp_path/f'{stem}.pdb').write_bytes(FIXTURE.read_bytes())
    (tmp_path/'broken.pdb').write_text('invalid')
    out = tmp_path/'surface.csv'
    result = CliRunner().invoke(app, ['surface', '-s', str(tmp_path), '-t', '1', '-o', str(out)])
    assert result.exit_code == 0, result.output
    assert calls == [('tcr', 2, 1), ('mhc', 2, 1)]
    table = pl.read_csv(out)
    assert set(table['structure.id']) == {'one', 'two', 'broken'}
    assert table.filter(pl.col('structure.id')=='broken')['error'][0]
    result = CliRunner().invoke(app, ['surface', '-s', str(tmp_path/'broken.pdb'), '-o', str(out)])
    assert result.exit_code == 1, result.output
    assert pl.read_csv(out).height == 1
    empty = tmp_path/'empty'; empty.mkdir()
    result = CliRunner().invoke(app, ['surface', '-s', str(empty), '-o', str(out)])
    assert result.exit_code == 1, result.output
    assert pl.read_csv(out).height == 0


def test_cpu_modes_and_initialized_native_limits(monkeypatch):
    import os
    import pytest
    from threadpoolctl import threadpool_info, threadpool_limits
    from tcren import runtime

    monkeypatch.setattr(runtime, 'available_cpus', lambda: 2)
    assert [runtime.resolve_threads(n) for n in (0, 1, 2)] == [2, 1, 2]
    for invalid in (-1, 3, True, 1.5):
        with pytest.raises(ValueError):
            runtime.resolve_threads(invalid)
    monkeypatch.setenv('OMP_NUM_THREADS', '7')
    before = dict(os.environ)
    with runtime.native_threads():
        assert all(os.environ[k] == ('4' if k == 'POLARS_MAX_BLOCKING_THREAD_COUNT' else '1')
                   for k in runtime._NATIVE_ENV)
        assert all(p['num_threads'] == 1 for p in threadpool_info())
        assert pl.thread_pool_size() == 1
    assert dict(os.environ) == before
    if threadpool_info():
        with threadpool_limits(limits=2):
            assert all(p['num_threads'] == 2 for p in threadpool_info())
            with runtime.native_threads():
                assert all(p['num_threads'] == 1 for p in threadpool_info())
            assert all(p['num_threads'] == 2 for p in threadpool_info())
    with pytest.raises(RuntimeError, match='probe'):
        with runtime.native_threads():
            raise RuntimeError('probe')
    assert dict(os.environ) == before


def test_preinitialized_polars_is_rejected_with_actionable_error():
    import os
    import subprocess
    import sys

    code = '''import polars as pl
assert pl.thread_pool_size() == 4
from tcren.runtime import native_threads
try:
    with native_threads(): pass
except RuntimeError as e:
    assert 'fresh process' in str(e)
else:
    raise AssertionError('unbounded pool accepted')
'''
    result = subprocess.run([sys.executable, '-c', code], text=True, capture_output=True,
                            env={**os.environ, 'POLARS_MAX_THREADS': '4'}, timeout=30)
    assert result.returncode == 0, result.stdout+result.stderr


def test_surface_annotation_failure_is_not_success(tmp_path, monkeypatch):
    import tcren.annotation.batch as batch
    def failed(*args, **kwargs):
        raise RuntimeError('reference unavailable')
    monkeypatch.setattr(batch, 'annotate_batch', failed)
    out = tmp_path/'errors.csv'
    result = CliRunner().invoke(app, ['surface', '-s', str(FIXTURE), '-o', str(out)])
    assert result.exit_code == 1, result.output
    assert pl.read_csv(out)['error'][0] == 'annotation: RuntimeError: reference unavailable'


def test_shared_cli_annotation_passes_one_budget(monkeypatch):
    import tcren.cli as cli
    import tcren.mhc as mhc
    calls=[]
    monkeypatch.setattr(cli, 'resolve_threads', lambda n: 2)
    monkeypatch.setattr(cli, '_iter_typed', lambda sources, organism, **kw:
                        calls.append(('tcr', sources, kw['threads'])) or iter([object(), object()]))
    monkeypatch.setattr(mhc, 'annotate_mhc_batch', lambda items, **kw:
                        calls.append(('mhc', len(items), kw['threads'])))
    assert len(cli._typed_batch(['one.pdb', 'two.pdb'], 'human', threads=2)) == 2
    assert calls == [('tcr', ['one.pdb', 'two.pdb'], 2), ('mhc', 2, 2)]


def test_all_mode_honors_scheduler_and_affinity(monkeypatch):
    from tcren import runtime
    monkeypatch.setattr(runtime.sys, 'platform', 'linux')
    monkeypatch.setattr(runtime.os, 'cpu_count', lambda: 16)
    monkeypatch.setattr(runtime.os, 'sched_getaffinity', lambda pid: {0, 1, 2, 3}, raising=False)
    monkeypatch.setattr(runtime.Path, 'is_file', lambda p: False)
    monkeypatch.setenv('SLURM_CPUS_PER_TASK', '2')
    assert runtime.resolve_threads(0) == 2
    monkeypatch.setenv('SLURM_CPUS_PER_TASK', '1')
    assert runtime.resolve_threads(0) == 1


def _native_worker_probe(_):
    import os
    import time
    from threadpoolctl import threadpool_info
    time.sleep(.1)  # Keep both coarse workers occupied while collecting their identities.
    return os.getpid(), pl.thread_pool_size(), [p['num_threads'] for p in threadpool_info()]


def test_coarse_workers_are_distinct_and_native_serial():
    import pytest
    from tcren.runtime import available_cpus, parallel_map
    if available_cpus() < 2:
        pytest.skip('two allocated CPUs required for the actual spawn check')
    rows = parallel_map(_native_worker_probe, list(range(8)), 2)
    assert len({row[0] for row in rows}) == 2
    assert all(row[1] == 1 and all(n == 1 for n in row[2]) for row in rows)


def test_nested_cgroup_quota_limits_all_mode(monkeypatch):
    from tcren import runtime
    import pytest
    monkeypatch.setattr(runtime.sys, 'platform', 'linux')
    monkeypatch.setattr(runtime.os, 'cpu_count', lambda: 16)
    monkeypatch.setattr(runtime.os, 'sched_getaffinity', lambda pid: set(range(4)), raising=False)
    files = {'/proc/self/cgroup': '0::/job/child\n', '/sys/fs/cgroup/job/cpu.max': '150000 100000'}
    monkeypatch.setattr(runtime.Path, 'is_file', lambda p: str(p) in files)
    def read(path):
        if str(path) not in files:
            raise FileNotFoundError(str(path))
        return files[str(path)]
    monkeypatch.setattr(runtime.Path, 'read_text', read)
    assert runtime.resolve_threads(0) == 1
    with pytest.raises(ValueError, match='allocation'):
        runtime.resolve_threads(2)
    files.clear()
    files.update({'/proc/self/cgroup':'2:cpu,cpuacct:/job/child\n',
                  '/sys/fs/cgroup/cpu,cpuacct/job/cpu.cfs_quota_us':'200000',
                  '/sys/fs/cgroup/cpu,cpuacct/job/cpu.cfs_period_us':'100000'})
    assert runtime.resolve_threads(0) == 2


def test_wide_modes_fail_before_the_live_thread_ceiling(monkeypatch):
    import pytest
    from tcren import runtime
    monkeypatch.setattr(runtime, 'available_cpus', lambda: 1000)
    with pytest.raises(ValueError, match='live-thread ceiling'):
        runtime.resolve_threads(0)
    with pytest.raises(ValueError, match='live-thread ceiling'):
        runtime.check_parallel_entrypoint(122, 122)


def test_search_allowance_does_not_force_requested_parallelism_serial(monkeypatch):
    import os
    from tcren import runtime
    monkeypatch.setattr(runtime, 'available_cpus', lambda: 2)
    previous = dict(os.environ)
    with runtime.native_threads():
        with runtime.search_threads(2):
            assert os.environ['OMP_THREAD_LIMIT'] == '2'
            assert os.environ['OMP_NUM_THREADS'] == '2'
            assert os.environ['POLARS_MAX_THREADS'] == '1'
        assert os.environ['OMP_THREAD_LIMIT'] == '1'
    assert dict(os.environ) == previous


def test_stdin_parallel_map_preserves_results_and_errors():
    import subprocess
    import sys
    import pytest
    from tcren.runtime import available_cpus
    if available_cpus() < 2:
        pytest.skip('two allocated CPUs required')
    script = """
import os
import tcren
from tcren.runtime import parallel_map
previous = dict(os.environ)
assert parallel_map(abs, [-2, -1, 0, 1], 2) == [2, 1, 0, 1]
assert dict(os.environ) == previous
try:
    parallel_map(int, ['1', 'invalid'], 2)
except RuntimeError as exc:
    assert 'invalid literal' in str(exc)
else:
    raise AssertionError('worker failure was hidden')
assert dict(os.environ) == previous
print('interactive batch OK')
"""
    result = subprocess.run([sys.executable, '-'], input=script, text=True,
                            capture_output=True, timeout=30)
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == 'interactive batch OK'


def test_stdin_parallel_cancel_cleans_transport():
    import os
    import subprocess
    import sys
    import pytest
    from tcren.runtime import available_cpus
    if os.name != 'posix' or available_cpus() < 2:
        pytest.skip('POSIX and two allocated CPUs required')
    script = """
import signal
import tempfile
from pathlib import Path
import time
import tcren
from tcren.runtime import parallel_map
previous = set(Path(tempfile.gettempdir()).glob('tcren-workers-*'))
def interrupt(signum, frame):
    raise KeyboardInterrupt
signal.signal(signal.SIGALRM, interrupt)
signal.setitimer(signal.ITIMER_REAL, 1.5)
try:
    parallel_map(time.sleep, [30, 30], 2)
except KeyboardInterrupt:
    pass
else:
    raise AssertionError('cancellation was hidden')
finally:
    signal.setitimer(signal.ITIMER_REAL, 0)
assert set(Path(tempfile.gettempdir()).glob('tcren-workers-*')) == previous
print('cancelled cleanly')
"""
    result = subprocess.run([sys.executable, '-'], input=script, text=True,
                            capture_output=True, timeout=15)
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == 'cancelled cleanly'
