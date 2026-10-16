"""One allocation-aware CPU budget; serial native kernels inside coarse workers."""
from __future__ import annotations

from contextlib import contextmanager
from functools import wraps
import operator
import os
from pathlib import Path
import subprocess
import sys
from threading import RLock

_NATIVE_ENV = ('POLARS_MAX_THREADS', 'POLARS_ASYNC_THREAD_COUNT',
               'POLARS_MAX_BLOCKING_THREAD_COUNT', 'RAYON_NUM_THREADS', 'OPENBLAS_NUM_THREADS',
               'MKL_NUM_THREADS', 'OMP_NUM_THREADS', 'OMP_THREAD_LIMIT',
               'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS')
# ponytail: process-global native limits serialize concurrent API calls; separate processes
# own independent budgets when simultaneous independent jobs are required.
_LOCK = RLock()


def available_cpus() -> int:
    """Cores available under affinity, scheduler, container and physical-core limits."""
    counts = [os.cpu_count() or 1]
    if hasattr(os, 'sched_getaffinity'):
        affinity = os.sched_getaffinity(0)
        counts.append(len(affinity))
        try:
            cores = set()
            for cpu in affinity:
                topology = Path(f'/sys/devices/system/cpu/cpu{cpu}/topology')
                cores.add((topology.joinpath('physical_package_id').read_text(),
                           topology.joinpath('core_id').read_text()))
            counts.append(len(cores))
        except OSError:
            pass  # CPU topology can be hidden by a container; retain affinity limits.
    if sys.platform == 'darwin':
        try:
            result = subprocess.run(['sysctl', '-n', 'hw.physicalcpu'], capture_output=True,
                                    text=True, timeout=5, check=False)
            if result.returncode == 0:
                counts.append(int(result.stdout))
        except (OSError, subprocess.TimeoutExpired):
            pass  # Optional topology is unavailable; affinity/allocation caps still apply.
    for name in ('SLURM_CPUS_PER_TASK', 'PBS_NP', 'NSLOTS'):
        if name in os.environ:
            try:
                value = int(os.environ[name])
            except ValueError as exc:
                raise ValueError(f'{name} must be a positive integer') from exc
            if value < 1:
                raise ValueError(f'{name} must be a positive integer')
            counts.append(value)
    # cgroup v2 paths may be nested below the mounted root.
    paths = [Path('/sys/fs/cgroup'), Path('/sys/fs/cgroup/cpu'),
             Path('/sys/fs/cgroup/cpu,cpuacct')]
    try:
        for line in Path('/proc/self/cgroup').read_text().splitlines():
            _, controllers, group = line.split(':', 2)
            if not controllers:
                paths.append(paths[0] / group.lstrip('/'))
            elif 'cpu' in controllers.split(','):
                paths.extend(root / group.lstrip('/') for root in paths[1:3])
    except FileNotFoundError:
        pass
    for root in paths:
        for parent in (root, *root.parents):
            if parent == Path('/sys/fs'):
                break
            quota = parent / 'cpu.max'
            if quota.is_file():
                amount, period = quota.read_text().split()
                if amount != 'max':
                    counts.append(max(1, int(amount) // int(period)))
            quota = parent / 'cpu.cfs_quota_us'
            if quota.is_file():
                amount = int(quota.read_text())
                if amount > 0:
                    period = int(quota.with_name('cpu.cfs_period_us').read_text())
                    counts.append(max(1, amount // period))
    return max(1, min(counts))


def resolve_threads(threads: int = 1) -> int:
    """1 is serial, 0 is all allocated cores; reject invalid/excessive counts."""
    if isinstance(threads, bool):
        raise ValueError('threads must be an integer: 0 (all allocated) or a positive count')
    try:
        count = operator.index(threads)
    except TypeError as exc:
        raise ValueError('threads must be an integer') from exc
    if count < 0:
        raise ValueError('threads must be nonnegative (0 means all allocated cores)')
    if count == 1:
        return 1
    available = available_cpus()
    if count > available:
        raise ValueError(f'requested {count} CPUs exceeds the available allocation ({available})')
    selected = count or available
    if selected + 16 >= 1000:
        raise ValueError("CPU request exceeds the live-thread ceiling; select fewer than 984 CPUs")
    return selected


@contextmanager
def native_threads(*, check_polars=True):
    """Limit new and already-loaded numerical kernels; restore caller settings."""
    from threadpoolctl import threadpool_limits

    with _LOCK:
        previous = {key: os.environ.get(key) for key in _NATIVE_ENV}
        try:
            os.environ.update(dict.fromkeys(_NATIVE_ENV, '1'))
            # Tokio needs spare blocking I/O workers to avoid dependency deadlocks.
            os.environ['POLARS_MAX_BLOCKING_THREAD_COUNT'] = '4'
            with threadpool_limits(limits=1):
                if check_polars:
                    import polars as pl
                    if pl.thread_pool_size() != 1:
                        raise RuntimeError('Polars was initialized with multiple threads before tcren; '
                                           'use the tcren CLI in a fresh process, or import tcren '
                                           'before initializing Polars')
                yield
        finally:
            for key, value in previous.items():
                if value is None:
                    os.environ.pop(key, None)
                else:
                    os.environ[key] = value


@contextmanager
def search_threads(threads):
    """Give the single external search its CPU budget, then restore serial worker limits."""
    threads = resolve_threads(threads)
    with _LOCK:
        keys = ('OMP_NUM_THREADS', 'OMP_THREAD_LIMIT')
        previous = {key: os.environ.get(key) for key in keys}
        try:
            os.environ.update(dict.fromkeys(keys, str(threads)))
            yield
        finally:
            for key, value in previous.items():
                if value is None:
                    os.environ.pop(key, None)
                else:
                    os.environ[key] = value


def initialize_native_threads():
    """Initialize immutable native pools before eager public-API imports."""
    with native_threads(check_polars=False):
        import numpy  # noqa: F401
        import scipy.linalg  # noqa: F401
        import polars as pl
        pl.thread_pool_size()


def single_threaded(function):
    """Apply native limits to a public batch call, including initialized interpreters."""
    @wraps(function)
    def call(*args, **kwargs):
        with native_threads():
            return function(*args, **kwargs)
    return call


def check_parallel_entrypoint(count, n_items):
    # Include bounded native/I/O pools plus parent, interactive broker and tracker/queue overhead.
    if 8 * min(count, n_items) + 24 >= 1000:
        raise ValueError("worker request exceeds the live-thread ceiling; use at most 121 workers")


def parallel_map(function, work, threads):
    """Ordered map over at most B coarse contiguous chunks; never parallel searches."""
    from multiprocessing import get_context

    threads = resolve_threads(threads)
    count = min(threads, len(work))
    check_parallel_entrypoint(count, len(work))
    if count <= 1:
        return [function(item) for item in work]
    with native_threads():
        entrypoint = getattr(sys.modules.get('__main__'), '__file__', None)
        if 'ipykernel' in sys.modules or not entrypoint or not Path(entrypoint).is_file():
            return _interactive_map(function, work, count)
        with get_context('spawn').Pool(count) as pool:
            return pool.map(function, work, chunksize=(len(work) + count - 1) // count)


def _interactive_map(function, work, count):
    """Run the existing spawn pool from an importable module, without replaying a notebook."""
    import pickle
    import signal
    from tempfile import TemporaryDirectory

    # Private, single-use transport, removed on success/failure; never a result cache.
    with TemporaryDirectory(prefix='tcren-workers-') as folder:
        root = Path(folder)
        with root.joinpath('input').open('wb') as stream:
            pickle.dump((function, work, count), stream, protocol=pickle.HIGHEST_PROTOCOL)
        with root.joinpath('stderr').open('w+b') as log:
            process = subprocess.Popen([sys.executable, '-m', 'tcren._parallel', folder],
                                       stderr=log, start_new_session=os.name == 'posix')
            try:
                code = process.wait()
            except BaseException:
                try:
                    if os.name == 'posix':
                        os.killpg(process.pid, signal.SIGTERM)
                    else:
                        process.terminate()
                except ProcessLookupError:
                    pass  # The broker may have completed as cancellation arrived.
                try:
                    process.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    if os.name == 'posix':
                        os.killpg(process.pid, signal.SIGKILL)
                    else:
                        process.kill()
                    process.wait()
                raise
            if code:
                log.seek(0, os.SEEK_END)
                log.seek(max(0, log.tell() - 8192))
                raise RuntimeError('interactive descriptor workers failed (no serial fallback):\n'
                                   + log.read().decode(errors='replace'))
        with root.joinpath('output').open('rb') as stream:
            return pickle.load(stream)
