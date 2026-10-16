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

_NATIVE_ENV = ('POLARS_MAX_THREADS', 'RAYON_NUM_THREADS', 'OPENBLAS_NUM_THREADS',
               'MKL_NUM_THREADS', 'OMP_NUM_THREADS', 'OMP_THREAD_LIMIT',
               'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS')
# ponytail: process-global native limits serialize concurrent API calls; separate processes
# own independent budgets when simultaneous independent jobs are required.
_LOCK = RLock()


def available_cpus() -> int:
    """Cores available under affinity, scheduler, container and physical-core limits."""
    counts = [os.cpu_count() or 1]
    if hasattr(os, 'sched_getaffinity'):
        counts.append(len(os.sched_getaffinity(0)))
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
    paths = [Path('/sys/fs/cgroup')]
    try:
        for line in Path('/proc/self/cgroup').read_text().splitlines():
            if line.startswith('0::'):
                paths.append(paths[0] / line[3:].lstrip('/'))
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
    quota = Path('/sys/fs/cgroup/cpu/cpu.cfs_quota_us')
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
    return count or available


@contextmanager
def native_threads(*, check_polars=True):
    """Limit new and already-loaded numerical kernels; restore caller settings."""
    from threadpoolctl import threadpool_limits

    with _LOCK:
        previous = {key: os.environ.get(key) for key in _NATIVE_ENV}
        try:
            os.environ.update(dict.fromkeys(_NATIVE_ENV, '1'))
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
    if count > 1 and n_items > 1:
        entrypoint = getattr(sys.modules.get('__main__'), '__file__', None)
        if not entrypoint or not Path(entrypoint).is_file():
            raise ValueError("parallel descriptors require an importable script with an "
                             "if __name__ == '__main__' guard; use the CLI for batches "
                             "or threads=1 in notebooks/stdin")


def parallel_map(function, work, threads):
    """Ordered map over at most B coarse contiguous chunks; never parallel searches."""
    from multiprocessing import get_context

    count = min(threads, len(work))
    check_parallel_entrypoint(count, len(work))
    if count <= 1:
        return [function(item) for item in work]
    with native_threads():
        with get_context('spawn').Pool(count) as pool:
            return pool.map(function, work, chunksize=(len(work) + count - 1) // count)
