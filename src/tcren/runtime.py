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


def _windows_cpu_counts():
    """Return process-affinity logical/physical counts using Windows topology.

    Multiple processor groups cannot be represented by GetProcessAffinityMask;
    reject them instead of silently counting only the primary group.
    The current job rate is honored; Windows does not expose parent-job rates here.
    """
    import ctypes as C
    from ctypes import wintypes as W
    kernel = C.WinDLL("kernel32", use_last_error=True)
    kernel.GetCurrentProcess.restype = W.HANDLE
    process = kernel.GetCurrentProcess()
    groups = W.USHORT(64)
    group_ids = (W.USHORT * 64)()
    kernel.GetProcessGroupAffinity.argtypes = [W.HANDLE, C.POINTER(W.USHORT), C.POINTER(W.USHORT)]
    if not kernel.GetProcessGroupAffinity(process, C.byref(groups), group_ids):
        raise RuntimeError("Windows CPU topology or affinity query failed")
    if groups.value != 1:
        raise RuntimeError("CPU budgeting requires a single Windows processor group")
    mask, system = C.c_size_t(), C.c_size_t()
    kernel.GetProcessAffinityMask.argtypes = [W.HANDLE, C.POINTER(C.c_size_t), C.POINTER(C.c_size_t)]
    if not kernel.GetProcessAffinityMask(process, C.byref(mask), C.byref(system)):
        raise RuntimeError("Windows CPU topology or affinity query failed")
    size = W.DWORD()
    query = kernel.GetLogicalProcessorInformationEx
    query.argtypes = [W.DWORD, C.c_void_p, C.POINTER(W.DWORD)]
    query(0, None, C.byref(size))  # RelationProcessorCore
    if not size.value:
        raise RuntimeError("Windows CPU topology or affinity query failed")
    buffer = C.create_string_buffer(size.value)
    if not query(0, buffer, C.byref(size)):
        raise RuntimeError("Windows CPU topology or affinity query failed")
    data = buffer.raw
    offset, physical = 0, 0
    width = C.sizeof(C.c_size_t)
    while offset < len(data):
        length = int.from_bytes(data[offset + 4:offset + 8], "little")
        if length < 32 or offset + length > len(data):
            raise RuntimeError("Windows CPU topology or affinity query failed")
        count = int.from_bytes(data[offset + 30:offset + 32], "little")
        for index in range(count):
            start = offset + 32 + index * (width + 8)
            if start + width + 8 > offset + length:
                raise RuntimeError("Windows CPU topology or affinity query failed")
            bits = int.from_bytes(data[start:start + width], "little")
            group = int.from_bytes(data[start + width:start + width + 2], "little")
            if group == group_ids[0] and bits & mask.value:
                physical += 1
                break
        offset += length
    if not mask.value or not physical:
        raise RuntimeError("Windows returned empty CPU affinity or physical topology")
    logical = mask.value.bit_count()
    # Job CPU rate is a fraction of the host capacity, not the affinity mask.
    class Rate(C.Structure):
        _fields_ = [("flags", W.DWORD), ("rate", W.DWORD)]
    rate = Rate()
    query_job = kernel.QueryInformationJobObject
    query_job.argtypes = [W.HANDLE, C.c_int, C.c_void_p, W.DWORD, C.c_void_p]
    in_job = W.BOOL()
    kernel.IsProcessInJob.argtypes = [W.HANDLE, W.HANDLE, C.POINTER(W.BOOL)]
    if not kernel.IsProcessInJob(process, None, C.byref(in_job)):
        raise RuntimeError("Windows job membership query failed")
    if in_job.value and not query_job(None, 15, C.byref(rate), C.sizeof(rate), None):
        raise RuntimeError("Windows job CPU rate query failed")
    if rate.flags & 1:
        maximum = (rate.rate >> 16) if rate.flags & 16 else rate.rate
        if rate.flags & (4 | 16) and maximum:
            kernel.GetActiveProcessorCount.argtypes = [W.WORD]
            kernel.GetActiveProcessorCount.restype = W.DWORD
            host = kernel.GetActiveProcessorCount(0xffff) or os.cpu_count() or logical
            logical = min(logical, max(1, host * maximum // 10000))
    return logical, min(physical, logical)


def available_cpus(*, physical: bool = True) -> int:
    """Cores available under affinity, scheduler, container and physical-core limits."""
    counts = [os.cpu_count() or 1]
    if sys.platform == "win32":
        logical, cores = _windows_cpu_counts()
        counts.append(cores if physical else logical)
    if hasattr(os, 'sched_getaffinity'):
        affinity = os.sched_getaffinity(0)
        counts.append(len(affinity))
        try:
            cores = set()
            for cpu in affinity:
                topology = Path(f'/sys/devices/system/cpu/cpu{cpu}/topology')
                cores.add((topology.joinpath('physical_package_id').read_text(),
                           topology.joinpath('core_id').read_text()))
            if physical:
                counts.append(len(cores))
        except OSError:
            pass  # CPU topology can be hidden by a container; retain affinity limits.
    if physical and sys.platform == 'darwin':
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


def resolve_threads(threads: int = 0) -> int:
    """1 is serial, 0 selects all allocated physical cores; reject invalid/excessive counts."""
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
    available = available_cpus() if count == 0 else available_cpus(physical=False)
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
                                       stderr=log, start_new_session=os.name == 'posix',
                                       env={**os.environ, 'PYTHONPATH': os.pathsep.join(
                                           os.path.abspath(p) for p in sys.path)})
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
