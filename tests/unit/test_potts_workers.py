"""Worker payload isolation and bounded process startup, without scientific fits."""
import os
import threading

import numpy as np
import pytest

from tcren.potts import score


def test_simultaneous_serial_calls_keep_their_own_payload(monkeypatch):
    barrier = threading.Barrier(2)

    def body(s):
        barrier.wait(timeout=5)
        return score._unpack(s)[0]

    monkeypatch.setitem(score._BODIES, 'probe', body)

    def run(name):
        payload = score._payload(np.array([name]), np.array([0, 1]), np.zeros(1),
                                 np.zeros(1), np.zeros((0, 3), dtype=int),
                                 np.array([0, 0]), np.zeros(0))
        return score._map_structures('probe', 1, payload, {}, workers=1)

    results = {}
    def work(name):
        results[name] = run(name)
        assert score._SHARED.get(None) is None

    workers = [threading.Thread(target=work, args=(name,)) for name in ('first', 'second')]
    for worker in workers:
        worker.start()
    for worker in workers:
        worker.join(10)
        assert not worker.is_alive()
    assert results == {'first': ['first'], 'second': ['second']}


def test_default_workers_respect_allocation_and_reuse_batched_runtime(monkeypatch):
    seen = {}
    monkeypatch.setattr(os, 'cpu_count', lambda: 64)
    monkeypatch.setattr('tcren.runtime.available_cpus', lambda **kwargs: 2)
    monkeypatch.setitem(score._BODIES, 'probe', lambda s: s)

    def parallel(function, chunks, threads):
        seen.update(threads=threads, chunk_sizes=[len(chunk[1]) for chunk in chunks])
        return [function(chunk) for chunk in chunks]

    monkeypatch.setattr(score, 'parallel_map', parallel)
    assert score._map_structures('probe', 8, {}, {}, workers=None) == list(range(8))
    assert seen['threads'] == 2
    assert seen['chunk_sizes'] == [4, 4]
    assert score._map_structures('probe', 1, {}, {}, workers=3) == [0]


def test_real_spawned_workers_preserve_serial_results():
    from tcren.runtime import available_cpus
    if available_cpus() < 2:
        pytest.skip('two allocated CPUs are required to exercise a real process split')
    payload = score._payload(np.array([f'p{i}' for i in range(8)]), np.arange(9),
                             np.zeros(8), np.ones(8), np.zeros((0, 3), dtype=int),
                             np.zeros(9, dtype=int), np.zeros(0))
    options = dict(particles=2, steps=2, seed=0)
    serial = score._map_structures('score_sites', 8, payload, options, workers=1)
    parallel = score._map_structures('score_sites', 8, payload, options, workers=2)
    assert parallel == serial
