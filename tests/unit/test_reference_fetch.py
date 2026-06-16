"""The canonical reference structures fetch from HF for an installed library/CLI.

When tcren orients a new, non-canonical structure but no local ``data/Native2026`` exists
(a pip-installed package), :func:`tcren.paths.reference_structure_path` must lazily download
the reference (1ao7/1fyt) from the HF dataset. This verifies that path works and is fast
(network fetch bounded; the cached re-lookup near-instant). Skipped if HF is unreachable.
"""

from __future__ import annotations

import time

import pytest

pytest.importorskip("huggingface_hub")

from tcren.paths import reference_structure_path
from tcren.structure import parse_structure


def test_reference_fetched_from_hf_when_absent(tmp_path, monkeypatch):
    # Point the data dir at an empty location so the local lookup misses and HF is used.
    monkeypatch.setenv("TCREN_DATA_DIR", str(tmp_path))

    t0 = time.perf_counter()
    try:
        path = reference_structure_path("1ao7")
    except FileNotFoundError as exc:  # offline / HF down — not a code failure
        pytest.skip(f"HF reference not fetchable here: {exc}")
    first = time.perf_counter() - t0

    assert path.exists()
    # Timing: the network fetch is bounded; a small reference must not take minutes.
    assert first < 120.0, f"reference fetch too slow: {first:.1f}s"

    # Cached re-lookup hits the HF cache — no network, near-instant.
    t1 = time.perf_counter()
    again = reference_structure_path("1ao7")
    cached = time.perf_counter() - t1
    assert again == path
    assert cached < 5.0, f"cached reference lookup too slow: {cached:.2f}s"

    # The fetched asset is a usable structure (1ao7 = MHCa/B2M/peptide/TRA/TRB).
    s = parse_structure(path, pdb_id="1ao7")
    assert sum(len(c.residues) for c in s.chains) > 100
