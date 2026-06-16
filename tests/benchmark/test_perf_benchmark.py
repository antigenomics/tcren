"""Speed + memory benchmark for the core pipeline.

Guarded by ``RUN_BENCHMARK=1`` (excluded from the normal suite). Reports wall time per stage
and peak RSS; needs arda + a local Native2026 structure. Run with::

    RUN_BENCHMARK=1 pytest -k benchmark -s
"""

from __future__ import annotations

import os
import platform
import resource
import time

import pytest

pytestmark = pytest.mark.skipif(not os.getenv("RUN_BENCHMARK"), reason="set RUN_BENCHMARK=1 to run")
pytest.importorskip("arda")

from tcren.contactmap import ContactMap
from tcren.paths import reference_structure_path
from tcren.potential import tcren
from tcren.scoring import score_peptides
from tcren.structure import parse_structure


def _best(fn, n=5):
    return min((lambda: (t := time.perf_counter(), fn(), time.perf_counter() - t)[-1])() for _ in range(n))


def _rss_mb():
    div = 1024 * 1024 if platform.system() == "Darwin" else 1024
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / div


def test_pipeline_benchmark():
    from tcren.annotation import classify_chains

    path = reference_structure_path("1ao7")
    t_parse = _best(lambda: parse_structure(path, pdb_id="1ao7"))
    s = parse_structure(path, pdb_id="1ao7")
    classify_chains(s, organism="human")
    t_cm = _best(lambda: ContactMap.from_structure(s))
    cm = ContactMap.from_structure(s)
    cands = ["GILGFVFTL", "NLVPMVATV", "KQWLVWLFL"] * 333
    t_score = _best(lambda: score_peptides(cm, cands, tcren(), interface="tcr_peptide"))

    print(f"\nplatform: {platform.platform()}")
    print(f"parse_structure          : {t_parse*1000:7.1f} ms")
    print(f"ContactMap.from_structure: {t_cm*1000:7.1f} ms")
    print(f"score {len(cands)} peptides    : {t_score*1000:7.1f} ms ({t_score*1e6/len(cands):.1f} us/peptide)")
    print(f"peak RSS                 : {_rss_mb():7.0f} MB")
    assert t_parse < 5.0 and t_cm < 2.0 and t_score < 5.0  # generous sanity ceilings
