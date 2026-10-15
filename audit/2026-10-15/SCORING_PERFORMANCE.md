# Scoring performance audit

The timed workload is the public 1ao7 fixture, annotated once, no orientation,
512 directed native contact rows. Both executions used one Polars/BLAS thread.
The reproduction command (from the library root) is:

```sh
TCREN_HOME=$PWD POLARS_MAX_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  PYTHONPATH=src .venv/bin/python audit/2026-10-15/profile_table_scoring.py
```

Before: installed structural-reference wheel from master b42b8e3.
After: current source with a single reference-map build shared across interfaces.
Times are the minimum of three repeated batches, not a whole-cohort throughput claim.
A concurrent manuscript run was active, so absolute timings are host-load dependent.

| Operation | Before (ms/structure) | After (ms/structure) |
|---|---:|---:|
| Three preselected contact-table energy sums (1000 repetitions/batch) | 0.03347 | 0.03184 |
| Typed pipeline including native/structural-reference contacts (30 repetitions/batch) | 20.468 | 13.645 |

The typed pipeline is about 1.50 times faster on this fixture. Profiling counted
90 contact-map builds for 30 structures before and 60 afterward: the native map
plus one poly-alanine map, rather than two identical reference builds. Potential
matrices and interface selections already reuse in-memory representations; no
new persistent intermediate cache, native extension or dependency is needed.

The energy sum itself was already negligible. Coordinate construction, annotation,
SASA and repeated command-level passes over structures must not be attributed to
contact-table arithmetic. SASA cannot be recovered from a collapsed residue-contact
table. Old/new potentials can use the same freshly generated contact maps within
one execution; the present manuscript driver still launches separate CLI passes,
so this change alone does not eliminate that orchestration cost.

Validation: ten pipeline tests passed, including exact equality to independently
computed per-interface structural references for all/cdr/cdr+fr region selections
and residue/atomic contact weighting. The test also requires exactly one reference
map build. No scores, thresholds, populations or fitted readouts were changed.

Built-wheel checks: 18 pipeline/potential tests passed in 77.83 seconds from a
separate temporary installation; Ruff and warning-free Sphinx build passed.
Wheel SHA256: `84ab1ef149502746e4a8ca65c3a27b70609d1a4877b9a941a3231490061d74a8`.
The active manuscript reproduction kept its original pinned wheel throughout.
