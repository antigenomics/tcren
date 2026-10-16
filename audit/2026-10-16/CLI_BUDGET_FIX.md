# CLI batching and CPU budget fix — issue 43

Implementation checkpoints: `9514ec4` and `26e92e7`; final tests/evidence accompany
this report. No potential table, fitted parameter, score formula or scientific
benchmark population was changed. The two pre-existing dirty research outputs
remain outside the commits.

## One complete MMseqs batch, using its own multicore implementation

There is **no Python thread/process pool around MMseqs** and no per-structure
sequence search. `annotation.batch.annotate_batch` owns one complete query list
per species reference. MHC annotation receives the complete candidate-chain set
in one call. These searches run sequentially; descriptor/scoring workers start
only after annotation finishes.

An independent observation wrapper around the real ARDA MMseqs APIs, using
1AO7 + 1FYT in one manifest, recorded exactly:

| Reference / API | Query sequences in the one batch | MMseqs CPU allowance | Python caller |
|---|---:|---:|---|
| Human TCR / `search` | 10 | 2 | Main caller thread |
| Mouse TCR / `search` | 10 | 2 | Same main caller thread |
| MHC / `easy_search` | 4 | 2 | Same main caller thread |

The runtime process trace also observed **two actual MMseqs threads**, not merely
`--threads 2` in an argument list. Internal MMseqs createdb/search/conversion child
processes are its own pipeline, not additional independently scheduled query batches.
No search/descriptor-worker overlap was observed in any of the 23 final probes.
The regression test checks complete query membership and a single caller thread.

## Shared implementation

- `surface` validates positive grid dimensions before annotation; batches TCR/MHC
  searches; preserves parsing, typing, search and map errors; and returns nonzero
  for empty/all-failed output. A mixed batch keeps its failed rows and succeeds
  only if at least one map is produced. Malformed/empty structures remain explicit.
- `_typed_batch` and `iter_typed` combine all supplied sources before one annotation
  batch. The requested CPU count reaches both TCR and MHC searches. Exact exports,
  prediction inspection, features, recognition, scoring and mechanics share the
  repaired resource controls. Annotation query/record cardinality is checked.
- Default 1 is serial, 0 resolves all allocated cores, and a positive request must
  fit the allocation. Detection accounts for affinity, available physical topology,
  scheduler counts and nested Linux cgroup v1/v2 quotas. Excessive live-thread
  projections fail before worker creation; no silent serial fallback.
- Native numerical pools use one compute thread per coarse spawn worker. Searches
  temporarily receive their own selected OpenMP allowance, then restore serial
  numerical-worker settings. Worker pools never wrap sequence searches.
- Already-loaded BLAS/OpenMP limits and caller environment settings are restored.
  An existing multithreaded Polars pool cannot be resized and fails explicitly;
  ordinary CLI processes initialize their own native pools without user exports.
- Polars async I/O uses one thread and blocking I/O is capped at four auxiliary
  workers per process. An initial diagnostic cap of one produced Polars' deadlock
  warning and was rejected. The accepted configuration completes without it.
  Bounded auxiliary I/O threads are distinct from numerical compute threads.

The vectorized/table calculations and existing native kernels sufficed. No new
scientific control, replacement score, fitting, derived-result cache or worker-per-query
implementation was introduced. `threadpoolctl`, already present through sklearn,
is now an explicit dependency because initialized native pools require runtime caps.

## Fresh examples, scaling and manual inspection

Authoritative local deliverable: `/tmp/tcren43-examples-accepted-20261016/`.
All 23 example/probe commands have their expected exit status. The original
15-command receipt is extended with a two-structure surface batch, two fresh runs
at each of one/two CPUs on eight public structures, restricted-all/excessive-count
checks and two-worker scoring. Every run rebuilds from raw structures. There is
no resume/skip-if-exists path. CLI logs and annotation/potential sidecars remain
with the outputs; portable summary receipts are committed beside this report.

The eight inputs are 1AO7, 1FYT, 2AK4, 3HG1, 3UTS, 5MEN, 5NME and 7NMG. These are
software examples, not a newly labelled scientific benchmark. Both classes remain
identifiable; no pooled class-performance claim is made.

| Eight-structure feature run | Fresh times, s | Median, s | Descriptor workers |
|---|---|---:|---:|
| One CPU | 20.886; 21.059 | 20.9725 | 0 (serial caller) |
| Two CPUs | 13.332; 13.517 | 13.4245 | 2 distinct processes |

The median ratio is 1.56 on this workload (two repetitions per mode), not a universal
speedup estimate. Peak simultaneous process-tree RSS across all probes is
1,223,294,976 bytes (1.14 GiB); maximum live OS threads is 23, including parent,
MMseqs internals, tracking/queue and bounded I/O threads. Each descriptor worker's
Polars/loaded BLAS compute cap is independently checked as one. Restricted-all mode
resolves two CPUs under a two-CPU allocation; requesting three fails before search.

Manual review opened both newly rendered 520 × 760 class-I/class-II previews
individually, and inspected the scores, pMHC features, annotation, contacts,
position/residue contributions and CPL table. No clipping or peptide-identity
change was found. Both SVGs are byte-identical to the preceding examples. The
native diagnostic maps' missing colour-axis/outline legend and individually scaled
colour ranges remain the earlier documented presentation limitations. No manuscript
heatmap was added, and no publication figure was replaced.

On physical keys, raw scores, structural references, SASA/hydropathy, CPL values,
contacts and decompositions match the earlier examples exactly. All finite/missing
masks and row sets match. Position contributions still sum to TP energy; residue
endpoints count each contact twice. The 180 CPL cells retain zero WT effects and
`phi = phi_tcr + phi_mhc`.

Known unordered-reduction variability remains visible: 27 changed floating cells
across the 16 matched table comparisons, maximum absolute difference
5.551115123125783e-17. Only `co_pep`, `co_mhc`, `partcoef_tcr` and `partcoef_pmhc`
change; every signed difference is retained in `cli_budget_numeric_differences.tsv`.
No rounding suppresses these differences. They do not establish metric/rank
nonregression for a labelled scientific cohort; that prior limitation remains open.

## Verification and identity

Installed immutable wheel SHA256:
`670d0145b5cd62563990e40e31aaf2e1c35407f5446ca8d5b11c9a6c6548d73b`.
All **181** package payload files match the wheel, and packaged Python matches source.
Versions and all comparison shapes are in `cli_budget_identity.json`.

- Final affected installed-wheel checks: **54 passed**, 18 slow cases deselected,
  7.36 s; export/prediction, native budgets, surface and annotation checks.
- Explicit whole-query-set/caller-thread and budget checks after the author's
  clarification: **14 passed**, 2.75 s. No sequence-search pool was added.
- Final installed-wheel pipeline checks with the documented repository data root:
  **10 passed**, 58.91 s, including canonical orientation.
- Earlier broader checkpoint: **1,101 passed**, 129 deselected, one canonical-data
  lookup failure, 327.45 s. Its installed-wheel default root lacked the canonical
  structures. The configured final pipeline result above resolves that test setup;
  this is not a claimed complete final scientific reproduction or clean bootstrap.
- Ruff, diff checks and local warning-as-error Sphinx pass. The first docs CI attempt
  exposed a missing `threadpoolctl` in its manually enumerated lean toolchain;
  its recipe now installs the explicit runtime dependency rather than suppressing
  import warnings. Exact-head CI is checked before integration.

Linux allocation edge cases have synthetic checks and the normal Linux CI legs;
this report does not certify target-cluster deployment or every optional engine.
The stopped/earlier diagnostic attempts remain local and are not promoted as final
outputs. Historical scores and recorded metric regressions are untouched.

## Issue reconciliation

This resolves #43's software acceptance criteria. #21 retains the separate
interactive-safe parallel-backend request; guarded scripts/CLI and serial notebooks
remain supported. #26's long MD reconstruction, #32's historical Potts CLI gap,
#33/#35's publication/archive closure and #34's per-render CLI lineage audit remain
separate work. No replacement fitted model or scientific rerun closes them here.
The release containing these post-4.0.0 fixes still needs its own artifact/publication
checks and human approval gate. Local 4.0.1 metadata is not a public release.
