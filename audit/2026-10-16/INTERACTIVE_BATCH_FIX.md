# Interactive batched descriptor workers — issue 21

Implementation: `b6c0ce9`, on master base `8f860a2`. The shared descriptor path now
supports stdin and real Jupyter notebooks with multiple workers. Ordinary scripts
retain their main-guard requirement. No score, descriptor equation, potential,
model parameter, benchmark population or publication asset was changed.

## Existing APIs and CPU ownership

`inspect_predictions` keeps its existing whole-set `annotate_batch` and
`annotate_mhc_batch` calls, then reuses `recognition_table(..., annotate=False)`.
The shared `runtime.parallel_map` still owns an ordered, coarse spawn pool.
Interactive callers now launch that same pool from `python -m tcren._parallel`,
rather than asking spawn to reload a notebook or stdin entrypoint. The caller
waits; the broker does not score alongside its workers. CLI and guarded-script
calls retain their existing direct pool route. No fork, executor, outer search
pool, new dependency or serial fallback was introduced. The stdlib rung sufficed.

The private temporary transport directory is created with restricted permissions,
used once and removed on success, failure or handled cancellation. This is IPC
transport, not a reusable calculation cache. Worker errors include the child
traceback. POSIX cancellation signals the broker/process group, and the existing
Pool context terminates/joins its workers. The live-thread guard now includes
broker overhead. Parent native settings are restored after each call.

Actual real-kernel searches on the complete 1AO7 + 1FYT manifest:

| Reference | Queries per batch | Serial budget | Explicit-two / all-allocated budget |
|---|---:|---:|---:|
| Human TCR | 10 | 1 | 2 |
| MHC | 4 | 1 | 2 |

Each search was called once per reference from the kernel's same main caller
thread, before the broker/workers started. `inspect_predictions` specifies human
annotation here; no mouse search was requested. Actual MMseqs processes reached
two threads, matching their explicit allowance. A request for three CPUs under a
two-CPU allocation failed before any search. No search/worker overlap was observed.

## Fresh examples and manual inspection

Local outputs: `/tmp/tcren21-examples-final-20261016/`. Every example rebuilt from
raw public 1AO7/1FYT structures; no prior feature table was used to skip work.
`interactive_identity.json` pins raw-input and immutable-wheel SHA256 identities.
All native thread-limit environment overrides were unset at startup. Only the
scheduler allocation was set to two CPUs for these budget checks.

| Fresh model-inspection call, two input structures | Seconds | Descriptor workers |
|---|---:|---:|
| Real Jupyter kernel, serial | 9.195 | 0 |
| Real Jupyter kernel, two CPUs | 7.609 | 2 |
| Real Jupyter kernel, all allocated (two) | 7.277 | 2 |
| stdin, two CPUs | 7.809 | 2 |

These are single-run execution/resource diagnostics, not a scaling estimate.
Peak sampled simultaneous process-tree RSS was 1,196,720,128 bytes (1.11 GiB);
peak live threads was 39, including the kernel's own I/O, broker, tracking and
bounded Polars I/O threads. Both workers did work. A separate eight-item native
probe used two distinct workers, each reporting a one-thread Polars compute pool.
Loaded supported BLAS/OpenMP pools, when present, are checked by the regression
tests; the local Apple native probe reported no threadpoolctl-discoverable pools.
No user-exported native-limit workaround was needed.

Fresh CLI features, structural-reference scoring and mechanics also completed
with two workers. Fresh serial mechanics matches its parallel output exactly.
Scoring matches the preceding issue-43 artifact exactly on both inputs; features
retain the same 128 columns and physical IDs against the preceding eight-input
record. All seven inspection tables were compared against the fresh serial call,
for explicit-two, all-allocated and stdin arms: **21 table comparisons**, plus
three CLI comparisons (**24 total**). Keys, memberships, categorical values and
finite/missing masks match. All five changed floating cells across these checks
are retained, unrounded, in `interactive_numeric_differences.tsv`; maximum absolute
difference is 5.551115123125783e-17, confined to `partcoef_tcr`/`partcoef_pmhc`.
The known four-column unordered-topology limitation remains on the roadmap;
these software examples do not establish labelled metric/rank nonregression.

Manual inspection covered every output type: model descriptors/statuses, chains,
residue regions, directed contacts, angles, confidence interfaces, assessment,
CLI features, structural scores and mechanics. Class I and II remain identifiable:

| Public crystal | Class | Peptide | Residue rows | Directed contact rows | TP / TM / PM contacts |
|---|---|---|---:|---:|---|
| 1AO7 | I | LLFGYPVYV | 605 | 512 | 29 / 41 / 55 |
| 1FYT | II | PKYVKQNTLKLAT | 416 | 648 | 24 / 48 / 87 |

Each has five chain rows, one model row and one angle row. 1FYT retains both groove
chains (DRA/DRB), floor and helix regions; neither chain is replaced by class-I
markup. The two models report descriptors available and historical assessment not
requested. Assessment and confidence-interface tables are correctly empty without
a fitted opt-in or supplied predictor confidence. TP raw energies remain
−1.020363318891668 (1AO7) and −0.7241452366597441 (1FYT), dimensionless potential
sums, not physical binding energies. No fitted readout or publication figure was
created or removed. This change produces tables; no affected rendered figure
required redrawing.

Independent cancellation observation saw the broker, tracker and two workers,
then **zero remaining descendants** after cancellation (2.529 s including startup
and observation). Success/failure cleanup and caller-environment restoration are
also regression-tested. No kernel or worker is left running.

## Reproduce the examples

Use the local installed wheel/declared environment. Create a fresh output directory
and a manifest with unique `complex.id` and absolute raw `structure` paths for
1AO7 and 1FYT. Run the following in a real notebook and stdin, selecting each
supported CPU mode; all-allocated was tested under an allocation of two:

```python
import tcren  # initialize bounded native pools before Polars work
import polars as pl
from tcren.prediction import inspect_predictions

manifest = pl.read_csv("manifest.tsv", separator="\t")
tables = inspect_predictions(manifest, threads=2)  # also 1 and 0
for name, table in tables.items():
    table.write_csv(f"fresh-output/{name}.tsv", separator="\t")
```

`structures.txt` contains the two raw paths, one per line. The affected CLI probes:

```sh
tcren features -s structures.txt -i placement,interface,topology --no-metadata -t 2 -o fresh-features.tsv
tcren scoring -s structures.txt --delta --reference-mode structural --no-superimpose -t 2 -o fresh-scoring.csv
tcren mechanics -s structures.txt -t 2 -o fresh-mechanics.csv
tcren mechanics -s structures.txt -t 1 -o fresh-mechanics-serial.csv
```

## Verification and remaining gates

Final installed immutable wheel SHA256:
`d24f599861773b9d3ff69c174c6059e14aa245b46e64878e57954611197dceba`.
All 182 package payload files match the wheel; packaged Python matches source.
Ruff, diff checks and Sphinx with `SPHINXOPTS=-W` pass. The final installed-wheel
affected suite passed **58 tests in 166.48 s** (168.07 s including interpreter startup),
within its 180-second bound:

```sh
TCREN_HOME=$REPO .venv/bin/python -m pytest tests/unit/test_cli_budget.py tests/unit/test_prediction.py tests/unit/test_annotation_batch.py tests/unit/test_features_chunking.py tests/unit/test_exports.py tests/unit/test_pipeline.py -m 'not slow' -q --tb=short
```

`REPO` denotes the library checkout with configured public reference data. Remote
supported-Python CI and integration are verified separately on the associated PR.

Caveats, collected here: the first broad restricted-shell unit attempt encountered
ARDA index-lock permissions (994 passed, 41 failed, 112 errors, one network skip).
The normal-access full reference-backed attempt printed no failures before its
360-second cap at roughly 70%; it is **incomplete**, not a full passing suite.
The narrower affected suite has its own 180-second bound. Initial notebook
monitoring stopped on a short-lived-process permission race; that diagnostic
attempt is retained separately. Accepted real-kernel cells have no unreadable
process samples; four later CLI/stdin probes report five such transient samples
in total. Peaks are sampled observations, not continuous proofs. Preinitialized
multithreaded Polars pools remain unsupported; import tcren before Polars work.
No scientific benchmark, manuscript reproduction, MD campaign or release was rerun.
