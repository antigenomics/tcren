# Deterministic literature topology — 16 October 2026

The roadmap reduction gate is resolved. `co_pep`, `co_mhc`, `partcoef_tcr`
and `partcoef_pmhc` retain their original dimensionless definitions, finite
support and raw increasing direction. Their floating means now consume a fixed
physical-key order. No labels, fitted weights, rounding, tie merging or thresholds
enter the fix. No Stratos command was run and the MD campaign remains independent.

## Cause and implementation

`topology.literature.participation_coefficient` previously reduced floating
participation values in the unspecified order of Polars hash groups. Integer
counts and squared counts are unchanged; final residue rows are now sorted by
chain/index, including target identity on the pMHC side.

`contact_order` already sorts contacted sequence indices inside each group.
It now also traverses target-chain/CDR-loop groups in sorted order before
averaging their original length-normalised sequence separations. Class-II chains
retain distinct sequence axes. Existing callers in the descriptor/recognition
table and installed `features -i topology` CLI all use this shared implementation.

An exact row-permutation test failed before the fix, with different last bits in
contact order and participation. Both class-I and class-II edge-set permutations
now pass. Sorting uses the existing native Polars operations; no new parallel
layer, native kernel or dependency was needed. The added bound is O(E log E)
sorting time and O(E) storage for E contact edges; original pairwise separation
work is unchanged. Guarantees apply within the pinned numerical environment,
not across arbitrary architectures or dependency versions.

## Fresh verification and retained regressions

The immutable wheel, built from numerical source commit `2cb08f7`, has SHA256
`c2d09c0440da04f7b58e51543d79fb06cc04b76bd9779d1de2e1ec1ee943a846`.
All 183 installed payload files match it. The existing raw-reference inputs and
actual searched ARDA/MHC payload identities are recorded in the owning benchmark.

Two complete fresh CLI topology passes, under different process hash seeds,
rebuild 2,102 CPL, 618 raw TCRvdb and 1,089 VDJdb structures. **All 71 exported
columns agree exactly** across passes. Publication comparisons retain the existing
617 eligible TCRvdb rows and preserve the excluded raw row separately. Among 70
comparable topology quantities, only the intended four change from the prior
publication baseline; the other 66 agree exactly in each cohort.

The four descriptors produce 15,232 eligible value comparisons. There are 3,603
changed values, at most 1.1102230246251565e-16 in absolute dimensionless value,
and 1,888 changed within-case ranks. Finite masks, labels, confidence, class,
template strata and the eligible populations are unchanged. R uses the existing
matched metric definitions and exact-double inputs, without sign optimisation.

| Baseline value set | Metric entries | Decreased | Increased | Unchanged |
|---|---:|---:|---:|---:|
| Prior publication rows | 568 | 54 | 83 | 431 |
| Earlier independent A3 values, raw cohort | 400 | 22 | 34 | 344 |
| Earlier independent B values, raw cohort | 400 | 23 | 32 | 345 |

These are repeated metric/endpoint views, not independent experiments or a pooled
performance score. All losses, gains, signed values, matched row membership and
original tables remain. Raw historical values are re-evaluated under explicit
assay grouping; original historical metric exports are not rewritten. No fitted
model, replacement score or manuscript claim is promoted by this audit.

The private benchmark owns the full `DETERMINISTIC_TOPOLOGY.md` record and
`results/deterministic-topology-20261016/`: raw-source/command/wheel identities,
both CLI passes, exact before/after values and ranks, class/template-stratified
metrics and membership, resource receipts and the deliverable checker.

## Budgets and checks

The six cohort CLI calls complete in a monitored **592.623 s** end-to-end run,
at maximum sampled process-tree RSS **11,023,007,744 bytes** and **33 live
threads**, including auxiliary threads. The compute budget is four CPUs;
observed MMseqs `--threads` values and peaks are four. Annotation consumes each
whole cohort in one batch per reference, before four coarse single-kernel workers;
there is no per-model search or outer search pool. Sampling is 0.1 s, with 15
transient process-read misses. Neither the 1,700 s nor 16 GiB cap was reached.
The R summary uses one thread and completes in 2.301 s.

Separate public 1OGA/class-I and 1FYT/class-II CLI examples also agree exactly at
budgets 1, 2 and 0 (allocation capped at two). The three calls take 18.900 s,
peak at 1,190,559,744 sampled RSS bytes and 19 live threads. Observed search
threads are 1/2/2. These are crystal software checks, not modelled class-II ranking.

Library/family checks: 23 passed, two slow tests deselected; benchmark accounting:
14 passed. Ruff, zero-warning Sphinx docs and exact-head Python 3.10/3.12 CI pass.
A legacy slow test initially encountered a sandbox reference-lock permission
failure; the bounded installed CLI examples supply the reference-backed checks.
No scientific failure was hidden or relabelled as a passing test.

The audit additionally caught default R TSV formatting collapsing adjacent doubles.
The existing benchmark TSV helper now provides an opt-in 17-significant-digit
round-trip writer. Its runnable adjacent-double test passes, and every final
before/after cell matches the CLI/baseline input bit-for-bit. The initial-format
diagnostic is retained separately. This fixes representation, not the calculations
or tie policy. All 57 pre-existing dirty files across the two repositories are
preserved. Release/PyPI approval and MD endpoint acceptance remain separate gates.
