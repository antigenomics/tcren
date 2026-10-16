# MD diagnostics and issue 26 — 16 October 2026

The installed CLI now exposes the existing trajectory, close-atom and numerical
energy/force APIs. No simulation, scientific score, fitted model or recorded
benchmark/manuscript result was rerun or replaced. Issue 26 remains open for its
long-run acceptance gate; this change resolves its library/CLI gap.

## Ownership and exact quantities

| CLI | Existing calculation owner | Export |
|---|---|---|
| `md inspect` | `validation.md.inspect_md_trajectory`, `component_relative_rmsd`, `minimum_potential_frame` | Partner-relative backbone RMSD in Å/ps, exact atom/component identities, whole-system minimum saved-frame potential in kJ/mol, first/minimum/final coordinates, optional contacts/occupancy/docking |
| `md quality` | `validation.structures.coordinate_diagnostics`, `close_atom_pairs` | Coordinate diagnostics and named nonlocal heavy-atom pairs, including author residue/insertion/altloc identities, possible-disulfide flags and explicitly unknown bond classification |
| `md compare-forces` | `validation.md.compare_same_frame_energy_forces` | Same-coordinate, equal-time numerical energy/force errors in kJ/mol and kJ/mol/nm; no simulation or binding-energy interpretation |

Each invocation records the caller-declared MHC class, original input hashes,
output hashes and units. No class pooling, fitting, new controls or confidence
assignment. The original component-fit and minimum-frame calculations are retained.
The shared trajectory reader now rejects ambiguous/overlapping component mappings
and author-residue/insertion/altloc correspondence mismatches instead of silently
accepting sequence-position-only matches. It owns allocation-aware native CPU
limits for all callers. Exporting uses existing table/array/native APIs; no new
native kernel, worker pool or dependency was necessary.

`scripts/inspect_gromacs_trajectory.py` retains its original GROMACS preparation
interface and outputs. `md inspect` accepts its prepared continuous PDB, not raw
XTC. The installed CLI does not modify active campaign orchestration.

## Verification

Implementation base: `4f6d1f1`; first implementation commit: `b04ab33`.
The immutable task wheel's SHA256 and all 183 byte-matched installed payload files
are recorded in [MD_DIAGNOSTICS_CHECKS.json](MD_DIAGNOSTICS_CHECKS.json). Source
Python files match that wheel. The existing two unrelated dirty files retain their
original SHA256 values. No shared Stratos checkout was modified.

The installed wheel passed **64 focused tests in 21.62 s**, covering MD CLI/API
equality, required acknowledgements, failed provenance, ambiguous mappings,
author identities, saved-time selection, partner displacement, force dimensions/
units, alternate locations, CPU allocation and both CLI entry points. Ruff passes
for `src`, `tests` and the example checker. `make -C docs html SPHINXOPTS=-W`
passes. Local checks do not establish release publication.

Reproduce the public software examples into a new output directory:

```sh
.venv/bin/python audit/2026-10-16/check_md_cli_examples.py \
  --structures data/Canonical2026 --output /NEW/OUTPUT
```

These fixtures comprise two rigidly translated public crystal frames with synthetic
potential values. They test software and are **not MD evidence**. Raw structures
remain in the owning mirror; the input hashes and exact commands are in the receipt.
Every derived fixture/export is freshly built; there is no resume/reused result path.

| Example | Class | TCR atoms | Peptide atoms | MHCα atoms | MHCβ atoms | β2m atoms (whole fit only) | RMSD rows | Chain-pair count rows | Occupancy rows | Docking rows |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1OGA | I | 876 | 36 | 1,104 | 0 | 396 | 24 | 20 | 259 | 2 |
| 1FYT | II | 900 | 52 | 720 | 716 | 0 | 24 | 20 | 391 | 2 |

For class I, `-t 1`, `-t 2` and `-t 0` (allocation capped at two CPUs) each match
the API exactly for RMSD, contact counts, occupancy and docking; the selected
minimum-frame metadata and coordinates also match exactly. Class II passes the
same exact checks independently at `-t 1`. Both quality exports match the API
and retain typed headers with zero pairs below 2 Å. Nonempty atom-pair outputs,
alternate locations and possible-disulfide flags are covered by focused fixtures.
The force CLI fixture reproduces the existing one-atom API check: a 1 kJ/mol
potential difference and 0.2 relative RMS force error; these are synthetic units
checks, not a PME scientific result.

Manual inspection covered the physical atom/component selections, partner-relative
rows, all chain-pair count identities including zero-contact pairs, occupancy,
docking rows, selected timestamp/coordinate files, empty quality headers and
source/annotation provenance. Class II retains both groove chains. Both examples
select the second saved timestamp, 10 ps, at synthetic potential −20 kJ/mol.
No figure or render was changed by this software task.

Six installed CLI invocations took **17.7478 s combined**, with maximum sampled
process-tree RSS **1,175,601,152 bytes** and **10 live threads**, including auxiliary
threads. The 20 ms monitor records per-process search budgets; observed MMseqs
processes use at most the selected one/two compute threads. Searches annotate the
single reference as one complete batch per reference database, never per frame or
inside an outer pool. Numerical/frame phases remain serial. All native pools in the
checker report one thread. These timings exclude the separate API comparisons.
Every CLI call is bounded at 90 s and 4 GiB; neither bound was reached.

## Remaining campaign gate and limitations

Read-only Stratos inspection verified the original production controller and active
GROMACS job. The PME-order 4 arm has a terminal 50,000,000-step, 100,000 ps record.
The PME-order 8 arm was still running, last recorded at 10,418,000 steps / 20,836 ps
of its 100,000 ps target. Final paired analysis was waiting for MD. The observed
job uses one MPI rank and eight OpenMP threads on a 16-physical-core host; its
snapshot has 14 total live threads and 670.9% process CPU. The latter is a snapshot,
not a full-campaign resource bound or proof of numerical convergence. No simulation,
controller, source checkout, setting or output on Stratos was changed or restarted.

Issue 26 must stay open until both requested trajectories and final paired analysis
pass their existing terminal/input/time/CPU-rerun checks. Missing original MDP/TPR,
repaired inputs and trajectories remain historical provenance gaps; reconstructed
work is not exact historical reproduction. Only issue 26 was open at task start;
other completed issues were not reopened or closed without new acceptance evidence.

PBC preparation is declared, not inferred or performed by `md inspect`; no-unit-cell
status remains explicit. Whole-system potential is not binding energy. Proximity
exclusions remain sequence heuristics and covalent bond classification stays unknown.
`compare-forces --same-system` is a caller declaration of common topology,
coordinates and atom order; original rerun commands/MDPs must substantiate it.
The class argument records the caller's assignment, not a guessed annotation.

The first example attempt's 1e-12 Å rigid-translation assertion failed because the
native PDB parser stores float32 coordinates. The observed maximum is 6.66e-7 Å
for 1OGA and 5.67e-7 Å for 1FYT. The checker now records a coordinate-scale bound
of `4 * eps(float32) * max(atom coordinate norm)` (3.03e-5 / 2.83e-5 Å).
**CLI/API comparisons still use exact equality**, without rounding or tolerance.
All fixtures/exports were rebuilt in a new directory after this checker correction;
the failed attempt remains local. A documentation underline warning was corrected
before the zero-warning acceptance build. These corrections change no scientific
quantity or recorded result.
