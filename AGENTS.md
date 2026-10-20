# AGENTS.md — tcren

`tcren` owns structure → descriptors → scores and frozen-model fitting for alphabeta
TCR:peptide–MHC complexes. Read `skills/tcren/SKILL.md` for APIs, `AUDIT.md` for
verified findings and unresolved mathematical claims, and `ROADMAP.md` for publication gates.

## Final audit directive — 15 October 2026

**Follow-up authorization:** the author subsequently authorized a separate fix branch,
fresh reproduction, PR to master, Copilot review and release preparation to the human
PyPI approval gate. Fixed empirical potential tables are explicitly permitted with
provenance; learned Q/T weights and cohort normalization remain excluded from final
models. This supersedes the audit-only stopping point below. Current implementation
and verification status: `audit/2026-10-15/RELEASE_FIXES.md`.

The author now requires **fit-free final models**, archival of fitted readouts, and CLI-owned
scientific calculations. Gaussian readouts, confidence GLMs, fitted transition links and fitted
Potts final scorers are historical evidence, not final-model candidates. The existing library APIs
and coefficient files remain inspectable for historical reproduction; their presence does not
approve them for the final paper. Native-reference covariance/correlation weights and empirical
potential tables must expose their estimated dependencies rather than claiming coefficient-free
status by default. See the manuscript's `FIT_FREE_DISPOSITION.md`.

Metric regressions are allowed when explicitly marked alongside gains on identical populations;
never hide losses, change floors, mix populations, or substitute a new fitted scorer. Audit direct
Python producers when the CLI lacks an equivalent quantity; report the gap/error as a repository
issue instead of swapping in a different calculation. The current assignment ends after reconciling
and recording the audit evidence. Stop before follow-up library fixes or replacement-paper work.
This directive supersedes older frozen-final-model and metric-nonregression publication rules below.

## Ownership and evidence

- Library code, equations, tests and API documentation belong here.
- `~/vcs/projects/2026-tcren2-code` owns datasets, evaluation and reproduction.
- `~/vcs/manuscripts/2026-tcren2-ms` owns prose and generated publication assets;
  `FINDINGS.md` owns standing results and `FIGURES.md` owns figure dependencies.
- Read recorded outputs for an audit. Do not rerun completed analyses without explicit
  authorization. An authorized fresh run rebuilds every derived block from raw inputs;
  no computed caches, resume or skip-if-exists paths. Raw references and frozen parameters
  are inputs. Never overwrite a frozen model as a side effect of evaluation.
- Preserve dirty edits and measured results. Report any contradiction that weakens a standing
  claim before changing publication prose. Do not introduce controls, nulls, sweeps or thresholds.
- Keep private provenance and identifiers out of public artifacts. `SOURCES` owns library
  parameter provenance; the benchmark owns dataset provenance.

## Scientific contract

Five blocks: CPL peptide ranking, validated receptor ranking, receptor ranking under template
scarcity, dynamics/kinetics, and confidence diagnostics. A capability must serve a named block.
Report VDJdb cohorts by template coverage, never a pooled strata average. Keep the two TCRvdb
epitopes separate and compared arms on identical rows. Do not confuse deposited-model ranking
with fixed-template threading. Keep the fit-free predecessor scores alongside the frozen readouts.

`features` → `assess` is the frozen-model route; `recognize` owns the predecessor Q/T/S route.
Nothing is fitted at scoring time. Training dependencies, statistical assumptions and physical
interpretations must describe the implementation: see the unresolved tier/confidence/ergodicity
questions in `AUDIT.md`. A dimensionless contact log-odds is not automatically a binding free energy.

- Scope is TRA+TRB, a standard-amino-acid peptide and an MHC-I/II groove. Drop out-of-scope cases.
- IMGT CDR3 excludes junction anchors. Recover Cys104 and the actual F/W118 from annotated
  FR3/FR4; never construct a junction with unconditional `"C" + cdr3 + "F"`.
- `pitch_angle` is generator-confidence leakage and must not enter the structural score set.
- Define every equation's symbols, domain and units. Distinguish exact identities, approximations
  and measured associations. Preserve values when correcting a misleading mathematical label;
  catalogue changes also require a deliberate model/provenance migration.
- `m_erank_*` is the order-2 Hill effective rank, not the order-1 entropy effective rank.
- Class-II crystal evidence does not establish modelled class-II benchmark performance.

## Implementation and verification

Use existing APIs first. Batch annotation/search once per set (or memory-bounded chunk), never a
search per item or an outer pool around an internally threaded search. Inspect callers before edits.
Validate shapes, dimensions and indices at native boundaries before dereferencing pointers.
All data roots go through `paths.tcren_home()`; all-failed commands exit nonzero with diagnostics.

`./setup.sh` creates the local venv and editable install; `./setup.sh --tests` runs fast tests.
Use `.venv/bin/python -m pytest`, focused checks first, then `-m 'not slow'` with a bounded runtime.
Do not run scientific refit/regression stages as incidental audit checks. Rebuild after native edits.
Lint gate: `.venv/bin/ruff check src tests`. Docs gate: activated venv, `make -C docs html SPHINXOPTS=-W`.
New public capabilities update the API skill and docs together. Update CHANGELOG before release.
Do not edit `appendix/generated/*`: `scripts/gen_appendix.py` owns those tables, and
`gen_descriptor_table.py` / `gen_family_graph.py` own their generated catalogues.

`bootstrap.sh` prepares environment/raw inputs; `reproduce.sh` coordinates the three repositories.
An explicit subset is diagnostic, never complete reproduction. No `--skip-slow` reproduction flag.
A local test does not prove CI, publishing or remote installation.

## Retired work

`attic/2026-10-audit/` preserves original instructions/status, the retired geometry/chain/cross
experiment, resources and tests, with checksums. It is historical evidence, not another instruction
source or supported package. Do not restore discarded experiments merely because they exist there.
Retain the contact model, self-consistent regeneration and peptide–MHC MD model: current consumers
use them. `OBSOLETE.md` records the dependency-based disposition; absence from one descriptor family
alone is not evidence that code is unused.

## Contact energy follow-up — 20 October 2026

The author authorized the complete energy-score audit across all tasks, including
poorly performing tasks and legacy 2022 benchmarks, plus sparse atom-potential
research with antibody transfer and carefully scoped exchange symmetry.
`CONTACT_ENERGY_SUBROADMAP.md` owns this subgoal and its completion gates.
Retain original TCRen, identical populations and all gains/losses; fit no assay
weights. The earlier 3808-model ranking comparison does not certify all-task or
legacy performance. Use all allocated CPUs for parallelizable work and fix thread
explosions rather than retreating to serial; keep one bounded compute layer.
