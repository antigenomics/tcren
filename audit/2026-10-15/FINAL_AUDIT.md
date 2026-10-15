# Final audit reconciliation — 15 October 2026

The author selected fit-free final models, archival of fitted readouts, current packages, CLI-owned calculations, and explicit reporting of gains and regressions. The current task stops at this audit record and filed issues, before follow-up library fixes or a replacement paper.

## Evidence and publication status

No complete successful raw-input or CLI-only reproduction is claimed. Attempts 8 and 9 stopped before publication. Attempt 8 completed fresh binder/pose/contact joins but found an obsolete raw-reference alias; that exact local input path was repaired and guarded. Attempt 9 was stopped when a matched record comparison exposed regressions; the author subsequently permitted marked regressions. All 810 calculation/figure sources stayed unchanged during that run. Source/log/runtime/RSS records are under `reproduce/attempt-9-*`.

The current main/SI are visibly marked historical audit drafts. Gaussian readouts, confidence GLMs, fitted transition products and fitted Potts final scorers are retired from the final specification. Their numerical evidence is preserved. The manuscript's `FIT_FREE_DISPOSITION.md` maps all 27 historical figures and every estimated dependency. Its checked archive contains 504 files / 36,476,442 bytes; all hashes match, and numerical producers/outputs were unchanged by disposition. Library fitted bundles and sources have a separate checksummed snapshot in `attic/2026-10-fit-free-disposition/`.

Earlier cleanup moved 40 unused/duplicate assets with exact hashes. Trees and the withdrawn HPV posterior panel were archived. The 12 selected structure scenes passed a temporary runtime/visual check, but direct Python rendering is not a CLI-only reproduction certificate. Known height and pixel-provenance differences remain explicit in the figure ledger.

## Matched recorded gains and regressions

These are predecessor/historical-arm comparisons from attempt 9, using the repaired local tcren wheel with its then-locked runtime. They are not results from a new final fit-free model, nor a full cohort rerun under the later package update. ROC-AUC and PR-AUC are dimensionless. Binder rows refer to the benchmark candidate `QxdPhi_TCR (Eq7)`; TCRvdb rows refer to `Q x dPhi_tcr (coupled)`. Macro ROC is the arithmetic mean of within-cohort ROC values in the named stratum, with no pooled-stratum average.

| Population / quantity | Observations (positive / negative) | Before | After | Delta |
|---|---:|---:|---:|---:|
| VDJdb templated S macro ROC-AUC, 6 cohorts | 297 (142 / 155) | 0.734546 | **0.742799** | +0.008253 |
| VDJdb template-free S macro ROC-AUC, 16 cohorts | 792 (381 / 411) | **0.510021** | 0.506157 | −0.003864 |
| CLGGLLTMV S ROC-AUC | 46 (23 / 23) | **0.512287** | 0.463138 | −0.049149 |
| LLYDANYFL S ROC-AUC | 43 (19 / 24) | **0.425439** | 0.394737 | −0.030702 |
| YLQPRTFLL coupled S ROC-AUC | 423 (174 / 249) | 0.8544 | **0.8557** | +0.0013 |
| YLQPRTFLL coupled S PR-AUC | 423 (174 / 249) | 0.7746 | **0.7760** | +0.0014 |
| GLCTLVAML coupled S ROC-AUC | 195 (135 / 60) | **0.7495** | **0.7495** | 0 |
| GLCTLVAML coupled S PR-AUC | 195 (135 / 60) | 0.8590 | **0.8591** | +0.0001 |

The compared key sets, counts, schemas and finite status were retained. Q and confidence ranking records were unchanged. Complete recorded deltas are in `/tmp/tcren-attempt9-completed-record-comparison/`; the retained manifest is `reproduce/attempt-9-matched-record-manifest.json`. Component and region-label diagnostics remain in `/tmp/tcren-binder-component-deltas.csv` and `/tmp/tcren-binder-contact-region-deltas.csv`. These diagnostics were not used as cached inputs to a subsequent scientific run.

The benchmark historical-report driver now records all keyed numeric values, finite flags, original floor-test exit/failure details and regression flags. Metric losses do not silently disappear or abort solely because they are losses; schema/key/count/support and finite-availability errors remain hard failures. Original thresholds/floors were not changed. This reporting path is historical evidence, not final-model selection.

## CLI and mathematical discrepancies filed for follow-up

| Issue | Verified discrepancy / required reconciliation |
|---|---|
| [30](https://github.com/antigenomics/tcren/issues/30) | Historical annotation/potential lineage is incomplete. Region labels changed for 643 of 26,868 matched contacted receptor residues; archive artifacts span different potential generations. Human-only annotation is forced in one mixed-species native producer. |
| [31](https://github.com/antigenomics/tcren/issues/31) | Exact virtual receptor-side CDR3/CDR12/all poly-Ala referencing is implemented outside the CLI. Peptide `--delta` and structural single-residue scans are different quantities. The direct helper also silently drops exceptions without complete-key accounting. |
| [32](https://github.com/antigenomics/tcren/issues/32) | A direct CPL Potts producer ignores custom model radius/cutoff. Uncoupled log Z₀ and first-order fixed-marginal effects need accurate approximation labels and equivalent CLI semantics. |
| [33](https://github.com/antigenomics/tcren/issues/33) | Retire fitted final readouts and make the final simple formula/parameter boundary explicit. Native-Q covariance weights, cohort correlation/standardization and empirical potentials require honest dependency descriptions. |
| [34](https://github.com/antigenomics/tcren/issues/34) | Exact position/residue decompositions and selected render inputs have no equivalent complete CLI export. Leave-one-out `explain` is not the same quantity. |
| [35](https://github.com/antigenomics/tcren/issues/35) | Benchmark-local grouped-binomial Potts and Gaussian fits are independent models, not equivalent frontends to library CLI algorithms. Preserve historical evidence; do not silently migrate by command-name substitution. |
| [36](https://github.com/antigenomics/tcren/issues/36) | Current Polars 2 rejects the Float64-label / Int64-list `is_in` validation expression before the intended library error. One existing fast test fails. No dependency downgrade or fix was made in this stopping phase. |

The observed smoothed-reference invalid-arithmetic warning is already tracked in [17](https://github.com/antigenomics/tcren/issues/17). Existing assessment-availability and Gaussian numerical issues are [18](https://github.com/antigenomics/tcren/issues/18) and [19](https://github.com/antigenomics/tcren/issues/19); they were not duplicated.

Earlier local audit repairs cover Potts physical-chain identity and final-state warm starts, orientation coefficient accounting, IRLS/domain checks, compact-Cα mapping, native buffer validation, constrained shear geometry, and deliberate metadata/unit migration. Their earlier tests and exact array-preservation records remain in AUDIT.md. The mathematical descriptions now distinguish reduced spring units, truncated Betti difference, linearized contact-loss surrogates, uncoupled marginal approximations, and heuristic confidence combinations. Fitted transition widths are score-scale parameters, not calibrated free energies or absolute occupancies. These are local checked changes/definitions, not a release or proof that every possible input is correct.

## Current-package verification

Both isolated runtimes now use ARDA 2.36.1, NumPy 2.5.4, Polars 2.0.0, SciPy 1.18.1 and current compatible runtime packages. Exact version records are `current-packages-library.json`, `current-packages-benchmark.json`, and the two `reproduce/current-package-audit-*.lock` files. The immutable local tcren 3.10.0 wheel remains SHA256 `17768f46d92a247d5864290c34470918edc6a2a0a37b1ff6d0fa957741e7c6fb`; all 169 installed payload files still match it. Historical frozen arrays were not overwritten.

- Current-package CLI `scoring --delta` processed all 374 raw Native2026 structures and wrote `/tmp/tcren-cli-audit-current/native-scoring.tsv`. This validates that route, not every publication producer.
- Existing bounded fast unit audit: **989 passed, 1 failed, 129 deselected**, 198.84 s, 33 warnings. The failure is issue 36 above. Command: `TCREN_HOME=<library> TCREN_DATA_DIR=<raw mirror> python -m pytest tests/unit -m 'not slow' -k 'not refitting and not reproduc' -q --tb=short -p no:cacheprovider`.
- Current-package Potts ensemble checks: **18 passed**, 7.61 s.
- Source Ruff and diff whitespace checks passed. Benchmark focused checks passed **24 tests**, 3.49 s, with warnings treated as errors. Manuscript read-only gates passed after notice-only wording was made compatible with the unchanged corpus gate: 18 main labels, 31 note sections, 45 supplementary floats, 53 plot scripts and 28 existing PDFs; no numerical body or vocabulary-gate rule changed. A read-only gate never certifies fresh PDFs.

No follow-up tcren fixes, new fitted models, new scientific controls, threshold changes, replacement paper or release follow from this record. The next implementation phase starts from the filed discrepancies and the explicit fit-free specification.

## Visible integration checkpoints

Library audit/disposition: `6d5526e`; benchmark historical-report source: `c7143ad`; manuscript archive/status: `ee071f8`. These are local commits in the visible integration checkouts. No remote push or release was performed. The committed blob contents of all 9 library archive files and all 504 manuscript snapshot files were independently SHA-256 verified; see `committed-archive-verification.json`. Historical archived originals intentionally retain their original line endings and whitespace; active source/status changes pass whitespace checks.
