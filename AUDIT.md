# TCRen audit and reconciliation — 2026-10-14

## Final author-directed audit disposition — 2026-10-15

The author permits explicitly marked metric regressions, requires current packages and CLI-owned
calculations, and selected **fit-free final models with fitted readouts archived**. That audit assignment
ended at the reconciled issue record. The subsequent author request authorized library fixes,
fresh fit-free reproduction and release preparation; see [release fixes](audit/2026-10-15/RELEASE_FIXES.md).
Replacement-paper work remains outside this follow-up.
This supersedes the earlier frozen-final-model/full-publication promotion plan; historical results
and their coefficients remain preserved as evidence. See
[the final audit record](audit/2026-10-15/FINAL_AUDIT.md),
[the historical fitted-readout snapshot](attic/2026-10-fit-free-disposition/MANIFEST.json), and the
manuscript's `FIT_FREE_DISPOSITION.md` and 504-file verified publication archive.

New issues are filed in antigenomics/tcren: [annotation/potential lineage](https://github.com/antigenomics/tcren/issues/30),
[receptor-reference CLI/accounting](https://github.com/antigenomics/tcren/issues/31),
[Potts geometry/approximation semantics](https://github.com/antigenomics/tcren/issues/32),
[fit-free final-model disposition](https://github.com/antigenomics/tcren/issues/33),
[exact residue/render CLI exports](https://github.com/antigenomics/tcren/issues/34), and
[independent local fitted algorithms](https://github.com/antigenomics/tcren/issues/35).
The observed smoothed-reference arithmetic warning is already tracked by
[issue 17](https://github.com/antigenomics/tcren/issues/17).

Attempt 9 was stopped before publication when a matched comparison first exposed binder losses.
The author's later directive explicitly permits these marked losses. Templated predecessor S macro
ROC increased 0.734546→0.742799 (6 cohorts;297 structures,142 positive/155 negative); template-free
macro ROC decreased 0.510021→0.506157 (16 cohorts;792 structures,381 positive/411 negative).
CLGGLLTMV ROC decreased 0.512287→0.463138 (46 structures,23/23); LLYDANYFL decreased 0.425439→0.394737
(43 structures,19/24). YLQPRTFLL coupled S ROC increased 0.8544→0.8557 and PR-AUC0.7746→0.7760
(423 structures,174/249). GLCTLVAML ROC stayed 0.7495, PR-AUC0.8590→0.8591 (195 structures,135/60).
Q and generator-confidence ranking records were unchanged. These are matched **historical-arm**
comparisons from attempt 9, not new fit-free final-model claims or a complete reproduction.

The historical region partition changed 643 of 26,868 matched contacted receptor residues; the
archive also mixes potential generations. An actual raw-input diagnostic returned identical
receptor-reference values under current extant ARDA 2.5.7 and 2.20.0, ruling out that version change
as a demonstrated explanation. No guessed region masks or old coefficients were restored.
The latest package audit and CLI/software verification are recorded in the final audit record.

## Resolution checkpoint — 2026-10-15

The author authorized blocker repairs and a full fresh reproduction, then explicitly approved
retiring tree figures/claims and correcting ergodicity and naive-Bayes interpretations while
preserving empirical values. The earlier inventory below records the entry state, not the current
implementation status. Reproduction and publication certification remain open until the complete
driver and matched-result gates finish.

- Potts sites/kernel/marginal joins retain physical receptor and partner chain identity. MHC tables
  lacking partner identity require raw regeneration; single-peptide legacy tables retain compatibility.
  Threaded candidates require the peptide arm and exact full length, and carry all final Gibbs states.
  Orientation coefficient accounting, hybrid-imputation convergence, probe provenance, sampler
  domains, IRLS convergence/backtracking and final Hessian are repaired. A fresh run additionally
  exposed compact-Cα orientation indexing; missing geometry is preserved and rejected by oriented fits.
- Definitions now distinguish reduced spring units and finite-grid work, the truncated Betti
  difference, linearized nearest-pair contact loss, heuristic confidence composition, transformed-axis
  conditional means and shared preprocessing fitted on both classes. Future Gaussian/Schur-complement
  and mechanics inputs have explicit domain checks. No SI calibration or likelihood-EM theorem is claimed.
- The two-step peptide link uses binder-oriented referenced scores and a conditional ranking fit.
  Its fitted sigmoid transition widths are not calibrated free-energy conversions or absolute
  occupancies; one width changes a factor's log-odds by one nat. Pooled parameter estimates,
  structure-grouped held-out performance and likelihood-ratio versus cluster-robust Wald diagnostics
  are now distinguished. All standing values and result macros were preserved. The manuscript and
  `FINDINGS.md` originals and numeric-preservation proofs are archived in their follow-up directory.
- Frozen metadata was deliberately migrated; all five coefficient arrays and original numerical fit
  metadata are unchanged. [Migration evidence](audit/2026-10-14/frozen-metadata-migration.json)
  retains original/new bundle hashes and catalogue digests. Corrected-feature compatibility still
  requires fresh comparison. Metadata migration is not a refit or numerical equivalence certificate.
- Publication producers, dependency order, removal of cached/resumed calculations and tree retirement
  are implemented and source-checked across the benchmark and manuscript; end-to-end reproduction
  remains pending. The recovered topology producer comes from the
  recorded repository history. Historical paragraphs/data and dirty edits were archived before edits.
- The immutable wheel and both reproduction environments share one payload; preflight verified
  all 169 installed library payload files and all 12 required raw sets. The wheel SHA-256 is
  `17768f46d92a247d5864290c34470918edc6a2a0a37b1ff6d0fa957741e7c6fb`.
  This local artifact retains version 3.10.0 and is not a release or remote CI assertion.
- Applicable fast unit run: 981 passed, 1 network-dependent skip, 129 deselected, 182.87 s.
  After the missing-Cα repair, all 123 affected Potts checks passed and installed-wheel boundary
  checks passed 22 tests. The final installed-wheel fast run passed 987 tests, with one skip and
  one environment failure: the pipeline fixture looked in the installed package's default cache
  instead of the supplied canonical mirror. That exact test passed with the reproduction's
  explicit data roots (14.74 s). A subsequent connected-core selection regression failed before
  repair; its affected source/pipeline checks passed 12 tests, and the final installed wheel passed
  all 10 selected boundary checks. Ruff, diff and Sphinx `-W` checks passed. Native extensions
  built and imported.
- Full-run attempts exposed and repaired the multiword recipe command bug and compact-Cα indexing.
  The next attempt passed all six potential checks, including the sampled recipe (246.80 s), and
  freshly generated 7,816 derivation contacts over 358 contributing structures and scores for all
  374 native inputs. The benchmark was explicitly terminated before obsolete parsed driver
  functions could run cached calculations or retired tree fits. A full restart after those driver
  repairs was started from raw inputs with the final immutable wheel. Its six potential checks
  passed again (241.52 s), then stopped at an obsolete binder `recognize --scores` call.
  Earlier TCRvdb processing had reused a legacy descriptor table and therefore does not establish
  fresh reproduction. That path, withdrawn-score dependencies and two additional physical-chain/order
  joins in CPL consumers were repaired before subsequent complete restarts.
  A reporting edit to the running shell driver also caused a terminal parser error; the next
  driver's source will remain fixed during execution. No publication outputs were promoted.
- The subsequent full attempt passed native scoring, dynamics inputs, ATLAS, potency,
  neoepitopes, legacy, fresh TCRvdb and CPL processing, then stopped on three auxiliary CDR1/2
  reference fields absent from the binder catalogue. All 1,089 binder rows were retained.
  These fields are now explicitly retained outside the established scan/design; no new arms
  were introduced. Parent-shell raw-path exports and clean-output creation are also repaired.
- The author additionally authorized unused panel cleanup and a reproducible live figure set.
  Forty legacy/duplicate assets were moved with verified SHA-256 preservation. All 27 live
  figures have source/build entries; four caption maps and statistical-unit/marker labels were
  reconciled without changing values. The selected 12 structure scenes passed temporary runtime
  verification in 89.6 s. Four geometry-height diagnostics changed under the corrected library;
  the [complete deltas](audit/2026-10-14/render-verification-deltas.json) were reported and remain
  flagged before canonical promotion.
  One attempt was stopped explicitly before further computation while these missing render
  edges were repaired. The existing late CPL heatmap producer was confirmed and given a complete
  260-cell output check. Thirteen focused producer/identity checks and the read-only manuscript
  gates pass. Attempt 7 passed the fresh binder recognition/mechanics/reference/scan exports,
  then stopped at a pose-table join retaining only the 523 positives. All 1,089 structures had fresh
  recognition and pose outputs. Negative metadata IDs retained suffixes absent from the descriptor
  key; the same metadata also carried historical computed fields that could precede fresh join
  columns. Exact key coverage and explicit raw-only metadata projection are repaired and checked;
  diagnostic joins retain all 566 negatives and 523 positives without rerunning descriptors.
  The stopped run took 3,405 s; sampled summed process RSS peaked at 18,660,065,280 bytes (shared
  memory may be counted in more than one process). Full fresh reproduction, matched-result gates
  and rendered publication QA remain open.

- The HPV producer still called a withdrawn cohort posterior and an old energy column. The author’s
  legacy-panel disposition applies to that already-retired panel: original three-panel source/PDF/data
  and benchmark outputs are archived; the live route retains only the existing Q and referenced-energy
  arms. All 32 retained display macros are exact in a temporary layout check. The 4 Native2026 crystals,
  56 Shcherbinin mutants and 21 original archive decoys keep their original coordinate sources. Raw
  metadata supplies the exact label/path mapping; a small source-owned identifier-order input preserves
  seeded report/bootstrap ordering without storing computed features. Matched point/statistic gates
  protect both retained arms. The nested all-positive/crystal-core contrast is labelled descriptive.
- The existing native TCR-reference producer is now scheduled before QTF consumes it. Nineteen combined
  producer/identity/preservation checks pass, selected changed sources pass Ruff, and wheel/raw/render
  preflight passes. The temporary two-panel HPV plot renders correctly, and manuscript read-only gates
  pass. These are source/layout checks; another whole raw-input run remains necessary.

- Attempt 8 completed the fresh binder pose join with all 1,089 rows (566 negatives/523 positives)
  and both contact-geometry cohorts (1,089/618), then stopped after 3,980 s at an obsolete VDJdb
  reference alias. The documented local release is byte-identical to the pinned archive member;
  its existing path and exact SHA-256 guard now replace the absent alias, with unchanged counting.
  No download or release substitution was needed. All 810 calculation/figure source files stayed
  unchanged during the run. Sampled summed process RSS peaked at 18,648,170,496 bytes (shared
  pages may be counted more than once). A read-only audit of the remaining scheduled raw inputs
  found no further missing mandatory references. Another complete fresh restart is required.

The [publication reconciliation](audit/2026-10-14/publication-resolution.md) inventories all 34
live reshapers, recovered producers, exact archived tree routes, fresh-input/cache repairs and
matched-result gates. The full driver's finite-score report now counts both null and NaN values;
null counts alone cannot describe missing numeric scores. Native panels may retain diagnostic
rows with missing geometry; raw input counts do not imply every score is finite.

`scripts/check_fresh_holdout.py` now requires all 8,292 manifest structures, fits only a temporary
output, records numerical differences and retains the original exact full-scale identity contract.
It is a mandatory gate after fresh harmonization and before publication production.

## Original audit outcome and evidence boundary

This audit found and repaired concrete numerical, data-selection, native-buffer and reproduction
failures. The library, benchmark and manuscript now have canonical instructions and an ordered
publication roadmap. Retired geometry experiments were moved as a dependency group into attics.
**This is not a certification that every algorithm or equation is correct, nor a successful
scientific reproduction.** The remaining publication blockers below prevent that claim.

All three checkouts were already dirty. Existing scientific outputs, potential arrays and frozen
model arrays were preserved. Only synthetic numerical checks and software tests were run; there
was no corpus refit, benchmark rerun, new scientific control or release. Some repaired descriptors
can change future outputs: unchanged arrays are preservation, not proof of compatibility.

`ROADMAP.md` contains the exit criteria. The dated review notes below give precise coverage,
counterexamples, tests and unresolved domains. Their chronological observations are superseded by
the final dispositions here; they are not independent release approvals.

## Original publication blocker inventory, in priority order

| Priority | Finding and owning source | Required disposition |
|---|---|---|
| P1 | `score/fit.py` fits nonlinear preprocessing on both classes before estimating binder moments. The tier-1 “no negatives/no label” claim is not true of the training pipeline. | Describe the actual dependency, or approve a separate transform/model revision and matched reproduction. No silent refit. |
| P1 | Corrected class-II chain identity and constrained shear direction change descriptor semantics. Catalogue units/labels also need correction. | Deliberate model/catalogue provenance migration and assessment of affected outputs before release. Current digest equality cannot certify unchanged numerical semantics. |
| P1 | Benchmark `requirements.txt` pins 3.0.1 while producers use later APIs; confidence GLM/correction/operating-point producers are absent from the driver. | One immutable tested library identity, complete producer edges, and explicit ownership of benchmark GLM versus library conditional residual. |
| P1 | Supplementary panel macros still pool template strata; GLM gain is described as against the better term but is computed against the structural term. Fits are in-sample in that route. | Reconcile intended comparison/population and prose with the existing stratified route. Preserve all recorded results; do not silently replace macros. |
| P1 | Mechanics catalogue assigns N/m, N and J to uncalibrated spring surrogates. Default stiffness is multiplicity divided by distance squared in Å. | Correct to reduced units through a deliberate catalogue/model migration; do not invent an SI conversion. |
| P1 | Cross-complex energy correspondence is described as a direct ergodicity test; kinetics correlation is given a stronger mechanistic interpretation than the computation establishes. | Author decision on inference scope. Recorded R²=0.460 over 184 complexes and R=−0.5398 over 23 complexes remain unchanged. |
| P2 | `fp_chi` computes β₀−β₁, not the Euler characteristic of an arbitrary flag complex. Six octahedron vertices give β₀=1, β₁=0, β₂=1: output 1, Euler characteristic 2. | Name the truncated Betti difference or explicitly assume vanishing higher homology; preserve historical values and migrate provenance deliberately. |
| P2 | Contact-loss score is a linearized nearest-atom-pair surrogate, not an exact residue-contact probability. | Source doc corrected; reconcile catalogue and manuscript definitions without silently replacing values. |
| P2 | `binder_score + logit(ipTM)` has no established naive-Bayes interpretation; ipTM is not demonstrated binder likelihood-ratio evidence. | Retain any measured performance as a heuristic result; settle probability language and independence assumptions. |
| P2 | Confidence back-transformation uses sigmoid of a conditional mean, which is not generally the conditional mean of sigmoid. | Distinguish transformed-axis residual, plug-in original-scale prediction and conditional expectation. |
| P1 | Active Potts sites independently drop partner-chain identity; class-II sites can conflate chains. Threaded ensemble restarts after candidate one with abbreviated burn-in rather than its advertised warm start. | Explicit site/model migration and sampling-policy review before affected readouts are promoted. See the Potts audit; neither was silently changed. |
| P2 | Potts parameter counts/convergence omit orientation coefficients; some EM edge cases and weighted-diagnostic labels are inconsistent. | Close the named implementation/definition gaps before relying on these fit diagnostics. Full list and boundaries in the Potts audit. |
| P2 | Coarse Potts `p_any` applies an independence union formula to coupled marginals. | Documentation now identifies the approximation; retain existing values and decide whether publication uses it as a heuristic or requires joint-sample occupancy. |
| P2 | Rupture's finite displacement horizon does not guarantee all springs break; direction/scalar domains also need tighter validation. | Define the actual finite-grid observable or approve a corrected horizon with model/results review. No unnoticed numerical change. |
| P2 | Confidence cross-covariance plus shrunk descriptor covariance does not guarantee a PSD joint block for arbitrary future fits. | Validate the Schur complement for future fitting. The shipped conditional variance is positive, 0.0477092387; no shipped failure was found. |

The author must resolve the scientific interpretation choices before publication prose or accepted
results are changed. This audit does not propose additional controls, sweeps or metrics as a default.

## Confirmed repairs

| Area | Defect repaired | Evidence |
|---|---|---|
| Missing scores | All-missing Q/T blocks became zero; constant standardization erased missingness. | Synthetic missing rows now remain NaN. |
| Marginal scores | Missing unrelated descriptors suppressed receptor/channel scores. | Coordinate selection now precedes missing-row selection; excluded-column null/drop agree. |
| Score contracts | Legacy signs ignored; catalogue digest unchecked; missing confidence fit returned a spurious residual. | Sign, stale-model and explicit confidence-availability regressions. |
| Fitting inputs | Duplicate/null IDs and invalid labels could silently alter the joined population. | Unique keys, binary labels, 1:1 join and omitted/incomplete counts; synthetic fitter cases only. |
| Transform | Direct Transformer use converted infinities to finite clipped values. | All nine transform families preserve nonfinite input as missing. |
| Attribution | Selecting peptide residues silently switched scoring mode. | Separate `receptor` keyword; selected residues retain the same full-score baseline. |
| Potts numerics | Coupled marginals were applied in the wrong row order; uncentred pins reloaded centred; tiny log-partitions cancelled to zero. | Identity joins, exact pinned-array serialization and stable `expm1`; 13 focused tests. |
| Contact typing/batching | Mouse autodetection repeated mouse searches; archives fell back to per-item searches; v2 stacking orientation and other flags were inconsistent. | Shared batch API, canonical pair keys and synthetic regression checks; frozen v1 classifier unchanged. |
| Native alignment | Affine traceback lost its gap state even when the optimum score was right. | Score-34 synthetic alignment previously traced score 30; exhaustive 961 short pairs and 500 seeded longer pairs checked. |
| Native buffers | Five extension families had unchecked array dimensions, indices or domains before raw access. | Focused malformed-input regressions; finite valid-input checks; extensions rebuilt. |
| MHC reference | Marker files made a missing ~915 MiB prefilter data file look ready. | Component/length checks, staged rebuild, local backup and real 1AO7 annotation smoke. Raw FASTA unchanged. |
| Mechanics | Projecting the leading 3D eigenvector did not maximize stiffness in the shear plane. | Solve the 2D constrained eigenproblem; diagonal tensor counterexample repaired. |
| Class II | MHC chains with equal residue numbers were merged. | Preserve target chain in deduplication/participation and normalize contact order per chain. |
| Surfaces | Missing distances became zero similarity; incompatible hydropathy scales could be compared. | Reject undefined distance and invalid/mixed chemistry scales. |
| Features | All-error commands succeeded; computed chunk intermediates could be reused. | Nonzero failure with diagnostics; resume removed and batches computed afresh in memory. |
| Paired metrics | Bootstrap arms filtered missing observations independently. | Shared finite rows before resampling, including group labels; durable synthetic regression in benchmark tests. |
| Reproduction | Relative interpreter broke after chdir; consumers preceded producers; caches/skips and benchmark-side library model overwrites violated the contract. | Absolute interpreter, repaired known dependency edges, fail-fast fresh features, removed fitting/install side effects and skip-slow flag; static checks. |

## Mathematical validation and qualifications

The Gaussian class log-odds, marginalization and conditional mean algebra are consistent on the
valid positive-definite domain. Reliability-weight covariance agrees with `numpy.cov(aweights=...)`
on a synthetic example. Yeo–Johnson branches agree with SciPy over seven signed values and six
lambdas (maximum absolute difference 1.82×10⁻¹¹). These are focused checks, not universal proofs of
numerical stability. Direct `Joint` validation now guards the identified weights, class support,
shapes, descriptor subsets, widths and positive-definite covariance/Schur domains. These checks do
not establish stability for every possible floating-point input.

For unit contact directions uᵢ and nonnegative spring constants kᵢ, K=Σᵢ kᵢuᵢuᵢᵀ is symmetric PSD
and tr(K)=Σᵢkᵢ. In-plane maximum stiffness is the largest eigenvalue of UᵀKU, where U has two
orthonormal columns perpendicular to the pulling axis. The corrected implementation follows this
constrained Rayleigh quotient. It does not supply an energy calibration or physical rupture law.

Participation 1−Σₘpₘ² and degree evenness (Σᵢdᵢ)²/(nΣᵢdᵢ²) match their definitions for module
fractions pₘ and degrees dᵢ of n selected residues. The literature contact-order helper is an adapted
interchain statistic normalized by contacted span, not the original full-chain intrachain statistic.
The footprint H₁ computation needs triangle-boundary rank and is correct at that order; this does
not make β₀−β₁ a full Euler characteristic. Order-2 Hill effective rank remains distinct from the
order-1 entropy effective rank.

For one atom pair at initial distance d>0, cutoff R and isotropic displacement of fixed magnitude
δ>0, the exact loss probability is clip(((d+δ)²−R²)/(4dδ),0,1). The implemented linearization is
clip((δ−(R−d))/(2δ),0,1). At d=4.5 Å, R=5 Å and δ=1 Å these are 0.2916667 and 0.25. Neither
single-pair formula alone is the loss probability of a multiatom residue contact. No values changed.

Fixed pseudocounts prevent invariance to global scaling of weighted contact counts: two synthetic
contact types changed the resulting matrix by up to 1.0352 when weights were multiplied by ten.
Scale invariance requires zero pseudocount or scaling it too. Documentation now states this.
Independent sidechain repacking minimizes each frozen-neighbour objective; combining the choices
does not guarantee a globally improved assembled pose. Proline construction is an unconstrained
open-chain approximation, not a validated closed ring. These limitations were clarified, not hidden.

The annotation backend now propagates thread-cap failures, primary-complex selection retains only
the connected TCR:pMHC core, and degenerate PCA frame inputs fail explicitly. Alternate-conformer
union versus writer deduplication remains a documented qualification. Frozen hydrogen-bond labels
remain distance proxies, not full donor/acceptor-angle or DSSP-energy determinations. The original
source review did not infer corpus effects; the subsequent selected-render check measured the four
height deltas recorded in the resolution checkpoint above.

## Cleanup and retained scope

The library attic contains the geometry prior, chain generator, cross-product module, their three
resources and three tests (nine moved files), plus original instruction/status/CLI snapshots.
The benchmark attic contains eleven dependent geometry exploration scripts and snapshots of
changed drivers/instructions. Manuscript instruction history is also retained. Checksummed path
maps make the moves reversible. The library attic is excluded from distribution builds.

Live contact-map fitting, self-consistent potential regeneration, peptide–MHC MD models, shared
geometry primitives and historical oracle paths remain: identified consumers still use them.
Scientific tables and unused-looking research outputs were not deleted. No module was archived
merely because it was absent from one descriptor family. Broad removal of every heuristic would
break published producers; the roadmap instead freezes expansion and names the remaining owners.

## Review coverage and remaining limits

Detailed notes: [scores](audit/2026-10-14/score-audit.md),
[native/topology/potentials](audit/2026-10-14/native-topology-audit.md),
[additional mathematics](audit/2026-10-14/additional-math-audit.md),
[benchmark/manuscript evidence](audit/2026-10-14/evidence-audit.md), and
[all 35 publication reshapers](audit/2026-10-14/publication-dag.md),
[Potts mathematics](audit/2026-10-14/potts-math-audit.md), and
[structure/annotation/contact flow](audit/2026-10-14/structure-audit.md).

Complete-file reviews cover score/model/transform/fit/explain, cohort and reliability, all five native
kernels, potential derivation/redundancy, topology graph/footprint/surface/literature/pose and spring
mechanics, the active Potts model/fit/kernel/sample/sites/score/ensemble/EM/selfcon stack, and
contact/annotation/docking-frame construction; related descriptor, CLI, MHC, dynamics and benchmark callers were traced. The publication
review traced important macro producers and contradictions; it did not independently verify every
bibliographic entry, chemistry constant, dataset label, raw-source availability or manuscript equation.
All native caller domains, every annotation/refinement module and all scientific dependency edges
are not exhaustively certified. Remaining direct-input gaps are listed in the detailed notes.

## Original 14 October verification checkpoint

- Applicable unit suite after the final source changes: **948 passed, 1 failed, 1 skipped
  (RCSB unreachable), 129 deselected**, 175.03 s. The failure was an obsolete storage-contract
  assertion requiring an uncentred named pin to reload through the centred-name representation.
  The test now explicitly checks centred named pins; the new raw-pin regression checks exact
  matrix/JSON preservation. The complete affected Potts/self-consistency suites then passed
  **78 tests**, 19.22 s. No numerical floor or scientific result was changed to pass this check.
  The preceding broad run passed 938 tests before the final extensions to coverage.
  Command: `.venv/bin/python -m pytest tests/unit -m 'not slow' -k 'not refitting and not reproduc'
  -q --tb=short -p no:cacheprovider`. The suite emitted 33 warnings, including known unannotated-MHC
  fixture warnings and Biopython deprecations; it was not a warnings-as-errors run.
- Follow-up score/surface checks after final guards: **36 passed, 20 deselected**, 4.57 s.
- `ruff check src tests` and `git diff --check`: passed.
- Sphinx HTML with `SPHINXOPTS=-W`: passed.
- All five native extensions rebuilt; wheel built and installed in a temporary venv without
  dependency installation, using this checkout's existing dependency environment for the smoke test.
  Imported the final installed package and native modules, loaded the frozen model and 8,292-row manifest,
  and verified retired geometry/attic exclusion. The source distribution also built and contains active source/resources while excluding audit/attic.
  This is an artifact smoke, not a clean dependency solve.
- Benchmark paired-row regression and shell/static driver checks: passed. No scientific evaluation.
- Existing manuscript check-only gates: **passed**, 3.78 s; 55 figure scripts and 29 existing PDFs,
  zero gate violations. Command `make MAKE='make -o all' gate` deliberately suppresses the ordinary
  figure-rebuild prerequisite. No PDF compilation or visual page review is claimed.
- Archive verification: **40 archived files** matched SHA-256 across the original and supplementary
  metric/orphan manifests. Orphan geometry reads were removed from the live reshaper and
  panel-stat driver after checking they had no live output use. Frozen-object metadata export now
  reads public APIs and matches the existing columns without fitting.
The first applicable unit run failed because the MHC index data file was missing (817 passed,
18 failed, 67 errors); the root cause was repaired before the final run. Model-refitting and
scientific-reproduction tests are intentionally excluded, as are tests marked slow.

Original 14 October wheel SHA-256: `142f38cdcab2192c454ca157c796c6294684e12de2372df1ee960fee40c2f9fb`.
It is an unreleased macOS arm64/Python 3.12 build carrying the checkout’s existing 3.10.0 metadata;
it does not establish PyPI or CI state.

Integration scope: audited library source and supporting pre-existing sidechain implementation were
reconciled together; unrelated generated appendix edits remain in the working tree. Benchmark and
manuscript commits contain the selected drivers, source cleanup and instructions; their pre-existing
scientific outputs and other ongoing edits remain uncommitted. No remote was pushed.
