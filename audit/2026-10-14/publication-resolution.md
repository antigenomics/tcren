# Publication producer reconciliation — 2026-10-14

This is a source/graph reconciliation record. It does not certify the running fresh reproduction,
new benchmark values, rendered PDFs or publication readiness. The parent orchestration owns the
fresh-run status and final wheel digest. All changes below preserve recorded scientific outputs;
exact dirty source/PDF originals are retained under each repository's
`attic/2026-10-publication-reconciliation/`, with source SHA-256 guards before editing. No commit or
scientific analysis was run by this producer worker.

## Identity and promotion gates

The historical benchmark `requirements.txt` release pin remains an explicit baseline recipe.
Fresh reproduction requires `TCREN_REPRO_PYTHON` (absolute interpreter), `TCREN_REPRO_WHEEL`,
`TCREN_REPRO_WHEEL_SHA256`, and `TCREN_BASELINE_ARCHIVE`. The preflight verifies wheel archive
SHA-256, noneditable installed direct-URL path and version, and every installed `tcren/` Python,
native and resource payload byte. `uv` may omit its archive hash in `direct_url.json`; this is
handled by actual payload comparison, never by a version-only assumption. The shell drivers and
active CPL/receptor/contact/kinetics subprocesses use the same interpreter.

The preflight requires 12 raw structure sets and the supplied ATLAS/MD measurement tables. The
read-only inspection observed: Canonical2026 374, Native2026 374, VDJdb models 10,884, VDJdb native
369, balanced receptor panel 1,089, free pool 2,000, TCRvdb 618, IMMREP negative reference 832, CPL
2,103, ATLAS substitutions 378, HPV 81 and Shcherbinin23 3,240. These are raw-source counts, not
fresh scored/eligible population counts. Declared historical minimum counts guard named principal
sets; all sets must be present and nonempty. CPL, balanced-panel and IMMREP archives are extracted
into fresh unique raw directories; all harmonization raw sets are freshly linked/extracted too.
Raw original sources are never deleted or overwritten. Transient feature chunks are removed after
concatenation and on process exit.

Before publication reshaping, the driver (1) checks the exact 8,292 shipped holdout IDs against the
fresh full harmonized corpus through the library's `check_fresh_holdout.py`, without installing any
fit, and (2) runs unchanged historical floors plus `compare_established_results.py`. The latter
compares PCA identities/labels, per-cohort counts and every recorded ROC/AP value, and legacy
per-unit keys/counts/ROC/AP. Every delta is retained in a requested report. Any lost population,
missing value or decrease beyond floating roundoff stops promotion and requires author review.
The baseline archive contains 3,809 per-structure PCA rows, 31 cohort summary rows and 9,217 legacy
ranking rows. These are baseline records, not successful new reproduction counts. The guard does
not establish that an arbitrary unexported arm used identical finite support: absent per-arm
support records remain a limitation of the earlier summary format.

## Ordered graph repairs

- Balanced-panel scoring now precedes analysis that consumes it. `binder_recognize.py` (recognize
  and mechanics), `binder_tcr_reference.py`, `pose_binder_bm_score.py`, contact geometry and balance
  metadata run before merge/ranking, receptor feature joins and QTF designs. The existing merge,
  rank, split, differential and export routes remain. The historical anchor forest/tree and its
  500 permutation refits are excluded from the supported full route.
- `qtf_design.py` runs first in the Potts stage: both `potts_fit.py` and the constrained scorer need
  its fresh label joins. Both peptide and MHC site tables are rebuilt before joint fitting.
  The legacy site adapter preserves `chain.id.rec` and `chain.id.par`, closing a second loss of
  chain identity outside the repaired library. Joint-model cache reads are removed.
- Confidence GLM, rescale, gain and operating-point producers are explicitly ordered; the GLM
  belongs to the benchmark, distinct from the installed library's frozen conditional residual.
  Curves, paired within-epitope comparisons, cohort dispersion and template ambiguity are wired.
- All existing two-step descriptor/baseline/product/link producers run after their fresh template,
  marginal, legacy-ranking and harmonized-corpus inputs. Surface and self-consistency display
  producers are also wired. Metadata is exported from the shipped object without fitting.
- The missing `contact_topo.py` source was recovered from commit `ac84314`. Its exact former source
  is archived. The supported version reuses the installed topology helpers and the existing
  contact-partition/Hill-number helpers, and rebuilds raw structures. Four documented combinations
  are restored from that commit's `results/notes/contact_entropy.md`: standardized D2 plus signed
  standardized patch fraction, H0 persistence entropy, graph cycle rank or truncated Betti
  difference. No new combination, radius, threshold or control was added. The old topology CSV
  exists; an earlier worker statement that the CSV itself was absent was incorrect.
- `atlas_feats.py:rate_span` already produces `atlas_rates.csv` and `atlas_kd_split.csv` from supplied
  ATLAS measurements. The earlier audit's missing-producer concern is resolved by this writer.
- Historical quadrature of two percentile-CI half-widths is now computed by
  `receptor_panels.py:template_drop_summary`, with its approximation explicitly labelled. The
  manuscript reads the exact unchanged difference/half-width values, computes no interval, and
  does not label this as a bootstrap CI of the difference.

## Complete live reshaper inventory (34)

Paths in the producer column are benchmark-relative unless marked L (library). The source counts
above and fresh-run output counts travel through the authoritative producers; the exact wheel
identity is recorded in `results/recompute/input_identity.json`. Publication files are reshapes of
these outputs. This inventory establishes scheduled edges, not scientific freshness by itself.

| Live manuscript reshaper | Authoritative producer/input route | Stage |
|---|---|---|
| `make_arms_dat.py` | `pca/curves.py`, curve/summary CSVs from fresh per-structure scores | publication |
| `make_correction_dat.py` | `pca/confidence_correction.py`, benchmark GLM coefficients/cohort corrections | publication |
| `make_cpl_dat.py` | `eda/cpl_figure.py`, `cpl_predict.py`, CPL scores/models and response matrices | cpl/analysis |
| `make_cplfreq_dat.py` | `eda/cpl_contact_freq_model.py`, measured CPL MD contacts | analysis |
| `make_cplfusion_table.py` | `pca/cpl_fusion.py`, fresh CPL descriptor pack/legacy ranking | intervals |
| `make_datasets_table.py` | maintained `data/datasets.csv`, supplied provenance metadata | metadata input |
| `make_diag_dat.py` | `eda/af_bands.py`, active confidence-band/contact-family diagnostics | compose/Potts |
| `make_dispersion_dat.py` | `pca/cohort_dispersion.py`, `template_ambiguity.py`, fresh receptor/native features | publication |
| `make_erg_dat.py` | `producers/ergodic_inputs.py`, `eda/ergodic_bridge.py`, supplied MD measurements | ergodic/analysis |
| `make_footprint_dat.py` | `scripts/contact_subset.py`, recovered `contact_topo.py`; raw contact sets, exact labels | binder/publication |
| `make_hpv_dat.py` | `models/hpv_sanity_check.py`, two retained fit-free arms; original raw sources and explicit report order | analysis |
| `make_inversion_dat.py` | receptor/PCA per-cohort read-outs, raw two-epitope labels | native/harmonize |
| `make_legacy_dat.py` | `scripts/legacy_peptide_ranking.py`, `pca/cpl_fusion.py` | harmonize/intervals |
| `make_legacy_table.py` | same per-unit legacy ranking and accepted candidate pools | harmonize |
| `make_matrix_dat.py` | `eda/tasks_extra.py`, `matrix.py`, `verify_reframe.py`, named benchmark outputs | matrix |
| `make_md_dat.py` | `eda/cpl_md_boltzmann.py`, deposited MD occupancy/calibration records | analysis |
| `make_native_dat.py` | `scripts/receptor_panels.py`, fresh endpoints, strata, display approximation | native |
| `make_potential_dat.py` | fig2 table shell, `dphi_polyala.py`, `panel_stats.py`; raw crystals/accepted pools | analysis |
| `make_potts_dat.py` | `eda/potts_md_contacts.py`, `potts_md_baselines.py`; supplied MD contact maps | Potts |
| `make_qpca_dat.py` | harmonization/PCA invariant, screen, benchmark and figure writers | harmonize |
| `make_raqa_dat.py` | `rebuttal/raqapppsw_diagnostics.py`, `raqapppsw_contacts.py` | rebuttal |
| `make_rates_dat.py` | `producers/atlas_feats.py:rate_span`, kinetics uncertainty producer | atlas/rebuttal |
| `make_render_dat.py` | three scheduled selected-render recipes, fresh manifests and supplied generator confidence | figures |
| `make_rescale_dat.py` | confidence GLM → rescale → gain/operating-point CSVs | publication |
| `make_rescue_dat.py` | `pca/paired_within_cohort.py`, fresh per-epitope/paired-support outputs | publication |
| `make_resid_dat.py` | PCA conditional residual/per-structure diagnostic outputs | harmonize/publication |
| `make_score_dat.py` | PCA scores/intervals; safe L frozen-object metadata export | harmonize/intervals/figures |
| `make_selfcon_dat.py` | `eda/selfcon_compare.py`, frozen parameters plus fresh contact counts | publication |
| `make_tcren_dat.py` | potential/ablation/sidechain/epitope-weighted producers | analysis/Potts |
| `make_tcrvdb_dat.py` | `eda/tcrvdb_figure.py`; fresh validated-receptor scores/labels | tcrvdb |
| `make_tiers_dat.py` | legacy ranking plus existing two-step descriptor/product/baseline/link family | harmonize/publication |
| `make_tmpl_dat.py` | balanced-panel outputs, template confidence ANOVA, PCA template merge | binder/analysis/intervals |
| `make_tmpl_table.py` | same template-stratified authoritative tables | binder/intervals |
| `make_topo_dat.py` | `scripts/interface_topography.py`, fresh corpus/block/literature outputs | publication |

## Author-approved scientific wording and retirement

Exact recorded association values/macros are preserved. The static-score versus time-averaged MD
regression is described as a measured correspondence across complexes, not a direct ergodicity
proof. The kinetics association is retained without assigning an uncalibrated transition-state
energy. Mechanical stiffness/force/work use reduced units and a fixed displacement grid; contact
loss is the linearized nearest-atom-pair surrogate. `fp_chi` is a truncated Betti difference.
Binder log-odds plus logit ipTM is described as a heuristic: ipTM is not an established binder
likelihood ratio. Original-scale confidence prediction is sigmoid of a transformed conditional
mean, a plug-in prediction/conditional median under the Gaussian, not the conditional expectation
of sigmoid. Pose/confidence label-use claims are scoped to scoring time; shared nonlinear
preprocessing is fitted on both hold-out classes and binder-selected moments define these readouts.

Both conditional-inference tree publication routes are retired: main decision tree and validated
receptor tree panels. Their R sources, reshaper/plotting sources, tables, data and pre-existing dirty
PDFs are preserved in attics. Forced/clean receptor scatter findings remain. Pooled panel
family/GLM prose is archived; the live manuscript refers to the existing stratified and
per-epitope routes. Historical GLM gains were in sample against the structural term, not against
the better term. No historical numerical value was silently swapped for a new one.

The old diagnostic correction/per-structure table route is no longer invoked by `make_diag_dat.py`:
those producer-less withdrawn correction artifacts remain preserved. Active frozen GLM display
routes are separate. `selfcon_trace.csv` carries explicit **historical frozen-model fitting trace**
provenance. It is retained parameter-fit evidence, not a claim that evaluation freshly reproduced
that convergence trace.

## Source-level cache and skip review

The driver names 131 benchmark Python producer/checker invocations. An automated source inventory
checked their direct `exists`/`is_file` branches, followed by repair of demonstrated computed reuse:
affinity terms, ergodic interfaces, contact-probability TSVs, joint Potts fit, legacy ensemble
checkpoint/resume, kinetics annotation/scoring, feature chunks and raw extraction refreshes.
Legacy ensemble scores all eligible prepared units afresh, reports each unit, and lets its native
batched sampler own parallelism (`--workers=1`); an unsupported outer-worker count fails explicitly.
The optimization rung is existing batched native sampling, not a new parallel scoring engine.
Frozen potentials/models and accepted candidate/reference inputs remain inputs. The IEDB pool
helper preserves a single accepted decoy draw across compared arms; removing that would change the
baseline population and violate pairing.

Missing required feature/legacy arms, raw sets, frozen trace, confidence-correction and receptor
fixtures now fail explicitly. Some older matrix/checker and exploratory display code still has
optional-existing-input branches; their mere presence is not proof of cached computation, but a
static driver inventory alone cannot certify freshness of every such historical optional row.
Likewise metadata/render manifests require their own producer-or-input provenance. The fresh-run
acceptance and publication gates must evaluate actual outputs and stop if a required route or
matched standing result cannot be reproduced. This report does not turn a passed source check into
a green scientific claim.

## Verification performed

- Three focused synthetic/driver checks passed in 1.92 s: ordering/no-tree contract, exact historical
  display approximation, and preservation of distinct partner/receptor chain IDs in the legacy
  adapter. No corpus refit or scientific evaluation was run by this worker.
- Three active zsh drivers passed syntax checks. All changed existing benchmark Python sources and
  new/recovered producers passed AST parsing; recovered topology producer imports in the immutable
  wheel environment and its `--help` succeeds.
- Read-only manuscript gate passed: 18 main labels, 31 note sections, 45 supplementary floats,
  53 figure scripts, 28 existing PDFs, no structure/xref/vocabulary/figure violations.
  This used `make MAKE='make -o all' gate`; it did not build the modified scatter PDF, compile TeX
  or inspect freshly rendered pages. The parent fresh pipeline owns those remaining checks.
- Both cross-repository `git diff --check` checks passed before the final source freeze. All original
  dirty scientific outputs were excluded from this worker's source edits/retirement operations.

## Preservation-gate extension: dynamics and kinetics

The late acceptance checker now also protects the two named standing associations, without fitting
or recalculating them:

- Energy correspondence: compare the original `erg_fit.csv` R² = 0.45989207147101174 (printed 0.460)
  and its reported 184 fitting observations, using the same term keys. `erg_points.csv` exposes the
  exact structure/heldout identities and observed MD endpoint. The checker requires that complete
  multiset to match before comparing the reported R²; it comprises 184 fitting structures plus
  three heldout structures in the archived record. A lower R², lost endpoint/identity or changed
  count stops publication promotion.
- Dissociation rate: compare the original signed Pearson R = −0.5397988480206806 (printed −0.5398)
  for `log k_off` versus `W_tensile  [incumbent]`, with n = 23. Its actual support can be recovered
  from the archived/current mechanical rows and native feature table: class I, inner joined to
  successful native features, finite `W_tensile` and `logkoff`. The checker matches the exact
  identity/observed-log-rate multiset, preserving duplicate multiplicities. The archived support
  was verified as 23 rows over 23 distinct structures, agreeing with the original summary.

Every raw signed metric delta is exported unchanged. `strength_delta` additionally records
Δ|R| for the signed kinetic association (ordinary ΔR²/ΔROC/ΔAP elsewhere). Lower |R| or a sign
reversal is a review condition; a more negative R of the established direction is correctly
recognized as a stronger association. This does not reinterpret R as an unsigned reported result.
No new statistical test, control, threshold, refit or producer execution was introduced. The only
source changes during the running scientific pipeline were its late comparison checker and focused
synthetic checks. Five tests in the producer-regression file passed in 0.90 s, covering the earlier
display/identity repairs and the association-direction/support guards. Archive support inspection
read recorded inputs only.


### 15 October API/freshness follow-up before restart

The aborted fresh binder stage exposed the retired `recognize --scores` option. All executable literal callers now use the descriptor route, with the established `Phi`/`dPhi` → `F`/`dF` schema adapter where downstream readers require it. TCRvdb descriptors regenerate from the full raw set before its energy table; retired per-structure probability/strain fields are no longer required. Binder uses one batch per class, retains stderr, and its dynamic options are tuples. Pose recognition uses the current keyword contract. Contact subset geometry imports the canonical docking module and counts uppercase PEPTIDE chains. The optional Potts observed-map check uses the existing ContactMap and chain-balance implementation, without the removed interface_features import.

The Kit p_bind producer/export calls have no live manuscript reader; they are retired without substitution, with source snapshots retained in the follow-up archive and original outputs preserved by the established-results archive. The two active QC diagnostics remain. Harmonization reads exactly the eleven completed cohort tables, rejects any missing table, excludes hidden interrupted chunks, and featurization deletes earlier partial chunks before rebuilding them.

Installed-wheel static checks covered 132 directly scheduled scripts, 403 named tcren imports and 450 resolvable API calls: zero missing symbols or unsupported named keywords. CLI inspection resolved 21 literal calls with zero unsupported commands/options. Two dynamic wrappers remain: binder is exercised by a mocked subprocess argv test; kinetics delegates the reviewed existing mechanics/scoring calls. Static checks establish API syntax, not numerical correctness or full runtime schema. The fresh full run and preservation gates remain mandatory.

Focused verification: 9 synthetic preservation, dynamic argv/schema, completed-cohort/hidden-chunk and CPL identity checks passed in 0.90 s with warnings treated as errors and pytest cache disabled. zsh syntax and git diff whitespace checks passed. No scientific producer was rerun in this follow-up. The first test invocation's nine bodies passed but pytest cache creation was sandbox-denied; disabling that test cache yielded the successful exit above.


### Stopped-run binder catalogue and raw dispatch fixes

Three fresh `dFT_*_cdr12` decomposition columns are retained in the binder CSV but explicitly auxiliary: they do not enter the accepted split/design or add scan arms. Existing TCR-reference features remain catalogued. Fresh CPL/binder path assignments and exports now execute in the parent zsh before piped stage dispatch; later producers therefore retain the same raw extraction. The two-step CPL reader uses the shared paths.CPL_PDB override. Binder CLI outputs create their parent directory. Ten focused producer/dispatch tests passed in 0.92–0.93 s; the actual zsh pipeline regression verifies identical raw paths inside the piped function and in subsequent parent scope. Existing tests also verify a previously absent nested output directory and preservation/exclusion of auxiliary values. No scientific computation was run.


### Active figure asset closure after legacy-panel disposition

All 27 manuscript includegraphics roots have existing figure producers: Figure 1 self-contained TikZ; Figure 2 TikZ plus potential-bc gnuplot through its wrapper; 26 numerical gnuplot outputs include the potential-bc intermediate. No numerical figure driver is an unused duplicate. The parent agent owns the archived legacy assets and FIGURES status ledger.

The missing SB27 predicted heatmap is an existing late edge, not a missing scientific algorithm: fresh cpl_predict_cells.csv has all 260 SB27 cells; cpl_sidebyside_prep.py rank-normalizes its recorded pred values within each position, emitting sb27_tcren.dat together with experimental/cognate matrices. That producer already ran late in the figure stage; it now explicitly rejects missing/duplicate SB27 cells rather than leave an old PDF. Its focused synthetic check pins prediction ranking without computing a new score.

The twelve active structure PNGs now have a driver edge before manuscript reshaping: publication_panels.py runs the selected MEL5 cpl_importance_panels, f4_groove_side, and f5_tilted_cutaway recipes with their established structures, cameras, palettes and default parameters. All scenes write into a fresh temporary directory; only complete selected PNGs/manifests are copied to canonical names. Retired candidate names are not recreated live. Each selected asset gets recipe, raw-byte and PNG hashes in publication_panels.provenance.json. The early preflight checks rendering tools and exact selected raw structures (archive members in a clean checkout); derived table checks occur when the render stage starts.

Provenance limits: four canonical tcrvdb PNGs are byte-identical aliases of the settled f4 groove-side images. Both RAQ canonical PNGs exactly match archived f5 tilted-cutaway images after the established 13-pixel edge shave (ImageMagick absolute error zero). The archived FLR trial pixels differ from canonical FLR; current settled camera/default source parameters are retained, but that historical pixel derivation is not fully recovered. This must remain marked for visual comparison, never asserted as bit-identical. Unused withdrawn Q_geom/combo/p_bind render-manifest dependencies are removed, preserving historical manifests in the archive; live readers use the actual energy/contact/height quantities.

Source checks: 12 focused tests pass, including canonical-only render coordination and recorded SB27 reshape. The installed-wheel render API scan resolves 16 named imports and 19 calls across three exact producers with no missing symbols or incompatible keywords. Both current-tree and deliberately absent-extraction archive preflights pass. Actual bounded render verification is recorded separately below when complete; none of these checks reruns benchmark scoring.

Actual selected-render verification completed in 89.6 s, writing only /tmp/tcren-selected-render-verification: all twelve PNGs and three scene manifests were generated/cropped, PyMOL traceback/error scanning passed, and temporary intermediate scenes were removed. Existing MEL5 numerical manifest fields are identical; Fig4 fields are identical except the corrected colour unit kT→dimensionless. Fig5 keeps contact counts 20/0/7/18 and energy/label/confidence fields; RAQ binder CDR height changes 43.42→43.23 Å and signed height -40.36→-40.00 Å; FLR binder CDR height 11.68→13.06 Å and signed height +11.74→+13.14 Å. All matched changes are recorded in /tmp/tcren-render-verification-deltas.json and reported to the author before canonical promotion. Full cohort scoring was not rerun. SB27 recorded-table reshaping and gnuplot rendering in /tmp succeeded with 260 predicted cells and a 28,702-byte PDF. The manuscript render-macro reshaper now rejects missing manifests or selected panels.


### Attempt 7 pose-join root cause and source repair

The stopped run produced 1,089 recognition rows and 1,089 pose rows but retained only the 523 positives after metadata joining: negatives in deposited metadata carry full hash_epitope_allele stems, while the two fresh producer passes encode bare hash plus class. Normalizing the metadata key yields exactly 1,089 unique (hash,class) keys; 21 bare hashes occur in both classes and must remain distinct. The repaired pose producer validates exact key-set equality and uniqueness against the raw structure manifest for metadata, recognition, pose and TCR-reference inputs; it rejects collisions rather than unique/first or permissive row thresholds. Raw class counts are fixed at the accepted 566 negative and 523 positive structures, and unknown source directories fail.

A second defect was found at the same join: deposited metadata also includes historical Q/energy/probability/geometry fields. Taking the complete metadata table kept those values in unsuffixed columns while fresh descriptors acquired _right suffixes. The producer now projects only raw labels, generator confidences, sequence and bookkeeping metadata; all descriptors come from fresh recognition/pose/reference tables. Current Phi/dPhi names use the existing F/dF schema adapter. Unexpected overlapping non-key columns fail, and missing TCR-reference input fails. Fresh diagnostic CSVs remain output records and are never inputs to the reproduction.

Verification: 13 focused synthetic tests passed in 1.33 s, including a repeated hash across classes, full negative filename suffixes, stale metadata sentinel values, missing keys, duplicate normalized keys and unexpected descriptor/metadata collisions. A read-only join of the stopped-run diagnostic tables confirmed all 1,089 keys (566/523), fresh recognition burial and no suffix collision, without writing outputs or running scientific scoring. The next full run must regenerate both producer passes; these diagnostic rows are not a resume route.

Narrow sibling audit found no additional required patch: binder_benchmark_final already normalizes the deposited key and projects labels/confidence only; contact_subset uses stem-normalized (hash,class) labels from the repaired pose table; receptor_panels selects confidence/bookkeeping/reference columns and aliases fresh feature energies; harmonize.datasets drops stale membership descriptor columns before joining fresh features; build_feature_table is unscheduled and projects confidence metadata. No scientific metrics or prose conclusions changed in this repair.


### Retained HPV panels and native-reference edge

The live HPV producer still called the retired cohort-fitted `p_native` API. Its panel c had been relabelled as a composition while plotting posterior values; no active result paragraph used that panel's macros. Under the author's legacy-panel disposition, its original producer, four CSVs and three-panel publication source/PDF/data are preserved in the two repository attics. The live route retains Q and the binder-oriented referenced energy only, with no replacement score. The manuscript labels the nested all-positive/core contrast descriptive and identifies the actual populations; existing tests, intervals and values are unchanged.

Owning mut2023 raw metadata exactly reproduces the old positive identifier/origin/kind mapping. The original coordinate choices remain: four Native2026 crystals, 56 Shcherbinin23 mutants and 21 original archive decoys. The four crystal variants in the mut2023 archive differ and are not substituted. `data/hpv_report_order.csv` freezes the original positive identifier order from the existing provenance inventory; labels and paths come from owning raw metadata. Its exact positive-key set is validated, negatives retain raw metadata order, and all fresh feature rows are explicitly restored to report order before seeded statistics. `SOURCES.md` provides a regeneration command using the archived identifier-only inventory. No old feature or score row is an input. Preflight fingerprints this small metadata input.

The preservation gate compares exact retained point support/values and every existing reported statistic/count in both directions; the retired posterior is excluded. The existing `native_tcr_reference.py` is now called before `qtf_design.py`; its scoring recipe is unchanged, paths portable and import inert. Nineteen combined synthetic producer checks passed in 1.16 s; changed producer sources pass Ruff. Temporary reshaping of recorded HPV values retained all 32 live macros exactly and produced a visually inspected two-panel PDF. The twelve temporary selected structure scenes were also visually inspected. Full scientific freshness and publication certification remain pending the complete run.
