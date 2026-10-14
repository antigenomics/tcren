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
| `make_hpv_dat.py` | `models/hpv_sanity_check.py`, raw HPV models/measurements | analysis |
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
| `make_render_dat.py` | render manifests and existing supplied render metadata | metadata input |
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
