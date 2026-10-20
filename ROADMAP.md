## 18 October: multiplicity-calibrated weighted TCRen

Author-requested extension: calibrate native/external compatibility ranks conditional
on atom-distance and distinct orientation opportunity counts, freeze before evaluation,
and compare weighted/original TCRen and native ipTM on unchanged benchmark rows.
Use bounded Stratos scoring and R diagnostic figures. Completed on all 3808 rows,
with exact baseline/graph/membership checks and six R figure pages. Keep all 22 VDJdb cohorts by
coverage, both TCRvdb epitopes and seven CPL clones separate; expose every gain/loss,
missing/transfer support and signed mechanism. No default promotion or paper-prose
rewrite. Design and verification: [CONTACT_WEIGHTED.md](audit/2026-10-18/CONTACT_WEIGHTED.md).

## 18 October: completed contact-penalty mechanism investigation

The author requests a meticulous case-by-case explanation of the completed penalty
comparison on unchanged rows: signed energy attribution in GLC/YLQ, chemical and
continuous geometry contact classifications, binder/nonbinder differences in every
VDJdb case and template stratum, matrix preservation under filtering/weighting,
and CPL extremes where ipTM succeeds or fails relative to TCRen. The author approved
both chemistry-supported/proximity-only counts and continuous equivalents Σ(1−u)
and Σu; no binary geometry cutoff is introduced. New diagnostic geometry exports
trace the worst atom/marginal and compare closest-atom with any-atom typing. They
reuse frozen parameters and do not recalculate or replace accepted scores. Literature
and exact contact decompositions must support the mechanism; AF output patterns do
not by themselves establish why its network chose a pose. No label-fitted score,
threshold search, default change or publication-prose rewrite is authorized.
Completed evidence: [CONTACT_MECHANISM.md](audit/2026-10-18/CONTACT_MECHANISM.md).
All 3808 rows and 85669 assessed edges are verified; 31 separate case records
identify contact-burden and signed-entry mechanisms. No new score is promoted.
Calibrated residue-edge reliability, full angular/steric chemistry and the broader
reference-prior goals remain separate follow-ups.

## 18 October: directed context correction and CPL contact-count audit

The author requires binder → target direction throughout residue statistics,
potentials, reference marginals and orientation priors. Remove residue folding
and triangular reductions; only individual atom-pair chemistry can be exchangeable.
The earlier unordered residue interpretation is withdrawn. The shipped original
scores remain fixed comparisons; no historical coefficients or moments are refitted.

CPL audit separates the exact count term from mean geometry and uses existing
source-matched receptor-free pMHC SASA/surface descriptors on the same rows.
Scalar SASA is not a residue/atom opportunity distribution. Test a fixed
reference-partition correction g=u+log(mean_reference(exp(-u))) within exact
native class-I geometry-opportunity strata, frozen before labels. It cancels a
constant background penalty and normalizes the empirical relative factor, without
a fitted score coefficient. It is not a physical binding-free-energy claim.
Retain original, unnormalized and corrected scores, all rows, all gains/losses,
zero contacts and explicit unsupported strata. No new default is promoted.
Private evidence: private-evaluation-record

## 18 October: active contact-plausibility and positive-penalty goal

Follow the author's manuscript `plan181026.txt` to summarize directed binder-to-target
native residue contacts and exchangeable individual atom-pair chemistry, verify quantitative literature priors, and derive frozen
chemical/geometric contact-plausibility rules. The author explicitly added a
candidate **positive-energy penalty for unnatural contacts**, beyond deletion
or down-weighting. GLC's recorded shared-motif inversion remains part of the
mechanism to explain; rarity or assay-negative membership is not physical truth.
[CONTACT_PLAUSIBILITY_PLAN.md](CONTACT_PLAUSIBILITY_PLAN.md) owns the design,
source/API reconnaissance, phase verification and current review status.
Compare frozen candidates with original energy and matched ipTM,
preserving populations, class/template strata, gains and losses. No label-fitted
classifier weights, cohort normalization, scoring-time fitting or altered
historical default. The first implementation stage is complete:
[CONTACT_OBSERVATIONS.md](audit/2026-10-18/CONTACT_OBSERVATIONS.md) records the
exact CLI export and verified fresh native summary (357 complexes, 339921 atom
pairs), including two explicit AMN-peptide exclusions versus historical parsing.
The first frozen transformed-marginal resource and four-arm matched comparison are
complete: [CONTACT_PENALTIES.md](audit/2026-10-18/CONTACT_PENALTIES.md) records
3808 freshly scored models, zero matched-row exclusions, 1.358 GiB sampled peak
RSS, exact baseline reproduction and every cohort-specific gain/loss. No default
is promoted. Broader chemistry/PPI priors and class-II model evaluation remain open.
The earlier bounded broad unit gate timed out; completed focused tests, docs and
raw-coordinate checks are recorded in the stage reports.

The author's follow-up authorizes a finite comparison of fixed additive,
matrix-maximum-scaled additive, literal sign flip and favourable-only sign flip.
The literal flip is a signed diagnostic: it can reward already positive entries.
Keep the same frozen unnaturalness assignment across arms; retain every matched
gain/loss. A common +matrix-maximum replacement is recorded as a proposal.

The [quantitative source dossier](audit/2026-10-18/CONTACT_PRIOR_SOURCES.md)
now resolves the large antibody survey to its 1833-complex publication, DTU
deposit and pinned contact-table producer, and records inspectable PPI/Arpeggio
criteria. The deposit is retrieved and the approved first marginal posterior is frozen;
broader quantitative PPI/atlas priors remain open.

## 17 October: YLQ/GLC contact-set difference investigation

The author-requested complete eligible-cohort/native comparison is recorded in
the benchmark's `CONTACT_DIFFERENCE.md` and `results/contact-difference-20261017-a/`.
Fresh exact CLI exports cover 617 models and five matching native complexes;
independent checks cover physical graphs, potential entries, set membership and
native preservation. Counts, peptide-site occupancy, receptor-region/amino-acid
motifs and local chemistry/context remain distinct. Native-reference dependencies
and all prior ranking losses are explicit; no label-derived filter or fitted
scorer is promoted. Existing chemical typing was audited through its direct API;
[issue 59](https://github.com/antigenomics/tcren/issues/59) records the exact-export
CLI gap. The next mechanism scope is region-specific native-compatible contact
preservation/replacement. Stratos MD and release approvals remain untouched.

## 17 October: contact-local topology/topography goal

The author clarified that contact filtering concerns local graph and surface
properties, explicitly not chemical `--drop-untyped` typing. The new
`export --quantity contact-context` reuses the contact map, pMHC surface and
potential; it retains every original contact and adds endpoint degrees,
four-cycle support, local height/hydropathy and exact pair contributions.
Benchmark `CONTACT_CONTEXT_GOAL.md` owns the current manual discordance audit
and matched mask evaluations. Keep original scores, zero-contact models, all
losses and class/template strata. No fitted classifier or threshold-selected
ipTM switching rule is promoted. The matched 2719-model audit checks all 67416
contacts and 504 accuracy values, with 44 views of 20 models inspected.
`contact-summary` owns canonical table-only reductions; all original baselines
and the differences from the first sums remain recorded. Stratos MD remains untouched.

## 17 October: full model-only assay extension

The complete 2719-model CPL and separate YLQ/GLC backbone export and fixed-direction
CLI evaluation are verified. Existing APIs suffice: 34 model-only diagnostics,
1476 accuracy values, 111479 memberships and 1008 label-conditioned spatial
associations. CPL chemistry gains and the MEL8 pose-to-assay transfer failure
are retained together; 126 prior baseline comparisons are unchanged. See benchmark
`ASSAY_DEFORMATION_VALIDATION.md`. All class/template populations remain separate.
No fitted score, native-reference classifier or Stratos MD operation. Additional
independent nonbinder structural truth remains the next evidence scope.

## 17 October: local conformation and spatial-chemistry validation

The next matched post-analysis is complete: 384 associations on 20 cognate saved
models, four targets kept separate, all class I/covered. MEL8 maximum CDR3alpha
virtual angle tracks recognition DockQ without a native reference; covariance
associations retain target-specific signs and MEL8's constant-zero ridge readout.
Existing CLI surface and backbone exports suffice; no new scorer, library module
or Stratos operation. The benchmark `POSE_CONTEXT_VALIDATION.md` owns exact
quantities, matched memberships, physical windows and all regressions. Class-II,
scarce-template and independent nonbinder structural truth remain separate scopes.

## 17 October: scientific and dependency reconciliation

The existing CLI exports already cover separated contact graphs, intrinsic peptide/CDR
geometry and receptor-relative placement. No replacement scorer or new annotation
framework was needed. Cognate bound references can be used with recorded framework
variants; framework identity is not an automatic pose-validation exclusion. Peptide,
MHC and relevant loop correspondence stay explicit. This supersedes treating exact
full-framework identity as a prerequisite for the cognate 868 comparison.

Historical DockQ artifact authentication remains an archival task, not a blocker to
fresh cognate structural validation. Class-II and scarce-template extensions remain
separate evidence scopes rather than prerequisites for the existing class-I result.
No fitted final score or Stratos MD operation.

Dependency requirements are implemented: Arda's batch region markup and content-bound
index provenance, and seqtree's bounded batch/alignment/grouped APIs. Release follow-up
is explicit in [seqtree #18](https://github.com/antigenomics/seqtree/issues/18) and
[Arda #149](https://github.com/antigenomics/arda/issues/149). The author owns publication.

## 17 October: interface shape definitions and exports

Implemented and checked separate CDR:peptide / CDR:MHC / peptide:MHC contact graphs
and peptide/CDR intrinsic C-alpha geometry through the batched CLI; existing scores
unchanged. See [SHAPE_EXPORTS](audit/2026-10-15/SHAPE_EXPORTS.md). Publication follow-up:
restore fixed-template residue reconstruction, define normalized effective rank before
its first result, retain all historical differences, compare sequence-identical loop
windows without treating grafted 868 as an identical full construct. Independent
forced-deformation validation remains open. Stratos MD untouched.

# Publication and software roadmap

Reconciled 16 October 2026 against master base `68a1521`, the deterministic topology fix, live GitHub issues and
recorded benchmark/manuscript sources and the linked issue-34 audit.
[AUDIT.md](AUDIT.md) owns scientific findings; this file owns remaining work.
[The example audit](audit/2026-10-16/EXAMPLE_AUDIT.md) records current execution,
manual inspection and limitations. Earlier roadmap text remains in Git history;
[the release record](audit/2026-10-15/RELEASE_FIXES.md) preserves the attempts,
metric regressions and annotation-environment corrections.

## Scientific boundary

Final scoring uses fixed empirical potential tables with explicit derivation and
identity, structural alanine references and individual structural descriptors.
Gaussian/GLM, fitted Potts/transition and reference-estimated Q/T/S final readouts
are retired; their APIs, parameters and results remain historical evidence.
A contact log-odds sum is not calibrated affinity. No replacement fitted scorer,
new control, threshold scan or cohort normalization is on this roadmap.

CPL peptide ranking and receptor ranking remain separate. Compare arms on identical
rows. Keep TCRvdb epitopes, VDJdb epitope/template strata and MHC classes separate.
Record gains and regressions; do not change labels, populations or precision to
conceal a loss. Class-II crystal coverage does not establish modelled class-II
ranking performance.

## Completed, with recorded evidence

- Correctness/domain fixes, exact CLI exports and searched-index provenance were
  integrated. Subsequent PRs #40–42 add structural alanine references, shared
  reference contact construction and co-best MHC alleles. Master tests and docs
  workflows passed at `c04cb94`. Subsequent issue-43 and issue-21 example receipts
  pin their distinct byte-matched installed 4.0.1 wheels; none establishes a public release.
- The five-block fit-free reproduction completed. The obsolete “attempt 8 stopped”
  status is superseded by `RELEASE_FIXES.md`. Later source/annotation corrections
  remain separately recorded rather than being called the same reproduction.
- Coordinate/assay pairing audit retains all 618 raw TCRvdb models and explicitly
  excludes one incompatible pMHC pairing from the eligible 617-model comparison.
  The existing matching GLC model is retained. Raw sources were not overwritten.
- The original published fixed-template procedure was reproduced with unchanged
  scientific Java/R calculations and modern CLI on 2,102 later CPL candidates.
  These are later data, not an original published CPL benchmark. Original versus
  modern published-potential scores agree within 4e-14. Six clones improve ROC/AP
  with the current potential and MEL5 ties; its graded-response regression remains.
- Historical energy-migration rows were recovered. All 272 energy/cohort comparisons
  remain; 42 violate at least one old universal bound. The claim and exceptions
  coexist in the historical record. Class II and template strata remain separate.
- The active main paper/SI restores valid fit-free findings and archives fitted-only
  claims. The recorded clean Git-worktree build has 19 main + 69 SI pages; all 39
  active PDF assets pass font/bounds gates. Manual inspection and exact hashes are
  in the manuscript's `MANUAL_FIGURE_QA.tsv`; the complete historical claim ledger
  remains in `QUANTITATIVE_CLAIM_AUDIT.tsv`.

Benchmark evidence: `STRUCTURE_IDENTITY_AUDIT.md`,
`results/published-cpl-final-20261016/`,
`results/historical-energy-migration-stratified-20261016/`.
Manuscript evidence: `PUBLISHED_CPL_REPRODUCTION.md`, `FIT_FREE_DISPOSITION.md`,
`FIGURES.md`, `HISTORICAL_AUDIT.md`. These live in their owning repositories, not
in the software distribution. Completed analyses are not rerun to update this list.

## Remaining work, in order

| Priority / owner | Work and acceptance evidence |
|---|---|
| Completed / library; issue [43](https://github.com/antigenomics/tcren/issues/43) | Fixed empty/all-failed surface exits, query batching and shared CPU-budget propagation. Actual MMseqs traces verify one complete multicore query batch per reference and no outer search pool; 23 fresh probes and per-process resource/equality checks are recorded in [CLI_BUDGET_FIX.md](audit/2026-10-16/CLI_BUDGET_FIX.md). Release publication remains separate. |
| P1 / library release | Prepare a release containing post-4.0.0 changes and ARDA dependency requirements. Build/install immutable artifacts, run supported-Python CI/docs and verify packaged resources and ordinary CLI behavior. The latest observed GitHub release is 4.0.0; local 4.0.1 metadata is not publication. Preserve the human PyPI approval gate. |
| Completed / benchmark + manuscript; issue [34](https://github.com/antigenomics/tcren/issues/34) | All 15 active molecular panels of 13 structures consume exact CLI tables. Physical residue/region identities, potentials, actual colours, cameras and geometry hashes are verified; all panels and composed affected pages were manually inspected. Class II stays separate. [MOLECULAR_CLI_LINEAGE.md](audit/2026-10-16/MOLECULAR_CLI_LINEAGE.md) records acceptance, archives and resource evidence. |
| Completed / manuscript; issues [33](https://github.com/antigenomics/tcren/issues/33), [35](https://github.com/antigenomics/tcren/issues/35) | The final scientific include/producer graph excludes fitted final readouts; all ten method families and estimated parameters have explicit dispositions. All 540 archive files, 242 local model source/output identities and 17 accepted result tables are preserved. Retained descriptive regressions remain. [FINAL_PRODUCER_DISPOSITION.md](audit/2026-10-16/FINAL_PRODUCER_DISPOSITION.md) records the exact reference-identity correction, historical execution guard and future structural CLI semantics; no scientific rerun or replacement fitted score. |
| Completed / library + benchmark | Fixed canonical physical-key reduction order for all four literature-topology columns. Two complete fresh raw-coordinate topology passes reproduce exactly across 2,102 CPL, 618 raw TCRvdb and 1,089 VDJdb models; the established 617-row TCRvdb eligibility remains. All original quantities, rank changes and 54 metric decreases versus the publication baseline are preserved alongside gains. Class/template strata remain separate; no rounding, refit or replacement score. [DETERMINISTIC_TOPOLOGY.md](audit/2026-10-16/DETERMINISTIC_TOPOLOGY.md) records wheel, CPU, exact-double TSV and verification evidence. |
| Completed / library; issue [21](https://github.com/antigenomics/tcren/issues/21) | Notebook/stdin calls now use an importable subprocess entrypoint for the existing bounded spawn pool, after whole-set annotation. Real Jupyter, stdin, serial/all/explicit budgets, worker failure/cancellation, physical-key comparisons and manual class-I/II table inspection are recorded in [INTERACTIVE_BATCH_FIX.md](audit/2026-10-16/INTERACTIVE_BATCH_FIX.md). Final affected suite: 58 passed; the full local reference-backed suite hit its 360-second cap and remains incomplete. Preinitialized multithreaded Polars remains unsupported. |
| Closed by historical disposition; issue [32](https://github.com/antigenomics/tcren/issues/32) | Exact scan overlap is verified for centred full-field and fixed-marginal arms; all other decomposition/normalized-weight routes remain explicitly unsupported, not relabelled as scan. Original quantities and regressions are preserved. Custom geometry, physical residue joins, units/signs and the frozen-marginal surrogate are audited in [HISTORICAL_POTTS_DISPOSITION.md](audit/2026-10-16/HISTORICAL_POTTS_DISPOSITION.md). The feature gap is retired with the fitted final model, not claimed implemented. |
| CLI fixed; campaign pending; issue [26](https://github.com/antigenomics/tcren/issues/26) | Installed `md inspect`, `md quality` and `md compare-forces` expose existing exact readouts, atom/component provenance, class-separated exports and CPU limits. 64 focused tests, installed-wheel class-I/II examples and manual table inspection pass; [MD_DIAGNOSTICS.md](audit/2026-10-16/MD_DIAGNOSTICS.md) records acceptance. The last recorded original eight-thread Stratos campaign had one 100 ns arm complete and the other running; its live state was not rechecked; final paired acceptance and missing historical campaign inputs stay unresolved. Do not restart the simulation or relabel reconstruction as exact reproduction. |

Six existing issues (#21, #26, #32–35) were open and no PR was open at the
initial checkpoint. This audit added issue [43](https://github.com/antigenomics/tcren/issues/43)
for the demonstrated CLI gaps. The subsequent [fix](audit/2026-10-16/CLI_BUDGET_FIX.md)
resolves #43; the [interactive fix](audit/2026-10-16/INTERACTIVE_BATCH_FIX.md) resolves #21.
The molecular-panel closure resolves #34. The final producer/disposition audit resolves #33 and #35. The historical Potts disposition retires #32 without claiming complete CLI parity. Only issue #26 remains open: its CLI gap is fixed, while the active Stratos campaign's final acceptance remains pending. Open/closed GitHub state is distinct from the implementation and
publication evidence above. The current example audit does not certify every CLI,
reproduce the long MD campaign or authorize paid review, fitting or publication.

## Submission exit

Freeze one reviewed source/wheel/environment identity; ensure each active figure
and number has its raw-input → CLI → summary → render lineage. Run structural,
xref, font and bounds gates; inspect the actual final PDFs. Verify data/code
availability and deposits within their access boundaries. Finish release approvals
and verify the public install before claiming a release. There is no date promise.

Joint scalar pMHC context is complete: repeated CLI feature products, exact matched
score/ipTM preservation, separate template/class strata, component/product support
checks and published negative changes. See
[JOINT_FEATURE_PRODUCTS.md](audit/2026-10-16/JOINT_FEATURE_PRODUCTS.md).
Spatial hydropathy co-localization, independent interaction and a validated ranking
boost remain distinct scientific gates. Issue 26's Stratos campaign is untouched.

Spatial pMHC co-localization is now CLI-owned and validated: exact cell TSVs,
receptor-independent frames, population height/ridge-excess hydropathy covariances,
component identities and missingness. The native/model atlas and matched score
associations preserve class/template strata and all weaker results; see
[SPATIAL_PMHC.md](audit/2026-10-16/SPATIAL_PMHC.md).
Independent causal interaction, a validated ranking boost and forced-pose
sensitivity/specificity remain evidence gates. Assay negatives are not pose truth.
Issue 26 and Stratos MD are untouched.

## Current structural-validation goal — 17 October 2026

Existing batched library APIs cover the required measurements; no replacement
score or fitted model is justified. The complete updated evidence and review
record is [CPL_LINEAGE_REVIEW.md](audit/2026-10-17/CPL_LINEAGE_REVIEW.md).

Resolved: the original1108-model interface definition;15 saved predictions of
three sequence-verified targets; MEL8 reference-copy configuration; all-copy
868/SB27 comparison with504 matched diagnostic records; complete546-model
reference-deviation reconstruction. Every finite original868 total (215) matches
5NME copy2, and106 original missing totals reproduce empty peptide alignments.
ForSB27,133/225 original totals match at least one physical copy; the later
automatic-mapping replay resolves the other92 through biologically invalid assignments
(see the follow-up below). The new records preserve all old scores, class/template groups,
missing values and adverse results. Copies/ranks are not independent targets.

Remaining objective, in order:

1. Recover original868/SB27 interface files, native-coordinate/executable identities
   and authentication of their producing lineage. The92 SB27 numerical discrepancies
are explained by incorrect automatic chain roles, including cross-copy assignments. Newly reconstructed values
   do not fill missing historical records; original missingness stays explicit.
2. Validate on independent identical-complex targets within each intended
   class/template stratum before claiming general forced-pose detector accuracy.
   Current evidence establishes within-target confidence selection errors and
   receptor-relative associations; no fitted detector or causal forcing claim.
3. Retain mechanistic bulge/hydropathy interaction and a validated fit-free ranking
   boost as separate scientific gates. Existing context associations and every
   score regression remain; no optimized sign or threshold.
4. Preserve separate modelled classII and full five-block locked-reproduction
   gates. Verify code/data availability and public installation through the
   existing release process, leaving the human PyPI approval untouched.
5. Keep issue26's MD campaign acceptance separate. Stratos was not accessed or
   changed; no live campaign-state claim is made by this follow-up.

The simulated structural-immunology review is retained in the benchmark, without
attributing opinions to named scientists. Its concrete Methods/provenance changes
are integrated into main/SI, retaining the task distinction and valid findings.
GitHub checked17October2026: only issue26 remains open; latest published release
is v4.0.0. Source/CLI implementation, scientific validation and public release
remain distinct gates.

## Current author presentation goal — 17 October 2026

The author requests concise, intuitive main/SI following the actual August draft,
with the roadmap and goals updated accordingly. The manuscript's
`READING_COPY_ROADMAP.md` now owns that presentation goal. The August comparator
is commit92942d8, compiled and visually inspected; its task-first organization is
reused without restoring fitted claims. Seven main figures, twelve SI figures
and three task tables replace the accumulated reading atlas. Full values,
regressions, old figures and lineage are preserved in the electronic index and
checksummed expanded draft. Existing library APIs suffice; no library calculation
or scientific rerun is needed for this edit. The independent truth/source-lineage,
class-II, locked reproduction and release gates above remain open. Stratos MD
remains untouched.


## Structural class-II registers — 17 October 2026

`export --register-reference-id/--register-p1` now emits all forward/reverse
core correspondences after the existing MHC-groove alignment. Exact physical
anchor identities, residuals, competing fits and the rigid transform are retained.
One annotation batch serves the full input set and reference. Source-verified
`--peptide-selection` supports MHC- and TCR-fused antigens without removing annotated
receptor regions. Optional modified-residue parsing retains chemical geometry as X;
no standard potential is assigned to citrulline. Historical sequence prediction
and ordinary parsing defaults remain available for reproduction.

The benchmark owns source-copy reconciliation, all 94 native class-II register
checks (including 42 sequence ties and the 6BGA auxiliary-peptide error), the
standard-AA calibration and matched regression tables. Manuscript/SI use the
recomputed calibration only; previous values remain electronic evidence.
Stratos MD is not part of this change.

## Automatic DockQ mapping follow-up —17October2026

Full546-model replay and supplied-reference checks: CPL_AUTOMATIC_MAPPING_AUDIT.md. SB27 automatic assignment reproduces 225/225 finite recorded totals; 92/225 assignments have incorrect biological roles. Every old value, unmatched case and descriptor regression remains. 5NMD is unbound868;4MNQ is engineeredILA1, not SB27. Strict pose truth still requires matching constructs and physical roles. Original artifact authentication and independent class/template validation remain open. No fitted score or Stratos operation.

Details: audit/2026-10-17/CPL_AUTOMATIC_MAPPING.md. This updates the numerical-discrepancy gate above; original artifact authentication remains open.


### 17 October: figure interpretation and matching 868 follow-up

The current primary-record search returns unbound 5NMD and bound 5NME/5NMF/5NMG; none matches the complete saved CPL 868 construct. All five saved predictions share the same receptor sequence; the observed receptor substitutions persist across physical native copies. The benchmark records the exact query and 28 physical-role/copy comparisons. Pose accuracy and forced peptide/CDR deformation are now explicitly separate. Independent construct-matched deformation validation remains open; existing diagnostic APIs suffice pending appropriate structural truth. No replacement score or Stratos MD operation.

### Directed antibody context correction — 18 October 2026

Completed binder→target cleanup across count derivation, potential moments,
angular priors and active benchmark producers. No residue reversal pooling or
transpose averaging remains in the supported routes. Only individual atom-pair
chemistry is exchangeable. Intrinsically symmetric published MJ tables are
unchanged; they do not impose symmetry on TCRen.

Both fresh directional-v2 runs retain 3808 models,7616 model–potential rows and
88073 physical edges. Original energies, labels, membership and ipTM reproduce
exactly. The frozen native-opportunity partition removes the constant geometry
background without assay-fitted coefficients. Peptide-prior/TCRen2 SB27 AUROC:
0.861510 original,0.176339 additive,0.820167 normalized. Gains and losses across
all 31 cases are retained; no default promotion. Full pocket/exposure/water
normalization and independent assay validation remain open. Private visible
evidence: private-evaluation-record and
private-evaluation-record.

### Directed physical-context goal reinstated — 20 October 2026

Active goal: improve contact definition and TCR ranking, including complementary
information to geometry/topology/ipTM and hard CPLs. Continue from the preserved
directed-v2 results; the existing multiplicity partition is not exposure/pocket
normalization. No assay-fitted weights or default promotion. New arms must use
the same 3808 raw models, retain original scores and all 31-case gains/losses.

Ongoing dependency order: atom/residue exposure with intact target context →
local pocket opportunity → area-conserving competing partners → frozen directed
reference/partition → fresh same-row evaluation and independent transfer checks →
R figures, tables and explicit examples. Physical definitions and invariants
precede benchmark inspection; scalar SASA or observed degree alone is insufficient.

Separate private independent-audit handoff and living state:
`private-evaluation-record`.
It records provenance, exact completed outcomes, implementation limits, pending
work and reproducible audit questions. Full physical-context evaluation remains open.

20 October implementation checkpoint: 3088cdf adds CLI-owned atom-resolved
`contact-exposure`, retaining intact pMHC occlusion and separate peptide-only MHC
shielding. Seven focused checks pass; existing scores are untouched. A 4848-atom
implementation diagnostic took 0.74 s locally. The subsequent surface-reference
checkpoint below supplies whole-cohort extraction, competing patch ownership and
a first exposed-composition partition; complete physical-context matching remains open.

20 October surface-reference checkpoint: d6d45a4 adds native conserved patch
allocation and `contact-penalties --surface-context`. Exposed noncontact CDR and
intact-pMHC peptide area define separate directed reference marginals; each
structure's partition is explicit, with no assay fit. Original arms remain.
Twenty-five installed-wheel checks and four visible-integration checks pass;
Ruff and warning-fatal Sphinx pass. Independent six-model verification preserves
all 12 original/prior rows and reconstructs areas and sums. The fresh 3808-model
Stratos job `ab-surface-full-20261020` completed successfully in 1049.78 seconds.
Full verification preserves all 7616 original/prior score rows and every benchmark
model; independent area/partition reconstruction passes. Final per-case changes
are mixed: the surface-context TCRen2 arm gains on three CPL cases and loses on
four; receptor-ranking gains and losses remain separated by epitope and template
coverage. The private `private-evaluation-record` separates the change
of contact measure from the partition correction for all 248 case/potential/metric
combinations. This candidate is not promoted. Frozen definition, provenance and audit commands are in the private
handoff and `private-evaluation-record`. Local pocket/hydration matching,
antibody atom-prior transfer and independent predictive validation remain open.

Final-product audits now include 496 pair-ranking comparisons against original
scores, shape complementarity, graph connectivity and ipTM. All model identities
remain; unavailable descriptor pairs are explicit. Discordant correct pairs exist,
but do not establish an improved combined predictor. The atom audit finds 11.572%
of peptide-assigned binder area attributed to target atoms with zero sampled
isolated-pMHC exposure. The current deepest-sphere ownership conserves area but
can select an atom buried inside the target union; pooling exposure by amino-acid
type masks most such assignments. Before another scoring run, establish an
accessible-boundary ownership definition and its geometric invariants. Preserve
the completed run, its losses, and its full atom tables. Private evidence:
`private-evaluation-record` and `private-evaluation-record`;
the R figure retains every case and template stratum.

Accessible-boundary correction implemented as an explicit opt-in:
`contact-penalties --surface-context --surface-ownership accessible_boundary`.
It reuses native exposed samples and SciPy batched nearest-boundary queries;
exact nearest atom-owner ties split area. Twenty-seven installed-wheel focused
checks pass, including unchanged default output, enclosed-atom exclusion,
conservation and end-to-end provenance. Ruff and warning-fatal Sphinx pass.
The subsequent full 3808-model comparison and independent final-product checks
completed successfully. Original scores and atom exposures are unchanged; all
boundary owners are exposed. Predictive gains remain case-specific and losses
are retained. The follow-up attribution on completed tables identifies changes
in raw residue-table chemistry as the main mechanism behind the selected rank
losses, with separate unchanged-logZ × changed-area partition effects.

20 October author subgoal: use contact information beyond a closest-heavy-atom
5 Å residue edge to distinguish modeled supported and mock receptor–pMHC pairs.
Reconcile generator/label/native lineage before attributing wrong contacts; keep
the original score and distinguish an atomic observation layer from the existing
residue-matrix surface score. Private subroadmap and completed-output mechanisms:
`private-evaluation-record`.
No new structure extraction or fitting was performed. Atom-reference opportunity,
context transfer and independent classification validation remain separate gates.

## 20 October: sparse contact energy and complete task audit

Follow [CONTACT_ENERGY_SUBROADMAP.md](CONTACT_ENERGY_SUBROADMAP.md). The author
requires every current energy task and legacy 2022 benchmark, retaining bad
outcomes and original potential identities. The completed residue+atom candidate
improves ILA1/YLQ and loses SB27; it does not restore all CPL classifications.
Chemical pooling, antibody transfer and exchange symmetry must address support
before introducing more atom/residue/context cells. No assay-fit or default change.

Contact-export scheduling now defaults to all allocated CPUs, parent-only batch
annotation, coarse contiguous scoring workers and ordered reduction. An installed
immutable-wheel128-model diagnostic measured28.331s serial versus8.428s on16CPUs
(3.361× cold speedup),76 peak live OS threads and4319936512bytes sampled peak
process-tree RSS. All256 score rows and complete contact/residue outputs match
serial and frozen originals exactly. Seventeen focused checks, lint and
warning-fatal docs pass. This is a measured subset, not a whole-cohort timing claim.
