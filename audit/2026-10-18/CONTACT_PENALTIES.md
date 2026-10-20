# Frozen contact penalties — 18 October 2026

Current state: the reference resource is frozen and the complete matched benchmark
finished and passed independent contribution/accounting checks on Stratos. The author's approved model uses normal transformed
marginals, actual external observation counts as prior strength, declared pooled
chemistry backoff and the worst supported atom/marginal per physical residue edge.
The four penalty formulas were fixed before accuracy inspection. Definitions and
CLI commands: `docs/contact-penalties.rst`; source task: manuscript `plan181026.txt`.

## Reference and implementation evidence

- External deposit: 856751 rows, 529 exact repeated atom keys removed, 1275 rows
  excluded through the complete 374-ID native inventory. The selected reference
  contains 854947 atom pairs from 1831 deposited complexes. No overlap was silently
  treated as independent evidence.
- Current raw mmCIFs were already downloaded by the preceding authorized Stratos
  job. Exact float32 amino-acid/atom/coordinate correspondence recovered both
  endpoints for 775760 atom pairs; 79187 retain missing recovered geometry.
  IMGT and current author numbering were never equated. Raw files are processed
  one at a time; no high-dimensional angle model, guessed atom or numerical floor.
- The corrected full external export completed in 228.43 s, peak sampled process
  tree RSS 886013952 bytes. The fresh native export accounted for 370 inputs:
  357 successes, 13 explicit eligibility errors, 339921 atom pairs. It completed
  in 75.51 s, peak sampled tree RSS 2505129984 bytes. Existing raw inputs and old
  observation results were preserved.
- The frozen resource contains 15464 exact marginal groups and 31 chemistry
  pools; 777 groups use pooled external prior support. It records 14552 unordered
  atom-chemistry categories separately from geometry. External/native marginal
  observation counts are 1255607/506524; these are correlated observations, not
  independent complexes. The descriptive reference check is not held-out
  calibration. Unsupported categorical support remains unavailable.
- Resource SHA256:
  `63c141f0a8be2e090cda543193148dbd8b06f937dcaf5f4f428124c77d013f16`.
  Visible local artifact: `data/contact-penalty-run-20261018/frozen-contact-prior-c.json`.
  No scoring-time updating, assay-fitted coefficient, cohort normalization,
  confidence mixing or retired fitted readout enters scoring.
- Scoring wheel SHA256:
  `2c8736fd735fc33e813f64178bb373fd32ce03e0ea94929144d0aaea31aa654d`.
  The driver verifies every installed package file against this noneditable wheel.
  Both the 2022 TCRen and TCRen2 tables retain their original all-TCR peptide
  baseline. Only CDR1/2/3:peptide residue edges receive an increment; framework
  contacts remain explicitly unassessed. Missing potential entries stay null and
  retain the existing scorer's omitted-entry total, with exact counts exposed.
- The six-model complete raw-source diagnostic passed; original TCRen2 energy
  differences from the recorded producer were exactly zero on those six models.
  Its sampled tree RSS peaked at 1339981824 bytes and it completed in 9.97 s.

## Verification and corrected defects

The focused gate passed 81 tests (3 slow tests deselected), 2.41 s including
startup. Scope: new prior/penalty tests plus contact observations/types, exact
exports, associations, atomic scoring, original scoring and contact-map geometry.
Ruff `check src tests scripts/run_bounded.py`, `git diff --check` and warning-fatal
Sphinx HTML passed. This is the completed focused gate, not a full-suite claim.
The independent review found two boundary issues, both repaired and rechecked:
coordinate/label MHC-class disagreements now fail, and a stopped/exited root cannot
leave its process-group descendants running. Synthetic SIGTERM-ignoring-child
and early-root-exit checks terminated all owned children within 5.47 s.

Other failures were diagnosed before full accuracy evaluation: the isolated
checkout initially lacked the pinned MHC FASTA; Polars bypassed a gzip handle;
the first selected atom could hide a recoverable edge orientation; empty integer
index arrays needed an explicit dtype; and one provenance hash received a string
instead of a Path. Regression checks cover these code paths. The corrected
external producer reread all selected raw structures. No invalid export, failed
fit, incomplete diagnostic or earlier resource is an accepted benchmark result.

## Completed matched run

The immutable input ledger pins 3808 unchanged models: 2102 CPL, 617 TCRvdb and
1089 VDJdb. Each selected archive member/direct mirror file is SHA256-verified
before scoring. Native ipTM and existing labels are retained. All four arms,
original energy and ipTM share one finite row set within each potential/cohort/
MHC-class/template stratum. Directions and established ROC/AP/precision-at-recall
metrics are unchanged; no additional controls, sweeps or tests of significance.

Remote project: `private-evaluation-record`.
Job: `tcren-contact-penalties-full-20261018`. Output: `full/`; guard:
`full-guard.json`. One compute thread, 32-model annotation batches, sampled 8 GiB
process-tree guard and 1800 s time cap. The guard is sampled, not an OS hard
allocation limit. Each annotation call contains at most 32 models; the preceding
typed batch can coexist during handoff (at most 64 structures plus the current
per-model tables), and buffers do not accumulate across the cohort.
The complete run exited zero in 1069.85 s (17.83 minutes), with peak sampled
process-tree RSS 1458225152 bytes (1.358 GiB). All 3808 models were retained,
including 139 zero-contact models; eight contain omitted potential entries and
none has missing supported assessed geometry. All evaluated models are MHC-I.
Independent verification checked 88073 unique physical pairs, 176146
potential/contact rows and 7616 model/potential rows. Every original/four-arm
sum conserved its contributions and every finite-row membership matched exactly.
Fresh original TCRen2 energy differs from the accepted producer by exactly zero
on all 3808 models. The verifier reread outputs without rerunning scientific scores
or accuracy metrics. Local retrieval passed a remote checksum dry run.

Visible benchmark report:
`private-evaluation-record`.
It retains all 31 cohorts, both potentials, all four existing metric endpoints,
all arm/baseline/ipTM deltas and detailed contributions, with exact counts.

Measured AUROC examples (0–1 scale): original 2022 TCRen literal flip raises GLC
0.360013→0.642498 (194 models, 135 binders), but lowers YLQ 0.593062→0.452938
(423 models, 174 binders). ipTM remains 0.816698/0.776993 respectively. All three
nonnegative arms lower AUROC on both epitopes under both tables. TCRen2 favourable-only
flip raises scarce-template RAQAPPPSW/HLA-B*57:01 0.454167→0.866667 (44 models,
20 binders; ipTM 0.047917), and CPL 868 0.599728→0.665373 (321 models,
161 positives). These are cohort-specific gains alongside the fully retained losses;
no winner or default is promoted. Literal flip can reward original positive energies,
so its GLC gain does not show that a nonnegative penalty resolves the inversion.

The implementation used vectorized tables/arrays and existing native searches;
no new native kernel, outer worker pool or computed-part cache was introduced.

## Limits and retained goals

Predictive compatibility is not P(natural), nonbinding probability or binding free
energy. Approximate normal marginals discard signed torsion handedness; the worst
supported reduction gives larger residues more penalty opportunities. Current
raw correspondence is incomplete and reference contacts are correlated. The
reference includes both native MHC classes; modelled class-II ranking is not
established by it. Full angular hydrogen-bond chemistry, line-of-sight and the
remaining PPI/atlas quantitative priors remain separate open goals. These choices
were declared before this comparison, rather than selected by assay accuracy.
Labels already informed the broader investigation, so model/arm selection on
this benchmark is development evidence. No standing manuscript claim or default
score is promoted or rewritten by this run. MD and release gates remain separate.
VDJdb covered/scarce retains the existing epitope-level
`native_structure_exists` metadata; it does not assert exact-template inclusion
for an individual prediction.
