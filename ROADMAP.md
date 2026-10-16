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
| CLI fixed; campaign pending; issue [26](https://github.com/antigenomics/tcren/issues/26) | Installed `md inspect`, `md quality` and `md compare-forces` expose existing exact readouts, atom/component provenance, class-separated exports and CPU limits. 64 focused tests, installed-wheel class-I/II examples and manual table inspection pass; [MD_DIAGNOSTICS.md](audit/2026-10-16/MD_DIAGNOSTICS.md) records acceptance. Stratos has an active original eight-thread MD job: one 100 ns arm is complete, the other remains running; final paired acceptance and missing historical campaign inputs stay unresolved. Do not restart the simulation or relabel reconstruction as exact reproduction. |

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
