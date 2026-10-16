# Publication and software roadmap

Reconciled 16 October 2026 against master base `8f860a2`, the issue-21 fix, live GitHub issues and
recorded benchmark/manuscript commits `8657b89` / `26b912c`.
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
| P2 / benchmark + manuscript; issue [34](https://github.com/antigenomics/tcren/issues/34) | Exact position/residue/contact/placement CLI exports exist and pass API equality checks. Complete the per-panel audit showing every retained molecular rendering consumes authoritative exports, including physical residue IDs, scene/potential/geometry hashes and visual review. A successful PDF build does not close this lineage requirement. |
| P2 / manuscript; issues [33](https://github.com/antigenomics/tcren/issues/33), [35](https://github.com/antigenomics/tcren/issues/35) | Fit-free disposition, archives, reproduced comparisons and active-paper cleanup are recorded. Reconcile these issue acceptance criteria with the final producer/include graph; audit any retained historical methodological equation separately. Do not implement replacement fitted models merely to close retired CLI gaps. |
| P2 / library + benchmark | Four literature-topology columns have recorded differences up to 1.11e-16 from unordered reductions, which can alter tied ranks. Keep all deltas visible. A deterministic-reduction fix needs focused checks and explicit authorization for its affected scientific rerun; do not round scores to suppress differences. |
| Completed / library; issue [21](https://github.com/antigenomics/tcren/issues/21) | Notebook/stdin calls now use an importable subprocess entrypoint for the existing bounded spawn pool, after whole-set annotation. Real Jupyter, stdin, serial/all/explicit budgets, worker failure/cancellation, physical-key comparisons and manual class-I/II table inspection are recorded in [INTERACTIVE_BATCH_FIX.md](audit/2026-10-16/INTERACTIVE_BATCH_FIX.md). Final affected suite: 58 passed; the full local reference-backed suite hit its 360-second cap and remains incomplete. Preinitialized multithreaded Polars remains unsupported. |
| Historical/deferred; issue [32](https://github.com/antigenomics/tcren/issues/32) | Custom Potts geometry forwarding and approximation labels were repaired. Complete CLI equivalence for historical decompositions is absent; no final-paper dependency authorizes a new fitted scorer. Preserve the original quantities and identify the outstanding gap. |
| Separate campaign; issue [26](https://github.com/antigenomics/tcren/issues/26) | Component RMSD, frame selection and close-pair tools exist. The long MD reconstruction and missing original campaign inputs remain unresolved. Verify resources/input provenance before scheduling; reconstructed trajectories cannot become exact historical reproduction by relabelling. |

Six existing issues (#21, #26, #32–35) were open and no PR was open at the
initial checkpoint. This audit added issue [43](https://github.com/antigenomics/tcren/issues/43)
for the demonstrated CLI gaps. The subsequent [fix](audit/2026-10-16/CLI_BUDGET_FIX.md)
resolves #43; the [interactive fix](audit/2026-10-16/INTERACTIVE_BATCH_FIX.md) resolves #21.
Five other issues (#26, #32–35) retain their independent scope above. Open/closed GitHub state is distinct from the implementation and
publication evidence above. The current example audit does not certify every CLI,
reproduce the long MD campaign or authorize paid review, fitting or publication.

## Submission exit

Freeze one reviewed source/wheel/environment identity; ensure each active figure
and number has its raw-input → CLI → summary → render lineage. Run structural,
xref, font and bounds gates; inspect the actual final PDFs. Verify data/code
availability and deposits within their access boundaries. Finish release approvals
and verify the public install before claiming a release. There is no date promise.
