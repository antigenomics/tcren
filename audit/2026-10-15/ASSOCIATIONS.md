# pMHC feature–performance associations

Author authorization: implement and run associations using recorded scores, derive
pMHC features from coordinates, keep class and template-coverage strata separate,
and retain every matched loss. No fitted replacement score.

## Implementation → verification

1. Extend existing exact exports with pMHC-only surface/SASA quantities, using
   existing surface maps and Biopython SASA; strip receptor before defining the
   frame. Verify receptor independence, units and class-II handling.
2. Add an association API/CLI using installed sklearn metrics and SciPy ranks.
   Require explicit keys, strata, score directions and feature columns; preserve
   pair-specific finite support and undefined metrics. Verify ties, missingness,
   duplicate joins and strict stratum isolation with synthetic examples.
3. Use existing verified receptor/CPL result tables as authorized inputs. Derive
   only the new surface quantities from raw coordinates. Record hashes, lineage,
   missing support and matched counts in the private benchmark deliverable.
4. Document formulas, CLI examples and assumptions in Sphinx, README and API skill;
   run focused checks, bounded non-slow tests, lint and warning-as-error docs.
5. Reconcile old branch tips without restoring retired fitted scorers. Preserve
   unique historical commits with archive tags before deleting obsolete branches.

The independent analysis units are cohorts/pMHCs, not models. Receptor-cohort
features are label-blind medians over the supplied model pMHC conformations,
with within-cohort spread retained. These are modelled bound-conformation
descriptors, not measurements of unliganded pMHC or causal effects. CPL cohorts
contain different peptides and are analysed separately as peptide panels.

Continuous feature associations are descriptive Spearman correlations with each
score's performance and its paired difference from ipTM; no p-values, fitted
weights, score inversion selected from labels, or surface threshold scan.
Class-II structures without labelled ipTM comparisons contribute only to feature
coverage, never to inferred class-II discrimination.

Progress and branch reconciliation evidence will be appended after verification.

## Branch reconciliation

The library now has only `master` and `codex/pmhc-associations` locally.
Git-confirmed merged local branches `codex/fit-free-release-fixes`,
`manifest-alignment` and `selfconsistent-tcren` were removed. Exact historical
tips remain under these archive tags, all pushed and verified on origin:

| Tag below `archive/2026-10-15/` | Commit | Disposition |
|---|---|---|
| annotation-index-provenance | e089e48 | Tree exactly matches squash merge b593c42; current master adds verification. |
| verify-pypi-400 | 42a4132 | Release verification text is already retained in master. |
| md-validation | 54535ab | Newer library implementation retained; missing wrapper, docs and eight tests restored and passed. |
| docking-zscore | 83be6e4 | Historical fitted covariance scorer, excluded by final-model contract; preserved without revival. |

Remote PR inventory was empty. After automatic approval review required explicit
authorization, the author approved these exact 12 deletions. Archive tags were
pushed first; all 12 remote deletions then succeeded. `git ls-remote --heads origin`
confirmed only master remained immediately after cleanup:

- codex/annotation-index-provenance
- codex/fit-free-release-fixes
- codex/md-validation
- develop
- feature/appendix-methods-reference
- feature/clash-register
- feature/cpp-ranking-ddg
- feature/fix-alanine-scan
- feature/poly-ala-reference
- feature/review-aug26-surface-topology
- feature/shuffle-decoys
- feature/tcrdock-geometry

All except annotation-index-provenance and md-validation are Git ancestors of
master. Those two are preserved by exact archive tags and their content
dispositions above. The new association branch is the only active feature line.

The benchmark's visible feature branch already includes the receptor-reproduction
branch's changes plus later fixes. Its separate checkout has uncommitted changes;
that checkout is preserved. No dirty scientific result or generated figure was
included in the library commits.

## Verification and completed computation

The installed wheel completed new pMHC exports for all 3,809 ranking models:
618 TCRvdb, 1,089 VDJdb and 2,102 CPL. Existing score values were joined without
recomputation. Each receptor dataset has 143 arms; CPL has 131, including its
explicitly historical descriptor snapshot. Eleven pMHC features and four metrics
were evaluated with pair-specific membership. Complete local outputs, input
hashes and commands belong to the benchmark's `results/pmhc-associations-20261015`.

The broad bounded installed-wheel check had 1,105 passes, three skips and one
environment-dependent failure: canonical data were absent at the default cache
root. The failing test passed in 16.86 s with the configured `TCREN_HOME`.
The original broad run took 246.91 s. Scientific refitting/reproduction stages
were excluded. Focused association/export/CLI checks, restored MD tests, Ruff
and warning-as-error Sphinx builds passed.

Final review added an explicit heavy-atom SASA guard and a regression check:
parser defaults retain hydrogens, whereas the documented feature excludes them.
A raw-coordinate scan found hydrogens in 10 native structures and in none of the
3,809 ranking-model inputs. The complete native inventory was rebuilt in 62.1 s
with the final wheel (278 class I, 94 class II, two accounted validation failures); the ranking associations retain their original producing
wheel and are unaffected by this guard. No old feature rows are reused.

Final installed-wheel focused validation: 17 tests passed in 1.37 s. Native SASA
changed for eight successful hydrogen-containing structures; surface features and
ranking results are unaffected. All row-level scientific outputs remain local
because the TCRvdb endpoints are restricted. The benchmark report retains the
matched gains, losses, effect sizes, counts and explicit class-II evidence scope.
