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
