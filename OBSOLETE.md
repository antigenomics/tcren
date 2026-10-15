# Archival disposition — 2026-10-14

`attic/2026-10-audit/OBSOLETE.md` preserves the earlier measurements and removal proposals.
The following decisions follow actual imports and publication consumers, not catalogue membership.

| Material | Disposition | Reason |
|---|---|---|
| Potts geom/chain/cross, three geometry resources, their CLI paths and dedicated tests | moved to attic | retired geometry experiment; no current manuscript consumer after removing obsolete figure stage |
| Shared geometry helpers, contact map agreement | retained | active contact-map/site APIs and independent tests |
| selfcon, shipped self-consistent potential/contact model | retained | current figure and regeneration dependencies |
| peptide–MHC contact model | retained | current MD benchmark consumes it; old “unread” label was false |
| ensemble/EM paths | retained pending separate disposition | real CLI/benchmark callers exist; not deleted merely because no descriptor uses them |
| Historical potential and redundancy alternatives | retained | published oracle and explicit caller contracts; retirement requires matching caller cleanup |
| Research result tables and figures | preserved | archival of a method does not erase its measured record |

The potential-specific provenance registry remains unchanged. Its `not_in_tcren2` entries below
mean “not in the current potential estimator”, not “unused software”.

- `tcren.potential.derive.derive_tcren_loo` --- Diagnostic leave-one-out derivation.
- `tcren.potential.model.tcren` --- Published 2022 potential, retained as an oracle.
- `tcren.potential.redundancy.cluster_weights` --- Explicit alternative clustering weights.
- `tcren.potential.redundancy.epitope_weights` --- Explicit single-axis weighting alias.
- `tcren.potential.redundancy.nonredundant_ids` --- Explicit cluster-member exclusion.
