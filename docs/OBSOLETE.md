# Archival disposition — 2026-10-14

`research archive: OBSOLETE.md` preserves the earlier measurements and removal proposals.
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

## Disposition — 22 October 2026

Acted on: `integral.py` and `integral_glm.py` were **removed**, with their four CLI commands,
tests and `docs/integral.rst`. Their energy component carried no alanine reference, was
undefined wherever an interface had no TCR:peptide contact, and divided out the contact-burden
term; both fitted GLMs inherited that input. That was a wrong solution, not merely a retired
one, so it was withdrawn rather than archived as a supported path. Its measurements stand.

Checked and **not** removed: every module under `src/tcren/contacts/` has at least two callers
in `src` or `tests`, counted 22 October — the lowest are `class_audit`, `context_normalization`,
`prior_geometry`, `surface_reference` and `weighted_summary` at two each, and four of the set
back a live CLI command. Nothing there is orphaned, so the contact-plausibility campaign's
machinery is retained under the standing rule that retirement travels with caller cleanup rather
than preceding it. Their measured records are research evidence and are unchanged.

Release cleanup: removed six gitignored wheels/source archives for versions 2.15.0–2.20.0
from `dist/` (11,336,132 bytes). They were superseded build output; accepted scientific
records, historical source and required runtime/native assets are retained.
