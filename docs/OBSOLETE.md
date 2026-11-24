# Historical potential APIs

These explicit APIs are not part of the bundled TCRen2 potential estimator.
They remain available for reproducing their declared reference and redundancy
conventions; absence from the current estimator does not imply an unused API.

- `tcren.potential.derive.derive_tcren_loo` --- Diagnostic leave-one-out derivation.
- `tcren.potential.model.tcren` --- Published 2022 potential, retained as an oracle.
- `tcren.potential.redundancy.cluster_weights` --- Explicit alternative clustering weights.
- `tcren.potential.redundancy.epitope_weights` --- Explicit single-axis weighting alias.
- `tcren.potential.redundancy.nonredundant_ids` --- Explicit cluster-member exclusion.

# Deprecated paper bootstrap options

Remove extra keyword arguments passed to `tcren.paper.bootstrap`. Extra
keywords are ignored and now raise a `DeprecationWarning`. Use `structures`
to control the structure download and `canonical=True` to include the
canonical set. The legacy comparison data remain in the reproduction
workflow.
