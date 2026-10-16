# Label-independent joint pMHC context

The `associate` CLI now accepts repeated `--feature-product LEFT:RIGHT` options.
It uses vectorized Polars expressions to form raw Float64 products before cohort
aggregation, exports `joint_features.tsv`, and records each source pair in provenance.
The original scores, score direction, matched rows, metrics and template/class
strata are unchanged. No labels, fitted weights, centering, normalization or new
thresholds enter the feature definition. Missing/nonfinite inputs or overflow
produce a missing product; no observation is silently imputed.

The publication analysis declared relief × exposed SASA-weighted hydropathy as
primary and protrusion fraction × exposed hydropathy as secondary before outcomes.
Their units are Å × Kyte–Doolittle index and fraction × that index. Cohorts use
medians of per-model products, never products of marginal medians. Joint and
component associations share the same finite cohort support for each score/metric.
All weakening associations and sign changes remain recorded. Scalar products do
not establish spatial co-localization, an independent interaction or a ranking boost.
Class-II native evidence remains separate from modelled class-I classification.

Verification: all five association tests pass, including exact unchanged metrics
and membership with the CLI product option and synthetic class-II strata. Ruff
passes on src/tests; Sphinx builds with warnings as errors. The immutable wheel
was byte-verified across 183 installed payload files; SHA256:
`d1ae0e8f7ee9cf12af9632669ce8401a23412f6cd3fb2de5d4162fa0c81fe2a5`.
The recorded-table analysis and R summaries completed in 20.363 seconds with one
numerical CPU and peak sampled process-tree RSS 513,589,248 bytes. Fourteen live
threads included bounded auxiliary I/O/runtime threads; no search engine or worker
pool was launched. Vectorized table operations sufficed; no new native code.

No structure recomputation, sequence search, fitting or remote access occurred.
Stratos MD was neither accessed nor modified. Release publication and its human
approval gate remain separate.
