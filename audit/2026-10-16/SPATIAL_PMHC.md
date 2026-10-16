# Exact spatial pMHC exports

`export --quantity pmhc-cells` reuses the existing batched annotation, MHC region
markup and surface-map APIs. The receptor is removed before groove-frame fitting.
The CLI owns occupied-cell TSVs and the `.summary.tsv` sidecar; provenance hashes
that summary and records the grid, extent, probe, smoothing, KD scale and ridge.
The existing `pmhc` export and every scoring API remain unchanged.

`spatial_features` reduces common finite peptide-owned cells with grouped Polars
expressions. It reports population Cov(height, KD) and Cov(max(height-ridge,0), KD),
component means, conditional above-ridge chemistry and counts. The ridge remains
the existing 90th percentile of helix heights. Both covariances have units Å times
KD index. E(excess*KD)=E(excess)*E(KD)+Cov(excess,KD) is checked numerically.
No fitted coefficients, normalization, chemistry threshold or label-selected sign.
A zero ridge-excess covariance for a wholly below-ridge field is defined; an absent
ridge or above-ridge conditional chemistry remains missing. Duplicate cells,
missing keys and inconsistent metadata fail.

Thirteen focused tests passed, including an independent numerical covariance
oracle, constant/missing fields, class separation, receptor-translation invariance,
existing export behavior and exact unchanged association metrics/membership.
Ruff passes on src/tests; Sphinx builds with warnings treated as errors.
A fresh five-set application exported 4,187 structures and 7,610,879 occupied cells,
including separate native class-II evidence. All 18,140 original score/ipTM metric
entries and their exact per-score memberships were preserved. Every reduction was
independently reproduced from its authoritative exported cells, with 1e-12 absolute
floating-point tolerance. The original scalar-product findings and weaker spatial
associations remain inspectable; these descriptors are not replacement scores.

The measured coordinate/association/R-summary build used 237.649 seconds and
3,567,353,856 bytes peak sampled process-tree RSS under a four-CPU allocation.
All sampled MMseqs processes respected <=4 threads and the requested search budget;
queries ran through the existing whole-set batch APIs with no outer worker pool.
Twenty peak live threads include bounded auxiliary/I/O runtime threads, not twenty
concurrent numerical kernels. Numerical reductions are serial. Vectorized table
operations and the existing native surface API sufficed; no new native kernel.

Removing a TCR does not reverse induced fit. Smoothed cell chemistry can mix
neighbouring peptide/MHC values. Neither a single bound-state surface nor assay
negativity establishes apo strain, loop formation kinetics or a forced coordinate
pose. Continuous geometry/confidence/packing diagnostics remain separate from an
independently validated pose-quality detector. No new detector or ranking boost
is claimed. Native asymmetric-unit inputs require explicit single-complex selection;
the CLI deliberately rejects ambiguous multiple-peptide assemblies.

Stratos MD was not accessed or modified. Publication/release approval is separate.

Broader local verification remains explicitly partial. A sandbox-only fast-suite
attempt hit ARDA reference lock permissions (1,028 passed before reporting those
errors). The authorized reference-backed retry reached its 360-second cap during
historical Potts self-consistency tests, after 755 progress entries and one
observed failure. The isolated failure was a missing canonical database at the
installed wheel's default data root; rerunning that exact canonical-frame test
with `TCREN_HOME` set to the project root passed in 38.13 seconds. No broader-suite
pass is claimed, no historical test was removed, and no fitted scientific stage
was run as an analysis. Exact-head lean CI (without ARDA references), both Python
versions, lint/docs, the 13 focused tests and full spatial-export oracle passed.
