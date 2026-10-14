# Publication roadmap

Updated 2026-10-14. The goal is a reproducible, interpretable paper supported by the existing five
benchmark blocks, not a larger catalogue of scores. Prior plans and measured abandoned approaches
are preserved in `attic/2026-10-audit/ROADMAP.md`. `AUDIT.md` owns findings; this file owns exit criteria.

## 1. Correctness and bounded scope

Fix demonstrated scoring, native-boundary and pipeline defects with focused synthetic checks.
Archive only complete obsolete dependency groups, retaining source, resources and result records.
Keep the accepted potential and fitted arrays unchanged until a deliberate model revision is approved.
Exit: focused tests, full applicable unit suite, lint and docs pass; every remaining defect is named.
Current: Potts identity/sampling, Gaussian/mechanics domains and source-boundary repairs are
implemented and checked. The catalogue/frozen metadata migration preserves every coefficient
array; corrected-feature equivalence still requires gate 4. Passing software tests do not close
the scientific or publication gates.

## 2. Settle definitions before changing numbers

The author approved correcting the ergodicity and naive-Bayes interpretations and retiring tree
figures/claims. Tier-1 documentation now identifies shared preprocessing fitted on both classes;
the benchmark GLM is explicitly separate from the shipped conditional residual.
Correctly distinguish low-order Betti difference from full Euler characteristic and approximate
contact stability from an exact physical probability. Preserve every measured result and recorded
comparison. Do not add new metrics or fitted alternatives to address a wording problem.
Exit: one definition per readout, matched equations/code/prose, approved claim scope.
Current: source definitions reconciled without refitting or replacing empirical values. Final
rendered-PDF review remains in gate 5.

## 3. Close the publication dependency graph

For every main/SI number and figure, record raw inputs → library version/model → benchmark producer
→ reshaper → rendered artifact. Resolve the benchmark's old wheel pin and unwired confidence
producers. Fail on absent required input, stale model catalogue, failed stage or reduced corpus.
Evaluation never installs fitted models into the library. Separate restricted inputs from public
availability claims without exposing private provenance.
Exit: complete producer/input manifest and one tested, immutable software/environment identity.
Current: all 34 live reshapers have scheduled owning producers/inputs; tree routes are archived.
Preflight checks the actual wheel payload and 12 required raw sets. Historical frozen-model fit
traces are explicitly historical parameter evidence. See
`audit/2026-10-14/publication-resolution.md`; only execution can certify the complete graph.

## 4. Reproduce only after authorization

Run one complete fresh pass from raw sources after gates 1–3 and explicit authorization. Compare
matched units against the accepted recorded results, including row counts and missingness. Stop and
report material weakening; never silently replace a baseline, refloor a test or average away a loss.
No sweeps, additional controls or model-selection exercise is implied by this roadmap.
Exit: all five blocks and every publication artifact trace to that run; no resumed intermediates.
Current: full fresh reproduction authorized; the latest attempt stopped on an obsolete producer
CLI call. Producer/API and cached-table repairs precede another complete raw-input restart.
Both full-manifest frozen-parameter comparison and matched accepted-result gates precede
publication. Earlier interrupted diagnostic evidence is retained but does not complete this gate.

## 5. Submission and release

Run manuscript structural/xref/figure gates, render and inspect both PDFs, and verify equations,
units, denominators, citations and availability statements. Build/install the wheel independently,
check packaged resources and CLI behavior, then use the existing CI/release approval process.
Exit: clean build/tests, reviewed PDFs, verified deposits and approved release. No date promise or
claim of publication readiness until all gates pass.

## Frozen/deferred

No sequence-only geometry ladder revival, new Potts ranking composite, speculative descriptor
family, unmeasured compute-budget recommendation or broad framework rewrite. Existing MD contact
and self-consistent potential paths stay only for their identified current consumers. Historical
oracle paths remain available for the published work they reproduce.
