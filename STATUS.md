# tcren status — 2026-10-14

The checkout declares version 3.10.0. Release/PyPI status was not verified in this audit.
All three integration checkouts contained substantial existing uncommitted work; it was preserved.

The current audit is in `AUDIT.md`; its verification record distinguishes passing checks from
unverified scientific claims. `ROADMAP.md` defines publication gates. This checkout is not yet
certified publication-ready, and no benchmark rerun or fitted-array replacement was performed.

The sequence-only geometry prior, coordinate-chain and cross-product modules/resources/tests have
moved to `attic/2026-10-audit`. The corresponding benchmark experiments are in that repository's
attic. Recorded results remain intact. Contact probabilities on supplied coordinates, their fitting
and self-consistent regeneration, and MD consumers remain active.

Previous release descriptions, numerical records and open-loop history remain recoverable in
`attic/2026-10-audit/STATUS.md`; they are not assertions about today's live release state.

Verification: the broad final run passed 948 tests; its single obsolete serialization assertion
was corrected and the full affected Potts suites passed 78 tests. Lint, warning-fatal docs,
wheel/install and source-archive checks passed. See AUDIT.md for exclusions and unresolved defects.
