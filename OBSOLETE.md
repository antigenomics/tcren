# OBSOLETE — scheduled for deletion

Everything here is marked `@not_in_tcren2` in the source: correct, tested, and taking no part in
producing the shipped matrix. The recipe is
`tcren derive-potential --structure-dir Native2026 --balance both`, at the 5 Å heavy-atom
contact cutoff `contact_table` uses (the subcommand has no `--cutoff`), and nothing else. **Delete this list, and the code behind it, at a major version.** v3.0.0 came and went without it: that release was for the new public surface and the new package data, not for a removal. Keeping four
redundancy schemes and three reference states around is how three mutually inconsistent potentials
came to be in circulation at once.

Regenerate with `python -m tcren._provenance`; `tests/unit/test_provenance.py` fails if this file
and the markings disagree.

- `tcren.potential.derive.derive_tcren_loo` --- Leave-one-out derivation, for testing how much any single structure moves the matrix. Diagnostic, not a production path.
- `tcren.potential.model.tcren` --- The 2022 matrix, kept for reproducing published results. TCRen2 is tcren.potential.tcren2() and is the default since 2.11.0; the two correlate at r = 0.867 with a maximum absolute difference of 0.943 and are not interchangeable.
- `tcren.potential.redundancy.cluster_weights` --- Down-weights by sequence-distance clusters rather than exact identity. Needs a threshold, and conflates the epitope and receptor axes that --balance separates.
- `tcren.potential.redundancy.epitope_weights` --- TCRen2 balances the epitope AND receptor axes; this single-axis alias is what the manuscript's earlier matrix used. Receptor redundancy is the larger of the two on Native2026 (223 structures share a receptor against 212 an epitope).
- `tcren.potential.redundancy.nonredundant_ids` --- Excludes cluster members outright. TCRen2 down-weights instead, which keeps every structure's data.

## What deletion costs

| Group | Also delete | Kept by |
|---|---|---|
| Redundancy alternatives | `potential/redundancy.py` less `balanced_weights`, the `--nonred`/`--redundancy-t` CLI options | nothing |
| Leave-one-out | `derive_tcren_loo`, the `--loo` CLI option | nothing |
| The 2022 matrix | `potential.tcren()` and `data/TCRen_potential.csv` | published results predating TCRen2, and `-p karnaukhov2022` |

The DFIRE reference states and the substitution-matrix pseudocounts that were on this list are
gone: all five estimation improvements tried against TCRen2 were measured and none improved every
endpoint, so they were deleted rather than left to rot. See
`src/tcren/data/TCRen2_potential.NOTES.md` before re-attempting any of them.

## Already removed in 2.12.0 — the combiner zoo `P_native` replaced

Not scheduled: gone. Every one had zero callers in the library and zero in the benchmark repo's
reproduction path. Recorded here so a name that no longer imports can be traced to what took over.

| removed | what it was | what to use instead |
|---|---|---|
| `pose_sweep` module (605 lines) | pose-consistency experiment | nothing — `P_native` made it unnecessary |
| `pose.c_score` + `pose_af_reference.csv` + `pose_native_reference.csv` | manifold-referenced pose score; **492 KB off the wheel** | nothing. Its own docstring recorded why: scored against the crystal manifold it reads *provenance*, not model quality. `pose.pose_consistency` and the `POSE_FEATURES*` tuples are unchanged |
| `footprint.footprint_score` | the `fp_score` z-sum | `reliability.t_score`, which is what `tcren footprint --score` emits since 2.26.0, or the `shape` entry of `tcren.score.channel_scores`. 2.12.0 pointed here at `cohort.p_native`, itself removed in 2.26.0 |
| `cohort.q_iptm`, `q_f`, `q_f_iptm`, `f_invert_by_iptm`, `phi_bind`, `agreement` | hand-picked combination rules; `agreement` was the per-structure summand of `C*` | `tcren.score.binder_score`, and `binder_iptm` where a confidence is available. 2.12.0 pointed here at `cohort.p_native`, itself removed in 2.26.0 |
| `recognition.kit_score` | a z-sum of `p_bind` and ipTM | the `binder_iptm` column of `tcren assess` — `binder_score` plus `logit(ipTM)`, two log-odds added with no coefficient to fit. 2.12.0 pointed here at `cohort.p_native`, itself removed in 2.26.0 |
| `scripts/fit_pose_reference.py`, `scripts/fit_joint_reference.py` | regenerated the two deleted reference CSVs | nothing |

**Deprecated, not removed:** `cohort.coupling` and `cohort.q_coupled` remain importable, tested and
byte-identical in behaviour, so every number published under the 2.12.0-era `S` — which was
`q_coupled(Q, ΔΦ)`, not today's composite — still reproduces. The current composite is
`reliability.s_score`. `tcren.binder`'s fitted `p_bind` and `recognize --scores`, which this entry
once recorded as kept for v1 reproduction, were **removed in 2.26.0**; see below.

## Already removed in 2.26.0 — `P_native` and every frozen-coefficient composite

Not scheduled: gone, on the author's ruling to discard rather than deprecate and rebuild from a
recompute. `P_native` refitted a latent class on every call, raised when a cohort had fewer rows
than features, and its value depended on which rows the fit was anchored on — none of which
survives contact with a user holding one model. The v1 block's coefficients were frozen against
training sets that no longer exist, which made them the one part of the package a reader could not
reproduce.

| removed | what it was | what to use instead |
|---|---|---|
| `P_native`, `cohort.p_native`, `P_NATIVE_CHANNELS` / `_POOL` / `_ORIENT` / `_FEATURES` / `_BANNED`, and the per-channel posteriors `G` / `T` / `E` | a latent class over three channels, refitted on every call | `tcren.score.binder_score` and `channel_scores`, through `tcren assess`. Both are frozen on a hold-out that ships in the wheel, and `tcren fit-holdout` regenerates the shipped arrays from its manifest bit for bit |
| `p_bind`, `p_forced`, `q_bind`, `s_strain` as emitted columns, with `recognize --scores` and `--features-only` | the v1 score block | `tcren assess`. `cohort.q_score` and `cohort.strain_z` remain as functions; what went is the frozen-coefficient layer around them |
| `FORCED_POSE_MODEL`, `forced_pose_score` | the fitted forced-pose read-out | `reliability.inversion_flag`, the forced-pose detector, which fits nothing |
| `tcren.binder.binder_score`, `BINDER_MODEL`, the `tcren binder` command | the frozen 5-feature binder logistic | `tcren.score.binder_score` — the same name on a different object. `tcren.binder` now holds only `is_real_interface`, the pre-energy check that an interface is a plausible dock at all |
| `p_real`, `p_real_bn`, `real_probability`, `frozen_recognizers`, `encode_features`, `GaussianBNClassifier`, `BayesianLogisticRecognizer`, `_hill_climb`, and the shipped weights `shuffle_logistic.json.gz` / `shuffle_bn.json.gz` | the Bayes-net and Bayesian-logistic recognizers | `tcren assess` |
| the `score` descriptor family and the `score` invariance class | catalogue entries for fitted quantities | nothing — **the catalogue is descriptors only** |

Renamed in the same release, with no alias: `F_*` → `Phi_*`, `dF_*` → `dPhi_*`, `F_TERMS` →
`PHI_TERMS`, `f_score` → `phi_score`. The `d` in `dPhi` is the **reference difference**
ΔΦ = Φ(sequence) − Φ(reference), never a derivative; the only `dd` quantity in the package is
`ddG`, the change in binding free energy on mutation.

## Already removed in 2.27.0 — the `free` qualifier

`S_free` only ever meant *not `P_native`*, and `P_native` was discarded in 2.26.0, so the qualifier
distinguished the score from nothing. **There is no alias and no deprecation shim**: a caller on an
old name, or a table carrying the old column, fails loudly rather than reading a renamed quantity by
accident. The three blocks keep their symbols — `Q`, `T` and `Π` are unchanged, as is the
construction.

| removed name | what it was | what to use instead |
|---|---|---|
| `reliability.s_free` | the three-block composite | `reliability.s_score` |
| the `tcren recognize` output column `S_free` | the same quantity, emitted | the column `S` |
| the frozen calibration keys `<set>\|S_nat` | link names inside `data/reliability_moments.json` | `<set>\|S`; the coefficients are untouched and only the keys were renamed, in place |
| `correct_confidence`'s returned coefficient `b_s_free`, and `af_band`'s `s_free_roc_in_band` | returned field names | `b_s` and `s_roc_in_band`; `correct_confidence` itself went in 2.28.0 |

Careful with the letter `S`: in 2.12.0 it named `cohort.q_coupled(Q, ΔΦ)`, a different quantity from
the composite it names now. [BENCHMARKS.md](BENCHMARKS.md) records which measurement used which.

## Already removed in 2.28.0 — every out-of-fold-fitted read-out

The author's ruling: leave-one-epitope-out fitting was the mechanism behind the defects in the
discarded `P_native`, the implementations were not trusted, and the analysis would be rewritten
rather than patched. After this release nothing tcren returned was fitted against a binding label.
The v3.0.0 tier-2 read-outs are fitted, against a named hold-out that ships in the wheel and refits
bit for bit, which is the property this layer lacked.

| removed | what it was | what to use instead |
|---|---|---|
| `reliability.p_binder`, `reliability.available_links`, and the `calibration` section of `data/reliability_moments.json` | Platt links fitted out of fold and shipped as fold means | `tcren.score.binder_score` |
| `reliability.correct_confidence`, `reliability.available_corrections`, `CORRECTION_VALIDATED_ON`, and the `corrections` section of the same file | the confidence correction: four coefficients, fitted out of fold and frozen | `tcren.score.confidence_residual`, tier 1, which reads no binder label |
| the `tcren diagnose` command | existed only to run that correction | `tcren assess` |
| the `p_binder` column from `tcren recognize` and `tcren assess`, with `assess`'s `--link` / `--list-links` | — | `assess` keeps `--list-bands` |

Untouched by that release: `Q`, `T`, `Π` and `S`, plus `inversion_flag`, `screening_yield` and
`af_band`. `moments()` asserts in its own test that `calibration` and `corrections` are absent, so
neither can return quietly. The removed read-outs and the numbers they produced are recorded in the
manuscript repository's `LEGACY.md`.
