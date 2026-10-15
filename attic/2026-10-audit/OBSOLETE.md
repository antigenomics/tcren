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

## Measured and not adopted, 3.9.0 — a record, not a deletion schedule

Everything below is **still shipped and still works**. This section is the only one here that is
not about removal: it records geometry estimators that were measured against the crystals and did
not earn a default, so nobody re-runs them. Four of the six exist only in the benchmark repository
and have no library code at all. Written as a table because `tests/unit/test_provenance.py` asserts
that every `- \`name\`` bullet in this file matches a `not_in_tcren2` marking, and `not_in_tcren2`
is scoped to the TCRen2 **potential** derivation — applying it to a geometry would be a category
error.

| what | where | what was measured |
|---|---|---|
| mean-surface read-back (`geom.BUNDLED_PRIOR`, `predicted_sites`) | library, shipped | superseded by 3.7.0's ensemble: averaging surfaces places a contacting pair at 10.90 Å where it measures 7.94, expected contacts 0.47x observed |
| an adjoint graphon registration (`graphon_weights(convention="interior")`, `read_back(convention="endpoint")`) | library, shipped, opt-in | buys 0.19 Å on CDR2 and costs 0.49 Å on CDR3; overall +3.26/+3.28 Å against the shipped pairing's +2.91, over 7,732 contacts in 370 leave-one-out crystals. The shipped mismatch aligns loop termini across lengths, which is the operation an average performs |
| Boltzmann annealing over drawn geometries (`annealed_contact_map(beta > 0)`) | library, shipped, defaults to 0 | js 0.4342 → 0.3538 and E[contacts] 22.0 → 47.1 against 22.6 observed, from β = 0 to β = 1. Computed on a **union** site set, so the site-count confound is ruled out; `log Z` is a free energy and rises with occupancy |
| quenched energy average; Boltzmann-λ over the fixed crystal members; near-end quantile surfaces; Platt sharpening | benchmark repo | each worse on 366 held-out crystals; Platt fits γ = 1.009, i.e. the map was already calibrated |
| `epibal` / `countbal` / germline-V-family / species+length member selection | benchmark repo | germline reads 0.437 against count-centrality's 0.438 at 79.0 % availability |

**Not on this list, and not to be confused with it:** `geometry.reach_max`, `reachability_floor`,
`span_saturation`, `LoopInternalCoords` and `cdr3_internal_coords`. `reach_max` is in production —
it is the `feasible` column of `tcren potts predict` — and the rest are the inverse direction a
forward loop model still needs. `geom_tcr_peptide_ens` remains the shipped template-free fallback.

## Scheduled for removal — the Ca-geometry ladder and the Potts scoring surface (2026-09-08)

**A table, not bullets**, for the reason the section above gives: `not_in_tcren2` is scoped to the
TCRen2 potential derivation and `tests/unit/test_provenance.py` matches every `- \`name\`` bullet
here against a marking, so listing a geometry as a bullet would fail the test and misuse the
decorator. Everything below still works and still ships; this records the measurement that says it
should not, and what keeping it costs.

The decision (author, 2026-09-08): remove the Ca-geometry material from the manuscript and mark both
subsystems here. What replaces the Potts scoring claim is a narrower one it can carry -- MD contact
frequency -- and the ranking claim it cannot.

| what | where | why it is marked |
|---|---|---|
| the Ca geometry prior | `potts/geom.py`, `geom_tcr_peptide{,_ens}` , `tcren potts geom` / `predict` / `agreement` | reachable only from the CLI; imported by nothing in `descriptors`, `score`, `recognition` or `pipeline`, so it feeds no descriptor, no channel and no benchmark read-out. Its ceiling is measured: MAE **4.65** contacts and **41.3 %** within +/-3 from everything a sequence knows, against the standing ensemble's 4.82 unaided |
| the coordinate chain model | `potts/chain.py` (600 lines), `geom_tcr_peptide_chain`, `tcren potts chain` | not exported from `potts/__init__.py` at all. Held out it reads js 0.4342 against the crystal's own Ca at 0.7052 |
| the cross product | `potts/cross.py` (459 lines, 36.3 M rows), `tcren potts cross` | consumed by nothing in the score path; feeds only `chain.distance_agreement`. The cross product of two marginals carries no receptor:peptide dependence by construction, which that module's own docstring says |
| the Potts read-outs no descriptor reaches | `score.bound_unbound`, `count_profile`, `sample_maps`, `connected_correlations`; `potts/ensemble.py`, `em.py`, `selfcon.py` | none is in `recognition.DESCRIPTORS`; `reproduce.sh` runs `features -i placement,interface,topology,energetics`, so no published number reads any of them |
| the opt-in coupling classes | `g_orient` (3.9.0), `far` / `X` (3.8.0) | fitted and measured (+787.0 and +81.8 nats) and `[]` / `False` in every shipped model |
| the two unread models | `data/potts_tcr_mhc.json`, `data/potts_pep_mhc.json` | shipped and `-m`-addressable; no descriptor and no read-out reads either |

**What stays, and why.** `potts/model.py`, `kernel.py`, `sites.py`, `fit.py`, `sample.py`,
`score.score_sites` and `contact_probabilities`. The main text keeps the contact-probability
weighting -- worth **+0.0105 [+0.0040, +0.0185]** median per-unit ROC-AUC on 218 cognate-random
units -- and the independent marginals supply it. `available_pairs` is the site enumeration
everything else is defined on.

**The measurement behind the demotion**, all on the sets already deposited:

- The couplings buy pseudo-likelihood, not ranking: **+505.7 nats for 18 parameters**, and the
  paired difference against independent marginals is **+0.0000 [-0.0015, +0.0000]** over 218 units
  (67 better, 107 worse). Independent marginals read 0.9710 / 0.9745 against the coupled
  0.9715 / 0.9708.
- The Pi block costs on the functionally validated screen: `S = Q + T + Pi` reads macro ROC-AUC
  0.7886 against `Q + T` alone at **0.8128**, i.e. **-0.0243 macro ROC** over 590 models. On the
  VDJdb panel it is +0.0353 over 1,027. Opposite signs, net a wash.
- Its one usable descriptor duplicates one already present: `potts:n_contacts` reads 0.7206 on
  TCRvdb, the best of the family against `neg_energy`'s 0.5717, and correlates with
  `interface:n_contacts_tp` at **rho = +0.999**.
- Two of five carry theory-contradicting signs on TCRvdb: `log_lik` 0.4450 and `psi` 0.4747, where
  theory says both exceed 0.5.
- The benchmark's `potts_model_choice.csv` prefers **beta x TCRen2 (114 parameters)** over the free-J model
  (513) at **delta-BIC 3,982**, and `STATUS.md` describes the shipped default as pinned when it is
  free. Settle that before anything cites a parameter count.

**What the Potts model is repurposed to.** Predicting MD contact frequency, where it transfers with
no parameter adjusted. Measured on identical cells (`bench/eda/out/potts_md_cells.csv`, 682
TCR:peptide cells over 4 complexes), it **orders** better than a hard contact call and **levels**
worse: Pearson 0.7446 against the crystal's own indicator at 0.7760, Spearman 0.7199 against 0.6440.
Report both baselines beside it; a claim that omits them will not survive review.

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
