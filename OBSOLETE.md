# OBSOLETE — scheduled for deletion

Everything here is marked `@not_in_tcren2` in the source: correct, tested, and taking no part in
producing the shipped matrix. The recipe is
`tcren derive-potential --structure-dir Native2026 --balance both`, at the 5 Å heavy-atom
contact cutoff `contact_table` uses (the subcommand has no `--cutoff`), and nothing else. **Delete this list, and the code behind it, at the next major version.** Keeping four
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
| `footprint.footprint_score` | the `fp_score` z-sum | `cohort.p_native(t, channels=("topology",))` |
| `cohort.q_iptm`, `q_f`, `q_f_iptm`, `f_invert_by_iptm`, `phi_bind`, `agreement` | hand-picked combination rules; `agreement` was the per-structure summand of `C*` | `cohort.p_native`, which fits each channel's sign instead of asserting it |
| `recognition.kit_score` | a z-sum of `p_bind` and ipTM | `cohort.p_native`; join the generator's confidence yourself if you want it |
| `scripts/fit_pose_reference.py`, `scripts/fit_joint_reference.py` | regenerated the two deleted reference CSVs | nothing |

**Deprecated, not removed:** `cohort.coupling` and `cohort.q_coupled` remain importable, tested and
byte-identical in behaviour, so every published `S` reproduces. They are superseded by `p_native`,
which fits each channel's sign instead of measuring it. `tcren.binder`'s fitted `p_bind` and
`recognize --scores` are likewise kept for v1 reproduction only.
