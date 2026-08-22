# OBSOLETE — scheduled for deletion

Everything here is marked `@not_in_tcren2` in the source: correct, tested, and taking no part in
producing the shipped matrix. The recipe is
`tcren derive-potential --structure-dir Native2026 --balance both`, at a 6 Å contact cutoff, and
nothing else. **Delete this list, and the code behind it, at the next major version.** Keeping four
redundancy schemes and three reference states around is how three mutually inconsistent potentials
came to be in circulation at once.

Regenerate with `python -m tcren._provenance`; `tests/unit/test_provenance.py` fails if this file
and the markings disagree.

- `tcren.potential.derive.derive_tcren_loo` --- Leave-one-out derivation, for testing how much any single structure moves the matrix. Diagnostic, not a production path.
- `tcren.potential.model.tcren` --- The 2022 matrix, kept as the historical default of this loader. TCRen2 is tcren.potential.tcren2(); the two correlate at r = 0.875 with max |d| 0.846 and are not interchangeable.
- `tcren.potential.redundancy.cluster_weights` --- Down-weights by sequence-distance clusters rather than exact identity. Needs a threshold, and conflates the epitope and receptor axes that --balance separates.
- `tcren.potential.redundancy.epitope_weights` --- TCRen2 balances the epitope AND receptor axes; this single-axis alias is what the manuscript's earlier matrix used. Receptor redundancy is the larger of the two on Native2026 (223 structures share a receptor against 212 an epitope).
- `tcren.potential.redundancy.nonredundant_ids` --- Excludes cluster members outright. TCRen2 down-weights instead, which keeps every structure's data.

## What deletion costs

| Group | Also delete | Kept by |
|---|---|---|
| Redundancy alternatives | `potential/redundancy.py` less `balanced_weights`, the `--nonred`/`--redundancy-t` CLI options | nothing |
| Leave-one-out | `derive_tcren_loo`, the `--loo` CLI option | nothing |
| The 2022 matrix | `potential.tcren()` and `data/TCRen_potential.csv` | published results predating TCRen2 |

The DFIRE reference states and the substitution-matrix pseudocounts that were on this list are
gone: all five estimation improvements tried against TCRen2 were measured and none improved every
endpoint, so they were deleted rather than left to rot. See
`src/tcren/data/TCRen2_potential.NOTES.md` before re-attempting any of them.
