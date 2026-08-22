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
- `tcren.potential.dfire.apply_corrections` --- As corrections.
- `tcren.potential.dfire.corrections` --- The distance and rotation corrections are measured against TCRen2, not applied to the shipped matrix. Which interface they are estimated on decides their sign: transferred from peptide:MHC they improve CPL, pooled over all interfaces they harm it.
- `tcren.potential.dfire.geometry_set` --- As pair_geometry.
- `tcren.potential.dfire.pair_geometry` --- DFIRE reference states are an independent line of work, evaluated against TCRen2 rather than folded into it.
- `tcren.potential.dfire.radial_potential` --- As pair_geometry.
- `tcren.potential.dfire.select_scope` --- As corrections.
- `tcren.potential.model.dfire2` --- A physics-reference baseline TCRen2 is compared against, not a component of it.
- `tcren.potential.model.tcren` --- The 2022 matrix, kept as the historical default of this loader. TCRen2 is tcren.potential.tcren2(); the two correlate at r = 0.875 with max |d| 0.846 and are not interchangeable.
- `tcren.potential.model.tcren2_dfire` --- A DFIRE-corrected variant under evaluation, not the shipped TCRen2.
- `tcren.potential.redundancy.cluster_weights` --- Down-weights by sequence-distance clusters rather than exact identity. Needs a threshold, and conflates the epitope and receptor axes that --balance separates.
- `tcren.potential.redundancy.epitope_weights` --- TCRen2 balances the epitope AND receptor axes; this single-axis alias is what the manuscript's earlier matrix used. Receptor redundancy is the larger of the two on Native2026 (223 structures share a receptor against 212 an epitope).
- `tcren.potential.redundancy.nonredundant_ids` --- Excludes cluster members outright. TCRen2 down-weights instead, which keeps every structure's data.
- `tcren.potential.smoothing.blosum_background` --- Substitution-matrix pseudocounts are under evaluation against TCRen2, not part of it. See docs/potentials.rst for the measurements.
- `tcren.potential.smoothing.blosum_conditional` --- As blosum_background.
- `tcren.potential.smoothing.impute_thin_cells` --- As blosum_background.
- `tcren.potential.smoothing.smooth_counts` --- As blosum_background.

## What deletion costs

| Group | Also delete | Kept by |
|---|---|---|
| Redundancy alternatives | `potential/redundancy.py` less `balanced_weights`, the `--nonred`/`--redundancy-t` CLI options | nothing |
| Leave-one-out | `derive_tcren_loo`, the `--loo` CLI option | nothing |
| DFIRE | `potential/dfire.py`, the `derive-dfire` command, the `dfire2` and `tcren2_dfire` entries in `data/potentials.json` and their CSVs | the manuscript's reference-state comparison, which is measured and written up |
| Substitution pseudocounts | `potential/smoothing.py`, the `--smooth-beta` CLI option | under evaluation — see below |
| The 2022 matrix | `potential.tcren()` and `data/TCRen_potential.csv` | published results predating TCRen2 |

`potential/smoothing.py` is the one entry not to delete wholesale. Its two schemes measured
differently on the best-powered endpoint available (TCRvdb receptor ranking against DESeq2 log fold
change, clean cohort YLQPRTFLL, n = 423 poses, paired bootstrap over 5,000 resamples):

| | Delta Spearman vs the 6 A matrix | P(Delta > 0) |
|---|---|---|
| `smooth_counts`, beta = 20 | -0.009 [-0.049, +0.032] | 0.35 |
| `impute_thin_cells`, min_count = 10 | **+0.063 [+0.022, +0.109]** | **0.997** |

So **`smooth_counts` goes** with the rest of this list, and `impute_thin_cells` is a candidate for
the recipe rather than for deletion. It carries `blosum_conditional` and `blosum_background` with
it, which is why those stay listed but not condemned.
