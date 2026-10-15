# Publication dependency audit — 2026-10-14

Static inspection of all 35 `results/make_*.py` files, their named inputs, the current benchmark driver and selected owning writers. No scientific regeneration or figure/PDF rebuild. Repository-relative paths below use **B** for benchmark, **M** for manuscript, **L** for library. This inventory distinguishes a scheduled producer from a verified complete fresh-run dependency chain; it does not certify the latter.

## Producer families

| Publication reshapers (`make_…`) | Inputs / owning producers | Current driver status |
|---|---|---|
| `potential_dat`, `tcren_dat` | `fig2_*`, potential tables and ablations; `bench/figures/fig2/make_fig2_tables.sh`, `panel_stats.py`, `dphi_polyala.py`, `potential_figure.py`, `tcren_epitope_weighted.py`, `sidechain_reach.py` | Main producers scheduled in analysis; retired predicted geometry still has stale consumers below. |
| `cpl_dat`, `cplfreq_dat`, `md_dat` | `cpl_*`, `cplfreq_*`, `md_boltz_*`; CPL producers, `cpl_figure.py`, `cpl_predict.py`, `cpl_contact_freq_model.py`, `cpl_md_boltzmann.py` | Scheduled; input order needs complete-run verification, not a reshaper success alone. |
| `cplfusion_table`, `legacy_dat`, `legacy_table` | Per-unit legacy ranking + `pca/cpl_fusion*`; legacy harness, `legacy_peptide_ranking.py`, `cpl_fusion.py` | Scheduled. Legacy ranking was moved after harmonization datasets in this audit. |
| `tiers_dat` | Legacy ranking and `twostep_*` tables; `twostep_link_curves.py`, `twostep_link_scatter.py`, `twostep_product_afull.py`, `twostep_logsig_product.py`, `twostep_zsum_baseline.py` | Legacy input scheduled; no driver invocation of the two-step producer family found. Preserve its recorded results. |
| `tcrvdb_dat`, `native_dat`, `inversion_dat`, `tmpl_table` | `tcrvdb_*`, `native_bn_*`; `tcrvdb_figure.py`, `receptor_panels.py` | Scheduled; partial-support and pooled-macro publication issues remain in evidence audit. |
| `tmpl_dat` | `binder_*`, `tmpl_confidence_*`, `vdjdb_composition`, PCA template merge | Binder/model/template producers scheduled. |
| `erg_dat`, `rates_dat` | `erg_*`, `atlas_rate*`, kinetics uncertainty; `ergodic_inputs.py`, `ergodic_bridge.py`, `affinity_fits.py`, rebuttal kinetics producer | Energy/rebuttal paths scheduled. Rate-table source files require explicit raw-versus-derived provenance; no current driver edge regenerating `atlas_rates.csv` / `atlas_kd_split.csv` was identified. Do not infer absence of raw inputs. |
| `potts_dat` | `potts_md_contacts/cells/positions/saturation/baselines`; `potts_md_contacts.py`, `potts_md_baselines.py` | Scheduled. Baseline moved after cell-table producer. |
| `score_dat`, `resid_dat`, `qpca_dat` | PCA benchmark/interval/screen/kinetics/invariants, conditional residual, holdout metadata | PCA core scheduled; frozen-object metadata export gap remains below. |
| `arms_dat`, `rescue_dat`, `dispersion_dat` | `curves`, `curve_summary`, `paired_within_cohort`, `cohort_dispersion`, `template_ambiguity`; same-named `bench/pca/*.py` | Existing writers not invoked by current driver. |
| `correction_dat`, `rescale_dat` | `confidence_glm/corrected/rescale/operating/gain_by_stratum`; `bench/pca/confidence_*.py` | Existing writers not invoked by current driver; active publication consumers. |
| `decision_dat` | `pca/decision_{tree,leaves,edges,gate,holdout,by_problem}`; `bench/scripts/ctree_decision.R` | Not the `bench/eda/decisions.py` that the driver runs; tree publication input remains unwired. No new tree fit authorized. |
| `footprint_dat`, `topo_dat`, `selfcon_dat` | Contact subset/entropy, topography DATs, `selfcon_*`; `contact_subset.py`, `interface_topography.py`, `selfcon_compare.py` | Existing owning writers absent from driver invocations. |
| `diag_dat`, `matrix_dat` | AF bands, score-task matrix, legacy confidence correction/calibration files | AF-band/matrix writers scheduled; retired correction/calibration dependencies require explicit retention/retirement decisions. |
| `raqa_dat`, `hpv_dat`, `datasets_table`, `render_dat` | Rebuttal diagnostics, HPV results, dataset catalogue, render manifests | Rebuttal/HPV scheduled. Dataset catalogue/manifests are inputs; distinguish metadata from derived output before adding a producer. |

## Proven gaps and misleading gates

1. **Retired geometry still required by live reshaper.** M `results/make_potential_dat.py:281–295` loads `fig2_pred_<pdb>.tsv` and `fig2_pred_stats.csv`; `preds` is no longer consumed, and `ps_row` (`:338`) / `ps` (`:391`) are assigned but unused after the third plot row was removed. Those reads should be removed, preserving the underlying recorded data. B `bench/figures/fig2/panel_stats.py:193–240` retains a function computing predicted-geometry statistics, which requires the same retired TSVs; verify its invocation before removing it. A removed CLI dispatch alone does not close an artifact dependency.

2. **Fitted-object metadata has no safe evaluation-stage writer.** M `results/make_score_dat.py` requires `holdout_meta.csv`; B `bench/pca/freeze.py:79` is its identified writer, but that script also overwrites library model/manifest files and was correctly removed from evaluation. Export the metadata from the shipped frozen object without refitting/installing it. Do not restore the side effect just to make a reshaper pass.

3. **Confidence and display producers are unwired.** Verified writers include B `bench/pca/confidence_correction.py:166–190`, `confidence_rescale.py:140`, `confidence_operating.py:169`, `confidence_gain.py:130`, `curves.py:94`, `cohort_dispersion.py:67`, `template_ambiguity.py:91`, `paired_within_cohort.py:88`. Their output files are directly opened by the corresponding live reshapers, but current `recompute.sh` contains no invocation of those producers. Rewriting a macro file from their old CSV is not fresh reproduction.

4. **Old release pin remains.** B `requirements.txt:6` is `tcren==3.0.1`; the fig2 shell path documents using `tcren score --delta --reference structural`, introduced later. Pin selection requires agreement with the recorded results and the supported API; do not silently install a newer development version and claim the old release reproduced them.

5. **Reshaper import bans cannot establish “computes nothing.”** M `.github/check_structure.sh` checks forbidden imports, not mathematical operations. `results/make_native_dat.py:375–383` forms an error bar for a difference by combining percentile-CI half-widths in quadrature. That is a statistical approximation inside the publication layer; a lack of scipy imports does not validate it. Preserve existing figures until the author rules on moving/correcting that calculation.

6. **Normal `make gate` rebuilds figures.** M `results/figures/Makefile` declares `gate: all`. Therefore an ordinary gate invocation is not read-only and cannot be used under a no-regeneration constraint without suppressing that prerequisite.

## Checks performed

Inspected all top-level gate recipes and the checker scripts; no scientific computations are invoked by the check-only path. Ran:

```
make MAKE='make -o all' gate
```

The child make marks `all` old, so it checks existing figure PDFs without rebuilding them. Result: **exit 0 in 3.78 s**. Structure, vocabulary and cross-reference gates passed; 19 main labels, 32 note sections and 46 supplementary floats checked. Figure style: 55 scripts, 0 violations. Figure bounding box: 29 figures, 0 violations. Full log: `/tmp/tcren-manuscript-gate.log`.

These checks establish structural consistency of current sources/assets, not freshness, scientific correctness, label provenance, matching library version, unclipped rendered labels, or valid effect-size comparisons. No paper/figure compilation, bibliography verification or visual rendering review was performed in this task.

## Minimal remaining publication sequence

Finish orphan-read removal and safe shipped-model metadata export; reconcile the exact release and raw-input manifest; connect author-approved active producer families to the driver while retaining historical records; obtain rulings on the metric/interpretation inconsistencies in the evidence audit. Only then authorize one full fresh run, compare matched-unit results to standing records, and render/inspect the submission PDFs. Do not expand this gate into new scientific controls or rerun finished analyses during cleanup.

## Prepared follow-up closures

- `/tmp/tcren-make_potential_dat.py` removes the unused predicted-map loads and unused predicted-statistic assignments; retained crystal-map plotting and every live output expression are unchanged.
- `/tmp/tcren-panel_stats.py` removes `main()`'s active `pred_stats(a.out)` call and the orphan predicted-statistics helper block. Current panel-B computations are unchanged. Both originals must be archived before replacement. Original and prepared SHA256 values are in `/tmp/tcren-orphan-reads-plan.json`; AST parsing passed for both prepared files.
- `/tmp/tcren-export_holdout_metadata.py` is a proposed `bench/scripts/export_holdout_metadata.py`: it uses `holdout_model()` and `holdout_manifest()` to export the shipped object's metadata with stdlib CSV, never calls fitting or writes library files. Wire it before publication reshaping. A read-only model load/export to `/tmp/tcren-holdout-meta.csv` succeeded; every column shared with the existing `holdout_meta.csv` matched exactly. The export additionally records manifest row count 8,292, which was absent from that current table. No structural features or scores were computed.
