# tcren — achieved accuracy & performance

Validation of the Python pipeline against the legacy R/Java oracle and external references.
Reproduce with `pytest` (fast) and `RUN_BENCHMARK=1 pytest` (full-dataset sweeps).

## Accuracy (vs oracle / reference)

| Task | Metric | Result | Test |
|------|--------|--------|------|
| Classic TCRen potential | max\|Δ\| vs `TCRen_potential.csv` | **≤ 1e-9** (exact) | `test_potential_regression` |
| `am` (gap) potential | max\|Δ\| vs `tcren_am/tcren.txt` | **2.8e-17** (from matched data) | `test_potential_regression` |
| TCR↔peptide contacts | exact set vs `contact_maps_PDB.csv` | **312 / 312 structures** | `test_contacts_regression` (`RUN_BENCHMARK`) |
| Candidate scoring | max\|Δ\| vs `run_TCRen.R` on `example/` | **4.4e-15** (exact) | `test_score_regression` |
| TCR annotation sweep (mir set) | contacts reproduced / full-exact | **0 missing**, 278 full-exact, 31 region-label-only / 312 | `test_annotation_concordance_sweep` |
| MHC class + locus | sample concordance | **30 / 30**; 1ao7/5m01/4ozg exact | `test_mhc_regression` |
| MHC groove topology | TCR-on-helices / peptide-on-floor | satisfied (class I + II) | `test_mhc_groove` |
| TCR3D ground truth (60) | V-gene / CDR3 / class | **0.97 / 0.90 / 0.97** | — (annotation reproduction covered by `test_annotation_concordance_sweep`) |
| TCR3D epitope | concordance | 0.72 (CIF-content-bounded, see notes) | — |
| Canonical alignment | self / 1bd2→1ao7 groove RMSD | **0.000 / 0.44 Å** | `test_orient` |
| αβ/γδ from C-gene | 1ao7 / 1hxm | **ab (TRBC2) / gd (TRDC+TRGC1)** | `test_cgene` |
| Re-derived TCRen (analysis) | max\|Δ\| vs published | **< 1e-9** | `test_analysis` |
| v2 configurable potentials (default) | per-interface scores vs built-in families | **byte-identical** | `test_default_equals_explicit_equal_mapping` |
| v2 TCR regions (`tcr_regions="all"`) | region filter membership/ordering | matches definition | `test_real_asset_region_ordering_and_fr_membership` |
| v2 oracle facade `scores` | vs `pipeline.run` for same args | **byte-identical** | `test_scores_reproduce_run_byte_exact` |

Notes: J-gene and class-II MHC allele names differ between pipelines by design (arda locus
follows the J segment — TCR3D's 1bd2 `TRDJ1` is a mislabel; class-II TCR3D uses serotypes).
Epitope < 1.0 is driven by domain-split/multi-copy TCR3D CIFs lacking a separable peptide
chain plus ±1 unresolved terminal residues — not a tcren error.

## Performance (Apple M3, uv-managed CPython 3.12)

| Operation | Scale | Time |
|-----------|-------|------|
| Contact computation | 1 structure | ~40 ms |
| Full contact sweep | 312 structures | ~13 s |
| arda annotation | 1 TCR chain | ~1 s |
| MHC mapping (mmseqs `easy_search`) | 1 structure | ~7 s (per-call index build — TODO prebuild) |
| Fast test suite (`-m "not slow"`) | ~291 tests | **~75 s** |
| Slow test suite (`-m slow`) | ~46 tests | **~22 min** (arda/mmseqs per structure) |
| Annotation concordance sweep | 312 structures | ~20 min |
| Analysis notebook | full `contact_maps_PDB.csv` | < 30 s (no arda) |

## Example / analysis / benchmark tasks

| Artifact | What it shows |
|----------|---------------|
| `example/` | end-to-end scoring, reproduces `candidate_epitopes_TCRen.csv` |
| `notebooks/complementarity_map_2d.ipynb` | 2D interface map (SVG) + contact tables + polars summaries |
| `notebooks/pocket_cdr_3d.ipynb` | 3D groove + peptide + CDR overlay (py3Dmol) + matplotlib fallback |
| `notebooks/tcren_analysis.ipynb` | potential heatmaps (TCRen/MJ/Keskin), contact distributions per region & peptide/CDR3 position-vs-length |
| `tcren derive-potential` | re-derive TCRen from TCR3D native structures |

## Peptide modelling — open-source fold engines (draft, milestone S6 QC)

License-free replacement for FlexPepDock / MODELLER (`tcren.refine.model_peptide`). The benchmark
(`scripts/fold_benchmark.py`) is a **recovery** test, deliberately *not* native-in/native-out: it
threads the native peptide, applies a **rigid displacement** (default σ ≈ 1 Å translation, 15°
rotation) as the shared start for every engine, re-models, and measures peptide RMSD to the native
crystal pose (MHC-groove superposition). FlexPepDock is the optional **oracle** (accuracy ceiling) via
`$ROSETTA_BIN`; the open engines should approach it. Numbers below are a smoke subset — interpret with
the per-engine caveats, not as headline accuracy.

Smoke run (**n = 6 Native2026 class-I complexes**, rigid displacement σ = 1 Å / 15°, Apple M3, env
`tcren-fold`, all recovered, 0 failed/skipped); full guarded sweep: `RUN_BENCHMARK=1 python
scripts/fold_benchmark.py`. Engines installed via conda (OpenMM/OpenStructure+ProMod3) and
`pyrosetta-installer` (FlexPepDock).

| Engine | Backend | backbone RMSD (med) | anchor-Cα RMSD (med) | time (med) | Note |
|--------|---------|---------------------|----------------------|------------|------|
| `dope` | `tcren._refine` rigid-body MC (DOPE) | 0.35 Å | 0.43 Å | ~0.18 s | best here† |
| `ccd` | `tcren._fold` CCD Cα closure | 1.02 Å | 0.93 Å‡ | ~0.003 s | fastest |
| `openmm` | OpenMM AMBER, receptor frozen | 1.47 Å | 1.49 Å | ~6 s | local min§ |
| `promod3` | ProMod3 rotamer repack | 1.47 Å¶ | — | ~2.4 s | repack only¶ |
| `flexpep` (oracle) | PyRosetta FlexPepDock | 0.00 Å (native→native) | — | ~215 s | ceiling, opt-in |

All engines run; the oracle is validated (`flexpep_refine` on native 1ao7 → 0.00 Å in ~3.5 min) and is
opt-in via `--oracle` because it costs minutes/structure. Interpretation caveats (this is a diagnostic
harness — engine ranking depends on the displacement model, and a **rigid** displacement structurally
favours the rigid-body refiner):
† `dope` is a rigid-body MC refiner, i.e. the matched inverse of a rigid displacement, so it recovers
best *for this perturbation*; it is a **local** refiner (restrained to its input), not a global search.
‡ `ccd` is driven to the **native anchor Cα** (the only targets without de-novo pocket prediction), so
its anchor RMSD is a closure residual (input-driven), **not** an accuracy claim; its accuracy metric is
`bb`. Its output is a Cα-trace model (distorted peptide bonds) that must be energy-refined.
§ `openmm` freezes the receptor and does local gradient minimisation, so it settles in the basin near
the displaced start rather than searching back to native — a physics *relaxation*, not a docking search.
¶ `promod3` here does rotamer side-chain reconstruction only; it does **not** move the backbone, so its
backbone RMSD ≈ the displacement itself (it is a repack step, to be paired with a backbone engine).

The true accuracy ceiling is the FlexPepDock oracle; the native C++ engines (`CPP_REWRITE.md`) are
validated against it. These smoke numbers characterise pipeline behaviour, not final accuracy.

### Cross-peptide docking accuracy, scored against a held-out native

`scripts/fold_crossdock_benchmark.py` measures the real question, not self-reconstruction: take pMHC
structure **A**, replace its peptide with a *different* peptide **P_B** that binds the same MHC allele
and whose native complex **B** is known, model P_B into A's groove, and measure RMSD to P_B's **true**
pose in B (MHC-groove superposition A→B). Pairs are same-allele, same-length (class-I 9-mers dominate).

**40 pairs** (148 structures indexed, 2798 candidate pairs), peptide backbone RMSD to native P_B:

| Method | median | mean | p75 | vs baseline (per-pair) |
|--------|--------|------|-----|------------------------|
| baseline (P_B threaded on A, no refine) | **0.98 Å** | 0.97 | 1.25 | — |
| ccd | 0.98 Å | 0.97 | 1.25 | identical (0/40 differ) |
| promod3 | 0.98 Å | 0.97 | 1.25 | identical (repack only) |
| openmm | 0.98 Å | 0.97 | 1.25 | ≤ 0.001 Å (backbone unmoved) |
| dope | 1.03 Å | 1.01 | 1.22 | worse (+0.04 mean, 20/40 drift) |
| flexpep (oracle) | *excluded* | | | **no-op** — see `oracle_flexpep.py` |

Per-pair range 0.22 Å (near-identical peptides) to 1.18 Å. **Finding:** same-allele backbone transfer is
a ~1.0 Å baseline (MHC-I 9-mer backbones are groove-conserved), and none of the runnable engines beat
it — they are refiners, not pose predictors. `dope` is slightly worse (rigid MC drifts toward its own
energy optimum). The FlexPepDock oracle is currently a no-op on these 5-chain complexes (FoldTree/jump
setup needed), so the one method that might beat the baseline is not yet measured. The open problem for
the C++ rewrite is de-novo pocket/pose prediction, not refinement speed.

## Receptor ranking on AI-generated structures — the score set (v3.0.0)

Ranking candidate TCRs against a fixed pMHC on generated (AlphaFold/TCRmodel2) structures. The raw
TCR:peptide contact energy is at chance there — the forced-pose problem: the generator seats every
TCR in a plausible low-energy pose, binder or not. `tcren assess` reads the interface instead, and
returns the score set: `peptide_score` (tier 0, nothing estimated), `pose_score` and
`confidence_residual` (tier 1, a covariance over hold-out binders, no negative and no label),
`binder_score` and the five `channel_*` columns (tier 2, class means and covariances from hold-out
binder labels). Every one is a projection of one frozen object — a transform plus a Gaussian per
class over the transformed descriptor coordinates — so every one is **defined for a single
structure**.

**How to read every number in this section.** None of it was computed in this repository, and none
should be: computing an AUC belongs in the benchmark repo `~/vcs/projects/2026-tcren2-code`. Each
value is quoted from CHANGELOG `[3.0.0]` or `docs/assess.rst`, carrying the claim it was attached
to there. The metric is **ROC-AUC computed within epitope cohort**; where a source names the
aggregate across cohorts, the line below names it too. The two panels are the ones the superseded
section below defines: the functionally validated receptor screen (TCRvdb — n = 618 structures,
309 binders / 309 non-binders, 2 epitope cohorts) and the VDJdb real-versus-mock panel (n = 1,089
structures, 523 real / 566 mock, 22 epitope cohorts, of which 6 are template-covered and 16
template-free).

### The fit-free tier leads, and composes with the fitted one

| panel / stratum | score | ROC-AUC, within epitope cohort |
|---|---|--:|
| functionally validated receptor screen (2 cohorts, n = 618) | `S`, fit-free | **0.818** |
| functionally validated receptor screen (2 cohorts, n = 618) | generator ipTM | 0.795 |
| VDJdb template-covered stratum (6 of 22 cohorts) | `S` + `binder_score` | **0.783** |
| VDJdb template-covered stratum (6 of 22 cohorts) | `S` alone | 0.665 |
| VDJdb template-covered stratum (6 of 22 cohorts) | `binder_score` alone | 0.682 |

`S` leads the screen on its own, and on the template-covered stratum it **composes** with
`binder_score` rather than being replaced by it — 0.783 together against 0.665 and 0.682 apart.
That is why `cohort.q_score` (`Q`), `reliability.t_score` (`T`) and `reliability.s_score` (`S`) are
documented as the fit-free predecessor tier and were not retired with the rest. (Source: CHANGELOG
`[3.0.0]`, "Not removed".)

### The channels, and where one beats the whole model

A channel is the same posterior marginalized to one descriptor family, which is a sub-block of the
covariance — exact, closed form, no re-fit. The five channels do not sum to `binder_score` and
should not: the whole model also reads the correlations *between* channels. Sometimes a channel is
the better instrument, because the whole model dilutes it.

| panel | channel | channel ROC-AUC, within cohort | full posterior |
|---|---|--:|--:|
| VDJdb template-free stratum (16 of 22 cohorts) | `channel_shape` | **0.637** | 0.615 |
| combinatorial peptide library | `channel_energetics` | **0.700** | 0.542 |

(Source: `docs/assess.rst`.)

### Two things measured on the way, both reading no binder label

- **The transform is variance-stabilising and never rank-based, and that is measured rather than
  argued.** Mapping each marginal onto a uniform CDF took the per-cohort median ROC-AUC from
  **0.630 to 0.543** on the six template-covered cohorts, and from **0.613 to 0.507** on the sixteen
  template-free ones. The signal *is* the marginal scale, and flattening the distribution deletes
  it.
- **Some directions of the binder manifold belong to the structure generator, not to the
  interface.** `reliability.artefact_directions` is a label-free test for which. The tightest band
  is broken **4.55×** by native crystals against **3.18×** by decoys, and it scores
  **0.504** ROC-AUC over the sixteen template-free cohorts — a coin — where the loose
  band reads **0.606**. Ledoit-Wolf shrinkage floors the shipped model's smallest direction at
  s.d. **0.0797**, above that band, so `pose_score` cannot read it. (Source: CHANGELOG `[3.0.0]`.)

### Reproducing the frozen model

`holdout_manifest()` returns the **8,292 hold-out structures** the fit used, with dataset, epitope,
binder label and ipTM. The manifest (329 kB) and the model arrays (315 kB) ship inside the wheel;
the 19 MB descriptor table deliberately does not, and `tcren fetch-data` + `tcren features` +
`tcren fit-holdout` regenerates the shipped arrays **bit for bit**. That is the difference from
`P_native` below, whose coefficients were frozen against a training set nobody could reconstruct.

## Superseded — receptor ranking with `P_native` (shipped v2.12.0, discarded v2.26.0)

**`P_native` and `cohort.p_native` were removed in v2.26.0**, with `P_NATIVE_CHANNELS` / `_POOL` /
`_ORIENT` / `_FEATURES` / `_BANNED` and the per-channel posteriors `G` / `T` / `E`. Discarded, not
deprecated, and the changelog gives the reasons in one place: it refitted a latent class on every
call, raised when a cohort had fewer rows than features, and its value depended on which rows the
fit was anchored on — none of which survives contact with a user holding one model — and its
coefficients were frozen against training sets that no longer exist, which made it the one part of
the package a reader could not reproduce.

**The measurements stay, because they are the record of what was measured.** Nothing below is
reproducible from the current package and no value in it is a current tcren result. Three names
used below have since been reused, and must not be read as today's quantities:

| name below | what it meant here | what the name means now |
|---|---|---|
| `S` | `cohort.q_coupled(Q, ΔΦ)`, already deprecated when these rows were measured | `reliability.s_score`, the three-block composite `Q/sd_Q + T/sd_T + (Π − μ)/sd_Π`, renamed from `S_free` in v2.27.0 |
| `T` | `P_native`'s topology **channel** — one posterior of the same latent-class fit, so it carried every defect above | `reliability.t_score`, a fit-free directional score against the Native2026 crystals |
| `binder_score` | `tcren.binder.binder_score`, a frozen 5-feature logistic, removed in v2.26.0 | `tcren.score.binder_score`, the tier-2 log-odds of the frozen score set |

Ranking candidate TCRs against a fixed pMHC on generated (AlphaFold/TCRmodel2) structures.
`P_native` read the interface as a latent class over three channels (geometry, footprint topology,
contact energetics), each a conditional-linear-Gaussian Bayes network, their log-odds added.

**How to read every number below.** All are **macro** averages over epitope cohorts — the mean of
the per-cohort value — because a pooled AUC on these panels reads epitope composition rather than
recognition. `P_native` and its channels were fitted **leave-one-epitope-out**: a cohort was scored
by a model anchored on the *other* cohorts' rows, so no scored row contributed its own label to the
fit that scored it. Numbers are quoted to three decimals, as the generator printed them. Source: the
benchmark repo `~/vcs/projects/2026-tcren2-code`, files
`bench/eda/out/native_bn_{endpoint,channel_auc,template,glm,glm_gain}_*.csv` and `results/ledger.md`
— the producer `bench/scripts/native_bn.py` no longer exists in that repository.

### TCRvdb — n = 618 structures, 309 binders / 309 non-binders, 2 epitope cohorts

TCRmodel2 models of a validated receptor panel on HLA-A\*02:01; GLCTLVAML n = 195, YLQPRTFLL
n = 423. **Raw labels** (`padj < 1e-5`, no cleaning).

| score | macro ROC-AUC | macro PR-AUC | macro precision @ 10% recall |
|---|--:|--:|--:|
| **`P_native`** | **0.832** | **0.849** | **0.955** |
| `T` — topology channel alone | 0.815 | 0.828 | 0.933 |
| `S` = `q_coupled(Q, ΔΦ)` — **deprecated**, kept as the harness check | 0.802 | 0.817 | 0.935 |
| `D2_pep24` — single topology descriptor | 0.787 | 0.804 | 0.925 |
| AlphaFold/TCRmodel2 ipTM | 0.795 | 0.783 | 0.912 |
| AlphaFold/TCRmodel2 ranking confidence | 0.795 | 0.783 | 0.912 |
| AlphaFold/TCRmodel2 pLDDT | 0.776 | 0.800 | 0.916 |
| `Q` — interface geometry, fit-free | 0.779 | 0.764 | 0.827 |
| `P_native`, flat network over the union of the same features | 0.750 | 0.758 | 0.892 |
| ΔΦ TCR:peptide (inverse `dF`) | 0.557 | 0.622 | 0.729 |

`S` reproducing macro ROC-AUC 0.802 / macro PR-AUC 0.817 on these 618 structures is the harness
check that nothing upstream had moved. `cohort.coupling` and `cohort.q_coupled` are still
importable and still byte-identical, so those two values reproduce; the current composite carrying
the name `S` is `reliability.s_score`, which is a different quantity (see the table above).

### VDJdb real-versus-mock — n = 1,089 structures, 523 real / 566 mock, 22 epitope cohorts

TCRmodel2 models of motif-supported VDJdb binders against mock complexes built by the identical
unrelaxed procedure.

| score | macro ROC-AUC | macro PR-AUC | macro precision @ 10% recall |
|---|--:|--:|--:|
| **`P_native`** | **0.718** | 0.685 | 0.812 |
| `P_native`, flat network over the union | 0.689 | **0.692** | **0.872** |
| `T` — topology channel alone | 0.648 | 0.640 | 0.774 |
| AlphaFold/TCRmodel2 pLDDT | 0.605 | 0.613 | 0.779 |
| AlphaFold/TCRmodel2 ipTM | 0.592 | 0.606 | 0.770 |
| `S` = `q_coupled(Q, ΔΦ)` — deprecated | 0.576 | 0.577 | 0.720 |
| `Q` — interface geometry, fit-free | 0.560 | 0.564 | 0.698 |
| `D2_pep24` — single topology descriptor | 0.528 | 0.523 | 0.646 |
| ΔΦ TCR:peptide (inverse `dF`) | 0.488 | 0.506 | 0.649 |

Both `P_native` rules are published: the factored (log-odds sum) rule lifts the hard cohorts, the
flat network holds the top of the ranking on the easy ones.

### What each channel says on its own

| channel | TCRvdb ROC / PR (n = 618) | VDJdb ROC / PR (n = 1,089) |
|---|--:|--:|
| geometry | 0.782 / 0.763 | 0.633 / 0.630 |
| topology | 0.815 / 0.828 | 0.648 / 0.640 |
| energetics | 0.426 / 0.502 | 0.674 / 0.652 |

Topology is the only channel above chance on both. The energetics channel sits on **opposite sides
of chance** on the two panels — 0.426 against 0.674 — which is the forced-pose inversion read
directly, and the reason a hand-picked weighting does not transfer while a fitted sign does.

### Template coverage is the variable that matters

Split the VDJdb panel by whether *some* receptor has already been co-crystallized with that peptide.
Template-covered: 6 cohorts, n = 297 structures, 142 real. Template-free: 16 cohorts, n = 792
structures, 381 real. Macro ROC-AUC within each split:

| score | template-covered | template-free | lost |
|---|--:|--:|--:|
| **`P_native`** | 0.721 | **0.716** | **0.005** |
| GLM(`P_native` + ipTM + pLDDT), in sample | 0.735 | 0.727 | 0.008 |
| `T` — topology channel | 0.756 | 0.608 | 0.148 |
| pLDDT | 0.675 | 0.579 | 0.096 |
| ipTM | 0.692 | 0.555 | 0.136 |
| `Q` — interface geometry | 0.729 | 0.497 | 0.232 |
| ΔΦ TCR:peptide | 0.622 | 0.438 | 0.185 |

Every score that reads the placement collapses when the template goes; `P_native` does not, and the
generator's own confidence does **not** fall to warn you that it should be distrusted.

### Composing with the generator's confidence

A plain logistic on `P_native`, ipTM and pLDDT, **fitted and read in sample** as a demonstration of
complementarity rather than a ranking claim. Δ is the joint model minus `P_native` alone, with a
paired percentile 95% CI from 2,000 resamples macro-averaged over cohorts.

| panel | score | macro ROC-AUC | macro PR-AUC | macro P @ 10% recall |
|---|---|--:|--:|--:|
| TCRvdb (n = 618) | ipTM | 0.795 | 0.783 | 0.912 |
| TCRvdb (n = 618) | pLDDT | 0.776 | 0.800 | 0.916 |
| TCRvdb (n = 618) | GLM(ipTM + pLDDT) | 0.796 | 0.787 | 0.917 |
| TCRvdb (n = 618) | `P_native` | 0.832 | 0.849 | 0.955 |
| TCRvdb (n = 618) | GLM(`P_native` + ipTM + pLDDT) | **0.840** | **0.861** | **1.000** |
| VDJdb (n = 1,089) | ipTM | 0.592 | 0.606 | 0.770 |
| VDJdb (n = 1,089) | pLDDT | 0.605 | 0.613 | 0.779 |
| VDJdb (n = 1,089) | GLM(ipTM + pLDDT) | 0.597 | 0.609 | 0.789 |
| VDJdb (n = 1,089) | `P_native` | 0.718 | 0.685 | 0.812 |
| VDJdb (n = 1,089) | GLM(`P_native` + ipTM + pLDDT) | **0.729** | **0.719** | **0.830** |

| panel | ΔROC-AUC vs `P_native` | 95% CI | ΔPR-AUC vs `P_native` | 95% CI |
|---|--:|---|--:|---|
| TCRvdb (n = 618) | +0.008 | [−0.009, +0.027] | +0.012 | [−0.012, +0.036] |
| VDJdb (n = 1,089) | +0.011 | [−0.012, +0.035] | **+0.034** | **[+0.005, +0.054]** |

On TCRvdb the generator's confidences add **nothing resolvable** above `P_native`: both intervals
contain zero, at P(Δ>0) = 0.814 for ROC and 0.818 for PR. On VDJdb only the PR gain clears zero
(P(Δ>0) = 0.994). They rank well on their own; on these cohorts they carry little the structure does
not already say.

### Why `P_native` carried no fitted coefficient

> **Fitting to a cohort does not transfer, so nothing shipped was fitted to one.** The frozen
> 5-feature `p_bind` (`tcren.binder`, then kept for v1 reproduction and reached only via
> `recognize --scores`; the model, `tcren.binder.binder_score`, `BINDER_MODEL`, the `tcren binder`
> command and the `--scores` flag were all removed in v2.26.0) reads macro ROC 0.796 / pooled ROC
> 0.810 / macro PR 0.804 on the same 618
> TCRvdb structures — competitive in sample, and it does not carry over. Train on VDJdb, test on
> TCRvdb: macro 0.466, **below chance**. Train on TCRvdb, test on VDJdb: macro 0.537. Against a
> within-cohort TCRvdb cross-validation of 0.811. Holding whole epitopes out of the VDJdb panel
> drops a random-fold pooled 0.808 to 0.471, because class balance there tracks epitope almost
> perfectly. `P_native` was fitted per cohort **without any binding label**, so there was no
> coefficient to carry and nothing to transfer. (Source: `results/ledger.md`, entries C25 and the
> receptor-ranking table.)

> **Pick the baseline deliberately.** Which AlphaFold confidence is strongest depends on the metric.
> On TCRvdb macro ROC, ipTM (0.795) and ranking confidence (0.795) lead pLDDT (0.776); on macro PR
> the order flips and pLDDT (0.800) leads both (0.783). Quote the margin against the best one for
> the metric being reported.

> **Label denoising is a separate algorithm and is not benchmarked here.** Filtering TCRvdb by
> TCRNET motif-cluster consistency raises every method's score, tcren's and AlphaFold's alike; a
> number computed that way measures the two algorithms jointly and is not a tcren result. All rows
> above use raw labels. For the record, the legacy `p_bind` coefficients were *fit* on such a
> subset; the within-TCRvdb 5-fold CV once offered as evidence for that (macro 0.776 denoised
> against 0.761 raw, 20/20 seeds) is **not a split** — TCRvdb is two epitopes, so random folds
> interpolate inside two clouds — and under the proper cross-dataset split the advantage vanishes,
> denoising changing transfer by ≤ 0.006 in either direction (`results/ledger.md`, C24).
