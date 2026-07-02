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
| TCR3D ground truth (60) | V-gene / CDR3 / class | **0.97 / 0.90 / 0.97** | `test_native_concordance_sweep` |
| TCR3D epitope | concordance | 0.72 (CIF-content-bounded, see notes) | — |
| Canonical alignment | self / 1bd2→1ao7 groove RMSD | **0.000 / 0.44 Å** | `test_native_uses` |
| αβ/γδ from C-gene | 1ao7 / 1hxm | **ab (TRBC2) / gd (TRDC+TRGC1)** | `test_cgene` |
| Re-derived TCRen (analysis) | max\|Δ\| vs published | **< 1e-9** | `test_analysis` |
| v2 configurable potentials (default) | per-interface scores vs built-in families | **byte-identical** | `test_default_equals_explicit_equal_mapping` |
| v2 TCR regions (`tcr_regions="all"`) | region filter membership/ordering | matches definition | `test_real_asset_region_ordering_and_fr_membership` |
| v2 oracle facade `scores` | vs `pipeline.run` for same args | **byte-identical** | `test_scores_reproduce_run_byte_exact` |

Notes: J-gene and class-II MHC allele names differ between pipelines by design (arda locus
follows the J segment — TCR3D's 1bd2 `TRDJ1` is a mislabel; class-II TCR3D uses serotypes).
Epitope < 1.0 is driven by domain-split/multi-copy TCR3D CIFs lacking a separable peptide
chain plus ±1 unresolved terminal residues — not a tcren error.

## Performance (Apple M3, base anaconda Python 3.12)

| Operation | Scale | Time |
|-----------|-------|------|
| Contact computation | 1 structure | ~40 ms |
| Full contact sweep | 312 structures | ~13 s |
| arda annotation | 1 TCR chain | ~1 s |
| MHC mapping (mmseqs `easy_search`) | 1 structure | ~7 s (per-call index build — TODO prebuild) |
| Fast test suite (`-m "not slow"`) | 74 tests | **~4 s** |
| Slow test suite (`-m slow`) | 25 tests | **~22 min** (arda/mmseqs per structure) |
| Annotation concordance sweep | 312 structures | ~20 min |
| Analysis notebook | full `contact_maps_PDB.csv` | < 30 s (no arda) |

## Example / analysis / benchmark tasks

| Artifact | What it shows |
|----------|---------------|
| `example/` | end-to-end scoring, reproduces `candidate_epitopes_TCRen.csv` |
| `notebooks/complementarity_map_2d.ipynb` | 2D interface map (SVG) + contact tables + polars summaries |
| `notebooks/pocket_cdr_3d.ipynb` | 3D groove + peptide + CDR overlay (py3Dmol) + matplotlib fallback |
| `notebooks/tcren_analysis.ipynb` | potential heatmaps (TCRen/MJ/Keskin), contact distributions per region & peptide/CDR3 position-vs-length |
| `tcren native derive-potential` | re-derive TCRen from TCR3D native structures |

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

### Cross-peptide docking accuracy (the honest test)

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

## Binder identification (`tcren.binder`)

Ranking candidate TCRs against a fixed pMHC on generated (AlphaFold/TCRmodel2) structures. The raw
TCR:peptide contact energy is at chance there (ROC-AUC ≈ 0.44, the forced-pose problem: the generator
seats every TCR in a plausible pose). The shipped 5-feature model — AF-orthogonal interface geometry
(interface size, dual-chain balance, H-bonds, buried ΔSASA; native `tcren._geom` C kernel) plus the
CDR1/2−CDR3α TCRen potential term — recovers it:

| model | denoised ROC-AUC | note |
|-------|------------------|------|
| **tcren.binder (5-feature)** | **0.928** | native, no external tool |
| AlphaFold/TCRmodel2 ipTM (confidence) | 0.872 | the baseline to beat |
| raw TCR:peptide TCRen energy | ≈ 0.44 | forced pose (below chance) |

TCRvdb (2 epitopes, HLA-A\*02:01; sequence-cluster-denoised labels). The model is ~ipTM-independent and
uses no generator-reported metric. Caveat: coefficients frozen from a 2-epitope training set;
cross-allele/epitope generalization untested (re-fit via `scripts/binder_validate.py`).
