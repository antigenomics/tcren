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
