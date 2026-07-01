# tcren — datasets & sources

All structure sets live on HF **[isalgo/tcren_structures](https://huggingface.co/datasets/isalgo/tcren_structures)**
(rule: only `.gz`/`.tar.gz` LFS structures + `.txt`/`.md` descriptions). The manuscript repo
(`2026-tcren2`) mirrors this map in `scripts/bootstrap_data.py` + `DATASETS.md`.

## Structure sets on HF

| set | use |
|---|---|
| `Native2026`, `Canonical2026` | non-redundant TCR:pMHC structures — potential derivation, ergodicity, orientation refs |
| `Native2022`, `PolyV2022` | 2022-paper structure sets (reproduction/oracle) |
| `tcrvdb` | 618 TCRmodel2 structures — TCR-ranking / specificity benchmark |
| `cpl` | peptide-swap best/worst — peptide-ranking benchmark |
| `as_case` | ankylosing-spondylitis B\*27:05/:02 — MHC-ranking benchmark |
| `Bobisse`, `Bigot` | neoantigen cohorts (see below) |

## Legacy 2022 benchmarks (reproduced with the tcren pipeline)

`notebooks/natcompsci2022/` reproduces Karnaukhov et al. (Nat Comput Sci 2022) with **only** the
`tcren` Python pipeline (no `mir.jar`). Cohorts under `data_legacy/{Bobisse,Bigot,Birnbaum,vdjdb}`;
recomputed results in `results_new/`:

| notebook | task | output | tcren2 result |
|---|---|---|---|
| `02_benchmark_cognate_unrelated` | specificity (crystal) | `benchmark_cognate_ranks.csv` | AUC 0.88 (top-5%), n=218 |
| `03_benchmark_yeast_display` | peptide (Birnbaum) | `benchmark_birnbaum.csv` | AUC 0.89 |
| `04_benchmark_neoepitopes` | neoantigen (Bobisse/Bigot) | `benchmark_{bobisse,bigot}.csv` | Bobisse #1; Bigot median 13/43 |
| `05_benchmark_shuffle_structures` | shuffle control | `benchmark_shuffle_auc.csv` | TCRen 0.73 (MJ 0.52) |
| `06_benchmark_models` | generated vs crystal | `benchmark_models_ranks.csv` | — |

## New benchmarks (manuscript repo)

CPL / TCRvdb / AS scoring + the unified three-task table live in the manuscript
(`2026-tcren2/scripts/{benchmark_suite.md,tcren_binder_score.py,tcrvdb_physics.py}`). ATLAS ΔΔG is
pending upstream; MD trajectories (`md100ns_*.xtc`, ~57 GiB) pending upload.
