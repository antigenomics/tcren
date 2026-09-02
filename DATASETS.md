# tcren — datasets & sources

All structure sets live on HF **[isalgo/tcren_structures](https://huggingface.co/datasets/isalgo/tcren_structures)**
(rule: only `.gz`/`.tar.gz` LFS structures + `.txt`/`.md` descriptions). `tcren fetch-data` pulls
the reference sets into `paths.tcren_home()/data/`; `$TCREN_DATA_DIR` overrides that directory.

The **benchmark repo is `~/vcs/projects/2026-tcren2-code`** and it is the catalogue of record:
`bench/scripts/bootstrap_data.py` fetches the benchmark cohorts, `SOURCES.md` carries the
per-dataset provenance, and `data/datasets.csv` is the machine-readable table (one row per dataset:
task, structure and instance counts, positives/negatives, epitopes, alleles, MHC class, species,
provenance, aldan3 path, HF path).

## Structure sets on HF

| set | use |
|---|---|
| `Native2026`, `Canonical2026` | non-redundant TCR:pMHC structures — potential derivation, ergodicity, orientation refs |
| `Native2022`, `PolyV2022` | 2022-paper structure sets (reproduction/oracle) |
| `tcrvdb` | 618 TCRmodel2 structures — TCR-ranking / specificity benchmark |
| `cpl` | peptide-swap best/worst — peptide-ranking benchmark |
| `as_case` | native B\*27:05 cognate complexes — held-out ergodicity validation |
| `vdjdb_positives`, `vdjdb_negatives` | TCRmodel2 real-versus-mock models; the 1,089-structure receptor-ranking panel (523 real / 566 mock, 22 epitopes) is assembled from these, and `data/datasets.csv` records the consolidated tarball as `vdjdb_binder_benchmark/` |
| `garcia_b27` | HLA-B\*27:05 crystals with a measured EC50 substitution series |
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

## New benchmarks (benchmark repo)

Everything current lives in `~/vcs/projects/2026-tcren2-code`, not here:

| what | where |
|---|---|
| dataset catalogue, provenance | `SOURCES.md`, `data/datasets.csv` |
| every reported number, with its n | `results/ledger.md` |
| fetch the benchmark cohorts | `bench/scripts/bootstrap_data.py` |
| receptor ranking — the score set (`tcren assess`), the fit-free `Q` / `T` / `S`, the template split, the composition with AlphaFold confidence | `bench/scripts/tcrvdb_panel.py` and the rest of `recompute.sh`'s `tcrvdb` stage → `bench/eda/out/`. The `P_native` producer `bench/scripts/native_bn.py` and its `native_bn_*.csv` outputs no longer exist there; `P_native` itself was discarded in tcren 2.26.0 |
| VDJdb real-versus-mock panel assembly | `bench/scripts/binder_benchmark.py`, `bench/scripts/tcren_binder_score.py` |

Sets `bootstrap_data.py` fetches from HF: `cpl`, `tcrvdb`, `vdjdb_positives`, `vdjdb_negatives`,
`as_case`, `garcia_b27`. Pending upload, aldan3-only for now: `md`, `atlas`, `neoantigen`,
`immrep23`, `immrep25`.
