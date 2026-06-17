# TCRen — Nat Comput Sci 2022 reproduction

These notebooks reproduce the analyses of `code_paper/*.Rmd` (Karnaukhov et al., *Nature
Computational Science* 2022) using **only the `tcren` Python pipeline** — no `mir.jar`. The
2022 results serve as a reference/oracle; the headline analyses are recomputed on the 2026
structure set with the current pipeline.

## Environment

These notebooks use the shared `tcren-nb` conda env for all notebooks. Set it up once from the
parent `notebooks/` directory:

```fish
cd .. ; bash setup.sh        # creates the tcren-nb env + editable installs + Jupyter kernel
```

See [`../README.md`](../README.md) for details. Select the **"Python (tcren-nb)"** kernel in
VS Code / Jupyter.

## Data layout

- **`../data/`** — shared inputs (gitignored structures): HF structure sets under
  `structures/` (`Native2022`, `Native2026`, `Tcr3d2026`, …), `vdjdb/`, and the allowed
  external published inputs (MJ/Keskin matrix, IEDB, Birnbaum, Bobisse, Bigot). New results
  are computed **only** from here + the `tcren` pipeline.
- **`data_legacy/`** — the legacy `mir` outputs (old TCRen matrix, contacts, MHC annotation,
  paper source data). **Comparison oracle only — never a pipeline input.**
- **`results_new/`** — outputs produced by these notebooks (contacts, markup, TCRen
  matrices, benchmark tables).

The 2026 analyses use **`Native2026`** as the structure source (covers the full
non-redundant set). Annotation is always **batched**: all sequences are gathered and passed
to a single `mmseqs2` call (it parallelises internally) — see `tcren.paper.mhc_annotation`
and `annotate_structure_set`.

## Notebooks

| # | Notebook | Reproduces | Output |
|---|----------|------------|--------|
| 00 | `00_bootstrap.ipynb` | data bootstrap (HF fetch, external inputs) | populates `../data/` |
| 01 | `01_nonred_and_derivation.ipynb` | Rmd 1 + 2 — non-redundant set + TCRen derivation | `TCRen_2022/2026[_LOO].csv` |
| 02 | `02_benchmark_cognate_unrelated.ipynb` | Rmd 3 — cognate vs random (Fig 2a) and vs IEDB-matched (Fig 2b) decoys | `benchmark_cognate[_iedb]_ranks.csv`, `mhc_2026.csv` |
| 03 | `03_benchmark_yeast_display.ipynb` | Rmd 4 — Birnbaum yeast-display (2b4/226/5cc7) ROC + enrichment correlation | `benchmark_birnbaum.csv` |
| 07 | `07_compare_legacy.ipynb` | regression vs the `data_legacy` mir oracle | — |

*Remaining (planned, sequential):* Rmd 5 (Bobisse neoepitopes), Rmd 6 (shuffled structures),
Rmd 7 (modelled structures + Bigot).

## Running

Run top-to-bottom from a clean kernel. Notebook 01 must run before 02/03 (they read
`results_new/`). If executing headless, pass the kernel explicitly:

```fish
jupyter nbconvert --to notebook --execute --inplace \
  --ExecutePreprocessor.kernel_name=tcren-nb 03_benchmark_yeast_display.ipynb
```
