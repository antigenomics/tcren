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
| `as_case` | native B\*27:05 cognate complexes — held-out ergodicity validation |
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

CPL / TCRvdb / native B\*27:05 held-out validation + the unified three-task table live in the manuscript
(`2026-tcren2/scripts/{benchmark_suite.md,tcren_binder_score.py,tcrvdb_physics.py}`). ATLAS ΔΔG is
pending upstream; MD trajectories (`md100ns_*.xtc`, ~57 GiB) pending upload.

## Evaluated and rejected: VDJdb template-based models (Zenodo 8143087)

Shcherbinin DS, Karnaukhov VK, Zvyagin IV, Chudakov DM, Shugay M. *The database of TCR-peptide-MHC
modeled structures from the paper: "Large-scale template-based structural modeling of T-cell
receptors with known antigen specificity reveals complementarity features."* Zenodo, 2023-04-19.
doi:[10.5281/zenodo.8143087](https://doi.org/10.5281/zenodo.8143087) (concept DOI
10.5281/zenodo.7845843). CC BY 4.0.

- Fetch: `curl -L -o models_annotation.csv 'https://zenodo.org/records/8143087/files/VDJdb_Models_annotation_w_new_names.csv?download=1'`
  (structures: `Models_VDJdb.tar.gz`, 685 MB).
- Content: 3,213 TCR-pMHC models. Columns `vdjdb.pdb.id, cdr3.alpha, cdr3.beta, orig.alpha,
  orig.beta, vdjdb.pdb.id.from, mutation.signature, pdb.id`. **`orig.*` are booleans, not
  sequences**; `vdjdb.pdb.id.from` is the template; `pdb.id` is the model's descriptive name.
- Provenance: **derived/computed, not experimental.** Every model is exactly **one** point
  mutation of its template (`mutation.signature` holds a single `pos X->Y` in all 3,213 rows),
  applied in cascades of depth 1-3.

**Not usable for gap placement or for loop-conformation claims.** Walking the cascade to its
roots leaves **32 independent crystal backbones** carrying all 3,213 sequences (median 28
models per root, max 890). Two consequences, both measured:

- **0 of 32 roots have descendants that differ in CDR3-beta length.** There is not one modelled
  indel anywhere in the dataset, so it carries zero information about where a gap goes.
- Backbone conformation is inherited from the root through point mutations, and root choice is
  made by *sequence similarity* -- which is exactly what defines a sequence island. Testing
  "do distinct islands share a conformation" here would measure template sharing. The set is
  also smaller than what it would supplement: 32 independent backbones against the 199 unique
  junction sequences already in `Canonical2026`.

It remains the right dataset for its own purpose -- side-chain repacking and contact
complementarity at a fixed backbone -- and for measuring how far a single point mutation moves
a loop.

## Derived: junction loop geometry

`scripts/single_block_test.py` — **derived/computed, not experimental.** Where the gap goes when
two junctions differ in length, and whether one gap block suffices.

- Origin: `data/Canonical2026/*.pdb.gz` (374 crystal TCR-pMHC complexes) + `markup_2026.csv`.
- Regenerate: `python scripts/single_block_test.py --d-max 4 --per-cell 400`
- Method: `scripts/_harvest.py` types each junction by which curated CDR3 it reproduces and
  collapses crystal redundancy (372 junctions → **199 unique sequences**, 1.87×), so pairs are
  resampled by junction, never by pair. `tcren.loops.structural_align` supplies a correspondence
  with **no** single-block assumption.
- Results: the true correspondence is a single contiguous block in **95.2–100%** of pairs for
  both chains and every `d = 1..4`; forcing one block costs no median CA-RMSD. The block starts
  at Cys-offset 6 (TRA) and 6–7 (TRB) — the loop apex — and does **not** drift with `d`.
  `central_prior`'s argmin hits it exactly 42.4% (TRA) / 30.1% (TRB) of the time; a
  germline untemplated-span rule hits 0.4% / 19.8% and is rejected.

`scripts/shape_noise_floor.py` — the resolution limit of this structure set. Between two crystals
of the **same** junction sequence (n = 441 pairs): CA-RMSD median **0.239 Å**, (κ,τ) median
**4.447°**. Nothing smaller is measurable here, at any sample size.

**Removed:** `scripts/fit_gap_prior.py` and `data/gap_prior.tsv`. Superseded and wrong in three
ways — λ was chosen by in-sample argmax on the pairs it was evaluated on; the ~29k pairs it drew
from come from only ~200 independent junctions, so its intervals were anticonservative by several
fold; and its ground truth (`structural_block_position`) is an argmin over the single-gap-block
family, so it could rank block positions but could never test whether one block was enough.
`single_block_test.py` answers both questions with an oracle that can disagree.

**Removed:** `scripts/island_conformation.py` → `scripts/shape_noise_floor.py`. Its finding
("distinct sequence islands for one epitope share a (κ,τ) conformation", −3.470°, p = 0.004) is
**withdrawn**: it compared κ/τ over the shared *prefix* of two loops (the wrong correspondence)
while its CA-RMSD sibling used the gap-block one on a different subset; its islands came from a
fixed threshold at which random control junctions cluster *harder* than real ones; and both
claimed effects are smaller than the crystal noise floor above (0.78× and 0.31×).
