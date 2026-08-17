# tcren — project status & TODO

Status of the Python re-implementation of TCRen (`src/tcren/`). The legacy R/Java pipeline
is preserved (tag `legacy-r-1.0`) and serves as the numerical oracle. Current release: **v2.8.0**
(feature table + AF-orthogonal kit: `recognize --full --scores`, `kit_score`, `forced_pose_score`,
interface mechanics, binder identification, configurable potentials, fast ΔΔG; `arda-mapper >= 2.5.7`).
See **[CHANGELOG.md](CHANGELOG.md)** for the authoritative per-release record, [BENCHMARKS.md](BENCHMARKS.md)
for achieved accuracy, and `docs/` (`features.rst`, `kit.rst`) for the current API.

> Note: the detailed "Done"/"TODO" sections below are **historical** (they predate v2.1+) and describe an
> earlier module layout — the `native/` module is now `orient/`, FlexPepDock lives in
> `refine/oracle_flexpep.py`, and the standalone `tcren mhc` command was removed. Treat CHANGELOG.md +
> README.md + SKILL.md as the current source of truth.

## Done

| Area | Module(s) | Notes |
|------|-----------|-------|
| **Potentials** | `potential/` | classic + `am` (gap) variants, LOO; wide/long CSV loaders; MJ/Keskin bundled |
| **Structure I/O** | `structure/` | biopython parse; `import_structure` (C-gene trim by default, `keep_c_gene` for MD) |
| **TCR annotation** | `annotation/` | arda V(D)J → CDR/FR markup; αβ/γδ from C-gene (`cgene`) |
| **Contacts** | `contacts/`, `contactmap.py` | cKDTree 5 Å + Cα matrix; TCR/peptide/MHC interfaces |
| **Scoring** | `scoring.py` | substitution scoring; drop-in for `run_TCRen.R`; opt-in TCR framework regions (`cdr`/`cdr+fr`/`all`) |
| **Configurable potentials** | `pipeline.py` | per-interface potential override (`Potential`, bundled name, CSV, or None) on `pipeline.run`/CLI |
| **Percentile rank** | `scoring_rank.py` | native peptide energy vs. random pMHC background (`tcren rank`) |
| **Fast ΔΔG** | `ddg.py` | virtual-matrix point-mutation ΔΔG + alanine scan + neoantigen ΔΔG (`tcren ddg`) |
| **Oracle facade** | `oracle.py` | `summarize_structure` composes S1–S4 into one frame bundle for the paper notebooks |
| **MHC** | `mhc/` | IMGT/HLA + mouse H-2 reference, mmseqs mapping, groove partitioning, linker-peptide split |
| **Native DB** | `native/` | TCR3D download/version/manifest; ground-truth comparison; align-to-canonical; potential re-derivation |
| **2D maps** | `project2d/`, `viz/` | groove-plane projection, canonical tables, metadata-rich SVG, py3Dmol pocket+CDR |
| **Analysis** | `analysis.py` | potential heatmaps/compare, contact distributions (per-structure/region/position) |
| **CLI** | `cli.py` | `info/annotate/contacts/derive-potential/score/rank/ddg/pipeline/recognize/orient/superimpose …` |
| **Docs** | `docs/` | Sphinx + 3 tutorial notebooks (`notebooks/`); zero-warning build |

## In flight (branch `feature/review-aug26-surface-topology`, 2026-08-17)

Acting on `review/rev17aug26.md` PART 1 + the surface-topology ask. All six items landed; see
CHANGELOG `[Unreleased]`. Open loops out of it:

- [x] **`_relax.repack` (C++)** — done. Same input, same atoms: side-chain RMSD 4.131 → **2.364 Å in
  6 ms**, where OpenMM returns 4.133 Å (unchanged) in 3.1 s, because a local minimiser cannot cross a
  torsional barrier. 8/8 improved. `tcren repack` / `tcren refine --repack`.
- [ ] **Side-chain *construction*** — `repack` rotates the side chains a model has; it cannot rebuild
  ones `substitute_peptide` stripped, so that path still returns 44 of 77 heavy atoms. Needs ideal
  internal geometry per residue type (the `Full-atom loop build` row). AlphaFold output is full-atom,
  so `repack` already covers the main case.
- [x] **Flexible-backbone MC** (`_relax.relax_interface` / `tcren.dynamics`) — done, and used to
  test Sewell's intra-peptide-stabilisation hypothesis on the CPL set (2102 structures). Stability
  beats the contact energy in 4/4 clones where the contact model fails and 0/3 where it works.
  The mechanism is supported; the specific 4C6 guess is not (4c6's contact AUC is 0.955 here).
- [ ] **Follow up the stability result**: it is currently a *diagnostic*, not a shipped score. Open
  questions — does it survive on crystal rather than modelled structures; does it hold at more MC
  steps and more seeds; is the combined contact+stability score worth shipping (n = 7 clones is too
  few to tell, Wilcoxon p = 0.22); and does the per-position stability profile localise to the P3/P6
  pair Dolton et al. Fig. S4 names.
- [ ] **Repack inside the MC loop** — `relax_interface` currently samples backbone only, with χ
  fixed. A repack per cycle is affordable now (6 ms) and is what would make it a real FlexPepDock
  analogue rather than a stability probe.
- [ ] **Full-scale fold benchmark on aldan3** — `scripts/fold_benchmark.sbatch`, n ≈ 374 with all
  oracles. FlexPepDock burned 21 min of CPU on six structures locally without finishing.
- [ ] **PART 2 of the review** (deferred, agreed): a C++ kernel for *de novo* peptide placement into
  an empty groove, benchmarked against FlexPepDock. `_refine` + `_fold` already cover template-based
  placement; the rotamer machinery above is its first half.
- [ ] `2wbj` is the one class-II Canonical2026 structure whose β-sheet core still fails to map.
- [ ] Lawrence–Colman shape complementarity (`src/_geom/geom.cpp:13`) — the surface work makes it cheap.

## TODO / pending

- [ ] **AI-model refinement** (`refine/`): batch-refine predicted PDBs → canonical → score; QC (anchor RMSD, plDDT, completeness). Inputs in `data/TCRpMHCmodels/`, `data/Bigot/`, `data/Bobisse/`.
- [ ] **FlexPepDock** (`flexpep/`, optional): peptide substitution + Rosetta relaxation; gated on a discovered Rosetta binary. Needs `keep_c_gene=True`.
- [ ] **Standalone `orient/` module**: generalise `native/align.py` (multi-structure overlay, canonical chain renumbering, write oriented PDBs).
- [ ] **Regenerate stale `tcren_am/` outputs** from the current contact data (see the spawned task).
- [ ] **MHC mapper speed**: prebuild the mmseqs index (currently ~7 s/structure from per-call `easy_search`).
- [ ] **2D map polish**: optional "contacting residues only" mode for less cluttered overlays.
- [ ] Mouse class-II MHC reference is sparse (TRGC3/4 skipped); extend if needed.

## Known caveats

- All bundled structure sets (`data/PDB_structures/`, TCR3D CIFs) are **variable-domain-only**; the C-gene classifier and full-complex geometry need full RCSB inputs (fixtures in `tests/assets/cgene/`).
- TCR3D `tcr_complexes_data.tsv` mislabels some TRAV/DV J calls (e.g. 1bd2 `TRDJ1`); arda is correct (locus follows J). Locked by a test in `arda` dev.
- arda is a runtime dependency, published to PyPI as `arda-mapper>=2.5.7` (imports as `arda`); installed by `pip install -e .` / `pip install tcren`. It auto-fetches its reference and a static mmseqs binary on first use (no conda).
