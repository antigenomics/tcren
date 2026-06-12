# tcren — project status & TODO

Status of the Python re-implementation of TCRen (`src/tcren/`). The legacy R/Java pipeline
is preserved (tag `legacy-r-1.0`) and serves as the numerical oracle. Gitflow: `master` =
legacy, `dev` = integration, feature branches off `dev` (merged `--no-ff`). See
[BENCHMARKS.md](BENCHMARKS.md) for achieved accuracy/performance and the plan in
`docs/` for design detail.

## Done

| Area | Module(s) | Notes |
|------|-----------|-------|
| **Potentials** | `potential/` | classic + `am` (gap) variants, LOO; wide/long CSV loaders; MJ/Keskin bundled |
| **Structure I/O** | `structure/` | biopython parse; `import_structure` (C-gene trim by default, `keep_c_gene` for MD) |
| **TCR annotation** | `annotation/` | arda V(D)J → CDR/FR markup; αβ/γδ from C-gene (`cgene`) |
| **Contacts** | `contacts/`, `contactmap.py` | cKDTree 5 Å + Cα matrix; TCR/peptide/MHC interfaces |
| **Scoring** | `scoring.py` | substitution scoring; drop-in for `run_TCRen.R` |
| **MHC** | `mhc/` | IMGT/HLA + mouse H-2 reference, mmseqs mapping, groove partitioning, linker-peptide split |
| **Native DB** | `native/` | TCR3D download/version/manifest; ground-truth comparison; align-to-canonical; potential re-derivation |
| **2D maps** | `project2d/`, `viz/` | groove-plane projection, canonical tables, metadata-rich SVG, py3Dmol pocket+CDR |
| **Analysis** | `analysis.py` | potential heatmaps/compare, contact distributions (per-structure/region/position) |
| **CLI** | `cli.py` | `info/annotate/contacts/derive-potential/score/mhc/native …` |
| **Docs** | `docs/` | Sphinx + 3 tutorial notebooks (`notebooks/`); zero-warning build |

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
- arda is installed from `git+https://github.com/antigenomics/arda.git@dev` (or `ARDA_DIR` editable).
