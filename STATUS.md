# tcren — project status & TODO

Status of the Python re-implementation of TCRen (`src/tcren/`). The legacy R/Java pipeline
is preserved (tag `legacy-r-1.0`) and serves as the numerical oracle. Current release: **v2.9.0**
(MHC class II docking geometry, contact typing `v2`, surface topology, rotamer-averaged contacts, the
native side-chain packer and the flexible-backbone sampler, on top of the v2.8 feature table +
AF-orthogonal kit: `recognize --full --scores`, `kit_score`, interface mechanics, binder
identification, configurable potentials, fast ΔΔG; `arda-mapper >= 2.5.7`).
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

## Shipped in v2.9.0

The August 2026 review's PART 1 plus the surface-topology ask: MHC-II docking geometry, contact typing
`v2`, rotamer-averaged contacts, the type-conditioned potential, peptide-position weighting,
`tcren.surface`, the native side-chain packer and the flexible-backbone sampler. CHANGELOG
`[2.9.0]` is the record, with the measurement that says whether each one worked. Defaults are
unchanged throughout, so no existing number moves unless asked.

## Roadmap

Moved to **[ROADMAP.md](ROADMAP.md)** — the single place for forward plans.

## Known caveats

- All bundled structure sets (`data/PDB_structures/`, TCR3D CIFs) are **variable-domain-only**; the C-gene classifier and full-complex geometry need full RCSB inputs (fixtures in `tests/assets/cgene/`).
- TCR3D `tcr_complexes_data.tsv` mislabels some TRAV/DV J calls (e.g. 1bd2 `TRDJ1`); arda is correct (locus follows J). Locked by a test in `arda` dev.
- arda is a runtime dependency, published to PyPI as `arda-mapper>=2.5.7` (imports as `arda`); installed by `pip install -e .` / `pip install tcren`. It auto-fetches its reference and a static mmseqs binary on first use (no conda).
