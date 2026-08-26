# tcren — project status & TODO

Status of the Python re-implementation of TCRen (`src/tcren/`). The legacy R/Java pipeline
is preserved (tag `legacy-r-1.0`) and serves as the numerical oracle. Current release: **v2.12.1**.

The recommended score is **`P_native`**, and it is two commands:

```console
tcren features  -s models/ -i placement,interface,topology,energetics -o feats.tsv
tcren recognize --features feats.tsv -o scores.tsv
```

`scores.tsv` carries `complex.id`, `Q`, the three channel posteriors `G` / `T` / `E`, and
`P_native`. `P_native` is the posterior of a latent class over three channels — geometry, footprint
topology and contact energetics — each a conditional-linear-Gaussian Bayes network fitted by
expectation maximization with **no binding label**, the three combined by adding log-odds. `Q` is
the fit-free single-structure interface-quality score. Underneath sit MHC class II docking geometry,
contact typing `v2`, surface topology, rotamer-averaged contacts, the native side-chain packer and
the flexible-backbone sampler, TCRen2 as the default TCR:peptide potential, interface mechanics,
configurable potentials and fast ΔΔG (`arda-mapper >= 2.5.7`).

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
| **Recognition** | `recognition.py`, `cohort.py` | the descriptor table (`tcren features`) and the scores built on it (`tcren recognize`): `Q`, the channels `G`/`T`/`E`, `P_native` |
| **Scoring** | `scoring.py` | substitution scoring; drop-in for `run_TCRen.R`; opt-in TCR framework regions (`cdr`/`cdr+fr`/`all`) |
| **Configurable potentials** | `pipeline.py` | per-interface potential override (`Potential`, bundled name, CSV, or None) on `pipeline.run`/CLI |
| **Percentile rank** | `scoring_rank.py` | native peptide energy vs. random pMHC background (`tcren rank`) |
| **Fast ΔΔG** | `ddg.py` | virtual-matrix point-mutation ΔΔG + alanine scan + neoantigen ΔΔG (`tcren ddg`) |
| **Oracle facade** | `oracle.py` | `summarize_structure` composes S1–S4 into one frame bundle for the paper notebooks |
| **MHC** | `mhc/` | IMGT/HLA + mouse H-2 reference, mmseqs mapping, groove partitioning, linker-peptide split |
| **Native DB** | `native/` | TCR3D download/version/manifest; ground-truth comparison; align-to-canonical; potential re-derivation |
| **2D maps** | `project2d/`, `viz/` | groove-plane projection, canonical tables, metadata-rich SVG, py3Dmol pocket+CDR |
| **Analysis** | `analysis.py` | potential heatmaps/compare, contact distributions (per-structure/region/position) |
| **CLI** | `cli.py` | `info/annotate/contacts/derive-potential/score/rank/ddg/features/recognize/scoring/orient/superimpose …` (`pipeline` is a hidden alias of `scoring`) |
| **Docs** | `docs/` | Sphinx + 3 tutorial notebooks (`notebooks/`); zero-warning build |

## Shipped in v2.9.0

The August 2026 review's PART 1 plus the surface-topology ask: MHC-II docking geometry, contact typing
`v2`, rotamer-averaged contacts, the type-conditioned potential, peptide-position weighting,
`tcren.surface`, the native side-chain packer and the flexible-backbone sampler. CHANGELOG
`[2.9.0]` is the record, with the measurement that says whether each one worked. Defaults are
unchanged throughout, so no existing number moves unless asked.

## Shipped in v2.11.0 – v2.12.1

- **v2.11.0** — **TCRen2 is the default TCR:peptide potential** (it was `karnaukhov2022`), and the
  shipped matrix is re-derived on the **362** fully annotated αβ TCR:pMHC complexes of Native2026,
  down from all 374. Scores are not comparable across that boundary; pass `-p karnaukhov2022` for
  the old one.
- **v2.12.0** — **`P_native` and its three channels are the recommended score**, and the
  hand-picked combiners it replaces are gone. Removed, with zero callers in the library and none in
  the benchmark reproduction path: the `pose_sweep` module and `pose.c_score` with its two bundled
  reference manifolds (492 KB off the wheel), `footprint.footprint_score`, `cohort.q_iptm` / `q_f` /
  `q_f_iptm` / `f_invert_by_iptm` / `phi_bind` / `agreement`, `recognition.kit_score`, and the two
  fit scripts that regenerated the deleted CSVs. `cohort.coupling` and `cohort.q_coupled` are
  **deprecated but still importable**, byte-identical in behaviour, so every published `S`
  reproduces. See [OBSOLETE.md](OBSOLETE.md).
- **v2.12.1** — **`paths.tcren_home()`**, the on-disk root for reference data: `$TCREN_HOME`, else
  the source checkout recognised by its `pyproject.toml`, else `$XDG_CACHE_HOME/tcren`. Before it,
  every root was derived from the checkout layout, so an **installed wheel could not find the MHC
  allele reference at all** and `tcren annotate` failed on every structure.

## Roadmap

Moved to **[ROADMAP.md](ROADMAP.md)** — the single place for forward plans.

## Known caveats

- **Reference data is fetched or built, not bundled.** The MHC allele reference is built on demand
  from IMGT by `tcren build-mhc-ref`; the structure sets come from `$TCREN_DATA_DIR` or
  `tcren fetch-data`. Both resolve under `paths.tcren_home()`.
- All bundled structure sets (`data/PDB_structures/`, TCR3D CIFs) are **variable-domain-only**; the C-gene classifier and full-complex geometry need full RCSB inputs (fixtures in `tests/assets/cgene/`).
- TCR3D `tcr_complexes_data.tsv` mislabels some TRAV/DV J calls (e.g. 1bd2 `TRDJ1`); arda is correct (locus follows J). Locked by a test in `arda` dev.
- arda is a runtime dependency, published to PyPI as `arda-mapper>=2.5.7` (imports as `arda`); installed by `pip install -e .` / `pip install tcren`. It auto-fetches its reference and a static mmseqs binary on first use (no conda).
