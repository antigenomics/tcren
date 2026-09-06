# tcren — project status & TODO

Status of the Python re-implementation of TCRen (`src/tcren/`). The legacy R/Java pipeline
is preserved (tag `legacy-r-1.0`) and serves as the numerical oracle. This tree is **v3.6.0**
(2026-09-06) and **not published**. The newest release on PyPI is **v3.0.1**, uploaded 2026-09-02
as 12 wheels (cp310-cp313 x macOS arm64 / manylinux x86_64 / win_amd64) plus the sdist. Six
releases therefore exist in the checkout only: 3.1.0 (batched annotation), 3.2.0 (`tcren explain`),
3.3.0 (the TCRen2 matrix re-derived on 359 crystals), 3.4.0 (the self-consistent arm on the command
line), 3.5.0 (`tcren potential-table`) and 3.6.0 (the Ca geometry prior and `tcren potts predict`).

## The model lineage — six shipped objects, and how they relate

**One estimator, counted three ways, and one contact model that reads whichever matrix it is
pinned to.** This is the single place that relationship is stated; `CHANGELOG.md` is the
per-release record and `data/potentials.json` the reproduction manifest.

| name on the command line | what it is | how it is fitted | ships as |
|---|---|---|---|
| `karnaukhov2022` | the residue-pair contact log-odds, `-ln[n_ab N / (n_a. n_.b)]` | tallied over the crystals available in 2022, unweighted | `TCRen_potential.csv` |
| **`tcren2`** *(default)* | the same estimator | the same tally over 359 alpha-beta `Native2026` crystals, weighted per structure by `--balance both` | `TCRen2_potential.csv` |
| `tcren2_em` | the same estimator over every **available** pair, each weighted by its contact probability | the self-consistent loop of `potts selfcon`, converged at step 9 | `TCRen2_sc_potential.csv` |
| **`potts_tcr_peptide`** *(default)* | `P(sigma)` over the 15 A Ca site set: one-body field plus couplings between pairs sharing a residue | penalised pseudolikelihood, pair term pinned to `tcren2` | `potts_tcr_peptide.json` |
| `potts_tcr_peptide_sc` | the same model, pair term pinned to `tcren2_em` | refitted at every step of the self-consistent loop | `potts_tcr_peptide_sc.json` |
| `geom_tcr_peptide` | `P(d_ca \| loop, i, L, j, M, MHC class)` — the Ca distance matrix itself | graphon average over the same 359 crystals, superimposed; 12 surfaces at grid 12 | `geom_tcr_peptide.json` |

Two things not to relearn:

- **Nothing was switched.** `tcren2()` is still the default potential and `PottsModel.bundled()`
  still the default contact model: on one identical site table the pinned self-consistent model is
  a *worse* pair term than a pin on TCRen2 at the same parameter count (k = 114 both). `tcren2_em`
  ships beside, not instead.
- **`tcren2_em` is the only shipped matrix that does not reproduce at 1e-9**, because its E-step is
  a block-Gibbs sample. Its fixed point does reproduce: three runs at different sampling budgets
  agree at Pearson *r* = 0.999997. The manifest records a 5e-3 tolerance for it and checks the
  recipe only under `RUN_BENCHMARK=1`; every other entry still reproduces at 1e-9.

The geometry prior is what makes the site set itself reachable without a structure. `available_pairs`
reads `d_ca` off the coordinates; `potts predict` reads it off `geom_tcr_peptide` instead, so a
contact map can be placed and scored from loop and peptide sequences alone.


The recommended read-out is the **score set** (`tcren.score`), and it is two commands:

```console
tcren features -s models/ -i placement,interface,topology,energetics -o feats.tsv
tcren assess   --features feats.tsv -o scores.tsv
```

`scores.tsv` carries one row per structure, and **every read-out in it is defined for a single
structure**: the transform, the class means and the covariance are frozen on a hold-out that ships
in the wheel, so nothing is estimated from the rows being scored and a score does not move
depending on what was scored beside it.

| read-out | tier | what is estimated |
|---|---|---|
| `peptide_score` | 0 | nothing; the direction is fixed by the potential |
| `pose_score` | 1 | a covariance over hold-out binders — no negative, no label |
| `confidence_residual` | 1 | the same covariance, read as a conditional mean |
| `binder_score` | 2 | class means and covariances, from hold-out binder labels |
| `channel_scores` | 2 | the same object, marginalized to one descriptor family |

The five channels are `placement` (where the receptor sits in the groove frame), `interface` (how
much interface it makes, of what chemistry), `shape` (the footprint free of its size), `energetics`
(the contact chemistry in kT) and `mechanics` (the interface as a network of breakable springs).
Where the feature table carries the generator's confidence, `binder_iptm` — `binder_score` plus
`logit(ipTM)`, two log-odds added with no coefficient to fit — is the recommended read.

**Why this is not `P_native` again.** `P_native` was discarded in 2.26.0 because its coefficients
were frozen against a training set no reader could reconstruct. These are frozen against a named
one: `holdout_manifest()` returns the **8,292 hold-out structures** with their dataset, epitope,
binder label and ipTM, they ship inside the wheel, and `tcren fit-holdout` regenerates the shipped
arrays from them bit for bit.

The **fit-free predecessor tier** is still shipped and still reported: `cohort.q_score` (`Q`),
`reliability.t_score` (`T`) and `reliability.s_score` (`S`), which `tcren recognize --features`
emits. `S` **composes** with `binder_score` rather than being replaced by it; the numbers are in
[BENCHMARKS.md](BENCHMARKS.md).

Underneath sit the **164-descriptor catalogue** in six families (`placement` 31, `interface` 26,
`topology` 70, `energetics` 15, `potts` 5, `kinetics` 17), MHC class II docking geometry, contact
typing `v2`, surface topology, rotamer-averaged contacts, the native side-chain packer and the
flexible-backbone sampler, TCRen2 as the default TCR:peptide potential, interface mechanics,
configurable potentials and fast ΔΔG (`arda-mapper >= 2.5.7`).

See **[CHANGELOG.md](CHANGELOG.md)** for the authoritative per-release record, [BENCHMARKS.md](BENCHMARKS.md)
for achieved accuracy, and `docs/` (`assess.rst`, `features.rst`, `kit.rst`) for the current API.

> Note: the detailed "Done"/"TODO" sections below are **historical** (they predate v2.1+) and describe an
> earlier module layout. The `native/` module is gone — RCSB fetching is `recent.py`,
> align-to-canonical is `docking/align.py` and re-derivation is `potential/derive.py`; `orient/` was
> renamed `docking/` in 2.29.0 when the flat middle of the package became the four sub-packages
> `docking`, `topology`, `energetics` and `mechanics`, each old top-level name kept as a transparent
> re-export; FlexPepDock lives in `refine/oracle_flexpep.py`; and the standalone `tcren mhc` command
> was removed. Treat CHANGELOG.md + README.md + SKILL.md as the current source of truth.

## Done

| Area | Module(s) | Notes |
|------|-----------|-------|
| **Potentials** | `potential/` | classic + `am` (gap) variants, LOO; wide/long CSV loaders; MJ/Keskin bundled |
| **Structure I/O** | `structure/` | biopython parse; `import_structure` (C-gene trim by default, `keep_c_gene` for MD) |
| **TCR annotation** | `annotation/` | arda V(D)J → CDR/FR markup; αβ/γδ from C-gene (`cgene`) |
| **Contacts** | `contacts/`, `contactmap.py` | cKDTree 5 Å + Cα matrix; TCR/peptide/MHC interfaces |
| **Recognition** | `descriptors/`, `recognition.py` | the 164-descriptor catalogue in six families and the batched table behind `tcren features`; `RECOGNITION_FEATURES` is a 40-column subset of it, not the catalogue |
| **Scores** | `score/`, `reliability.py`, `cohort.py` | the frozen score set behind `tcren assess`, and the fit-free `Q`/`T`/`S` behind `tcren recognize --features` |
| **Scoring** | `scoring.py` | substitution scoring; drop-in for `run_TCRen.R`; opt-in TCR framework regions (`cdr`/`cdr+fr`/`all`) |
| **Configurable potentials** | `pipeline.py` | per-interface potential override (`Potential`, bundled name, CSV, or None) on `pipeline.run`/CLI |
| **Percentile rank** | `scoring_rank.py` | native peptide energy vs. random pMHC background (`tcren rank`) |
| **Fast ΔΔG** | `energetics/mutation.py` | virtual-matrix point-mutation ΔΔG + alanine scan + neoantigen ΔΔG (`tcren ddg`); `tcren.ddg` re-exports it |
| **Oracle facade** | `oracle.py` | `summarize_structure` composes S1–S4 into one frame bundle for the paper notebooks |
| **MHC** | `mhc/` | IMGT/HLA + mouse H-2 reference, mmseqs mapping, groove partitioning, linker-peptide split |
| **Native DB** | `recent.py`, `docking/align.py`, `potential/derive.py` | RCSB fetch/discover (`tcren fetch-recent`), align-to-canonical, potential re-derivation. The `native/` package that once held all three is gone |
| **2D maps** | `project2d/`, `viz/` | groove-plane projection, canonical tables, metadata-rich SVG, py3Dmol pocket+CDR |
| **Analysis** | `analysis.py` | potential heatmaps/compare, contact distributions (per-structure/region/position) |
| **CLI** | `cli.py` | `annotate`, `assess`, `build-mhc-ref`, `contacts`, `cpl`, `ddg`, `derive-potential`, `energy`, `explain`, `features`, `fetch-data`, `fetch-recent`, `fit-holdout`, `info`, `mechanics`, `orient`, `rank`, `recognize`, `refine`, `score`, `scoring`, `shuffle`, `substitute-tcr`, `superimpose`, `surface`, plus the `potts` and `paper` sub-apps — 27 in all, and `tests/unit/test_cli_entrypoints.py` asserts `python -m tcren.cli` and the console script expose the same 27. `footprint` and `pipeline` are hidden; `binder` and `diagnose` were removed in 2.26.0 and 2.28.0 |
| **Docs** | `docs/` | Sphinx + 11 tutorial notebooks (`notebooks/`, all in the `docs/index.rst` toctree); zero-warning build |

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
  reproduces. See [OBSOLETE.md](OBSOLETE.md). `P_native` was itself discarded in 2.26.0 — see below.
- **v2.12.1** — **`paths.tcren_home()`**, the on-disk root for reference data: `$TCREN_HOME`, else
  the source checkout recognised by its `pyproject.toml`, else `$XDG_CACHE_HOME/tcren`. Before it,
  every root was derived from the checkout layout, so an **installed wheel could not find the MHC
  allele reference at all** and `tcren annotate` failed on every structure.

## Shipped in v2.13.0 – v3.0.0

Condensed; [CHANGELOG.md](CHANGELOG.md) is the record. Two reversals run through this stretch —
v2.26.0 discarding `P_native`, and v2.28.0 removing every read-out fitted out of fold — and between
them they are why the recommendation at the top of this file changed twice: `P_native` in v2.12.0,
`S_free` in v2.15.0, the score set in v3.0.0.

- **v2.13.0 – v2.14.0** — `tcren.potts`: the contact map modelled as a random variable rather than
  read off a structure, `E(σ) = −Σ η_a σ_a − ½ Σ A_ab σ_a σ_b`, with the *non*-contact as an
  observable. `Π`, the interface energy read against the partition function, comes from here.
- **v2.15.0 – v2.16.0** — **`S_free`, a composite you can put on one structure**:
  `Q/sd_Q + T/sd_T + (Π − μ)/sd_Π`, three fit-free directional blocks over the 374 Native2026
  crystals, each divided by its own native spread. `reliability.t_score` (`T`) is the block that
  survives template scarcity — on the balanced VDJdb panel `T` loses 0.06 ROC-AUC when the epitope
  has no solved complex to template on, against `Q`'s 0.24. Plus `inversion_flag` (the forced-pose
  detector) and `screening_yield` (what a testing budget buys), neither of which fits anything.
- **v2.17.0 – v2.25.0** — three defects in `p_native` fixed before it was discarded; `af_band`, the
  observed non-binder fraction per confidence quantile with a Wilson interval; the contact map read
  as a per-residue-pair marginal and then per residue; the presentation interface given its own
  Hamiltonian with a second potential and the gauge (`centred_potential`) that makes a pinned model
  reproduce a referenced score; and an alanine scan that moves atoms on both sides of the interface.
- **v2.26.0** — **`P_native` is discarded, not deprecated**, together with the whole v1 score block
  and every other frozen-coefficient composite; `F_*` becomes `Phi_*`; `tcren.provenance` stamps
  every feature table with a catalogue digest and `tcren recognize --features` refuses a table
  written under a different one. The catalogue is descriptors only. See
  [OBSOLETE.md](OBSOLETE.md) for the full removal list.
- **v2.27.0** — the composite is **`S`**: `reliability.s_free` → `reliability.s_score`, and the
  emitted column `S_free` → `S`. **There is no alias** — a caller on the old name fails loudly.
- **v2.28.0** — **every out-of-fold-fitted read-out is removed**: `reliability.p_binder` and
  `available_links`, `correct_confidence` and `available_corrections`, and the `tcren diagnose`
  command. `Q`, `T`, `Π`, `S`, `inversion_flag`, `screening_yield` and `af_band` are untouched.
- **v2.29.0 – v2.30.0** — the catalogue grows 123 → 141 → **164** descriptors (contact-graph
  evenness, spectral readings of the Cα/Cβ maps, then 23 from the published interface literature,
  each measured against all 141 incumbents before adoption), and the flat middle of the package
  becomes the four named sub-packages `docking`, `topology`, `energetics` and `mechanics`.
- **v3.0.0** — **`tcren.score`**, the score set at the top of this file, with `tcren assess` and
  `tcren fit-holdout`. The major version is for the new public surface and the new package data,
  not for a removal: every 2.30.0 entry point still works and returns what it returned.

**`S` names two different quantities across this history.** In v2.12.0 it was
`cohort.q_coupled(Q, ΔΦ)`; from v2.15.0 it is the three-block composite above, shipped as `S_free`
and renamed `S` in v2.27.0. A number quoted as `S` must name which.

## Roadmap

Moved to **[ROADMAP.md](ROADMAP.md)** — the single place for forward plans.

## Known caveats

- **A feature table is only valid under the catalogue that wrote it.** `tcren features` writes a
  `<name>.provenance.json` carrying a SHA-256 digest of the descriptor catalogue, and both
  `tcren recognize --features` and `tcren assess` **refuse** a table written under a different one,
  naming the command that regenerates it. The digest moved at 2.29.0 and again at 2.30.0.
- **`RECOGNITION_FEATURES` is a 40-column subset, not the catalogue.** The catalogue is the 164
  descriptors of `recognition.DESCRIPTORS`; do not read the shorter tuple as the whole of it.
- **Reference data is fetched or built, not bundled.** The MHC allele reference is built on demand
  from IMGT by `tcren build-mhc-ref`; the structure sets come from `$TCREN_DATA_DIR` or
  `tcren fetch-data`. Both resolve under `paths.tcren_home()`.
- All bundled structure sets (`data/PDB_structures/`, TCR3D CIFs) are **variable-domain-only**; the C-gene classifier and full-complex geometry need full RCSB inputs (fixtures in `tests/assets/cgene/`).
- TCR3D `tcr_complexes_data.tsv` mislabels some TRAV/DV J calls (e.g. 1bd2 `TRDJ1`); arda is correct (locus follows J). Locked by a test in `arda` dev.
- arda is a runtime dependency, published to PyPI as `arda-mapper>=2.5.7` (imports as `arda`); installed by `pip install -e .` / `pip install tcren`. It auto-fetches its reference and a static mmseqs binary on first use (no conda).
