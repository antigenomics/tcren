# Changelog

All notable changes to `tcren` are recorded here. Format follows
[Keep a Changelog](https://keepachangelog.com/); versions follow semantic versioning.

## [Unreleased]

### Added
- **`tcren.viz.pymol`** — the PyMOL figure layer, promoted out of `notebooks/pymol_canonical_figures.ipynb`
  into the library, where it can be tested and reused. Three scene presets (`overlay_scene`,
  `groove_scene`, `interface_scene`), one `render()` that ray-traces them headless, and shared
  styling so panels of one figure are comparable rather than each lit by its own bounding box.
- **A labelled axis gizmo on every panel.** Thin, arrow-headed, in a corner, turning with the
  camera, and named for what the axes mean rather than `x/y/z`: `width` (groove width, α1↔α2),
  `N→C` (groove axis toward the peptide C-terminus) and `TCR` (docking normal, MHC floor → TCR).
  `CANONICAL_AXES` carries those names with the definitions from `tcren.orient.frame` and their
  equivalents in the docking-geometry literature, and a test ties them to that module so the labels
  cannot drift from what orientation actually does. An axis pointing at the viewer foreshortens to
  a dot and its label falls to the lower left of it — the convention for an axis normal to the page
  — instead of piling onto the origin with the other two.

  The gizmo is rendered in its own pass and composited at pixel coordinates rather than projected
  into the corner: PyMOL's orthoscopic viewport does not span the world height that
  `field_of_view` and the camera distance imply (measured on a real scene it is out by about a
  quarter), so the arithmetic route puts the gizmo off-frame. Compositing also means the molecule
  can never occlude it.
- **`residue_importance()` / `importance_scene()`** — colour the interface by which residues carry
  the score. Φ is a sum over residue–residue contacts, so it decomposes exactly: a residue's share
  is the sum of `φ(a_i, a_j)` over the contacts it makes. Two columns come back because they answer
  different questions — `phi` is the energy share, `n_contacts` the geometric share, and a residue
  can be large on one and small on the other. The φ ramp is centred on zero, so blue and red mean
  *favourable* and *unfavourable* rather than merely less and more; a range-fitted ramp would redden
  the least-favourable residue even in an interface where every contact is stabilising. Each
  contact is attributed to **both** residues it joins, so the per-residue values sum to twice Φ —
  an attribution, not a partition, and a test pins that factor of two.
- **`notebooks/pymol_interactive.py`** — a [marimo](https://marimo.io) app
  (`pip install "tcren[marimo]"`): pick a structure and scene, swing the camera and watch the gizmo
  follow, restyle it, colour by residue importance with the numbers beside the render, and rotate a
  live 3Dmol.js view with the mouse. Renders are content-addressed on the scene text *and* every
  render option, so a changed option can never serve a stale panel.
- **A figure gallery in the docs** (`docs/gallery.rst`) — every view family as a rendered example
  with the code beside it, and the axis-gizmo convention written out once where readers will find
  it.

## [2.4.0] — 2026-07-28

### Added
- **`tcren recognize --mechanics` / `recognition_table(…, mechanics=True)`** — the koff proxies
  (stiffness tensor, steered rupture, coupling residues) appended to the descriptor table instead of
  returned as a second one. This is the shape a cohort actually wants: the manuscript's task needed
  both commands, and `tcren mechanics` as a separate run repeats the parse and both mmseqs searches
  to produce a CSV keyed `pdb.id` against `recognize`'s TSV keyed `complex.id`, which then has to be
  joined across the rename. Inside `recognize` the structures are already annotated, so the flag
  costs only the mechanics arithmetic — on 12 crystals, 19.0 s → 19.5 s against 22.5 s for the two
  commands. Values are bit-identical to `tcren mechanics`, and no existing column changes.
- **`mechanics.interface_mechanics(structure, …)`** — the union of `stiffness_tensor`, `rupture` and
  `coupling_residues` under their shipped defaults, and now the single definition of "the mechanics
  row": both `tcren mechanics` and `tcren recognize --mechanics` call it, so the two agree by
  construction rather than by two call sites being kept in step. A test pins that identity.
- **A `lint` job in CI** (`uvx ruff check src tests`). Nothing enforced ruff before and the tree had
  drifted to 77 reports, which is the same as having no linter.

- **`run_pipeline(…, reference_aa="A")` / `tcren scoring --delta`** — the poly-alanine-referenced
  ΔΦ alongside Φ, with the same per-interface breakdown. `F_total` is Φ = Φ_TP + Φ_TM + Φ_PM; the new
  `dF_tcr_pep` / `dF_tcr_mhc` / `dF_pep_mhc` / `dF_total` columns are ΔΦ_TP, ΔΦ_TM (≡ 0 — the peptide
  is not in that interface), ΔΦ_PM and ΔΦ. One command now yields both scores and the whole
  decomposition; ΔΦ is the one to use when each candidate carries its own generated pose. Off by
  default, so the existing `scores` dict is unchanged.
- **`tcren scoring --geometry`** — appends the interface descriptors (`burial`, `n_pep_contacted`,
  `chain_balance`, `n_hbond`, `pitch`, `crossing`) and `Q`, the directional decorrelated
  interface-quality score, by calling `recognition_table` + `cohort.q_score` rather than
  reimplementing them. `tcren recognize` remains the full 35-descriptor catalogue + P(real).
- **`tcren scoring -s` takes many inputs.** A file, a directory, a `.tar.gz`, a quoted glob, a
  `.txt`/`.list`/`.lst` manifest (one path per line, `#` comments, relative paths resolved against
  the manifest), a comma-separated list, or a repeated `-s` — mixed freely. New
  `tcren.structure.io.resolve_sources`; `structure_paths` now handles globs and manifests.
  Also `--contact-weight`, `--skip-errors`.

### Fixed
- **`cohort.q_f_iptm` and `cohort.f_invert_by_iptm` are exported.** Both were imported into the
  package namespace but missing from `__all__`, so the ipTM-gated F path the CLI uses was not part
  of the public API it appears to be.
- **The four `tests/regression/test_orient.py` tests run again.** They read
  `notebooks/data/Native2022/{pid}.pdb`, which fails two ways: that path is gitignored, so a fresh
  checkout has no such directory at all; and where a developer does have it, it holds `.pdb.gz`,
  which the hardcoded `.pdb` does not match. Either way the tests raised `FileNotFoundError` rather
  than skipping. They now resolve against `data/Native2026` — what `setup.sh` fetches — try both
  extensions, and skip cleanly when no reference structure is present.

### Changed
- **`tcren pipeline` is now `tcren scoring`** (breaking). The command never ran the preparation
  pipeline — canonicalisation, region mapping and the Cα/contact/atom-distance matrices are
  `tcren annotate`, `tcren superimpose` and `tcren contacts`. It scores structures. The old name
  is kept as a hidden command that errors with a pointer.
- **`score_row` columns are renamed to match `tcren recognize`** (breaking):
  `tcr_peptide.tcren` → `F_tcr_pep`, `tcr_mhc.mj` → `F_tcr_mhc`, `peptide_mhc.mj` → `F_pep_mhc`,
  `total` → `F_total`, and the `d_*` columns → `dF_tcr_pep` / `dF_tcr_mhc` / `dF_pep_mhc` /
  `dF_total`. The two tables now share one vocabulary and join on `pdb.id`.
- **Ruff is configured rather than merely present** (`[tool.ruff.lint]`). `E702` and `E402` are
  ignored with their reasons: the first is the deliberate `setup; assert` idiom used throughout, the
  second is `pytest.importorskip` before a guarded import. What remains is worth acting on.

## [2.3.2] — 2026-07-24

### Added
- **`cohort.f_score(table)`** — the binder-oriented TCRen contact-energy channel `z(-(F_tcr_pep +
  F_tcr_mhc))`, on the same z-scale as `q_score`. Unlike `Q` (geometry), `F` reads contact chemistry —
  and is **pose-conditional**: it works on well-modelled poses and *inverts* on forced ones (ledger
  C27/C42).
- **`cohort.q_f(table, sign=+1)`** — the pure-tcren combiner `z(Q_geom) + sign·z(F)` (no deep-learning
  term). `sign=+1` (`z(Q)+z(F)`) on clean poses beats raw-AF ipTM on both ROC and PR on template-covered
  epitopes (macro 0.759/0.725 vs 0.692/0.693, ledger C42); `sign=-1` (`z(Q)-z(F)`) is the form that ranks
  on forced poses (GLCTLVAML: 0.71 vs 0.52). Exported as `tcren.f_score`, `tcren.q_f`, `tcren.F_TERMS`.
- **`cohort.q_f_iptm(table, iptm, threshold=0.5)` + `cohort.f_invert_by_iptm(iptm, threshold)`** — the
  **AlphaFold-synergy** path: use ipTM (AF's own pose confidence) to auto-invert F per structure —
  `+z(F)` on confident poses, `-z(F)` on forced (low-ipTM) ones — turning the pose-conditional inversion
  into a single ranking.
- **`tcren recognize --cohort`** now also emits `F_score`, `z(Q)+z(F)` and `z(Q)-z(F)`; with `--iptm` it
  additionally emits `z(ipTM)+z(Q)+z(F)`, the `F_invert` flag and `z(Q)+z(F|iptm)` (pose-adaptive,
  threshold `--invert-f-thresh`), and **prints an advisory** naming how many poses are forced (so the user
  knows F inversion is in play); without `--iptm` it tells the user F is trusted unconditionally and how to
  gate it — the full fit-free panel for AF post-analysis in one line.

## [2.3.1] — 2026-07-24

### Added
- **`cohort.q_iptm(table, iptm, features=Q_FEATURES)`** — the fit-free synergy score `z(ipTM) + z(Q)`
  as one call. `Q` (interface geometry) and the generator's ipTM are near-orthogonal, so the
  standardized sum out-ranks either alone (macro ROC 0.83 vs ipTM 0.79 on TCRvdb; beats raw-AF ipTM on
  both ROC and PR on well-modelled epitopes — benchmark ledger C42). Previously hand-rolled in the
  benchmark; now shipped and exported from `tcren`.
- **`cohort.Q_FEATURES_GEOM`** — the four geometry-only descriptors (`Q_FEATURES` minus the `pp_combo`
  energy contrast). This is `Q_geom`, the AF-orthogonal channel robust to the forced-pose energy
  inversion (C27); pass `features=Q_FEATURES_GEOM` to `q_score`/`q_iptm`.
- **`tcren recognize --iptm META`** — single-line path: reads a metadata TSV/CSV (key column matched to
  `complex.id` + an `iptm`/`tcr-pmhc_iptm` column) and appends `Q_geom` and `z(ipTM)+z(Q_geom)` to the
  recognition table for a directory/tarball/glob of structures.

## [2.3.0] — 2026-07-24

### Added
- **`tcren.cohort` is the recommended fit-free scoring layer.** `q_score` (interface-quality `Q`),
  `strain_z` (crystal-calibrated forced-pose strain), `zscore`, `Q_FEATURES`/`Q_FEATURES_CORE`. Prefer
  these over the fitted `binder.binder_score` (`p_bind`) and `recognition.forced_pose_score`
  (`p_forced`): they carry no training set, so they cannot leak or go stale, and — unlike the fitted
  `p_bind` — `Q` generalises across cohorts. `tcren recognize --scores` now also emits `q_bind` and
  `s_strain` (cohort-relative over the input batch).
- **`tcren.recognition_matrix`** — the per-position × amino-acid substitution-energy landscape, the
  CPL/motif-matrix generalisation of `score_peptides` (either interface side; decomposes the full
  interface score exactly).
- **Graphon / loop geometry featurisation** (`structure → descriptor`, not binder scores):
  `contactmap.registered_map`, `contactmap.binding_mode` (`ModeCentroid`), and `tcren.geometry`
  (`reach_max`, `reachability_floor`, `span_saturation`, `cdr3_internal_coords` / `LoopInternalCoords`).
- **`tcren.stability.contact_stability`** — TCR:peptide contact fragility read straight off the contact
  map: per-contact margin `5 − dmin` to the cutoff, `mean_margin`, `frac_robust`, and `exp_lost` (expected
  contacts lost under a 1 Å shift) — a coordinate-only interface positional-confidence readout.
- **Native `_geom` kernels for interface quality.** `_geom.interface_clashes` (heavy-atom vdW-overlap
  scan, now backing `tcren.clashes`; numpy kept as the reference) and `_geom.contact_stability`.
- **`tcren recognize` emits five interface-quality columns** — `n_clashes`, `clash_score`, `exp_lost`,
  `mean_margin`, `frac_robust` (extra output columns, not part of the 35 model features).

### Changed
- `README`, `BENCHMARKS.md`, `docs/features.rst`, `docs/kit.rst` lead with the fit-free `Q` and disclose
  the AlphaFold baseline choice (ipTM is the weakest of the three confidences on the receptor task).
- Documented that `p_bind` and `FORCED_POSE_MODEL` are fit on labels/rows that no longer fully exist and
  that the fit-free `cohort` scores are the reproducible, transfer-robust alternatives.

### Deprecated
- **`cohort.phi_bind`** now raises `DeprecationWarning`: extending `Q` with the docking-angle term
  degrades ranking (the `z(-pitch)` term is below chance and derived from an AlphaFold-contaminated
  angle). Use `q_score`.

## [2.2.3] — 2026-07-19

### Changed
- **Install is now `uv`-based, no conda.** `setup.sh` creates a repo-local `.venv` with `uv` and
  runs `uv pip install -e .`; `environment.yml` removed. The only host requirement is a C++
  compiler — `arda-mapper` auto-fetches a static `mmseqs2` binary on first use, so no bioconda.
- Bumped `arda-mapper` pin to `>=2.5.7`.

### Fixed
- **Concurrency (SLURM array / Nextflow per-sample).** The on-demand MHC-reference mmseqs index
  build (`tcren.mhc.reference.reference_db`) now serializes through `arda._locking.build_lock`, so
  parallel jobs against a shared `data/mhc_cache` no longer race into a half-written index.

### Internal
- Audit pass: removed duplicated superposition / potential-sum / sigmoid / model-persistence code;
  vectorized the per-interface energy sum (`Potential.as_matrix` gather) on the recognition/pipeline
  hot path; cached bundled potentials and frozen recognizers; assorted docstring/doc fixes.

## [2.2.2] — 2026-07-17

Two data-integrity fixes. Both change output: MJ-based scores and MHC pseudosequence lookups
that previously failed silently now resolve correctly.

### Fixed
- **A–N pair in the bundled MJ/Keskin potentials** (`tcren/data/MJ_Keskin_potentials.csv`). The
  4th lower-triangle slot and its mirror carried a literal `1` where `N`/`A` belong, so the A–N
  pair was absent and a phantom `1` entered the inferred alphabet: `mj()` and `keskin()` built a
  21×21 matrix with 41 `NaN` cells instead of a complete 20×20. Because `as_matrix()` pre-fills
  `NaN` and `scoring.py` sums with `np.nansum`, **every Ala–Asn contact silently contributed 0
  energy** rather than raising. MJ is the default `tcr_mhc`/`peptide_mhc` potential, so MJ-based
  scores shift for any structure with an A–N interface contact. TCRen (a separate file) is
  unaffected, so headline TCRen results do not change. A–N is now 0.15 (MJ) / −2.06 (Keskin);
  the Keskin value is corroborated by `tests/assets/oracle/data/source_data/fig3.csv`, and the
  MJ value matches seqtree 0.6.0's `MJ_CONTACT`. Also regenerates the tracked
  `notebooks/natcompsci2022/data_legacy/MJ_Keskin_potentials.csv.gz` snapshot, which carried the
  identical corruption.
- **Collapsed-allele index in `build_pseudo_fasta.py`** — alleles sharing a 34-mer groove
  pseudosequence were collapsed to `alleles[0]`, discarding the rest (68% of `MHC_pseudo.dat`,
  80% of `pseudosequence.2023.all.X.dat` headers lost), so non-representative alleles such as
  HLA-B\*14:02 and C\*03:04 were unresolvable. Headers are now `>ALLELE [ALLELE ...]|n=<count>`.
  Separately, `_pseudo_index` never split headers on `|`, so 100% of its keys carried the suffix
  and every exact lookup missed.

### Added
- `build_pseudo_fasta.py --imgt-alignments` — derives class-I pseudosequences directly from
  IPD-IMGT/HLA 3.65.0 for alleles NetMHCpan does not cover (it lags IMGT and omits HLA-F).

## [2.2.1] — 2026-07-15

### Changed
- Bumped the default `arda-mapper` pin to `>=2.5.6`.
- PyPI-safe PNG logos in the README (raw SVG does not render on the PyPI project page).

## [2.2.0] — 2026-07-13

Feature table + AlphaFold-orthogonal scoring kit for AI-generated TCR–pMHC structures.

### Added
- **`tcren recognize`** — one flat per-structure table for a set of structures
  (`tcren.recognition`). Default: 35 core interface descriptors + `p_real`/`p_real_bn` (the
  real-vs-shuffled recognizers). `--full`: +18 CDR3-frame (FramePose groove-frame projection) +12
  matrix-swap (TCRen−MJ contrast) descriptors → 65 features. Column reference: `docs/features.rst`.
- **`--scores`** — appends the frozen good-results scores `p_bind` (binder-ID) and `p_forced`
  (`forced_pose_score` — the crystal-natural vs AF-forced strain classifier, 5-fold AUC 0.762).
- **`kit_score`** (`tcren.recognition.kit_score`) — the synergistic `z(p_bind) + z(iptm)` combination
  of the intrinsic binder score with the AlphaFold ipTM; on TCRvdb it beats either alone at precision
  (macro-PR 0.847, P@10% recall 0.969). Decision procedure: `docs/kit.rst`.
- Batched **`recognition_table`** — one arda + one mmseqs MHC call for a whole structure set
  (dataset-scale, ~3× faster than per-structure; byte-exact vs the per-structure path).
- Docs: new `docs/features.rst` (every feature + score) and `docs/kit.rst` (the AI-structure kit).

### Changed
- `recognition_features` gains `full=` and `annotate=` (skip re-annotation in the batch path); reads
  `mhc_class_bin` from `chain_supertype`.

## [2.1.2] — 2026-07-02
CI-health fixes folded into a published release.

## [2.1.1] — 2026-07-02
Re-cut of 2.1.0 with Windows-wheel (MSVC `M_PI`) + rapidfuzz import fixes.

## [2.1.0] — 2026-07-02
Binder identification (5-feature model + `_geom`/`_relax` C++ kernels + CLI), ATLAS ΔΔG harness,
interface mechanics, potential rederivation, legacy 2022 reproduction.

## [2.0.1] — 2026-06-30
Fix `rank` CLI no-candidates default path.

## [2.0.0] — 2026-06-30
Configurable potentials, TCR framework regions, percentile rank, fast ΔΔG, oracle facade.

## [0.1.0] — 2026-06-17
Initial PyPI release setup (publish workflow, `arda-mapper` dependency, lean sdist).
