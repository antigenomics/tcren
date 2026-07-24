# Changelog

All notable changes to `tcren` are recorded here. Format follows
[Keep a Changelog](https://keepachangelog.com/); versions follow semantic versioning.

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
