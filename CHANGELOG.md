# Changelog

All notable changes to `tcren` are recorded here. Format follows
[Keep a Changelog](https://keepachangelog.com/); versions follow semantic versioning.

## [2.12.1] — 2026-08-27

**An installed wheel could not find its reference data.** Every on-disk root was derived from the
source-checkout layout (``Path(__file__).parents[2]`` / ``[3]``), which resolves to
``site-packages``' parent once the package is installed. A wheel therefore looked for the MHC allele
reference at ``<venv>/lib/python3.x/database/mhc/alleles.aa.fasta``, so ``tcren annotate`` and every
command built on it failed on every structure — reproducible only from a git checkout.

### Added
- **`paths.tcren_home()`** — the root for tcren's on-disk reference data. ``$TCREN_HOME`` when set;
  otherwise the source checkout, recognised by its ``pyproject.toml``; otherwise
  ``$XDG_CACHE_HOME/tcren`` (``~/.cache/tcren``), which is writable and stable across upgrades.

- **`notebooks/pnative_channels.py`** — a marimo app that runs the released scoring path end to end:
  featurise a directory, fit each channel, combine into `P_native`, and read the geometry-versus-
  energetics correlation whose sign marks a templated pose.

### Fixed
- **`tcren footprint --score` was broken on every call.** It imported `footprint_score`, deleted in
  2.12.0. `--score`/`--group` now emit `T`, the shape channel's posterior, fitted per group where one
  is given. This is the replacement 2.12.0's own changelog named.
- **`tcren scoring` exited 0 when every structure failed.** The caller then died several stages later
  on a missing energy column, nowhere near the cause. A run in which every row carries an error now
  exits non-zero.
- `mhc.reference.DATABASE_DIR` / `CACHE_DIR` and `paper.bootstrap._REPO` resolve through
  `tcren_home()`, so `tcren build-mhc-ref` writes where an installed `tcren annotate` reads.
- `mhc.reference`'s docstring said the allele reference was *committed*. It is gitignored and built
  on demand from IMGT; the docstring now says so.

## [2.12.0] — 2026-08-24

**`P_native` is the recommended score, and the combiner zoo around it is gone.** Three channels —
geometry, topology, energetics — each fitted as its own latent-class Bayes network by EM, their
log-odds added. Adding log-odds is the exact posterior only across channels that are conditionally
independent given the class, so `geometry` pools the `placement` and `interface` descriptor
*families* into one network: they are the most dependent pair measured (|ρ| = 0.244 between
principal components on the VDJdb benchmark, against 0.023 for topology vs interface). Summing them
as two terms counts that dependence twice, worth 0.817 → 0.832 macro ROC on TCRvdb and 0.668 →
0.718 on VDJdb real-vs-mock.

### Added
- **`cohort.p_native(rule=...)`** — `"sum"` (default) adds per-channel log-odds, `"flat"` pools
  every channel's features into one network. `return_model=True` returns a `{channel: model}`
  mapping under `"sum"`.
- **`cohort.P_NATIVE_POOL`** maps each combination channel to the descriptor families it draws on,
  and **`cohort.P_NATIVE_ORIENT`** gives each channel's default orientation feature.
- **`GaussianBNClassifier.fit_em(orient_by=...)` accepts a leading `"-"`** meaning *lower is
  native-like*. The energetics channel needs it: Φ is a contact-preference sum in which lower is
  favourable, so orienting on the raw `F_tcr_pep` labelled the unfavourable component native.
- `tcren recognize --features` emits `P_native` and its three channel posteriors (`G`, `T`, `E`).

### Changed
- **`cohort.P_NATIVE_CHANNELS` is now `("geometry", "topology", "energetics")`**, three names where
  2.11 had four. `P_NATIVE_FEATURES` stays keyed by descriptor *family* (four keys); resolve a
  channel through `P_NATIVE_POOL`, or call `cohort._channel_columns(channel)`.
- **`p_native(..., return_model=True)` returns a dict, not a model.** Callers unpacking a single
  model must either index the channel they want or pass `rule="flat"`.

### Removed
Every one of these had zero callers in the library and zero in the benchmark repo's reproduction
path. They are the superseded combiners `P_native` replaces and the pose-consistency experiment it
made unnecessary.

- **`pose_sweep` module** (605 lines) and **`pose.c_score`** with its two bundled reference
  manifolds, `pose_af_reference.csv` and `pose_native_reference.csv` — **492 KB off the wheel.**
  The AF reference's own docstring recorded the reason: scored against the crystal manifold the
  score reads *provenance*, not model quality. `pose.pose_consistency` and the `POSE_FEATURES*`
  tuples are unchanged.
- **`footprint.footprint_score`** (the `fp_score` z-sum) — use `p_native(t, channels=("topology",))`.
- **`cohort.q_iptm`, `q_f`, `q_f_iptm`, `f_invert_by_iptm`, `phi_bind`, `agreement`** — hand-picked
  combination rules superseded by a fitted one. `agreement` was the per-structure summand of `C*`,
  which the manuscript no longer uses.
- **`recognition.kit_score`** — a z-sum of `p_bind` and ipTM that every caller already wrote inline.
- `scripts/fit_pose_reference.py`, `scripts/fit_joint_reference.py`, which regenerated the two
  deleted CSVs.

### Deprecated
- `cohort.coupling` and `cohort.q_coupled` remain importable, tested, and byte-identical in
  behaviour, so every published `S` reproduces. They are superseded by `p_native`, which fits each
  channel's sign instead of measuring it.

## [2.11.0] — 2026-08-23

**TCRen2 is now the default TCR:peptide potential, and it is re-derived on the fully annotated
αβ subset of Native2026.** Both changes move numbers: any score produced by an earlier release
without an explicit `--tcr-peptide-potential` used the 2022 matrix, and the shipped TCRen2 matrix
itself is not the one 2.10.0 shipped. Re-run anything you are comparing across this boundary.

### Changed
- **Default TCR:peptide potential is `tcren2`, not `karnaukhov2022`.** `pipeline._INTERFACE_POTENTIAL`,
  `tcren recognize` (`recognition.py`) and the CLI's `-p/--potential` fallback all resolved to the
  2022 matrix, so the released default disagreed with the matrix the manuscript reports. Pass
  `--tcr-peptide-potential karnaukhov2022` (or `-p karnaukhov2022`) for the old behaviour;
  `tcren.potential.tcren()` is unchanged and still loads the 2022 matrix.
- **`TCRen2_potential.csv` re-derived on the 362 fully annotated αβ TCR:pMHC complexes** of
  Native2026, down from all 374. The 12 dropped are 3 pMHC-only files (3gjf, 3hae, 4wuu), 8 carrying
  a single αβ chain (3nfj, 5xot, 5xov, 6bj3, 6bj8, 8yiv, 8yj2, 3tf7) and one γδ receptor (4qrr).
  Every TCR chain present in those files *is* annotated with its CDR3 — the chains are absent from
  the crystal, not from the annotation — but `balanced_weights` skips a structure with a null on any
  axis and `derive_tcren` then defaults it to weight **1.0**, the maximum, so three near-duplicate
  pairs (3nfj/5xov, 5xot/6bj3, 8yiv/8yj2) were each counted twice at full weight. Measured against
  the 374-structure matrix: TCRvdb receptor ranking macro-r **+0.034 [+0.003, +0.070]**, with the
  ergodic bridge, the Garcia EC50 series and CPL all inside noise. r = +0.966 against the old matrix,
  max |d| 1.011 — scores are not comparable across the change.

- **`derive-potential` now derives from αβ TCR:pMHC only, unconditionally.** A structure missing
  either CDR3 or the peptide is dropped, and there is no flag to turn that off: `tcren` is for
  αβ TCR : peptide-MHC (class I or II, standard amino acids) and nothing else. This is what
  produces the shipped TCRen2; the recipe in `data/potentials.json` is unchanged and
  `tests/regression/test_shipped_potentials.py` still reproduces the file bit-for-bit.

## [Unreleased]

### Added
- **`scripts/relax_openmm.py`** — full-complex OpenMM minimization (amber14 + GBn2 implicit solvent,
  all atoms free), the physics relaxation `refine_peptide` deliberately is not, and 10-30x faster
  than a Rosetta FastRelax. Relieves the interface strain of an AlphaFold forced pose without moving
  the model off its pose, and puts a deposited structure into the state an all-atom MD run scores it
  in rather than the state it was deposited in. Takes an input directory (`.pdb` or `.pdb.gz`),
  resumes, and shards across cores. Needs `openmm` and `pdbfixer`, which are not tcren dependencies.

## [2.9.0] — 2026-08-18

Everything below acts on an August 2026 review that raised seven points about what a contact
potential cannot see, plus a PART 2 asking how likely a T cell is to recognise a given pMHC at all.
All of PART 1 is implemented, each with the measurement that says whether it worked. Every one is
opt-in: defaults are unchanged, so no existing number moves unless asked.

**Still open from PART 2**: a fast in-house kernel for *de novo* peptide placement into an empty
groove, to compare against FlexPepDock. Two of its three pieces now exist — the side-chain packer and
the backbone sampler below — but building side chains that were never there, and placing a peptide
with no template, are not done. See [`refine/CPP_REWRITE.md`](src/tcren/refine/CPP_REWRITE.md).

### Added
- **MHC class II docking geometry.** `docking_geometry` raised on every class-II complex, so six of
  the 34 recognition features were silently NaN for DR/DQ/DP. Class II is the same β-sheet floor with
  its two pseudo-symmetric halves on separate chains, so the same six within-domain strand offsets
  taken from the α1 (MHCa) and β1 (MHCb) canonical sequences name the corresponding positions.
  **93/94 class-II Canonical2026 structures now compute** (was 0/94), d = 31.5 Å mean against class
  I's 30.6. The class-I path is bit-identical, verified against the previous implementation.
- **`tcren.surface` + `tcren surface`** — pMHC surface topology: a height field over the groove with
  hydropathy and charge, following SURFMAP (Schweke 2022) and Protein Surface Topography (Berkut
  2019), plus Manhattan map distance and a hierarchical epitope tree. Makes "featureless" a number:
  validated on all 374 Canonical2026 complexes, `relief`/`peak_to_valley`/`frac_above_ridge` all rise
  with class-I peptide length (`frac_above_ridge` 0.054 for 8-mers → 0.569 for 13-mers; Spearman on
  relief +0.414, p = 5.5e-13), and maps of the same epitope are closer to each other than to a
  different one (P[within < between] = 0.917, p = 1.6e-94). Against the literature: the epitopes
  named as bulged rank **2nd, 5th and 8th of 230** on `frac_above_ridge`, and both named featureless
  ones sit at exactly **0.000** — no peptide surface clears the helix crest at all.
- **`tcren.rotamers`** — rotamer-averaged contact probabilities. Under a deliberately wrong χ1 the
  hard 5 Å contact set keeps a Jaccard of 0.66; the averaged map keeps 0.95, and mean |ΔΦ| falls
  from 0.524 to 0.054 against energies whose own magnitude is 0.4–2.2.
- **Peptide position** — `peptide_positions`, `position_weights`, `position_profile`,
  `central_strain`, answering the review's "contacts in the centre of the peptide matter more than
  at the edges". The position was always on the contact table and `refine.anchors` always predicted
  anchors; they were never joined. The per-position Φ sums exactly to the total, and the class-II
  register resolves 4ozg's gliadin core to the published P1/P4/P6/P9.
- **`derive_tcren_by_type`** — the type-conditioned potential, with the occupancy report that
  answers whether to trust it. On 8002 typed contacts from 370 structures **no type reaches 5% cell
  occupancy**, so the review's own sparsity concern is confirmed and the filter is the usable half.

- **`_relax.relax_interface` + `tcren.dynamics`** — flexible-backbone Metropolis MC of the peptide,
  reporting conformational stability (`rmsf`, `drift`, `energy_gap`) rather than a better pose, with
  the intra-peptide term as a switch. Built to test Sewell's hypothesis that intra-peptide
  interactions stabilise the productive conformation and explain where an additive contact model
  fails. On 2102 CPL structures across seven clones: **stability beats the contact energy in 4/4
  clones where the contact model fails and 0/3 where it works** (ila1 0.348 → 0.862; sb27 0.570 →
  0.934). Removing the intra-peptide term lets best binders wander further (Δrmsf +0.021 Å, 4.4σ) but
  not worst binders (+0.002 Å); best vs worst p = 0.042. `scripts/sewell_stability.py`.
- **`_relax.repack` + `tcren.repack` + `tcren refine --repack`** — the native side-chain packer.
  On the same wrong-rotamer input and the same atom set, it recovers peptide side-chain RMSD from
  **4.131 Å to 2.364 Å in 6 ms**, where OpenMM's anchor-restrained minimisation returns 4.133 Å
  (unchanged) in 3103 ms — a local minimiser cannot cross a torsional barrier. 8/8 structures
  improved (3.93 → 1.66 Å median). The kernel reproduces the Python prototype's per-residue energy
  exactly (0.0, not a tolerance), and a crystal in gives the crystal back.
- **`notebooks/surface_topology.py`** (marimo) + its rendered page in the docs gallery: elevation,
  charge and hydropathy maps with the featureless-vs-bulged epitope comparison.
- `tcren score --soft` scores over rotamer-averaged contact probabilities.

### Changed
- **Contact typing rewritten** (`scheme="v2"`, default). `other` falls from **72.3% to 13.9%** of
  TCR:peptide contacts: `polar`/`vdw`/`cation_pi`/`stacking` classes added, apolarity decided per
  atom rather than per residue (which had excluded Tyr entirely), the H-bond cutoff widened to 3.9 Å
  with donor/acceptor typing, and `stacking.ring_stacking` finally joined. A contact may carry
  several types. The old scheme is kept verbatim as `scheme="v1"`, which `recognition.py` pins,
  because the frozen classifiers were fitted on its counts. Two things the typing then measured:
  interface ring stacks are genuinely rare (**1 across 10 crystals** — the 22–39 stacks per structure
  are core packing, not recognition), and `--drop-untyped` removes **~17%** of TCR:peptide pairs on
  crystals, every one of them `vdw`/`other`.
- **Hydrogens are filtered from contacts.** `all_atom_contacts` documented heavy atoms and did not
  enforce it, so the same complex scored differently depending only on whether the depositor modelled
  H (5jhd gained 7 of 28 TCR:peptide contacts, −58.5% on F_tcr_pep). Breaks legacy-oracle parity on
  the two H-bearing fixtures, which the regression test now records as a subset relation.
- `score_peptides` / `_interface_energy` take an explicit per-contact `weights` array; `tcren score`
  gains `--drop-untyped` and `--position-weights`. All default to the previous behaviour.

### Fixed
- `binder.noise.is_real_interface` compared a **signed** incident angle against an unsigned envelope
  whose floor is 0°, rejecting every downward tilt — including class-I crystals (5xot, −9.9°).
- `recognition_features` swallowed a failed docking geometry in a bare `except Exception: pass`; it
  now warns and says which six features are NaN.

## [2.8.0] — 2026-08-11

### Added
- **The one-body / pair split of a contact potential.** A contact energy is not purely an
  interaction: burying a residue against *any* partner costs something that depends on that residue
  alone, and only the remainder is chemistry between two identities. Summing a contact matrix over
  pairs and calling the total an interaction credits it with an additive component a per-position
  model already has.
  - `Potential.decompose()` → `e(a, b) = mean + H(a) + H(b) + J(a, b)`, by double-centring. Exact and
    unique; `J` has zero marginals and is **the only part a sum over positions cannot express**.
  - `Potential.hydrophobicity_fit()` → `C0 + C1(q_a + q_b) + C2 q_a q_b`, for a matrix that ships no
    solvent reference. Li, Tang & Wingreen (*Phys Rev Lett* 79:765, 1997) showed the MJ matrix is
    nearly rank one, so the one-body term can be recovered from the matrix itself. R² = 0.85 on the
    bundled `mj`, 0.98 on `mj1996`. The consequence is worth stating: where a potential has that
    shape, the interaction term is only `C2·q_a·q_b`, so it **cannot prefer one pair of side chains
    over another of equal hydrophobicity**.
  - Both refuse a directed potential — TCRen is TCR→peptide and must not be split this way.
- **Two Miyazawa–Jernigan reference tables, with recorded provenance.**
  - `mj1996()` — the 1996 Table 3 contact energies (`e_ij`, RT units), transcribed from AAindex
    `MIYS960101` and cross-checked against a second independent copy (same alphabet order, same
    Ala–Ala, same range). The companion repulsive packing-density term is deliberately excluded: it
    is a function of coordination number, not of a residue pair.
  - `mj_partition_energy()` — the 1985 effective partition energies (AAindex `MIYS850101`), the
    one-body term a pairwise matrix cannot supply. Cross-check: it correlates at r = +0.98 with the
    hydrophobicity axis `hydrophobicity_fit()` recovers from `mj1996`, which was transcribed from an
    unrelated source — a transcription slip would break that agreement rather than hide in it.
  - **What this settles about the bundled `mj` matrix**: it is *not* MJ 1996 Table 3. Table 3 is
    attractive everywhere (Ala–Ala −2.72, range −7.37 to −0.12); the bundled one takes both signs
    (Ala–Ala −0.12, range −1.19 to +0.76). They correlate at r = 0.89, but the bundled matrix is not
    Table 3's double-centred pair part either (r = 0.51), so **what it is remains open**. The file is
    left untouched — every score in the package is built on it — and `mj()` now says so.
- **`ring_stacking()`** (`tcren.stacking`) — a contact potential scores a residue pair by identity
  alone, so it treats two rings face to face at 3.5 Å exactly like the same two residues brushing
  past edge-on. This measures the difference from coordinates: centroid separation, interplanar
  angle, and the split of the separation into the gap between the planes (`vertical`) and the
  sideways slide (`lateral`) — enough to separate a parallel-displaced stack from an edge-to-face
  contact from two rings that merely happen to be nearby. Proline is included despite not being
  aromatic: its pyrrolidine ring packs face-on against aromatics through CH–π contacts, and omitting
  it would miss the interaction the module exists to measure. **Nothing here returns an energy** — it
  says the rings are or are not arranged the way a stack is.
- **`SOURCES`** records the origin of every bundled potential table: upstream accession or paper, the
  transcription check, and whether each value is measured, published or derived. Two of the five had
  no recorded origin before; one of those (`MJ_Keskin_potentials.csv`) still does not, and is now
  labelled unresolved rather than left to be assumed.

## [2.7.0] — 2026-08-11

### Added
- **The intra-peptide term.** Every energy in the package sums over contacts between two *different*
  chains, so a peptide held in its bound conformation by its own side chains scores the same as one
  that is not. That omission was one unconditional line in `all_atom_contacts` — invisible rather
  than deliberate. It is now a term you can switch on:
  - `all_atom_contacts(..., scope=)` — `"inter"` (default, unchanged), `"intra"`, or `"all"`.
    `peptide_internal_contacts()` wraps the intra case with the sequence-separation filter such a
    term needs: neighbours touch because they are bonded, not because the peptide folded that way.
  - `intra_peptide_energy(contact_map, potential, peptide=None)` — the energy itself, for the
    structure's own peptide or a candidate threaded onto its pose. The potential is **symmetrised**
    (`(F + Fᵀ)/2`): an intra-chain pair has no `from`/`to` orientation, and which residue lands on
    which side is an artefact of the contact table's canonical ordering, not chemistry. It defaults
    to MJ, since TCRen is derived from TCR↔peptide contacts and says nothing about a chain's
    contacts with itself.
  - `score_peptides(..., intra_weight=w, intra_potential=)` and `tcren score --intra-weight` —
    `score = Φ_interface + w · E_intra`, with the candidate threaded onto both sides of each
    internal pair.
  - `pipeline.run(..., intra_weight=w)` and `tcren scoring --intra-weight` — reports the energy raw
    as `F_pep_int` and folds `w ·` it into `F_total`, so the term and the weight given to it stay
    separable in the output. Its potential is overridable via `potentials={"peptide_internal": …}`.
  - `tcren recognize --full` emits `F_pep_int` and `n_pep_int`. Both are catalogued in `DESCRIPTORS`
    with `involves_tcr=False` — the peptide's contacts with itself are a property of the epitope's
    bound conformation, shared by every TCR that reads it, so `descriptors(tcr_only=True)` excludes
    them like the other cohort-identity columns.

  **The term is sparse, by construction.** At the 5 Å / `|i−j| ≥ 3` defaults — the same contact
  definition the rest of the package uses — a canonical extended class-I 9-mer makes **zero to two**
  internal contacts: over the 17 deposited complexes in `tests/assets/pdb` the totals are 18
  contacts at `|i−j| ≥ 3` against 134 at `|i−j| ≥ 2`, and that sevenfold jump is entirely `i`/`i+2`
  pairs of an extended chain — covalent geometry, not folding, which is what the separation floor is
  for. So the term separates candidates only where the peptide is genuinely bulged or packed against
  itself, which is the case the interface sum cannot see at all.

  **Everything above is off by default and changes nothing when it is.** `scope="inter"` is pinned
  byte-identical to the previous output on a deposited structure, `ContactMap.from_structure(...,
  peptide_internal=True)` stores the internal pairs *beside* `contacts` rather than in it, and
  `intra_weight=0.0` computes nothing.

## [2.6.0] — 2026-08-09

### Fixed
- **An installed tcren can orient again.** `fetch-data` downloads `Canonical2026`'s 374 structures
  but never its `orient_metadata.json`, which is not on the Hub — it was git-tracked in the repo
  `data/`, so `superimpose` found it in a checkout and nowhere else. Every user path that orients a
  new complex therefore failed off a maintainer's machine: `run_pipeline` raised, `tcren scoring`
  hid the same `FileNotFoundError` in an `error` column, and `superimpose`/`shuffle` reported
  success over empty output. The file now ships in the wheel (`src/tcren/data/`) and
  `_metadata_path` falls back to it **for the shipped database only**, so a user-supplied `--db` can
  never be described by someone else's metadata; the `FileNotFoundError` names the command *and* the
  library call that build a metadata file.
- **Failures are no longer silent.** `tcren scoring` prints the first error to stderr instead of
  only counting it, `run_superimpose` raises when *every* input failed rather than reporting
  "0/N" and exiting 0, and `run_shuffle` raises on an input directory it could not parse instead of
  writing zero decoys.
- **`substitute_peptide` keeps the peptide's region markup**, re-pointed at the new residues. Without
  it the contact map's `pos.from`/`pos.to` are null and `score_peptides` died on
  `int(None)` — the reason scoring a substituted peptide crashed. `score_peptides` now names the
  missing markup instead of raising a `TypeError` out of numpy.
- **`binder_score` names the missing descriptor** and how to build the input; the cohort column
  errors give the library call (`recognition_table(items, full=True)`) beside the CLI one.
- **`Structure`, `Chain` and `ContactMap` print a summary.** The dataclass repr expanded every atom
  and its coordinate — 474,504 characters for one complex — which floods a notebook cell and makes
  any error message that interpolates a structure unreadable.

### Changed
- **`tcren orient` writes `<out>/orient_metadata.json` by default** (was `orient_metadata.csv` in the
  working directory), so a database built by `orient` describes itself in the format `superimpose`
  reads. `--metadata` still takes a path, and a `.csv` suffix still writes CSV.
- **`binding_mode`'s default `contact` is 8.0 Å, up from 5.0.** The cutoff is a Cα–Cα distance, not
  the closest-heavy-atom 5 Å of `contacts`/`score`, so the old default made almost no contacts and
  returned `None` on real complexes — 8 Å is the reference proxy the docstring already cited.
- **`annotate_batch`'s `arda` argument is optional**, resolved lazily like every single-structure
  annotation call. Passing an instance still reuses one mmseqs handle across a batch.

## [2.5.0] — 2026-07-28

### Fixed
- **A failed PyMOL scene no longer returns the previous picture.** PyMOL exits **0** when a script
  raises — the traceback is merely printed — so `check=True` never fired and the only guard was
  "did a file appear". Re-rendering an edited-but-broken scene to a path that already held a good
  render therefore returned that old image and reported success, which is how a figure silently
  stops tracking the data it claims to show. `_run` now scans the output for a traceback and
  raises, and `render()` clears the target first so the existence check means something. Pinned by
  a test over four ordinary breakages (bad path, misspelled command, `NameError`, explicit raise).
- **Pillow is declared.** `render()` defaults to drawing the gizmo, which composites through
  Pillow — and Pillow was in no dependency group, satisfied only transitively via matplotlib in
  `[viz]`. A plain `pip install tcren` followed by the README's own figure example crashed on
  `from PIL import Image`. It is now an explicit `[viz]`/`[marimo]` dependency, and the import
  failure names the extra to install.
- **The sdist is lean again**: 2.73 MB → **1.24 MB**, 153 → 113 entries. `/appendix` (the LaTeX
  derivation and its PDFs — 2.4 MB, half the payload, and inert: the wheel builds without it) is
  excluded, as are `.claude/` and `.DS_Store`. Those last two are untracked and CI builds from a
  clean checkout, so no published artifact ever carried them — but scikit-build-core does not read
  git's *global* excludesfile, so a local `uv build` packaged them, and a local build should match
  what CI ships.

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
