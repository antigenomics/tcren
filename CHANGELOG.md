# Changelog

All notable changes to `tcren` are recorded here. Format follows
[Keep a Changelog](https://keepachangelog.com/); versions follow semantic versioning.

## [2.28.0] — 2026-09-01

**Every out-of-fold-fitted read-out is removed.** The author's ruling: leave-one-epitope-out
fitting was the mechanism behind the defects in the discarded `P_native`, the current
implementations are not trusted, and the analysis will be rewritten rather than patched. Nothing
tcren returns is now fitted against a binding label.

### Removed
- **`reliability.p_binder` and `reliability.available_links`**, with the `calibration` section of
  `data/reliability_moments.json`. The Platt links were fitted out of fold — leave-one-epitope-out
  on the 22-cohort VDJdb panel, within-epitope 5-fold on TCRvdb — and shipped as fold means.
- **`reliability.correct_confidence` and `reliability.available_corrections`**, with the
  `corrections` section of the same file and `CORRECTION_VALIDATED_ON`. This was the one shipped
  read-out that read a label: four coefficients, fitted out of fold and frozen.
- **The `tcren diagnose` command**, which existed only to run that correction.
- The `p_binder` column from `tcren recognize` and `tcren assess`, and `assess`'s `--link` /
  `--list-links` options. `assess` keeps `--list-bands`.

### Unchanged
`Q`, `T`, `\Pi` and `S` are untouched, as are `inversion_flag`, `screening_yield` and `af_band`.
The `blocks` and `phi` moments are native-crystal spreads, not fits, and `af_bands` is a quantile
binning of the generator's confidence carrying the observed non-binder fraction with a Wilson
interval — a measurement of the benchmark, not a model of it. `moments()` now asserts in its own
test that `calibration` and `corrections` are absent, so neither can return quietly.

The removed read-outs and the numbers they produced are recorded in the manuscript repository's
`LEGACY.md`.

## [2.27.0] — 2026-09-01

**The composite score is `S`.** It had accumulated three names for one quantity: `S_free` in the
code and the emitted column, `S_nat` in the frozen calibration keys, and `S_{\mathrm{free}}` in the
rendered maths. The `free` qualifier only ever meant *not `P_native`*, and `P_native` was discarded
in 2.26.0, so the qualifier now distinguishes the score from nothing.

The three blocks keep their symbols. `Q` (placement and interface), `T` (topology) and `\Pi` (the
partition-function-referenced energy) are unchanged, as is the construction: each is a fit-free
directional score over the 374 Native2026 crystals, divided by that block's native spread.

### Changed
- `reliability.s_free` → **`reliability.s_score`**, parallel to `cohort.q_score` and
  `reliability.t_score`, which are the same construction over other descriptors.
- The `tcren recognize` output column `S_free` → **`S`**.
- Frozen calibration link names: `<set>|S_nat` → **`<set>|S`**, and likewise inside the composite
  link names `z(ipTM)+z(pLDDT)+z(S)` and `min rank%(ipTM, S)`. The coefficients are untouched; only
  the keys are renamed, and `data/reliability_moments.json` is renamed in place.
- `correct_confidence`'s returned coefficient `b_s_free` → **`b_s`**, and `af_band`'s
  `s_free_roc_in_band` → **`s_roc_in_band`**.

### Removed
- **The old names, outright — there is no alias and no deprecation shim.** A caller on `s_free` or
  a table carrying an `S_free` column fails loudly rather than reading a renamed quantity by
  accident. The retired names, what each was, and the results each produced are recorded in the
  manuscript repository's `LEGACY.md`.

## [2.26.0] — 2026-09-01

**Every fitted composite is gone, `F` is `Phi`, and a generated table now says which catalogue
produced it.** The author's ruling: discard, do not deprecate, and rebuild from a recompute.

### Removed
- **`P_native` and `cohort.p_native`**, with `P_NATIVE_CHANNELS` / `_POOL` / `_ORIENT` /
  `_FEATURES` / `_BANNED` and the per-channel posteriors `G` / `T` / `E`. It refitted a latent class
  on every call, raised when a cohort had fewer rows than features, and its value depended on which
  rows the fit was anchored on — none of which survives contact with a user holding one model.
  `tcren footprint --score` now emits the fit-free `reliability.t_score` instead.
- **The v1 score block** — `p_bind`, `p_forced`, `q_bind`, `s_strain` as emitted columns, and the
  `--scores` and `--features-only` flags on `tcren recognize`. `cohort.q_score` and
  `cohort.strain_z` remain as functions; what went is the frozen-coefficient layer around them.
- **`FORCED_POSE_MODEL` and `forced_pose_score`**; **`tcren.binder.binder_score` / `BINDER_MODEL`
  and the `tcren binder` command**. Their coefficients were frozen against training sets that no
  longer exist, which made them the one part of the package a reader could not reproduce.
- **`p_real` / `p_real_bn`**, with `real_probability`, `frozen_recognizers`, `encode_features`,
  `GaussianBNClassifier`, `BayesianLogisticRecognizer`, `_hill_climb` and the shipped weights
  `shuffle_logistic.json.gz` / `shuffle_bn.json.gz`.
- The `score` family and the `score` invariance class. **The catalogue is descriptors only.**

### Changed
- **`F_*` → `Phi_*`, `dF_*` → `dPhi_*`, `F_TERMS` → `PHI_TERMS`, `f_score` → `phi_score`.** The
  symbol in the manuscript is Φ; the code now uses it too. `d` in `dPhi` is the **reference
  difference** ΔΦ = Φ(sequence) − Φ(reference), not a derivative — the docs say so at every site.
- **`Phi_total` is the commensurate sum.** The three interfaces are scored with different
  potentials, and those are Boltzmann-inverted from different contact statistics, so their matrices
  are not on one scale (sd: TCRen2 0.4880, MJ 0.3270, **Keskin 1.3181**). An unweighted
  Φ_TP + Φ_TM + Φ_PM is 2.70× more sensitive to a presentation contact than to a recognition one
  when Keskin scores presentation. Each term now enters divided by the standard deviation of that
  interface energy over the **374 Native2026 crystals** (1.6390 / 1.8697 / 4.3013), read from
  `reliability_moments.json`; an untabulated potential falls back to `Potential.scale()`, the sd of
  its own matrix. No label and no fit enters either.
- **`tcren.__version__` reads pyproject when imported from a source checkout.** An editable
  install writes its dist-info once, so `importlib.metadata` had been reporting **2.13.0** against a
  source tree at 2.25.0 — twelve releases stale, and every provenance record in between wrong.

### Added
- **`ddg.smoothed_reference`** — ΔΦ against the **free energy of the residue background** instead of
  one arbitrary poly-alanine sequence, plus `varPhi`, the variance of the same local field.
  Both interfaces containing the varying chain are summed and the third cancels identically:
  varying the peptide cannot change Φ_TCR:MHC, varying the TCR cannot change Φ_pep:MHC — verified,
  exactly 0.000000, on both the frozen-map and structural paths. β → 0 recovers the equimolar mean
  field a combinatorial library realises; β → ∞ the distance from the best residue at each position.
  Emitted as `dPhi_{pep,tcr,tra,trb}_soft` and `varPhi_{pep,tcr}_soft`. The receptor direction is
  **split by chain** so a linear model can form the TRB − TRA contrast rather than being handed it.
- **`tcren.provenance`** — every table `tcren features` writes gets a `<name>.provenance.json`
  recording the version, the invocation and a **SHA-256 digest of the descriptor catalogue** (names,
  families, invariance classes, units and definitions). `tcren recognize --features` checks it and
  **refuses** a table written under a different catalogue, naming the command that would regenerate
  it. This is the guard against quoting a number from a stale table that a fresh run would not
  reproduce; it fired on its first real use.
- **`recognition.STATUS`** — the descriptors that need a second look, and why: `pitch` (reads the
  generator's confidence, never a feature), eleven columns fixed by an exact algebraic identity over
  others, and five computed without the receptor, which carry cohort identity. Rendered as its own
  table in `docs/descriptor_table.rst`.
- **`Potential.scale()` / `.offset()`** — the sd and mean over a potential's defined pairs.

### Fixed
- **macOS AppleDouble sidecars are no longer parsed as structures.** Tarring a structure set on HFS
  writes `._x.pdb` beside every `x.pdb`, carrying a binary resource fork under the same extension;
  `is_structure_file` accepted them and the parser died on a decode error several frames from the
  cause. The corpus recompute lost two sets of 2,000+ models to exactly this.

## [2.25.0] — 2026-09-01

**The alanine scan moves atoms, on both sides of the interface.** `tcren ddg --alanine-scan` had
always taken the virtual path — the CLI never passed a structure — and the structural path it could
have taken was measuring against the wrong baseline.

### Added
- **`ddg.tcr_alanine_scan`** — the receptor-side scan, one row per *contacted* CDR residue, each
  truncated to alanine in 3D and rescored on the rebuilt contact map. There is no virtual variant:
  truncating a receptor side chain without moving atoms would leave every contact it made in place,
  which is the failure the function exists to fix. Scored over `tcr_peptide` alone, because a
  receptor substitution cannot change the peptide:MHC energy.
- **`ddg.tcr_alanine_reference`** — the per-loop aggregates `dPhi_ala_cdr12`, `dPhi_ala_cdr3a`,
  `dPhi_ala_cdr3b` and their total `dPhi_ala_tcr`. Each is the **sum of per-residue** ΔΔGs of that
  loop, which partitions the total exactly. Deliberately not the energy of mutating a whole loop in
  one pass: those differ once atoms move, because truncating every side chain at once loses contacts
  each residue alone retains.
- **`refine.substitute.substitute_residues`** — the general 3D primitive. Any chain, any set of
  residues, keyed `(chain_id, seq_index) -> one-letter`; every other residue and chain comes back
  byte-identical, and an empty mapping is the identity. `substitute_peptide` now shares its residue
  rewriter and chain rebuilder.
- `tcren ddg` gains `--side peptide|tcr|both` and `--virtual`, and now passes the structure, so the
  scan is 3D by default.
- **`recognition.INVARIANCE`** — what each of the 117 descriptors is invariant under, so
  *geometry* and *topology* are different questions rather than two names for the contact set.
  Five classes: `geometric` (a length, area, angle or direction cosine — preserved by isometry,
  **this is the docking**), `topological` (Betti numbers, the Euler characteristic and their
  normalized forms — preserved by continuous deformation, **this is the interface surface**),
  `compositional` (a count over the *labelled* contact set, or a share, entropy or Hill number of
  such counts), `energetic` and `categorical`. `descriptors(invariance=...)` filters on it and
  composes with `family` and `tcr_only`.
  Two things it makes visible. **The topology family is mostly compositional** — 19 of its 29
  columns read the labelling rather than the shape, and only 8 are topological invariants. And
  **neither shipped block is what its name says**: `Q`, called interface geometry, carries one
  continuous quantity of four (`burial`, an area) and no angle, distance or height at all; `T`,
  called the shape block, carries one topological invariant of five. Seven of the nine terms
  across both are compositional, which is why the two correlate more than their names suggest.
  `h0_pers_ent` is filed `geometric`, not `topological`: the H0 barcode's bar lengths *are* the
  minimum spanning tree's edge lengths in angstroms, and persistent homology is a metric
  construction.
- **Peptide coverage, normalized so class I and class II compare** —
  `footprint.PEPTIDE_COVERAGE_FEATURES`: `pep_free_frac`, `pep_cov_frac`, `pep_cov_even`,
  `pep_cov_d2n`, `pep_cov_centre`, `pep_cov_spread`. Every column divides by the peptide's own
  length, and **no position index, band or cutoff enters any definition**. Anchors are found from
  the coordinates: position *i*'s accessibility `a_i = n_TCR_i / (n_TCR_i + n_MHC_i)` is its share
  of contacts facing the receptor rather than the groove, so a buried anchor is discounted however
  many of its atoms sit within the contact cutoff of a CDR loop. A binary anchor test cannot do
  this — at 5 A **every** residue of a class I nonamer contacts the MHC, so the non-anchor set
  comes out empty.
  Measured, it recovers the canonical registers from contacts alone: the three class I complexes
  (8-, 9- and 10-mer) bury **both termini** and peak in the middle, while class II `4ozg` buries a
  **gapped core** — receptor-facing positions interleaved with groove-held ones — and leaves the
  receptor the least peptide of the four (`pep_free_frac` 0.125 against 0.157-0.356).


### Fixed
- **The structural alanine scan read every position against a poly-stub baseline.** It threaded the
  whole peptide through `substitute_peptide`, which truncates *every* residue to backbone + Cβ, so a
  single substitution silently stripped its neighbours' side chains too. On 1ao7 the **native**
  sequence threaded back through it keeps **14 of 29** TCR:peptide contacts and moves the interface
  energy −0.9603 → −1.7177, and that offset appeared at every position — including positions with no
  contacts, which must read exactly 0. The scan now substitutes one residue at a time through
  `substitute_residues`. `substitute_peptide` is unchanged and remains correct for the poly-alanine
  *reference*, where every residue genuinely is mutated.
- `alanine_scan` validates the peptide length against the structure's own peptide chain before
  scoring, so a mismatch names the chain rather than surfacing as "peptide was not scored".

- **`recognition.DETAIL`** — units and a one-line definition for all 124 catalogue entries, and the
  single source the docs table is generated from, so a descriptor cannot reach a feature table
  undocumented. Units come from a closed vocabulary because they are what a transform has to
  respect: a `count` is variance-stabilized by a square root, a `fraction` by the arcsine — the
  classical angular transformation — and an unbounded continuous quantity by neither.
- **`docs/descriptor_table.rst`** — all 117 descriptors by family, with invariance class, units,
  whether the receptor enters the definition, and what each measures. Generated by
  `scripts/gen_descriptor_table.py`; `--check` fails if it has drifted, and a test runs it.
- **`docs/_static/descriptor_families.{dot,svg,pdf}`** — the family-against-invariance graph,
  generated by `scripts/gen_family_graph.py` from the catalogue, Okabe-Ito so it survives greyscale
  and colour blindness. It draws the mismatch the classification exposes: the thickest edge out of
  `topology` runs to `compositional`, not to `topological`.

### Changed
- `docs/features.rst` gains *Two views of the same descriptors* and *The alanine scan, on both
  sides*; SKILL.md gains the matching section.
  The SKILL.md claim that `reference_delta` equals the sum of `alanine_scan().ddG` is now qualified:
  it holds on the virtual path, where a fixed contact map makes the energy additive over contacts,
  and not once atoms move.

## [2.24.0] — 2026-08-29

**A second presentation potential, and the reference-state trap that made the first comparison
unsound.** The bundled `MJ_Keskin_potentials.csv` holds two matrices in *different reference states*
— `mj()` is mixed-sign with mean −0.079, `keskin()` is negative everywhere with mean −3.547 — so
swapping one for the other changes the reference state as well as the derivation. The like-for-like
pairs are `mj` ↔ `betancourt` and `keskin` ↔ `mj1996`.

### Added
- `potential.betancourt()` — the Betancourt–Thirumalai `B` matrix (RT units), Miyazawa–Jernigan
  re-referenced with **Thr as the reference solvent**, so every Thr entry is exactly `0.00`. Parsed
  from AAindex3 accession `BETM990101`, never retyped; the Thr row, symmetry and cell count are
  asserted at build time and again in `tests/unit/test_potential.py`. Reference: Betancourt MR,
  Thirumalai D. *Protein Sci.* 1999;8(2):361–369. doi:10.1110/ps.8.2.361.
- `data/BT1999_contact_energies.csv`, with a `SOURCES` entry carrying the re-fetch URL.

- **`potential.aaindex` — the whole of AAindex3, bundled verbatim** (`data/aaindex3.txt`, 47 records,
  80 kB), so adding a published matrix to a comparison costs a string rather than a transcription.
  `catalogue()` lists every entry with its kind, symmetry, mean and citation fields; `aaindex(acc)`
  returns one as a `Potential` and **refuses** the 2 contact-count and 3 side-chain-distance tables,
  because scoring a contact map with a count table is a silent category error; `entry(acc)` returns
  them deliberately. 42 of the 47 are usable energies.
- `potential.identify(pot)` — compares a matrix cell by cell against every AAindex3 entry and
  returns the accessions ordered by maximum absolute difference.
- **`Potential.components()`** — the exact `mean + H(a) + H(b) + J(a,b)` split of `decompose()`,
  returned as three *scorable* potentials. Because an interface score is a sum over contacts, the
  split carries through to it: `size` sums to `mean x (contact count)`, `comp` to a degree-weighted
  composition term, `pair` to the interaction proper. Scoring each part in turn says whether a
  result is reading interface *area*, *composition* or *chemistry* -- a matrix with no positive
  entries has a large negative mean, so its interface sum is dominated by the contact count.

### Fixed
- **The bundled `mj()` and `keskin()` matrices are identified.** `mj()` is AAindex3 `MIYS990106`,
  Miyazawa--Jernigan **1999** (Proteins 34:49-68) -- not 1985 and not 1996 -- and `keskin()` is
  `KESO980101`, the solvent-mediated interfacial form. Both match 400 of 400 cells exactly, with
  runners-up off by 0.65 and 2.77, so the identifications are unique. The "upstream table
  unrecorded / do not cite" warning that had stood on `mj()` since 2026-08-11 is withdrawn; the
  files are unchanged and can now be cited.
- **Five cells of `mj1996()` disagree with AAindex3 `MIYS960101`** by 0.04-0.28 (M-V, D-M, E-M, H-R,
  A-P; the two matrices correlate at 0.99978). Ours was hand-transcribed and four of the five
  involve Met, Arg or His, so it reads as a transcription slip on our side. Left byte-for-byte
  untouched under the same rule as the MJ/Keskin file, **pinned by a test** so it stays visible, and
  `aaindex("MIYS960101")` is the curated alternative. It feeds no reported number.

### Changed
- `keskin()` gains a docstring stating its reference state and its citation; both it and the
  `SOURCES` entry now name the like-for-like partner, so a future swap cannot cross reference
  states silently.
- `docs/potentials.rst` gains three sections: the AAindex3 resource, identifying an unlabelled
  matrix, and the component split.

## [2.23.0] — 2026-08-29

**The gauge that lets a pinned model reproduce a referenced score.** `centred_potential` double-centres
by design, which is right for *ranking* potentials against each other and wrong for *reproducing* one.

### Added

- **`centred_potential(..., centre=False)` and `fit_potts(..., pin_centred=False)`.** Defaults are
  unchanged, so no existing fit moves. `reference_delta` is a difference of one-body sums, so a
  double-centred pin re-injects a burial-scaled composition term `n_i·c(a)` — the position's contact
  count times the potential's partner-residue column mean — and the referenced-energy identity
  fails. Measured on TCRen2, that column mean has s.d. 0.0668 and runs −0.212 to +0.027 over the 19
  residues the potential observes, while `n_i` runs 1 to 54. Pinned uncentred the coupling **is**
  the potential (max |Δ| = 0 against the raw sign-flipped matrix), so any linear read-out of the
  field reduces to the potential's own score up to the fitted scale.
- **`docs/potts.rst` gains "Three limits of one free energy".** `∂Φ/∂η_a = p_a`, so
  `Φ = Σ softplus(η)` is an interaction sum weighted by contact probability and a fixed contact map
  is the `p ∈ {0,1}` case. Hard contact reproduces `reference_delta`; smoothed (`p` free) fits a
  plastic interface; saturated fits a frozen one, where a fixed map is already exact. Which limit an
  interface is in is measurable: over four deposited 100 ns trajectories, **36 of the 38 engaged
  peptide positions (94.7%) reach a maximum groove-pair contact frequency above 0.98, against 14 of
  38 (36.8%) on the receptor side**. (Corrected after release: the figure first published here,
  "37 of 39", came from the raw fixture rather than from the producer, which restricts to the pairs
  inside the availability radius. The producer's population is the reproducible one.)

## [2.22.0] — 2026-08-29

**The presentation interface gets its own Hamiltonian.** `available_pairs` enumerated
receptor:partner pairs only, so the groove's grip on the peptide — the term an activation read-out
cannot do without — had no model at all and was scored with Miyazawa-Jernigan throughout.

### Added

- **`available_pairs(..., receptor="mhc", partner="peptide")`** — the peptide:MHC arm. The groove
  residue takes the receptor slot and `region.rec` is a groove region, so `MHC_RECEPTOR_REGIONS`
  replaces the TCR loop set in the `g_region` block. `fit_potts` and `site_codes` take
  `regions=` to carry it; `PottsModel` already stored its own level set.
- **`tcren potts fit --receptor mhc`**. On the 362 Native2026 crystals: 195,674 available pairs,
  23,492 contacts, 507 parameters, pseudo-logLik −17,121.5. The per-position contact profile is the
  textbook groove — peaks of 10, 9 and 10 contacts at P1, P2 and P9 against 1 at the bulge.
- **The groove's chemistry is not the receptor's.** The fitted pair field correlates with the
  shipped TCRen2 at **r = +0.017** and with Miyazawa-Jernigan at **r = +0.042**, so the two
  interfaces are not one field read twice.

### Fixed

Three instances of one defect: an MHC-side interface selects **zero rows** until `annotate_mhc`
splits `chain_type == "MHC"` into `MHCa`/`MHCb`, and nothing said so.

- **`tcren contacts --interface peptide_mhc` wrote a header and no rows** (and `--interface all`
  silently omitted every MHC-side contact). It now annotates, and `peptide_mhc` over Native2026
  goes from 0 to **24,648** contacts.
- **`ContactMap.interface` now emits a `RuntimeWarning`** when the map's MHC chains are still
  unrefined, matching the flag `footprint.cell_counts` has carried for the same condition.
- **`python -m tcren.cli potts …` reported "No such command 'potts'".** The `__main__` guard sat
  two thirds up the module, so the app ran before the nine `potts` subcommands below were
  registered; the `tcren` console script, which imports first, saw them. The guard is now at the
  end of the file.

## [2.21.0] — 2026-08-29

**A peptide could only be scored one interface at a time.** A CPL response-matrix *cell* has summed
both peptide-bearing interfaces since `cpl.response_matrix` existed — TCRen over TCR:peptide plus
Miyazawa-Jernigan over peptide:MHC — because an activation read-out fires only if the peptide is
presented *and* the receptor engages. A whole *peptide* had no such path: `ddg`, `neoantigen_ddg`
and `reference_delta` took a single `interface`, so a combinatorial-library ranking silently scored
the receptor term alone and was blind to a destroyed MHC anchor. On 1ao7 the C-terminal anchor
substitution `LLFGYPVYV -> LLFGYPVYA` reads **ΔΔG = 0.0000** over `tcr_peptide` and **-0.9740** over
the complex.

### Added

- **`interface="complex"`** on `ddg`, `neoantigen_ddg` and `reference_delta`, with a
  `mhc_potential=` argument defaulting to Miyazawa-Jernigan — the same convention
  `cpl.response_matrix` uses, so a library ranking and a matrix cell are now the same score.
  Exercised on real coordinates: the complex equals the sum of its two interfaces to machine
  precision. `weights=` still reaches the receptor channel only, matching `tcr_weights`.
- **`tcren ddg --interface complex`** and **`--mhc-potential`**.

### Fixed

- **`tcren ddg` never annotated the MHC**, so `--interface peptide_mhc` returned `0.0` for every
  mutant — a silent zero, not an error. The command now runs `annotate_mhc`, as `tcren cpl` always
  has. Any earlier peptide:MHC ΔΔG from the CLI is void; the library API was unaffected.

### Note

The two effects are **not separable in a combinatorial library that varies every position**: a
peptide can be inactive because it does not present or because it is not recognised, and the assay
reads only activation. `complex` is the right score for ranking; the per-interface terms are
reported alongside it so the confound is visible rather than absorbed.

## [2.20.0] — 2026-08-29

**The residue, not just the position.** 2.19.0 could say how engaged a peptide position was expected
to be; it could not say what happened when the residue there changed. The machinery was one
composition away and nothing composed it.

### Added
- **`potts.peptide_free_energy` and `tcren potts scan`.** `contact_map --by position` reads how
  engaged a peptide position is expected to be *before any residue identity is scored*; this reads
  what happens when the identity changes. The partner residue enters the one-body field twice, so
  threading a residue through position `i` moves eta at every available pair carrying it and the
  interface free energy moves with it:

      Phi_Potts(x) = log Z0(eta(x)) = sum_s log(1 + exp(eta_s(x)))
      dF_i(a)      = Phi_Potts(x_{i->a}) - mean_b Phi_Potts(x_{i->b})

  Exact and closed form for the coupling-free model, no sampling. `coupled=True` takes the linear
  response about the observed sequence, `d log Z / d eta_s = <sigma_s>`, so one Gibbs pass and a dot
  product per cell. `log Z0` is a sum over independent sites, so `dF` is additive over positions and
  one `L x 20` table scores a single substitution and any whole partner sequence alike. The
  reference is the **equimolar** one, the mean over the twenty residues at a position, which is the
  null a positional-scanning library holds its other positions at.

  **This one IS an energy**, unlike `contact_map`'s frequencies: `log Z0` carries k_B T and belongs
  in an energy block. A partner position carrying two different residues has no sequence to
  substitute into and raises rather than being averaged.

### Fixed
- **`n_contacts` was two different counts under one name, and the second silently won.** The
  footprint pass wrote the CDR-loop contact tally; the Potts pass wrote the available pairs that
  engaged. They are not the same number on the same structure — 1ao7 66 against 29, 1bd2 58 against
  21, 1fo0 42 against 11 — and since `_featurise_families` runs `topology` before `potts`, the
  emitted column meant whichever family the caller asked for, with no error either way.

  This mattered because `reliability.correct_confidence` standardizes `n_contacts` against
  `reliability_moments["blocks"]["n_contacts"]`, which is the **Potts** population. `tcren diagnose`
  on the default `-i placement,interface,topology,energetics` table therefore standardized the
  footprint tally against Potts moments: 1ao7's contact term read z = +7.2 where the Potts count
  gives +1.2, and the corrected probability moved with no warning and no NaN to notice. No published
  number is affected — the benchmark reads the count from `potts_bound_scores.csv` and standardizes
  it against the same population — so this was a library-user-facing defect only.

  Three changes close it. The footprint total is **`n_loop_contacts`** now, beside its existing
  `n_pep_contacts` / `n_mhc_contacts` siblings and the `interface` family's `n_contacts_tp` /
  `n_contacts_tm`. Bare `n_contacts` is **catalogued under `potts`**, where it was previously
  catalogued nowhere and so was dropped by `-i potts` and `-i topology` alike and emitted only when
  a caller happened to ask for `interface` as well. And `correct_confidence` **raises** when handed
  a table carrying `n_contacts` but none of the columns only the Potts pass emits (`neg_energy`,
  `log_z`, `log_lik`, `psi`, `n_sites`, `mu_star`), which is what a pre-fix or topology-only table
  looks like, rather than standardizing the wrong quantity;
  `tcren diagnose` reports that as a parameter error naming the rebuild command.
- **Documentation, audited paragraph by paragraph against the source.** `P_native` was still called
  "the recommended score" in three places while two others correctly said `S_free` is — it is
  cohort-refit, `S_free` is fit-free and defined for one structure. `SKILL.md` claimed anchors are
  never scored rows (`fit_em` pins the caller's own rows and keeps them in the design), that the
  energetics channel needs a leading minus (no shipped channel does), and that `f_score` feeds
  `P_native` (that channel moved to the `potts` family at 2.17.0). The bundled potential names, the
  footprint feature count, and `potts score`'s emitted column list were each wrong.
- **CHANGELOG history.** An `[Unreleased]` block sat between `[2.11.0]` and `[2.9.0]` describing
  work that had already shipped, and there was no `[2.10.0]` section at all.

## [2.19.0] — 2026-08-29

**The contact map, read as a map.** `contact_probabilities` has emitted a per-residue-pair marginal
since 2.16.0, but the grid an experiment measures is coarser, so nothing in the package could be
compared against one.

### Added
- **`potts.contact_map` and `tcren potts map`.** Closes the per-pair marginals onto the two grids a
  caller actually reads. `--by loop` gives one row per (structure, CDR loop, peptide position) —
  the contact-frequency map, which is what a molecular-dynamics trajectory reports as the fraction
  of frames in which any residue of that loop touches that position. `--by position` collapses the
  loops into peptide residue importance: how engaged the model expects each position to be, before
  any residue identity is scored. `--by pair` is the ungrouped passthrough.

  The closure is the Poisson-binomial "at least one", `P(N >= 1) = 1 - prod(1 - p_j)`, accumulated
  in `log(1 - p)` so a twelve-residue loop does not underflow and a saturated pair returns 1 rather
  than `nan`. Columns are `p_any`, `p_expected` (the expected contact count), `n_pairs`,
  `n_observed` and `observed`.

  These are **frequencies, not energies** — dimensionless and in [0, 1]. They carry no k_B T and
  belong to the diagnostic and importance side of the model; `score`'s `neg_energy` is the quantity
  with units, and it is what `S_free`'s Pi block reads.
- `notebooks/potts_contact_map.py`, a marimo app over the released path: the predicted map beside
  the contacts the structure made, and the importance profile under it.

### Changed
- **`contact_probabilities` gained `workers=`**, and `tcren potts score` / `potts contacts` / `potts
  map` gained `--workers`. The library has had process-parallel per-structure mapping since 2.17.0;
  three of its five samplers and both CLI entry points never got it. Bit-identical however the work
  is split, because each structure's numbers are a function of `(seed, pdb.id)` alone.

## [2.18.0] — 2026-08-29

**The generator says it is confident. What should you believe instead?** `af_band` answered how
often a confidence band is wrong; nothing answered what to believe in its place. This does.

### Added
- **`reliability.correct_confidence` and `tcren diagnose`.** One logistic over three terms —
  the generator's confidence, `S_free`, and the observed contact count, the last two in native-sd
  units — returning the corrected probability **and its decomposition**: `p_confidence` (the
  confidence alone, through the same link), `delta_logit` (what the coordinates added, in nats) and
  `p_corrected`. The decomposition is exact by construction, so a caller can see whether a number
  moved because of the generator or because of the structure.

  Coefficients are frozen, fitted out of fold, and rounded to one decimal — rounding costs under
  0.003 macro ROC-AUC, well inside the fold-to-fold spread. **This is the one shipped read-out that
  is not fit-free**: `s_free` takes no label anywhere, this learns four numbers from labels, the
  same standing as `p_binder`'s Platt links.

  Validated where the epitope has structural precedent. Leave-one-epitope-out on the balanced VDJdb
  panel it adds **+0.051** macro ROC-AUC to ipTM and **+0.068** to pLDDT over the 6 cohorts whose
  epitope has a solved complex (n = 284), and *subtracts* about 0.04 over the 16 that do not
  (n = 743) — the same template covariate the receptor benchmarks divide under.
- `reliability.available_corrections()` and `CORRECTION_VALIDATED_ON`.
- `notebooks/confident_negatives.py`, a marimo app over the released path: move the confidence
  slider to the top of the range and watch the corrected probabilities stay spread while the
  confidence-only reading collapses to one number.

### Changed
- **`data/reliability_moments.json` now has a producer**, and gained a `corrections` block. Until
  now the file every shipped score divides by was hand-assembled and committed once with 2.16.0 —
  no script wrote it and no benchmark stage regenerated it. It is now written by the benchmark's
  `bench/scripts/build_reliability_moments.py`, which defaults to a dry run and prints what would
  change, because a frozen constant that moves is a finding rather than something to overwrite.

  Installing it moved the four Potts blocks (`neg_energy`, `log_z`, `log_lik`, `n_contacts`, n 361
  → 362), which is exactly the 2.17.0 AIS seeding fix reaching the frozen file. Effect on `S_free`:
  max shift 0.0059 against the score's own spread of 2.0057 (0.3%), old against new correlating at
  0.99999994. `Q` and `T` are unchanged and verified to reproduce to 1e-10.

## [2.17.0] — 2026-08-28

**Three defects in `p_native`, and a sampler whose seed did not do what it said.** All four are
reproducibility fixes: the same structure now gets the same number whatever else is scored beside
it, and no channel reads a descriptor that is really a generator confidence.

### Fixed
- **`tcren.potts` samplers were seeded once per table, not per structure.** Every entry point in
  `potts/score.py` built one `np.random.default_rng(seed)` *outside* its loop, so each structure's
  AIS run consumed whatever the preceding structures had left: the same PDB scored on its own, in a
  subset, or in a reordered frame returned a different `log Z`. The generator is now derived from
  `(seed, pdb.id)` via `blake2b`, and `_prepare` sorts sites on their own identity
  (`pdb.id, chain.rec, region.rec, pos.rec, pos.par`, which is unique) rather than on arrival
  order, so the colouring is canonical too. Affects `score_sites`, `contact_probabilities`,
  `bound_unbound`, `count_profile`, `connected_correlations` and `sample_maps`. **Values move by
  AIS sampling noise relative to 2.16.0** — that is the defect being removed, not a new one.
- **`P_NATIVE_FEATURES["placement"]` carried `pitch`**, which is
  `orient.docking_angles(s).incident_angle` — the same quantity `pipeline.py` calls `pitch_angle`
  and treats as AlphaFold-confidence leakage rather than interface geometry. Replaced by
  `dock_torsion`. `Q_FEATURES_GEOM` and `T_FEATURES_TOPO` never contained it, so `q_score`,
  `t_score` and `s_free` are unchanged.
- **`p_native`'s energetics channel had never migrated to the Potts energy.** It read `F_tcr_pep`,
  `F_tcr_mhc` and `dF_tcr_pep` — the poly-alanine-referenced contact energies, which are the right
  instrument for ranking *peptides* and are at or below chance for ranking *receptors*. The channel
  now reads `neg_energy`, `log_z` and `log_lik`, and `P_NATIVE_ORIENT["energetics"]` orients on
  `neg_energy` (higher is more favourable) instead of `-F_tcr_pep`.

- **`p_native` silently reversed its own labelling when the orientation feature was missing.**
  A finite mixture is identified only up to a swap of its components, so `orient_by` (defaulting to
  the channel's `P_NATIVE_ORIENT` entry) is the only thing that decides which component is called
  native. When that column was absent from the table or constant across it, the code fell back to
  the first surviving column, whose direction is arbitrary. The energetics migration above turned
  that latent hazard into a measured reversal: a caller passing the old `F_TERMS` explicitly got
  Spearman **−0.63** against `-F_tcr_all` where the same call had read **+0.63**. It now raises,
  naming the column and pointing at `orient_by`. Every in-tree caller already passes `orient_by`
  explicitly, so nothing that worked stops working.

### Added
- **A `potts` descriptor family in `tcren features`.** `-i potts` emits `neg_energy`, `log_z`,
  `log_lik` and `psi` per structure, so the energy the receptor task wants
  comes out of the same command as every other descriptor rather than out of a separate scoring
  pass that has to be joined back on.
- `cohort.P_NATIVE_BANNED`, checked by `_channel_columns`: naming `pitch`, `pitch_angle` or
  `incident_angle` in any channel now raises instead of quietly fitting on it.

### Added
- **`workers=` on `score_sites` and `bound_unbound`: the per-structure loop now runs in processes.**
  `None` (the default) takes every core, `1` runs serially. Structures are split into as many
  contiguous chunks as there are workers, and the prepared arrays cross the process boundary once
  via the pool initializer rather than once per task. Measured on 616 TCRvdb structures:
  `score_sites` **82.6 s -> 19.8 s (4.2x)**, `bound_unbound` **199.3 s -> 35.9 s (5.6x)**, both with
  **max |diff| exactly 0** against the serial run.

  This is only sound because of the seeding fix above: a structure's numbers depend on
  `(seed, pdb.id)` and on its own sites, never on its position in the frame, so splitting the work
  cannot move a value. Threads were measured and are *worse* than serial (0.33x) — the arrays are
  too small for numpy to release the GIL usefully. `contact_probabilities`, `count_profile`,
  `sample_maps` and `connected_correlations` still run serially.

### Removed
- **The shipped `pnative_anchors.csv` and the `anchors="auto"` path.** Nothing read them: no
  benchmark producer, no test and no CLI subcommand passed `anchors="auto"`, and the one anchored
  path in the TCRen2 benchmark builds its own `{row: label}` dict from held-out rows. The file was
  8.3% of the wheel. `anchors=` still takes an explicit `{row_index: 0|1}`, and `anchors=None`
  (unsupervised, oriented by `P_NATIVE_ORIENT`) is unchanged and remains the default.

## [2.16.0] — 2026-08-28

**The two questions a score does not answer.** `S_free` says how native-like an interface looks.
It does not say *why* to distrust a confident one, and it does not say what testing the top of the
ranking would actually return. Both are read off things already computed.

### Added

- `tcren.reliability.inversion_flag` — the forced-pose detector. A generator pushed into a
  confident but wrong pose does not produce a random interface: to seat the chains it selects
  residue pairs it believes are favourable, so the recognition energy comes out **good**. The
  energy inverts under forcing rather than degrading. Measured on a 24-cohort forced-pose panel
  (1,707 structures) the TCR:peptide energy reads macro ROC-AUC 0.4952 and is below 0.5 in **15 of
  24 cohorts**, where ipTM reads 0.6093. The flag is the energy block's native-sd score minus the
  mean of the geometry and topology blocks: large positive is an energy vouching for a structure
  the shape does not. Needs the `neg_energy` term; returns NaN without it rather than scoring.
- `tcren.reliability.screening_yield` — the cut a testing budget implies: how many structures, at
  what score threshold and what percentile, plus the hits a blind test of that slice would return
  under a stated prevalence. Enrichment over random is deliberately **not** returned: it needs
  labels this function does not have, and a NaN would read like a measurement.
- `tcren assess` emits `inversion_flag` whenever the input carries the energy term, and takes its
  triage cut from `screening_yield` rather than an inline rounding.
- `tcren.ddg.ddg` takes `weights=`, the per-contact multiplier `tcren.scoring.score_peptides`
  already accepted, so a substitution can be scored against a contact **probability** —
  `tcren.potts.contact_probabilities`' `p_model`, or a rotamer-averaged occupancy — instead of the
  map's hard 0/1 indicator. It is dropped on a rebuilt mutant map, whose rows no longer align with
  the native's, rather than silently mis-applied.

## [2.15.0] — 2026-08-28

**A score you can put on one structure, and a number you can read as a probability.** `P_native`
refits a latent-class model per call, raises when a cohort has fewer rows than features, and its
published numbers depend on which rows the fit was anchored on. `S_free` has none of those
properties, and `tcren assess` is the command that turns a folder of models into the three things a
caller decides on.

### Added

- `tcren.reliability` — `s_free`, the recommended single-structure binder score:
  `Q/sd_Q + T/sd_T + (Pi - mu)/sd_Pi`, three fit-free directional blocks `z(x)' C^-1 s` over the
  Native2026 crystals, each divided by its own native spread. **The outer transform is one divide,
  not a z**: a block score's native mean is 0 by construction, but its variance is `s' C^-1 s`
  (1.43 for `Q`, 1.61 for `T`, 14.13 for `Pi`), so without the division the energy would carry ten
  times the weight of the geometry. `Pi` is `neg_energy` from `tcren.potts` — the interface energy
  read against the partition function rather than a poly-alanine reference, and the least redundant
  with `Q` of the five ways of spending `-E = log Z + L` (native Pearson +0.33, against +0.75 for
  the contact count).
- `tcren.reliability.t_score` — the topology block, the SHAPE of the contact set free of its size.
  It is the block that survives where the geometry block does not: on the balanced VDJdb panel `T`
  loses 0.06 ROC-AUC when the epitope has no solved complex to template on, against `Q`'s 0.24.
- `tcren.reliability.p_binder` and `af_band` — **ten** frozen out-of-fold calibration links
  (Platt, leave-one-epitope-out on the 22-cohort panel and within-epitope 5-fold on TCRvdb): the
  five scores `S_nat`, `ipTM`, `pLDDT`, `z(ipTM)+z(pLDDT)+z(S_nat)` and `min rank%(ipTM, S_nat)`,
  each under both the `binder_bm|` and `tcrvdb|` benchmarks. **Four** confidence-band tables behind
  the generator diagnostic — ipTM *and* pLDDT, under both benchmarks. `available_links()` and
  `available_bands()` list them; a value outside a band table's range clamps rather than
  extrapolating.
- `tcren.reliability.moments` — the accessor for every frozen constant above.
- `tcren.reliability.reliability_reference` — the four geometry and five topology descriptors plus
  the Potts energies over the 369 complete-case Native2026 crystals, so a single user structure can
  be standardized against the crystal manifold for all three blocks at once.
- **`tcren assess`** — reliability (`S_free`, `p_binder`), ranking within the set (rank, percentile,
  expected precision at a recall budget), and the generator diagnostic (`af_band`,
  `p_nonbinder_af`, and what `S_free` still separates INSIDE that band). On the balanced VDJdb panel
  the top ipTM decile is 26.2% [18.7, 35.5] non-binders and is also where `S_free` reads highest.
- `tcren recognize --features` emits `S_free` and `p_binder` beside the existing columns. When the
  table carries no `neg_energy` the two-block form is emitted and the message says so, rather than
  imputing the energy silently.
- `tcren.cohort.q_score` takes `signs=`, the per-descriptor orientation replacing `1`. It is what
  lets a block carry a term that runs the other way — the topology block's footprint fraction.
- Two shipped data files behind all of the above: `src/tcren/data/reliability_reference.csv`
  (the 369 complete-case Native2026 crystals, nine block descriptors plus the Potts energies) and
  `src/tcren/data/reliability_moments.json` (the block moments, the ten Platt links and the four
  band tables). Neither is refitted at call time, so a score computed today means what it meant
  when the paper was written.
- `tcren.potts.score_sites` and `bound_unbound` now also emit **`neg_energy`**, the sign-corrected
  interface energy. It is the `Pi` block `s_free` reads, and it is the column
  `docs/reliability.rst` and `SKILL.md` tell you to join — before this, only `energy` (the opposite
  sign) was written, so the three-block `S_free` was unreachable from the shipped package and every
  caller silently fell back to the two-block form. `tests/regression/test_cli_smoke.py` now runs
  `potts score` → `assess` and asserts the sign, which is the check that was missing.

### Notes

`P_native` is still emitted and still documented, now as cohort-refit and not the recommended
score. Nothing about it changed; what changed is that there is an alternative defined at n = 1.

`pyproject.toml` goes `2.12.1` → `2.15.0` in one step: 2.13.0 and 2.14.0 are recorded here for the
work they contain but were never uploaded, so no `tcren info` ever reported them.

## [2.14.0] — 2026-08-28

**Bound versus unbound, for the whole interface.** `eta_a` is the free energy between a single
site's two states; the same contrast for the whole interface needs a macrostate, and the contact
count `N(sigma)` defines one. Because `E(empty) = 0` exactly, three readings of that contrast come
out of the model already fitted, and one Gibbs pass serves all of them — every tilt in `N` is an
exponential family, so the tilted expectation is a reweighted average over draws taken at zero
tilt.

### Added

- `tcren.potts.bound_unbound` — `df_empty` = `log(Z - 1)`, the exact two-state contrast against the
  empty configuration; `df_threshold` = `log[P(N >= x)/P(N < x)]`, in which `Z` cancels; and
  `mu_star`, the chemical potential at which the model's mean contact count matches the observed
  one. `nan` outside the sampled support rather than a silent extrapolation.
- `tcren.potts.count_profile` — the pooled free-energy profile `F(N) = -log p(N)` along the contact
  count, beside the observed counts, so a threshold is read off the landscape rather than assumed.
- `tcren.potts.tilt_mean`, `mu_star`, `count_free_energy`, `delta_f_empty`, `delta_f_threshold` —
  the underlying reweighting and histogram primitives.
- `tcren.potts.gibbs` takes an optional `observer` callback, invoked on each kept draw with the
  `(chains, n)` configuration matrix. It exists so that a statistic of *whole configurations* — the
  kind a Lagrange multiplier couples to in Jaynes' construction — can be accumulated during
  sampling without materialising every draw. Default `None`; no behaviour change when unset.

### Notes

- A **linear** tilt in `N` is exactly a constant added to every field, `E - mu N = -(eta + mu).sigma`,
  so the reweighting identity is checked against direct simulation in `tests/unit/test_potts.py`
  rather than assumed. `delta_f_empty` is checked against exact enumeration of all `2^12`
  configurations.
- `df_empty` and `df_threshold` are not two estimates of one quantity. The unbound basin of a
  *docked* pose is astronomically improbable — the model is conditioned on an available set that
  already holds the receptor against the peptide — so no sampler reaches `N = 0` and only the
  `log Z` route gives it.

## [2.13.0] — 2026-08-28

**The contact map itself becomes a model.** Every scoring path in `tcren` reads a contact map that
a structure *has*; `tcren.potts` models the map as a random variable. A **site** is a residue pair
whose Cα atoms lie within 15 Å — a pair that *could* have contacted — and the configuration σ says
which of them did:

```
E(σ) = -Σ_a η_a σ_a - ½ Σ_ab A_ab σ_a σ_b ,   P(σ) = exp(-E(σ))/Z
```

That reference state is the point. A TCRen potential is a Boltzmann inversion conditioned on a
contact *existing*, so a residue that could have reached the peptide and declined contributes
nothing to it. Here the non-event is the observable, which is what lets the one-body fields
separate reach from chemistry.

### Added

- `tcren.potts` — `PottsModel`, `available_pairs`, `fit_potts`, `score_sites`,
  `contact_probabilities`, `sample_maps`, `score_structure` (the one-shot wrapper: enumerate one
  structure's pairs and score them), `connected_correlations` (the two-point generative check) and
  `kernel_table` (the coupling coefficients with cluster-robust s.e., what `potts fit` prints),
  plus the numerics (`irls`, `gauge`, `gibbs`, `ais_log_z`, `exact_log_z`, `colour`,
  `centred_potential`). Docs: `docs/potts.rst`.
- `tcren potts fit` / `score` / `contacts` — fit a model from structures, get each structure's
  energy, partition function and likelihood, and read out per-residue-pair contact probabilities.
  `--partner peptide|mhc|both`, `--coupling-matrix`, `--balance`, `--no-couplings`.
- Two bundled models, shipped as `src/tcren/data/potts_tcr_peptide.json` and
  `potts_tcr_mhc.json`: `potts_tcr_peptide` (the default; 362 αβ Native2026 crystals, 64,622 sites,
  7,865 contacts) and `potts_tcr_mhc` (239,093 sites, 15,451 contacts). Both reproduce from the
  CLI: `tcren potts fit -s data/Native2026 --balance both -o …`.

### What it measures

- **Contacts are strongly dependent, and the sign flips off-axis.** On the crystals every axial
  coupling is positive and every off-axis one negative: `K(+1,0) = +0.792 ± 0.057`,
  `K(0,+1) = +0.656 ± 0.053` against `K(+1,+1) = -0.816 ± 0.064` and `K(+1,-1) = -0.812 ± 0.066`
  (log-odds per contacting neighbour, s.e. clustered on the structure). A made contact recruits its
  own sequence neighbours onto the *same* partner residue and suppresses the diagonal one. The
  couplings buy +505.7 nats of pseudo-log-likelihood for 18 parameters. Unconditionally the
  diagonal offsets look *positive* — the negative sign appears only once the axial terms and the
  Cα distance profile are held fixed.
- **P(contact) ∝ TCRen, with a scale.** Fixing J to one coefficient on the double-centred TCRen2
  matrix — 1 parameter against 400 — gives `β = +1.131 ± 0.062` and costs 103.9 nats. The shipped
  potential is already at very nearly the right temperature on the log-odds axis.
- **Which potential belongs on which interface.** `--coupling-matrix` gives every candidate an
  identical parameter count and design, and the double-centring means none of them can win on
  composition. The ranking **inverts**: TCRen2 beats MJ by 103.3 nats on TCR:peptide; MJ beats
  TCRen2 by 35.5 nats on the TCR:MHC groove, where TCRen2's scale falls 5.4-fold (+1.131 → +0.209)
  while MJ's barely moves (+0.803 → +0.974). This is the measurement behind scoring `F_tcr_mhc`
  with Miyazawa–Jernigan and reserving TCRen for TCR:peptide — that default was a judgement call
  and is now a number.

### Numerics

Fitting is penalised pseudolikelihood and needs **no partition function**: the conditional of one
site given the rest is a logistic regression whose extra covariates are counts of contacting
neighbours, so it is convex. Scoring needs `Z`, and gets it by annealed importance sampling that
anneals *only* the coupling term — at β = 0 the model is the uncoupled one, whose
`log Z₀ = Σ log(1 + e^η)` is exact, so the reference is a verified model, not an approximation.
Transitions are block Gibbs on a greedy colouring of the real coupling graph, with conditional
independence asserted against the edge lists rather than argued. Verified against exact enumeration
of all 2¹³ configurations (within 0.12 nat at couplings up to ±1.5) and against the closed form at
zero coupling (marginals to 0.02, `log Z` to 1e-12).

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

## [2.10.0] — 2026-08-22

*Partial: this section records only what the `[Unreleased]` block that stood here had documented.
Neither 2.10.0 nor 2.11.0 was ever published to PyPI, so the entry below first reached users in
2.12.0.*

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
