# Changelog

All notable changes to `tcren` are recorded here. Format follows
[Keep a Changelog](https://keepachangelog.com/); versions follow semantic versioning.

## [4.0.1] — unreleased

### Added

- Exact pMHC-only surface/SASA exports and descriptive feature–performance
  associations versus ipTM, with matched rows, class/template strata, fixed
  score directions, input hashes and complete missingness accounting.

### Fixed

- Require ARDA's content-validated reference indexes and correct source-root detection.
- Capture actual searched sequence/header/index SHA256 identities around batched
  receptor and MHC annotation; reject mid-search changes and mixed input generations.
- Retain searched-index provenance in scoring, feature, recognition and exact-export sidecars.

## [4.0.0] — 2026-10-15

### Added

- Exact CLI receptor-reference, position/residue/contact/placement exports with
  complete source accounting, physical identities and annotation/potential lineage.
- Native model inspection and confidence joins; fitted assessment is explicit opt-in.
- Component MD/reference diagnostics, initial periodic-image validation and native
  predictor confidence adapters with unsupported token mapping explicitly unavailable.

### Corrected

- Undefined smoothed-reference terms avoid invalid arithmetic. Historical Gaussian
  evaluation uses Cholesky solves; Polars 2 validates binary labels without coercion.
- Historical assessment exports availability and missing descriptors; parallel
  notebook calls fail with actionable instructions before annotation.
- Documentation extra includes IPython's notebook lexer. Final-model eligibility is
  explicit: fixed empirical potentials allowed, learned/fitted readouts historical.

### Removed

- Retired geometry prior, coordinate-chain and cross-product experiment, geometry resources,
  and `potts geom`, `chain`, `cross`, `predict` commands moved with tests to the source attic.
  This is a breaking API removal for the next release, not a published 3.10.0 update.
- Feature resume and computed part reuse. Chunking rebuilds every batch from raw structures.

### Fixed

- Batched species/archive annotation and v2 contact stacking/other labels.
- Coupled Potts marginal row alignment, raw pinned-matrix serialization and near-zero
  partition subtraction. Coarse contact union is documented as an independence approximation.
- Missing Q/T blocks remain missing; signs apply in the legacy mean path; marginal scores
  ignore missing coordinates outside their own marginal.
- Frozen model catalogue checks, holdout key/label validation and explicit unavailable
  confidence residuals when confidence training is absent.
- All-failed feature commands exit nonzero and retain diagnostics.
- Affine-gap alignment traceback preserves the scored path; native kernels validate input
  shapes, index ranges and scalar domains before buffer access.
- MHC reference readiness validates actual data/index components and stages rebuilds.
- Shear direction solves the constrained in-plane eigenproblem; class-II topology retains
  separate MHC chain identities. These descriptor corrections require model/reproduction review
  before release. Attribution residue selection no longer changes the scoring mode.
- Surface clustering rejects undefined pairwise distances instead of treating them as zero.

### Documentation

- Canonical AGENTS, current status, dependency-based archival disposition and publication
  gates replace duplicated history. Original records and checksums remain in the attic.
- No frozen model arrays or recorded scientific result tables were regenerated. See AUDIT.md
  for unresolved mathematical interpretations and the scope of verification.

## [3.10.0] — 2026-09-08

### Added

- **The amino-acid background taken in 3D — `tcren.energetics.background`.** `smoothed_reference`
  is a *virtual* scan: it holds the contact map fixed and re-indexes the potential, so a contact a
  long arginine reaches across is still charged to the alanine that notionally replaces it. Measured
  on 370 Native2026 crystals over **8,002** TCR:peptide contacts at 5 Å, only **37.77 %** survive
  truncation to backbone plus Cβ and **45.03 %** are made by side-chain atoms on both sides, so that
  is not a correction to the reference — it is most of it.

  `scan_structure` builds every one of the twenty candidates at every rotamer of a χ grid and
  resolves each `(position, candidate, partner)` to **three** outcomes rather than one: no contact,
  contact with a probability, or **collision** — every rotamer of that candidate overlapping the
  environment. A colliding candidate leaves the background average with its prior mass renormalised
  over the survivors, rather than entering it at its least-bad rotamer, which is what a plain
  `w / w.sum()` does and which lets an impossible substitution set the reference free energy.
  `reference_free_energy` takes `F_i` over the **admissible set** `A_i`, and `|A_i|` is itself a
  read-out.

  Over 60 crystals and 596 peptide positions the three states run **14.48 % collision, 81.48 % no
  contact, 4.04 % contact** over 585,720 cells, with `|A_i|` mean **17.17** of 20 (p05 12, p95 20);
  24.3 % of positions admit all twenty and 4.2 % admit ten or fewer. The anchor pockets are where it
  bites — on 1AO7 the A\*02:01 P2 pocket admits 3 of 20.

- **Side-chain construction — `tcren.refine.sidechain`.** ROADMAP item 4, open since the
  substitution path was written. `substitute_residues` truncates every target at Cβ, which is
  *exact* for alanine (its heavy atoms are precisely N, CA, C, O, CB) and a stub for the other
  eighteen; this places every heavy atom of any of the nineteen non-glycine types from a backbone
  and a χ vector, by natural extension reference frame.

  **Chemistry is declared, geometry is measured.** `TOPOLOGY` names the atoms, their frames and
  which χ drives each torsion, and carries no number; the bond lengths, angles and fixed torsion
  offsets are fitted from **374 Canonical2026 crystals over 176,475 residues** by
  `scripts/fit_sidechain_geometry.py` and ship as `data/sidechain_geometry.json`. Round-tripping a
  deposited residue through its own χ reads **0.0784 Å mean heavy-atom RMSD** (median 0.0466, p95
  0.2323) over 54,954 residues, worst type ARG at 0.2210 Å and best SER at 0.0207 Å.

  `nerf` and `dihedral` invert each other to 2.4e-14. A sign error in the latter leaves every bond
  and angle correct and every atom 180° out — it reads as a plausible 4 Å RMSD rather than an
  exception — so `_selfcheck` and `test_nerf_inverts_dihedral` exist for exactly that.

- **`_relax.substitution_scan`, the C++ kernel.** The scan is positions × 20 candidates ×
  `n_grid ** n_chi` rotamers, each needing a DOPE energy against the environment and a heavy-atom
  distance to every partner. One batched call per structure, GIL released, threaded over positions
  inside — never one call per residue and never a Python pool around it. Measured against the
  `scan_numpy` reference on 1AO7 it agrees to **2.2e-16** on the contact probability and exactly on
  the admissible flag, the clashing-rotamer count and the rotamer count, at a **100× speed-up**;
  **27.5 ms** per structure over 60 crystals (median 23.4).

  **Parallelism is over cells, not positions.** The peptide side has nine positions, so a
  per-position split caps the pool at nine whatever is free; the `20 P` cells are equally
  independent and are handed out through an atomic cursor, because a tryptophan is nine atoms and
  an alanine none. With a clashed rotamer abandoned the moment one overlap is found — its energy
  is never weighted and its contacts never counted — the peptide side runs **505 ms single-threaded
  to 23.1 ms on 16 cores** and the receptor side **6,372 to 234 ms**. Both still agree with
  `scan_numpy` to 4.4e-16 and exactly on every integer output. `_geom.interface_clashes` is reused for the collision test
  rather than reimplemented.

  **The grid and the tolerance are calibrated against one invariant** — a residue the crystal
  deposited must be admissible at its own position, because the crystal is the proof that it fits.
  Over 400 deposited peptide residues in 40 crystals, 60° and 0.8 Å read **0.9725**; halving the step
  again does *not* help (0.9700 at 30°, 15× the cost), so the residual is strained crystal geometry
  rather than a coarse grid — rejections are Tyr (7 of 38) and Trp (4 of 14) and nothing else. `scripts/calibrate_background_scan.py` is that sweep.

### Changed

- **`potts.available_pairs` measures `d_heavy` out to `radius`, not to `cutoff`.** The old default
  populated it *if and only if the pair contacted*, so a named distance column was null on **90.1 %**
  of rows (154 of 171 on 3UTS) and null conditional on the outcome — any analysis regressing an
  outcome on `d_heavy` was conditioning on that outcome. Filling it costs nothing measurable and
  `sigma` is identical either way; the contact call keys on `cutoff` and never on the probe. Pass
  `probe=cutoff` for the old behaviour, which `selfcon.observation_field` reads as a confident
  non-contact.

### Deprecated

- **The Cα-geometry ladder and the Potts scoring surface are marked for removal**, with the
  measurement behind each, in `OBSOLETE.md`. Nothing is deleted and nothing changes behaviour. In
  short: the geometry prior is reachable only from the CLI and feeds no descriptor, channel or
  benchmark read-out; the couplings buy **+505.7 nats** of pseudo-likelihood and **+0.0000
  [−0.0015, +0.0000]** of ranking against independent marginals; the Π block costs **−0.0243** macro
  ROC-AUC on TCRvdb; and `potts:n_contacts`, the family's one usable descriptor, correlates with
  `interface:n_contacts_tp` at **ρ = +0.999**. What the coupled model is repurposed to is MD contact
  frequency, where it orders better than the crystal's own contact call (Spearman 0.7199 against
  0.6440) and levels worse (Pearson 0.7446 against 0.7760) over 682 cells.

## [3.9.0] — 2026-09-07

### Added

- **Every TCR read against every pMHC — `tcren.potts.cross` / `tcren potts cross`.** `Canonical2026`
  is `Native2026` after canonical MHC superposition, so every complex already sits in one frame —
  which is why `shuffle.graft_tcr` splices a pMHC and a foreign TCR with no transform. Reading the
  TCR of one crystal against the pMHC of another is therefore a pseudo-complex, and within MHC class
  there are **265² + 94² = 79,061** of them against the 359 solved. Written as a partitioned
  parquet dataset — **36,269,320 rows, 246 MB, 31 s wall clock** — the library's first parquet
  output.

  There is no loop over the 79,061 pairs: the cross product over *structures* is the cross product
  over *residues*, so the two residue tables are grouped by `(class, chain, region, len.rec,
  par.len)` and each group is one `(n_s, n_t, L, M)` broadcast.

  **The clash filter is measured, not assumed.** Over the 175,358 receptor:peptide Cα pairs of the
  359 solved complexes the closest approach is **3.366 Å**, so `--min-dist` defaults to 3.35 and
  drops no row a real complex produced. A foreign TCR dropped onto a groove does clash — 0.062 % of
  cross rows sit under 3.0 Å where **no** native row does — and the filter is per pair, not per
  pairing.

  **The frame is checked, and fixed if it is wrong.** A pseudo-complex is exactly
  `shuffle.graft_tcr` — keep the host pMHC, drop its TCR, put the donor's in — and that graft
  applies *no* coordinate transform, because every structure is already in the canonical MHC frame.
  Verified: over 40 built chimeras and 15,963 Cα cells the table matches the graft to **4.5e-5 Å**,
  float32 storage precision. On a set left in its deposited frames every pairing would be a number
  with no meaning and **nothing about the table would look wrong**, so `residue_tables` measures the
  peptide-centroid spread and `--superimpose auto` (the default) orients the set if it needs it, at
  about 3 ms a structure; a set still not co-framed afterwards is refused, and the CLI exits 1.
  `Canonical2026` reads **3.61 Å** and is passed through untouched. The bundled test fixtures read
  **141.7 Å as deposited and 0.95 Å oriented** — they were never co-framed, which the guard caught
  the moment it existed.

  What it buys is **length coverage** (131 `(L, M)` combinations for TRA:CDR3 against the 62 the
  crystals solved; 112 against 59 for TRB:CDR3) and hundreds of samples per cell. What it does not
  buy is receptor:peptide dependence beyond the 359 solved pairings — crossing every TCR with every
  pMHC *is* the product of the two marginals. `is_native` marks the diagonal, and `pdb.id.rec` /
  `pdb.id.par` stay separate columns so a caller can hold out by structure on either axis.

- **A geometry drawn from a fitted joint — `tcren.potts.chain` / `tcren potts chain`.** A mean of
  distance matrices is not a distance matrix any chain can realise: consecutive Cα are **3.7985 Å**
  apart (s.d. 0.0668 over 642 CDR gaps, 0.47 % *cis* at ~2.9 Å) and no per-cell average enforces it.
  This models the **coordinates** instead — a joint Gaussian over the concatenated graphon-knot Cα
  of all six CDR loops and the peptide, low-rank plus diagonal — so a draw is a set of real traces,
  `geometry.enforce_bond` puts each exactly on the virtual bond, and the distance matrix is derived.

  Three consequences. Every geometry is realisable with no convex relaxation anywhere. One draw
  moves all six loops and the peptide together, which is what carries the docking angle — a
  whole-receptor rigid-body degree of freedom no per-loop model can express. And **length stops
  being a joint problem**: `L` conditions the loop block and `M` the peptide block separately, so a
  CDR3 of 16 against an 11-mer needs no crystal that solved that pair.

  **The draw set contains the right geometry.** Held out by structure over 540 observed contacts in
  24 crystals, the best of 512 draws places every contact to **0.04 Å**. Fit 0.3 s, 256 draws
  0.04 s.

- **`tcren.geometry` gains the primitives the rest of this rests on** — `graphon_weights`,
  `resample_graphon`, `read_back` (which collapse a quadruplicated bilinear resample onto one
  matrix form that also batches), `enforce_bond`, `virtual_cb` and `orientation_cosine`. Plus
  `paths.canonical_dir()`.

- **"Does the residue turn toward the peptide?" is now a term.** Until now the site table carried
  `d_heavy` and `d_ca` and no Cβ, so the model could read that a pair was close and not whether the
  residue faced it. `orientation_cosine` is the cosine between the Cα→Cβ direction and the direction
  to the partner Cα. Held at a fixed Cα distance it separates contact from non-contact at **ROC-AUC
  0.69–0.84** across 4–12 Å.

  Fitted as `g_orient` — eight bins, 369 crystals, `J` pinned to `tcren2` — it is worth **+787.0
  nats** of pseudo-log-likelihood for eight parameters, against the 3.8.0 exclusion class's +81.8
  for one, and the profile is monotone across 3.05 log-odds:

  | cos bin | coefficient | sites | | cos bin | coefficient | sites |
  |---|---|---|---|---|---|---|
  | [−1.00, −0.75) | −1.2922 | 5,332 | | [+0.00, +0.25) | +0.0290 | 8,302 |
  | [−0.75, −0.50) | −1.0063 | 6,491 | | [+0.25, +0.50) | +0.5930 | 9,584 |
  | [−0.50, −0.25) | −0.8301 | 7,068 | | [+0.50, +0.75) | +1.2907 | 10,866 |
  | [−0.25, +0.00) | −0.5433 | 7,688 | | [+0.75, +1.00) | +1.7592 | 9,781 |

  The Cβ is **reconstructed, not read**: a drawn geometry has no side chains, so `virtual_cb` builds
  it from three consecutive Cα (0.547 Å mean residual) and `available_pairs` uses the same
  reconstruction on a crystal, so fitting and scoring describe one quantity. The reconstructed
  cosine correlates with the real one at Pearson **0.81** and keeps ROC-AUC 0.68–0.79. Glycine is
  not a special case — this is where a Cβ *would* sit, so there is no missing value to encode.


- **Where `g_orient` pays is a crystal, not a drawn geometry.** Held out by structure over 24
  crystals, read on the complex's **own** Cα the map improves on both placement and discrimination —
  Jensen–Shannon similarity **0.7052 → 0.7290** and contact ROC-AUC **0.9599 → 0.9691**, with the
  level unmoved (ratio 1.002 → 1.026). For scale, the 3.8.0 exclusion class moved js 0.7133 → 0.7162
  there and left the ROC-AUC flat, so this is roughly eight times the placement gain and the first
  block to move the ordering as well.

  On a geometry **drawn** from `potts.chain` it is neutral — js 0.4200 → 0.4203, ROC-AUC
  0.9255 → 0.9259 over 12 held-out crystals. The reason is not that the term is wrong there but
  that a drawn backbone carries the *population's* typical orientation rather than this receptor's:
  `cos_cb` is a deterministic function of the Cα trace, so it can only be as receptor-specific as
  the trace is. Use it when you hold coordinates.


- **`tcren potts predict --chain-prior`**, plus `--geometry-draws` and `--beta`.

### Measured and recorded, so it is not re-tried

- **Annealing over drawn geometries costs rather than pays.** `--beta` weights each draw by its own
  `log Z`, computed on the **union** of every draw's reachable cells — a draw that put a cell out of
  reach gets η = −60 — so the site count cannot confound the comparison. Over 10 crystals held out
  by structure, 96 draws, 22.6 observed contacts each:

  | arm | E[contacts] | ratio | js | w_eff |
  |---|---|---|---|---|
  | step 2b (in-sample baseline) | 24.7 | 1.10 | 0.4078 | — |
  | **chain, β = 0** | **22.0** | **0.90** | **0.4342** | 96.0 |
  | chain, β = 0.02 | 26.8 | 1.10 | 0.4278 | 69.6 |
  | chain, β = 0.05 | 36.4 | 1.46 | 0.4042 | 18.8 |
  | chain, β = 0.2 | 46.5 | 1.82 | 0.3601 | 1.9 |
  | chain, β = 1 | 47.1 | 1.84 | 0.3538 | 1.1 |

  The common site set rules out the site-count explanation, so the mechanism is the **contact**
  count: `log Z` is a free energy and rises with occupancy. `beta` defaults to **0**. Note the
  chain arm is held out by structure while the step-2b baseline is the bundled prior, which
  contains these crystals.

- **The graphon registration and read-back are not adjoint, and that is deliberate to keep.**
  `potts.geom._resample` registers residue *i* at grid coordinate `i(G−1)/(L−1)` while
  `predicted_distances` reads it back from `(i+1)(G−1)/(L+1)`, so a round trip contracts a loop's
  footprint toward its midpoint by `(L−1)/(L+1)` — 1.29 Å at the termini of a 6-mer, and **a fixed
  bias no grid refinement removes**. Over 7,732 observed contacts in 370 leave-one-out crystals,
  predicted minus measured Cα distance in Å:

  | registration / read-back | CDR2 | CDR3 | all |
  |---|---|---|---|
  | endpoint / interior (shipped) | +5.76 | +2.29 | **+2.91** |
  | endpoint / endpoint (adjoint) | +5.57 | +2.78 | +3.26 |
  | interior / interior | +5.63 | +2.79 | +3.28 |

  An adjoint pairing buys 0.19 Å on CDR2 and costs 0.49 Å on CDR3. `"endpoint"` registration aligns
  loop **termini** across structures of different length, which is the operation an average over
  crystals performs, and the `"interior"` read-back samples away from the surface's extrapolated
  edge cells. **So the CDR2 deficit is not a registration artefact** — it survives every pairing at
  5.6–5.8 Å against CDR3's 2.3–2.8 Å. `graphon_weights(convention=...)` makes both reachable and
  records the numbers.

### Unchanged

- **Every shipped model and prior is bit-for-bit what it was.** `g_orient` is `[]` in all four
  bundled Potts models and an empty block contributes nothing; `site_codes` emits the orientation
  code only when the model carries a fitted `g_orient` **and** the frame carries `cos_cb`, so a
  table written before 3.9.0 scores identically. `SITE_COLUMNS` is unchanged — `cos_cb` rides
  beside it. `geom_tcr_peptide` and `geom_tcr_peptide_ens` are untouched, as are
  `fit_geom_prior`'s and `read_back`'s defaults.

## [3.8.0] — 2026-09-07

### Added

- **A receptor residue's grip on the peptide can now be told to stop.** `X(same receptor residue,
  |dj| > 2)` is a fourth coupling class, the complement of `K`'s reach on one receptor residue:
  `OFFSETS` excludes `(0, 0)` and stops at `|dj| = 2`, `_edges_cross` needs two different loops and
  `_edges_class` needs two different partner classes, so nothing previously coupled two sites on one
  residue more than two peptide positions apart. `K(+0,+1) = +0.637` and `K(+0,+2) = +0.385` build a
  contiguous run and nothing said it ended.

  Fitted on 367 `Native2026` crystals with `J` pinned to `tcren2`: **K = −0.457, s.e. 0.049,
  z = −9.34**, worth **+81.8 nats** of pseudo-log-likelihood for one parameter — comparable in
  significance to the strongest cross-loop class and three times its magnitude.

  `fit_potts(..., far=True)`, `PottsModel.far`, `edges(q, far=...)`, `kernel_names(joint, far)`.
  Goes through the existing colouring, block Gibbs and AIS untouched; `tests/unit/test_potts.py`
  pins it against exact enumeration of all 2¹² configurations.

- **Where it bites is a predicted geometry, not a crystal's own.** On the crystal's own Cα the field
  is already sharp and the couplings have little left to move: Jensen–Shannon similarity to the
  crystal's map goes 0.7133 → 0.7162 and the pairs above p = 0.8 go 3,447 → 3,474 at a conversion
  rate of 0.905 → 0.907. Under a Cα geometry averaged over reference crystals, with each receptor
  held out of its own prior, three of the four `Native2026` receptors improve on **both** placement
  and level:

  | receptor | js, `far=False` | js, `far=True` | E[contacts], `far=False` | `far=True` | observed |
  |---|---|---|---|---|---|
  | `3qib` / 2B4 | 0.4846 | **0.4957** | 36.80 | **34.11** | 32 |
  | `3qiu` / 226 | 0.5226 | **0.5350** | 33.81 | **31.40** | 33 |
  | `4p2r` / 5cc7 | 0.4834 | **0.4924** | 29.97 | **27.82** | 25 |
  | `3qdg` / DMF5 | **0.4538** | 0.4529 | 19.10 | 19.00 | 24 |

  The ROC-AUC does not pay for it (0.9406 → 0.9431, 0.9547 → 0.9580, 0.9639 → 0.9640). `3qdg` is
  flat and is the one receptor already under-counted, where suppression has nothing to remove.

### Unchanged

- **Every shipped model carries `far=False` and is bit-for-bit what it was.** `potts_tcr_peptide`,
  `potts_tcr_peptide_sc`, `potts_tcr_mhc` and `potts_pep_mhc` all still hold 18 kernel coefficients,
  and a model JSON written before 3.8.0 round-trips unchanged. The class does **not** deliver the
  Beta-shaped map on a predicted geometry on its own: no pair reaches p ≥ 0.8 in either arm there.

## [3.7.1] — 2026-09-07

### Added

- **A `contacts` stage in `reproduce.sh`** — the 3.7.x claim end to end through released commands
  in 33 s: a leave-one-out geometry prior, a TSV of six loop sequences and a peptide, `tcren potts
  predict --sequences --ensemble`, and `tcren potts agreement` against the crystal the prior was
  held out of. It exits non-zero if the predicted geometry does not read below the crystal's own.


- `tcren.potts.map_agreement` and `tcren potts agreement` — score a contact map against the one a
  structure actually made, reporting the **level** (`ratio` = Σp/Σσ) and the **placement** (`js`,
  Jensen–Shannon similarity of the two maps normalised to distributions over sites, in bits, so
  0–1 and scale-free in the level) separately, plus `auc` for the ordering alone.

  `js = 1` requires a deterministic map, so it is read against the same model on the structure's
  own Cα geometry and not against 1: on the four Figure 2c receptors the crystal's own geometry
  reads 0.75–0.81 and the predicted geometry 0.42–0.58.

  **R² is deliberately not offered.** On a binary outcome `1 - SS_res/SS_tot` *is* the Brier skill
  score — they agree to three decimals on every map checked — so it is one number reported twice;
  and it scales with prevalence, which makes maps of different site counts incomparable, exactly
  the comparison the predicted and measured geometries require.

## [3.7.0] — 2026-09-07

**The predicted contact map recovers the contact count.** 3.6.0 scored one map read off the mean
Cα surface, and that map carried a fifth of the contacts the crystal makes. The defect was the
order of operations, not the parameterisation: a contacting pair sits at the near end of the
distribution being averaged, so the mean places it **2.6–4.5 Å too far out** (against 0.2–2.3 Å
over all pairs), and the distances within one map are not independent — consecutive Cα atoms are a
fixed virtual bond apart — so no set of per-pair marginals recombines into a geometry a chain could
adopt.

### Added

- `tcren.potts.predict_contact_map` — scores the map under **each** reference crystal's geometry
  and averages the contact probabilities, rather than scoring one map under their average. Every
  member is a chain that was actually solved, so the virtual-bond constraint holds by construction.
  Held out by structure over 364 `Native2026` crystals, against 21.02 contacts observed per
  structure on CDR1–CDR3 (median 21, s.d. 6.08):

  | estimator | mean Σp | ratio to observed | median contact ROC-AUC |
  |---|---|---|---|
  | mean surface | 9.80 | 0.466 | 0.9242 |
  | **geometry ensemble** | **24.62** | **1.171** | **0.9325** |

  Both the level and the ordering improve. One batched `contact_probabilities` call scores every
  member — the sampler already splits by structure, and a member is exactly a structure to it.
- `geom_tcr_peptide_ens` ships in the wheel: 2,154 registered surfaces from the same 359 αβ
  `Native2026` crystals at grid 32, 3.5 MB, distances stored to 0.01 Å as `uint16`.
  `GeomPrior.bundled()` resolves `.npz` before `.json`, so an ensemble prior and a mean-only one
  are asked for the same way.
- `GeomPrior.surfaces`, `GeomPrior.ensemble`, `to_npz` / `from_npz`, `predicted_ensemble`, and
  `fit_geom_prior(..., keep_surfaces=True)`.
- `tcren potts geom --keep-surfaces` (requires an `.npz` `--out`) and
  `tcren potts predict --ensemble [--members N]`. The scored table gains `p_avail`, the fraction of
  reference geometries placing that pair within the radius at all.

### Unchanged

`geom_tcr_peptide` and `predicted_sites` are untouched and remain the default, so nothing fitted
against 3.6.0 moves. What the ensemble does **not** buy is receptor specificity: the predicted
count correlates with the observed one at Pearson +0.23 (Spearman +0.18) and the per-structure MAE
of 5.95 contacts is about the spread of the population itself. It places the typical interface at
the right level, which is what loop lengths alone can carry.

## [3.6.0] — 2026-09-06

**A contact map can be scored without a structure of the complex.** `tcren.potts.geom` modelled the
receptor:peptide Cα distance matrix from sequence and lengths alone but had no fitted prior and no
command; both now ship, so the site set a Potts model scores no longer has to be read off
coordinates the complex was solved at.

### Added
- **`geom_tcr_peptide` — a fitted Cα geometry prior, in the wheel**, reachable as
  `GeomPrior.bundled()`. Averaged in graphon coordinates over the **same 359 αβ `Native2026`
  crystals `tcren2` is derived from** (12 not αβ and 3 withdrawn, of 374), so a prediction and a
  potential rest on one structure set. Twelve surfaces at grid 12 — TRA/TRB × CDR1/CDR2/CDR3 ×
  MHCI/MHCII, 265 class I and 94 class II crystals per key.

- **`tcren potts geom -s <structures> -o <json> [--grid] [--exclude]`** — fits it. `--exclude`
  leaves named `pdb.id`s out, which is how a held-out prediction is built; the held-out *design*
  belongs to whoever is evaluating, so the command only offers the exclusion. Withdrawn PDB entries
  are dropped unconditionally, as they are by `derive-potential` and `potts selfcon`.

  **The MHC groove is annotated whether or not the caller asks.** The surface is keyed on MHC
  class, and without the groove pass every crystal reports class I: a prior fitted without it
  silently carries six surfaces where it should carry twelve and reads class II off a class I
  average.

- **`tcren potts predict`** — contact probability per residue pair with the Cα distances
  *predicted* rather than measured. `potts contacts` reads the distance matrix off the structure,
  so both its site set and its distance field come from solved coordinates; here the distances come
  from the prior's surface read back at the lengths asked for, and the 15 Å availability cut is
  applied to the **predicted** distance.

  Two inputs, and the difference is the whole point. `-s` opens the structure for its loop
  sequences, its peptide and its MHC class only, and joins the crystal's own `sigma` and `d_ca`
  back on so the table can be scored against what the crystal did. `--sequences` takes a TSV of
  `id`, `peptide`, optional `mhc_class` and any of the six loop columns `TRA:CDR1 … TRB:CDR3`, and
  reads **no structure at all**; `sigma` is null there, because nothing observed a contact. One
  column order either way, `d_ca_obs` appended only when a crystal was read.

  Emits the `potts contacts` columns plus `d_ca_sd`, the surface's own spread at that cell, and
  `feasible`, false where no loop of that length could span the distance predicted — checked
  against `tcren.geometry.reach_max`, the closed-form discrete-Frenet reach bound, rather than
  against the fit.

- **`tcren.potts.feasible` is exported.** `predicted_sites` has always called it and written its
  result into the emitted `feasible` column, but it was in neither `geom.__all__` nor
  `potts.__all__`, so a caller could read the column and not the rule behind it.

### Changed
- **`_potts_structures` is the one batched-annotation path every `potts` command shares.** The loop
  was inline in `_potts_pairs`; `potts geom` and `potts predict` need the same structures without
  the pair table, and duplicating it would have meant a second place for an mmseqs-per-structure
  regression to appear.

### Fixed
- **`predicted_sites` crashed on a non-standard residue.** `X` has no code in the model's alphabet
  and two `Native2026` crystals carry one, so a set-wide prediction raised `KeyError: 'X'` rather
  than dropping the cell. `available_pairs` has always filtered these; now both do.
- `tcren2()`'s docstring said the matrix is derived over **362** crystals. It has been **359**
  since 3.3.0; the manifest and the CHANGELOG were already right.

## [3.5.0] — 2026-09-06

**The shipped residue-pair tables become readable from the command line**, so a figure or a table
can be drawn from the matrix the scoring commands resolve rather than from a file path into the
installed wheel.

### Added
- **`tcren potential-table -p <name|path> -o <csv>`** — writes out the residue-pair energy table
  behind `-p`, long form (`residue.aa.from, residue.aa.to, value`). Nothing is computed: it is the
  shipped matrix, resolved by the same name `score`, `ddg`, `cpl` and `rank` resolve, which is what
  lets `bench/figures/fig2/make_fig2_tables.sh` draw main-text Figure 2a from `tcren2`, `tcren2_em`
  and `karnaukhov2022` without reaching into `src/tcren/data/`. All three are 380 pairs over 20
  residues: cysteine appears on the peptide axis and never on the receptor axis, so 20 of the 400
  cells are absent in every one of them.

## [3.4.0] — 2026-09-06

**The self-consistent potential and the model pinned to it become addressable from the command
line**, and `tcren score` gains the referenced energy the benchmark had been computing by hand.

### Added
- **`tcren score --delta`, with `--reference structural|virtual`** — ΔΦ, the poly-alanine-referenced
  energy, per candidate. It is what makes two ranking units commensurate, and no command emitted it
  for a *candidate pool*: `tcren scoring --delta` scores a structure's own peptide, and `tcren score`
  had no reference at all. Two poly-alanine baselines exist and they are **not** the same quantity.
  `--reference virtual` re-indexes poly-alanine on the template's own contact map, so the baseline
  is charged for every contact the real side chains make; `--reference structural` (the default)
  truncates the peptide and **recomputes the map**, so the baseline is what backbone plus Cβ alone
  reach. On 2b4/3QIB under `karnaukhov2022` the two references read **−1.3916** and **−0.0166**
  natural-log units against a scored energy of −2.3104, and poly-alanine retains 62–85 % of the
  crystal's TCR:peptide contacts across the templates measured.

  Φ(reference) does not depend on the candidate, so it is computed **once per template** and
  subtracted from every row: the cost is one rebuild per structure, not one per candidate. On a
  single template ΔΦ is Φ minus a constant and cannot reorder anything — it changes the scale and
  the zero, and it separates candidates only where each carries its own structure.

- **`-c/--candidates` accepts a per-template manifest** — a TSV or CSV whose header names
  `complex.id` and `peptide`. Each structure then scores only its own rows. A plain one-peptide-
  per-line file is unchanged and every structure still scores every peptide. A benchmark of 204
  units each with its own 1,001-peptide pool was 204 CLI calls and therefore **204 mmseqs
  chain-typing searches**; it is now one call and one batched pass. `tcren rank` reads the same
  manifest.

- **`tcren2_em` resolves as a bundled potential name**, in `-p/--potential`,
  `--tcr-peptide-potential` and `potts fit --coupling-matrix`. The matrix shipped in the previous
  release with a Python loader only: `_BUNDLED_POTENTIALS` in `pipeline.py` and the loader dict
  inside `potts.centred_potential` are what every option resolves against, and neither knew the
  name, so the third arm of any three-matrix comparison was reachable from Python and nowhere else.

- **`potts` `-m/--model` accepts a bundled name** — `potts_tcr_peptide` (the default),
  `potts_tcr_peptide_sc`, `potts_tcr_mhc`, `potts_pep_mhc` — as well as a JSON path, across
  `potts score`, `contacts`, `map`, `ensemble` and `scan`. `PottsModel.bundled(name)` always took a
  name; the CLI could only pass a path, so the shipped self-consistent model had to be located
  inside the wheel by hand before it could be used.

- **`TCRen2-EM` ships: `data/TCRen2_sc_potential.csv`, `data/potts_tcr_peptide_sc.json`,
  `tcren.potential.tcren2_em()`, manifest key `tcren2_em`.** The self-consistent matrix, derived
  by `potts selfcon` over the same 359 `Native2026` crystals as `tcren2` with the same estimator,
  reference state and pseudocount — what differs is that every AVAILABLE residue pair contributes,
  weighted by the probability the coupled model gives it, rather than only the contacts one crystal
  shows.

  It is a **smaller** step from TCRen2 than the 2022 matrix is: Pearson r = 0.9444 against TCRen2
  over the 380 cells (Spearman 0.9371, mean |delta| 0.1052, largest 0.8010 at K:T), against
  r = 0.8670 between TCRen2 and the 2022 matrix. The movement lands where the crystal tally was
  thin — mean |shift| 0.1599 in the quartile of cells whose median cell rests on 2 observed
  contacts, 0.0564 in the quartile at 47.

  **`tcren2()` remains the default potential and `PottsModel.bundled()` remains the default
  contact model.** On one identical site table the pinned self-consistent model is a *worse* pair
  term than a pin on TCRen2 at the same parameter count (k = 114 both), so nothing was switched;
  see `Changed` below for the measurement.

  Reproducibility: three independent runs at `draws` 100/400/400 and `steps` 15/25/9 agree at
  **r = 0.999997**, mean |delta| 0.0008 — the fixed point does not depend on the sampling budget.
  The shipped artifact is the converged one, `--tol 5e-3 --draws 400`, converging at step 9.

  The manifest gains three optional keys for it, the first entry to use any: `command`
  (default `derive-potential`), `structure_flag` (default `--structure-dir`; `potts selfcon` takes
  `-s`), and `tolerance` / `reproduce`. **`tcren2_em` is the only shipped matrix that does not
  reproduce at 1e-9**, because its E-step is sampled; its recipe is checked at 5e-3 under
  `RUN_BENCHMARK=1`. Every other entry still reproduces at 1e-9 and must keep doing so.

- **`tcren.potts.selfcon` — the potential re-derived from the contacts the model predicts, until
  it stops moving**, as `tcren potts selfcon`.

  `derive_tcren` tallies the amino-acid pairs of the contacts a crystal *shows*; `fit_potts` reads
  the finished matrix back as its pair term. Nothing carried the model's own view of which contacts
  formed into the tally, and the two halves never met. Each step now recomputes the contact
  posterior under the current model, re-derives the potential from those contacts **weighted by
  their probability**, damps it against the previous matrix, and refits with the new matrix pinned
  as the pair term.

  **Weighting whole configurations by their probability costs nothing.** The tally
  `sum_delta P(delta | obs) n_ab(delta)` **equals** `sum_a p_a` by linearity of expectation, so the
  configurations never have to be enumerated to get the counts. And the estimator is untouched:
  `derive_tcren` already took a per-contact fractional weight, so the log-odds, its
  contact-fraction reference state and its pseudocount are the shipped ones.

  **The observation enters as one additive field.** `P(sigma_a = 1 | d_a) = sigmoid((5 - d_a)/m)`
  is a per-site constant on `eta`, so the E-step is the existing block-Gibbs sampler called with a
  modified vector — no clamping machinery and no second code path. One knob spans the family:
  `margin = 0` pins the crystal map, `margin = 0.5 A` releases the threshold band, a large margin
  lets the crystal enter only through the fit.

  **The reduction gate is exact.** At `margin = 0` the probability-weighted tally reproduces the
  plain unweighted derivation bit for bit — `max |soft - hard| = 0` over all 380 cells on the
  bundled crystals, weighted mass 309.000000 against 309 observed contacts. Over Native2026 with
  `--margin 0 --balance both` it reaches its fixed point at step 2 with `max_abs_de = 0.0` exactly
  and reads **R = 0.9977 against the shipped TCRen2**, mean |delta| 0.0075 over 380 cells. That is
  the reduction run and not the production one: at `--margin 0.5` the loop moves the matrix to
  R = 0.9439 and stops on the sampler's own noise, never on an exact fixed point. The residual is
  the availability set: the 15 A Ca radius holds **7,799 of the deposited derivation's 7,816
  TCR:peptide contacts (99.78 %)**, the 17 it misses spread one or two apiece over 17 structures,
  and they land in the thinnest cells — the largest single disagreement is K:K, the rarest pair.

  **This is EM-shaped but it is not EM**, and the docs say so. The E-step is a genuine conditional
  expectation, but the potential update is a re-derivation rather than the maximiser of one
  objective, so the pseudo-log-likelihood need not rise and on the bundled crystals it does not.
  Convergence is declared on `max_abs_de`, which contracts at the damping rate.

- **`fit_potts(coupling_matrix=...)` accepts a `Potential`, not only one of five bundled names.**
  It resolved through a hardcoded dict of `lru_cache`d loaders, and `coupling_matrix_name` was what
  serialised, so a model could not be pinned to a matrix derived in the same session — which is
  exactly what the loop above does on every step. `PottsModel` gains `coupling_pinned`; a pin on an
  unnamed matrix stores the matrix itself and round-trips, and `n_parameters()` still counts it as
  the single scale it is (1, not 400).

- **`tcren potts pairs --probe`**, so the published pair table can carry a heavy-atom layer measured
  past the cutoff. Without it `d_heavy` is null for every non-contact and both `potts fit --em` and
  `potts selfcon` can only ever remove contacts, never add them.

- **`tcren.potts.ensemble` — the contact map read as an ensemble rather than as a set of
  marginals**, with `tcren potts ensemble` on the command line.

  `contact_probabilities` settles the *mean* interaction energy and nothing else. `Phi` is linear
  in the contact indicator, so `<Phi> = sum_a <sigma_a> phi_a` holds for any joint, coupled or not
  — but only when `<sigma_a>` is the coupled marginal `p_model` and is read **at the sequence being
  scored**. The shipped soft-contact path (`ddg(weights=...)`,
  `cpl.predict_response(tcr_weights=...)`) takes one weight vector computed on the template and
  applies it to every threaded candidate, so it computes `sum_a p_a(template) phi_a(x)` and not
  `sum_a p_a(x) phi_a(x)`. `threaded_ensemble` recomputes the marginals per candidate;
  `energy_ensemble` adds the read-outs no set of marginals determines — `phi_var` (needs the
  two-point function), `phi_free = -1/beta log <exp(-beta Phi)>`, `map_entropy = <E> + log Z`, and
  `hamming_obs`. Over the 14 alpha-beta complexes in `tests/assets/pdb`, Var Phi runs **0.341 to
  1.286 squared TCRen units**, a 3.8-fold range.

- **`map_ensemble` returns the contact configurations themselves, ranked, each with its
  probability — coarse-grained, because a site-level map has none.** Measured on 1ao7's
  TCR:peptide interface (174 available sites, 6,400 draws): the raw site-level map recurs never,
  6,400 distinct configurations out of 6,400 draws at `p_hat = 1.6e-4` each, the state space being
  `2^174`. Coarse-grained to which peptide positions are engaged, the same draws collapse to **41
  distinct patterns with the most likely at `p_hat` = 0.3219, and that modal pattern is the one the
  crystal made**. `by="position"` is therefore the default; `by="loop"` (39 groups, 3,308 patterns,
  top 0.0216) is at the edge of what 6,400 draws resolve and `by="pair"` is only meaningful where
  `exact_log_z` can enumerate instead. Across the 14 complexes the modal pattern is the observed
  one in 6 of 14.

- **`tcren.potts.em` — Monte-Carlo EM over the contact calls the 5 A threshold cannot resolve**,
  reachable as `tcren potts fit --em`. A pair within half an angstrom of the cutoff is a coin toss
  recorded as a fact; the E-step resamples those sites from their conditional with every confident
  site clamped, and the M-step is the existing penalised pseudolikelihood on `draws` stacked
  completions at weight `1/draws`, so the effective sample size and the meaning of the ridge are
  unchanged. `consistency(reference, other)` reports how far a refit moved per parameter block,
  scaled by the reference's own cluster-robust standard errors where it has them. **This adds a
  fit; the shipped pseudolikelihood models are unchanged and are the comparison.** `margin=0`
  reduces it to `fit_potts` exactly.

- **`tcren.potts.geom` — a length-generalising prior over the receptor:peptide Ca distance
  matrix**, so a pair can be placed from sequence and lengths alone with no template of the right
  shape. `P(d_ca | loop, i, L, j, M, MHC class)` is averaged in graphon coordinates
  `(u, v) = (i/(L+1), j/(M+1))` — the same normalisation `contactmap.registered_map` uses — so a
  surface fitted at one set of lengths reads back at any other, including lengths the fit never
  saw. `predicted_sites` emits the schema `available_pairs` emits, so every `tcren.potts` function
  accepts it unchanged. All three hypervariable loops carry a surface, not just CDR3: on the
  bundled crystals CDR3 is **69.3% of TCR:peptide contacts**, CDR1 a further **26.5%**, and
  CDR1+CDR2+CDR3 **98.4%**, so a CDR3-only prior would cap its own recall below 0.70 by
  construction. Feasibility is checked against `tcren.geometry.reach_max`, the closed-form maximum
  apex reach of a length-`L` loop pinned at its neck, rather than against the fit. `sigma` on a
  predicted table is a placeholder and is 0 throughout.

### Fixed
- **A batched structure set is identified per SET, not per file.** `tcren.annotation.batch`
  resolved every `pdb.id` through `structure_id_from_path`, which takes the part before the first
  underscore — right for `4x5w_renumbered.cif`, wrong for any cohort whose names carry metadata
  there. The seven Bobisse models `302TIL_A0206_9_TCR-pMHC` … all came back as `302TIL`, so seven
  distinct structures landed under one id and a per-template candidate pool could not be addressed
  at all. `tcren.structure.io.resolve_structure_ids` has always made this decision for the set and
  `iter_structures` has always used it; `iter_annotated_set` and `annotate_batch` now do too, which
  is why `tcren scoring` and `tcren score` used to disagree about what a structure is called.
  **A set whose prefixes are already unique is unaffected** — all 374 `Native2026` ids are
  byte-identical either way.

### Changed
- **`-s/--structures` is resolved the same way by every command.** `_typed_batch` now runs the
  spec through `resolve_sources`, so a repeated `-s`, a comma-separated list, a glob and a manifest
  mean the same thing everywhere. `tcren score` declared `-s` as a single path and silently kept
  only the last one.
- **`fit_selfconsistent` and `potts selfcon` report `rms_de` beside `max_abs_de`, and the
  Monte-Carlo floor is now measured rather than asserted.** `max_abs_de` is a maximum over 380
  cells, so its floor is far above the r.m.s. of the same 380 numbers and no draw budget drives it
  to zero: on `Native2026` the plateau falls 0.00391 (`draws=100`, steps 10-15) to 0.00217
  (`draws=400`, steps 13-25), a ratio of **0.553 against the 0.500 that quartering the variance
  predicts**. Reaching 1e-3 on it would take roughly 1,800 draws and would measure nothing —
  across those same steps `r_vs_initial` moves by 6.2e-5 and `beta_matrix` by 2.1e-4, so the matrix
  has long stopped. `rms_de` is the statistic to read for that question, and at the step the
  default `tol` of 5e-3 trips it already reads 8.79e-4. `tol` still applies to `max_abs_de`, so no
  existing run changes behaviour.
- `tcren.potts.available_pairs` takes an optional `probe`, the radius the **heavy-atom** layer is
  measured out to, defaulting to `cutoff` so every existing table is byte-identical. It is not the
  contact definition — `sigma` is still called at `cutoff` — it only records how far a
  non-contacting pair missed by. Without it `d_heavy` is null for every non-contact and a latent
  band around the cutoff reaches only inwards: on the 17 bundled test structures the band selects
  44 sites and every one is a contact, against 44 contacts and 71 near misses once `probe` is set.
- `tcren.potts.gibbs` takes an optional `init` configuration, so a chain can be started at the
  structure's own contact map instead of at a draw from the factorised model. `init=None` is the
  default and leaves the random stream, and therefore every existing result, untouched. This is
  **not** a Hamming-ball constraint — an unconstrained chain started at the observed map still
  samples `P(sigma)`; `energy_ensemble`'s `hamming_obs` reports how far it travels.

## [3.3.0] — 2026-09-05

**The shipped TCRen2 matrix is re-derived on the 359 `Native2026` crystals the RCSB still holds**,
`python -m tcren.cli explain` reaches a command that the module entry point could not see, and
`fetch-data` reaches the whole HuggingFace deposit rather than five named folders of it.

### Changed
- **The default TCR:peptide potential is derived on 359 of the 374 `Native2026` structures, where
  every release from 2.11.0 through 3.2.0 shipped a matrix derived on 362.** Two filters now run rather than one: the
  alpha-beta restriction drops 12 files (3 pMHC-only, 8 single-chain, 1 gamma-delta), and the
  withdrawn-entry rule drops 3 more against the new bundled `data/pdb_withdrawn.csv`. `1jtr` and
  `4zez` are superseded by `1mwa` and `5jzi`, both already in the set, so one complex was entering
  twice in a set built to be non-redundant and the redundancy weighting could not see it because the
  two carry different PDB ids. `4euq` was withdrawn with no successor: the RCSB records that its TCR
  alpha chain was modelled from the wrong gene at positions 93-98, which is the CDR3alpha region this
  potential reads.

  **This is a provenance correction, not a change of instrument.** Against the 362-structure matrix
  it replaces, the new matrix agrees at Pearson *r* = 0.9996 over all 380 cells (19 x 20, cysteine
  dropped from the `from` axis), with a largest absolute shift of 0.104 on a TCRen2 range of 2.95
  (-1.25 to +1.70) and exactly one cell moving by more than 0.10. Scores computed under the two are
  still not strictly comparable; `data/potentials.json` records the new recipe and
  `tests/regression/test_shipped_potentials.py` reproduces the file from it.

### Added
- **`potential.withdrawn_entries()` and `potential.withdrawn_ids()`** — the bundled withdrawn-entry
  table (4 rows) and the set of its PDB ids, lowercased. The rule lived inline in the
  `derive-potential` CLI body, so nothing outside the CLI could apply it: the benchmark's
  epitope-weighted evaluation ran on 370 structures while the matrix it evaluates was derived on 359,
  and no gate could see the gap. `--keep-withdrawn` restores all 4 for auditing that difference;
  3 of the 4 change the derivation set, the fourth (`3nfj`) being a gamma-delta receptor the
  unconditional alpha-beta filter already excludes.
- **`derive-potential --contacts-out`** — write the contact table the matrix was derived from,
  **after** every filter. It reports 7,816 TCR:peptide contacts over 358 structures, one fewer than
  the 359 the derivation keeps, because one kept structure makes no TCR:peptide contact at all. A
  downstream analysis that wants to count the structures or contacts behind the shipped matrix had to
  assemble its own table, and that table held 370 structures.
- **`fetch-data` reaches the whole deposit.** `HF_FOLDERS` was a constant naming 5 sets against the
  14 the Hub lists today, so every other set was fetched by ad-hoc `curl`. The listing is now read from the Hub:
  `--list` shows what exists, `--sets` takes names or `all`, `--tables` adds the derived-table
  deposit, and a set added to the deposit is fetchable the same day.
- **`features --chunk N --resume`** — chunked featurisation with resumable parts, replacing three
  hand-rolled chunk-split-concatenate loops in the benchmark repository. Chunking is for memory, not
  speed: a parsed complex is about 1 MB, so a 20,000-structure corpus does not fit in one pass. The
  property that must not be traded away is the batched annotation — one arda call and one MHC search
  per chunk, never one per structure — so `tests/unit/test_features_chunking.py` asserts the chunked
  table equals the whole-set table exactly, that `--resume` reuses its parts, and that a run without
  `--resume` leaves none behind.
- **`potts pairs`** — the availability enumeration promoted from a private helper to documented
  surface. Two benchmark modules already imported the private form, and a private import is not an
  interface. It is also the first and most expensive stage of `potts fit`, so emitting it once and
  reusing it turns a repeated corpus pass into a single one; it stamps provenance like every other
  emitter.

### Fixed
- **`python -m tcren.cli explain` reported "No such command".** The `if __name__ == "__main__":`
  guard sat above the `explain` registration, so the module form ran the app before that command was
  registered while the console script, which imports the module first, saw it. **This was the second
  occurrence** — the comment above the guard already documented the identical failure for `potts`,
  and that fix moved the guard only as far as `potts` rather than to the end of the file, so the next
  command added below it inherited the bug. The guard is now the last statement in the module, and
  `tests/unit/test_cli_entrypoints.py` asserts the two entry points expose identical command sets
  (27 today), so a third recurrence fails the build. That test resolves the console script beside
  `sys.executable` rather than by `shutil.which("tcren")`, which had resolved to a stale `uv tool`
  install several versions behind the checkout and read as a difference between entry points when it
  was a difference between builds.
- **The HuggingFace tree endpoint pages at 1,000 entries** and hands the next cursor back in a `Link`
  header; `_hf_files` read the first page only. `Native2026` reported 346 files of 374 through that
  path, and a folder of 3,212 would have reported 1,000, downloaded 1,000 and raised nothing.
  `_hf_tree` now follows the cursor to the end. The skip-if-present check had a matching defect: it
  required a loose `.pdb`, so the 6 sets that ship a `.tar.gz` rather than loose structures — `cpl`,
  `as_case`, `atlas_ala`, `immrep23_negatives` and the two vdjdb benchmarks — were re-downloaded on
  every run.
- **`viz.pymol.importance_scene` crashed on a non-finite value** — it was formatted straight into
  `cmd.alter` as `b=nan`, which PyMOL evaluates as Python and raises `NameError` on, and PyMOL's own
  error printer then swallows the reason so the caller sees an empty traceback. A residue can
  legitimately have no value, because masking it can leave a descriptor undefined; those rows are now
  dropped and drawn as pale cartoon like any other unscored residue. The zero-centred ramp was gated
  on `by == "phi"`, the only signed column when it was written, and is now gated on the values, so
  `residue_deltas` colours correctly and `phi` is unchanged. The loader took `.pdb.gz` only, where a
  deposited model set holds plain `.pdb`; one `_struct()` helper resolves either at all 5 call sites.
- **`docs/potentials.rst` documented a `--ab-only` flag for `derive-potential` that does not exist.**
  The alpha-beta restriction is unconditional. The fix is to delete the flag from the documentation
  rather than add it to the command: making it optional would be a trap, because `--balance` skips a
  structure carrying a null on any weighting axis and the derivation then defaults it to weight 1.0,
  the maximum, so the least usable files would have counted for the most.
- **`test_the_hot_limit_is_the_distance_from_the_best_residue` hard-coded a tolerance of 1e-4** where
  its own docstring says the residual is the runner-up residue's `exp(-beta * gap)`. That held only
  while the smallest gap stayed where it was, so it failed on the re-derived matrix — and it was the
  assertion that was wrong, not the quantity. Both tolerances are now derived from the matrix's own
  gaps: log-sum-exp above the minimum for the delta, and the runner-up weight times the squared gap
  for the curvature. It is now a test of the identity rather than of one matrix's spacing.

### Documentation and reproduction
- **`bootstrap.sh` and `reproduce.sh` at the repository root**, and the split between them is the
  point. `bootstrap.sh` is the environment and the data — it wraps `setup.sh`, fetches the whole
  HuggingFace deposit (about 620 MB across 14 sets) through `tcren fetch-data`, and builds the
  28,579-allele MHC reference — and is **deliberately outside any time budget**, because how long a
  download and an `mmseqs` index take is a property of the machine. `reproduce.sh` assumes both are
  on disk and runs four library stages — the version and catalogue digest, the potential re-derived
  and asserted against the shipped matrix, the frozen score model refitted from its committed
  slice, and the whole score set computed live over all 374 `Native2026` crystals — then hands off
  to `$TCREN2_BENCH/recompute.sh` and `$TCREN2_MS/make`, each located by an environment variable
  with a documented default. Measured 2026-09-05 on macOS arm64: `--library-only` is **106 s**
  (env 1 s, potential 24 s, model 1 s, scores 80 s).

  **Neither script has a `--skip-slow` flag and neither should get one.** Both fail loudly, by
  name, on a missing input: a bad set name, an absent `Native2026`, an unbuilt MHC reference or an
  absent sibling repository is an exit code with the fix in the message, never a warning.
- **`docs/reproduction.rst`** — the recipe walked end to end, from `fetch-data` through `features`
  and `assess` to the figure and table build, naming the released version and the two provenance
  stamps that decide whether any number is comparable. Until now it lived only in a CLI docstring.
- **`docs/rederivation.md` renders at all.** It is 238 lines in Markdown, was in no toctree, and
  `conf.py` had no Markdown parser, so Sphinx never saw it. `myst-parser` is registered (two lines,
  against retyping eight tables in reStructuredText) and the page carries a dated banner: the
  candidate matrix it names, `TCRen_potential_rederived.csv`, was removed on 2026-08-22, and the
  default potential it calls `classic` has been `tcren2` since 2.11.0.
- `SOURCES` carried an entry for that removed file and none for `TCRen2_potential.csv`, which is the
  default. Both fixed, and `pdb_withdrawn.csv` added.

## [3.2.0] — 2026-09-04

**Which residue carries a score.** The five channels say which *part of the structure* a read-out
comes from; `residue_deltas` says which **residue**, for any read-out and any channel:

```
tcren explain -s model.pdb --score binder -o deltas.tsv
```

Each row is one interface residue and `delta` is `L(x) - L(x_without_i)` — the score of the complex
minus the score of the same complex with that residue's atoms removed. Positive means the residue
carries the score. The table feeds `tcren.viz.pymol.importance_scene` unchanged, so one
implementation colours a structure by `binder_score`, `pose_score` or any single channel, where
before only the peptide read-out had a colouring at all.

Why it is leave-one-out and not a decomposition: `peptide_score` is a sum over contacts and
`energetics.scoring.position_profile` already splits it **exactly**. The other four are functions of
~149 whole-structure scalars through a full covariance per class, so no residue owns a share of one
and there is nothing to split. The leave-one-out difference is defined for all of them.

* Chain typing and the MHC call run **once** for the whole structure; every masked copy is
  featurised with `annotate=False` against that one annotation. No mmseqs call is made per residue.
* Only residues touching a requested interface are masked — a residue making no contact moves no
  descriptor, and a test asserts its delta is exactly zero rather than assuming it.
* A masked row whose descriptors go undefined is reported as a **null** delta, never imputed:
  `ScoreModel.coordinates` drops such a row on purpose, and that is surfaced rather than filled in.
* Cost is one descriptor pass per interface residue, about 2 s; a complex is a few minutes.

Additive only: no signature changed and nothing was removed.

## [3.1.0] — 2026-09-04

**Structure sets are annotated once, not once per structure.** Every CLI command that iterates a
structure set was calling `classify_chains` and/or `annotate_mhc` inside the loop, each spawning an
`mmseqs easy-search`; `scoring` did it *twice* per structure, because `pipeline.run` annotates and
then calls `superimpose(annotate=True)`, which annotates again. A four-structure `scoring` run made
**48 mmseqs subprocess calls and spent 91% of its wall time in `select.poll`**.

Measured on the 374-structure Native2026 set, `tcren scoring` end to end:

| | before | after |
|---|--:|--:|
| wall time | **1,032 s** | **145 s** |
| output | — | identical, max abs difference **1.1e-16** (on `rmsd`; every Φ column exact) |

Fixed:

* `pipeline.run` passes `annotate=False` to `superimpose` — the structure it hands over is already
  chain-typed and MHC-annotated.
* `Structure.mhc_calls` caches the `MhcCall` list both annotation paths already compute, so
  `typed=True` now means chain-typed **and** MHC-annotated and nothing downstream re-derives it.
* `scoring` calls `annotate_mhc_batch` once over the whole input set.
* New internal `_typed_batch()`: two mmseqs searches for a whole set, mmseqs threading internally.
  `annotate`, `contacts`, `score`, `ddg`, `cpl`, `energy`, `refine` and the score-background path
  all take it. No command types or MHC-annotates per structure any more.
* `-t/--threads` on `scoring` no longer exists to hide annotation cost behind concurrent mmseqs
  processes: it sets the thread count of the two batched searches and the worker count for the
  contact maps and energy sums. The docstring said the opposite.

No behavioural change and no API removal — this is a performance release. Output is byte-comparable
to 3.0.1 within float noise.

## [3.0.1] — 2026-09-02

A packaging fix and a test-tolerance correction. **The frozen model is untouched** and no score
moves; this is the release that made a clean `pip install tcren` able to refit it.

### Fixed
- **`scikit-learn` was an undeclared runtime dependency.** `score.model._shrink` imports
  `sklearn.covariance.ledoit_wolf_shrinkage` and it sits on `Joint.fit`'s path, but scikit-learn was
  listed only under the `notebooks` extra, so a clean install of the 3.0.0 wheel raised
  `ModuleNotFoundError` inside `Joint.fit`:
  `test_refitting_the_committed_slice_reproduces_the_frozen_arrays` — the reproducibility contract
  the manuscript's Code availability paragraph names — **failed** rather than skipped. Declared as
  `scikit-learn>=1.3` under `[project] dependencies`, which is what it is, and added to the tests
  workflow's hand-maintained lean install list. Scoring from the frozen model was never affected.
- **`tests/unit/test_smoothed_reference.py` errored instead of skipping without `arda`.** It reaches
  `classify_chains` from a module-scoped fixture, so all 7 of its tests **errored at setup** on a CI
  runner with no arda reference database. Guarded with the `pytest.importorskip("arda")` the other 5
  arda-backed test modules already use. Verified in a venv built to the workflow's exact recipe
  (Python 3.12, `pip install -e . --no-deps` plus the explicit list, no arda):
  `pytest -m "not slow" --cov-fail-under=50` goes from 1 failed and 7 errors to **662 passed, 82
  skipped**, exit 0.
- **The refit reproducibility test asserted a tolerance the optimiser never promised.** It checked
  the two class means and the ipTM covariance row at rtol 1e-9, on a docstring premise that they
  "refit bit-identically". `Transformer.fit` takes each of its 55 Yeo-Johnson lambdas from
  `scipy.stats.yeojohnson_normmax`, whose Brent search stops at its own default tolerance of about
  1.5e-8, so a lambda is only ever determined to that precision and every array fitted through the
  transform inherits it. Measured between macOS arm64 and Linux x86_64: one lambda differs by
  8.4e-8 absolute (1.3e-6 relative), the refit class means by at most 6.0e-9 absolute on 18 of the
  149 coordinates, and the covariances by at most 2.5e-8 absolute on 50 of the 22,201 entries. The
  test passed on the machine it was written on and failed on CI. All the float comparisons now share
  one `_TOL` of rtol 1e-5 / atol 1e-6, declared once beside the fixture paths with the reason. The
  absolute part is needed because those 50 covariance entries sit near zero, at a magnitude of about
  9e-5, where a relative bound reports 2.9e-4 for a difference smaller than the one the means
  survive.

### Changed
- `appendix/generated/counts.tex` defines `\tcrenVersion`, so a version bump makes the generated
  appendix layer stale and `test_the_generated_appendix_layer_is_current` fails. Regenerated with
  `scripts/gen_appendix.py`.

## [3.0.0] — 2026-09-02

**The score layer is rebuilt. One frozen object — a transform and a Gaussian per class over the
transformed descriptor coordinates — and five named read-outs that are all projections of it, with
no further fitting. Every one is defined for a single structure, so a user holding one AlphaFold
model can run the whole set on it.**

The major version is for the new public surface and the new package data, not for a removal: every
2.30.0 entry point still works and returns what it returned.

### Added
- **`tcren.score`**, the score set, with semantic names in the API and in the CLI:

  | name | tier | what is estimated |
  |---|---|---|
  | `peptide_score` | 0 | nothing; the direction is fixed by the potential |
  | `pose_score` | 1 | a covariance over hold-out binders — no negative, no label |
  | `confidence_residual` | 1 | the same covariance, read as a conditional mean |
  | `binder_score` | 2 | class means and covariances, from hold-out binder labels |
  | `channel_scores` | 2 | the same object, marginalized to one descriptor family |

- **The five channels, and they are what makes a number explainable.** A marginal of a Gaussian is
  a sub-block of its covariance — exact, closed form, no re-fit — so asking *which part of the
  structure says this complex is real* costs an index and nothing else: `placement` (where the
  receptor sits in the groove frame), `interface` (how much interface it makes, of what chemistry),
  `shape` (the footprint free of its size), `energetics` (the contact chemistry in kT), `mechanics`
  (the interface as a network of breakable springs). Those are the same subsets the fit-free `Q`,
  `T` and the energy block were built on by hand.
- **`tcren assess`** — every read-out for a `tcren features` table, one row per structure.
  `--peptide` when the peptide is what varies, which restores the presentation descriptors the
  receptor read-out marginalizes out.
- **`tcren fit-holdout`** — refit the frozen model from its manifest and a descriptor table, and
  get the shipped arrays back bit for bit. This is the answer to why `P_native` was withdrawn: its
  coefficients were frozen against a training set nobody could reconstruct, and these are frozen
  against one that is named, deposited and shipped.
- **Package data**: `holdout_model.npz` (315 kB — the transform parameters, the class moments and
  the confidence coordinate) and `holdout_manifest.csv.gz` (329 kB — the 8,292 structures the fit
  used, with dataset, epitope, label and ipTM). The 19 MB descriptor table deliberately does **not**
  ship; `tcren fetch-data` plus `tcren features` regenerates it from the ids the manifest names.
- **`tcren.score.transform`** — nine transform classes keyed on a descriptor's unit and operator.
  Monotone and variance-stabilising throughout, and **never rank-based**, which is measured rather
  than argued: mapping each marginal onto a uniform took the per-cohort median ROC-AUC from 0.630
  to 0.543 on six template-covered cohorts and from 0.613 to 0.507 on sixteen template-free ones.
  The signal *is* the marginal scale, and flattening a distribution deletes it.
- **`reliability.artefact_directions`** — a label-free test for which directions of a binder
  manifold belong to the structure generator rather than to the interface. A direction real
  crystals break harder than modelled binders do is the generator's own regularity: measured, the
  tightest band breaks 4.55x on crystals against 3.18x on decoys and scores 0.504 over sixteen
  template-free cohorts, a coin, where the loose band reads 0.606. Nothing in it reads a binder
  label. Ledoit-Wolf shrinkage floors the shipped model's smallest direction at s.d. 0.0797, above
  that band, so `pose_score` cannot read it.
- **`OPERATOR`, `OPERATOR_NAME` and `OPERATOR_OBJECT` move into `tcren.descriptors.catalogue`**,
  which is where catalogue data belongs; `scripts/gen_appendix.py` imports them from there.

### Changed
- **`recognition.STATUS` gains four exact identities** — `sc_gap_index` = (`sc_dh` +
  `sc_gap_mean`)/2, `dPhi_tcr_soft` = `dPhi_tra_soft` + `dPhi_trb_soft`, and `n_contacts_tm` /
  `n_contacts_tp` as the sums of their contact-type tallies. Found by a rank-revealing search over
  21,939 corpus structures and verified to machine precision. **With all fifteen removed the
  catalogue is full rank**, so this list is now closed rather than partial. `STATUS` only: the
  registry digest does not move and no recompute follows.
- `Joint.fit` takes `shrink=False`, which the artefact test needs and nothing else should use.

### Documentation
- **The score set is exported at package level**: `peptide_score`, `pose_score`,
  `confidence_residual`, `binder_score`, `channel_scores`, `score_table`, `holdout_model`,
  `holdout_manifest`, `ScoreModel` and `CHANNELS` are in `tcren.__all__`, so
  `from tcren import score_table` works. The import costs nothing measurable — `tcren.recognition`
  was already pulled in transitively.
- **Two end-to-end example notebooks**, both executed and both running from deposited coordinates
  to a ranking. `notebooks/score_vdjdb_panel.ipynb` fetches the 1,089-model balanced VDJdb
  benchmark, featurises it, scores it with `score_table` and reports the receptor ranking **one
  cohort at a time** over the 22 epitopes, template-covered (6) and template-free (16) apart, with
  the five channel scores per cohort. `notebooks/rank_peptides_cpl.ipynb` does the peptide task on
  the combinatorial-peptide-library set (7 clones, 2,103 models): per-clone ROC, the rank
  correlation against the assay's graded activation, and a whole response matrix from one template.
- **`docs/tcren.rst` is rebuilt against the package tree.** It had been `automodule`-ing 17
  deprecated shim paths as canonical (nine top-level modules plus the whole `tcren.orient` package,
  eight of whose submodules no longer exist and produced twelve autodoc import warnings) while
  omitting `tcren.score`, `tcren.descriptors`, `tcren.topology.*`, `tcren.energetics.*`,
  `tcren.mechanics.*` and `tcren.docking.*` entirely.
- **The five "what a contact potential cannot do" essays moved out of `README.md`** into
  `docs/beyond-the-contact-sum.rst`; the README is 864 -> 730 lines and now carries the score set,
  which was absent from it.
- Every reference to `tcren binder` and `tcren diagnose` is gone from `README.md`,
  `docs/`, `CLAUDE.md`, `skills/tcren/SKILL.md` and `notebooks/`. Neither command has existed since
  2.26.0 and 2.28.0 respectively.

### Fixed
- **`annotate_structure_set` failed on every structure**, because `contact_table` was left behind in
  `tcren.paper.helpers` when the batching moved to `tcren.annotation.batch` at 2.29.0 and no import
  followed it. `on_error="skip"` swallowed the `NameError` per structure, so
  `tcren derive-potential --structure-dir` returned an empty frame and died several frames later on
  a missing `cdr3a` column. The import is now function-local, since `paper.helpers` re-exports this
  module and a top-level import would close the cycle the split opened. With it,
  `test_recipe_reproduces_the_shipped_matrix[tcren2]` passes again: the derived matrix matches the
  shipped `TCRen2_potential.csv` to within 1e-9 over all 380 amino-acid pairs.
- `tcren fit-holdout --help` pointed at `tcren score --model`, which is the epitope-threading
  command and has no `--model`. It is `tcren assess --model`.

### Not removed
- `cohort.q_score`, `reliability.t_score` and `reliability.s_score` stay, documented as the
  fit-free predecessor tier. `S` leads the functionally validated receptor screen on its own
  (0.818 against ipTM's 0.795) and **composes** with `binder_score` rather than being replaced by
  it — together they read 0.783 on the template-covered stratum against 0.665 and 0.682 alone.

## [2.30.0] — 2026-09-02

**Twenty-three descriptors, 141 -> 164: the published interface literature, measured on the
19,213-structure benchmark corpus before being adopted. The gap between the two faces is the find,
and it was already computable — `topology.surface.surface_complementarity` had shipped twelve
quantities and not one of them was catalogued.**

### Added
- **`tcren.topology.literature`**, and with it 23 catalogued descriptors in the `topology` family.
  Every one was measured on the whole corpus first: the R^2 of an OLS of it on **all 141 incumbent
  descriptors** decides whether it is a channel this package could already reach.
- **Nineteen surface descriptors, `sc_*`.** Both faces are rasterised as height fields on one
  shared groove-frame grid, so the gap is their difference cell by cell, `h_tcr - h_pmhc`, and no
  new geometry is computed — `surface_map(s, side="pmhc")` and `side="tcr"` already existed.
  `sc_gap_mean` and `sc_gap_sd` carry the largest binder/non-binder contrasts of any published
  descriptor tested (Cohen's *d* -0.651 and -0.681 out-of-panel, 4,907 labelled structures) at
  R^2 **0.131** and **0.255** on the incumbents. `sc_shape` — Lawrence & Colman's Sc on a raster
  rather than a dot surface — is the more familiar and the *less* novel: R^2 0.445, nearest
  incumbent `m_erank_tm` at rho 0.414.
- **The gap integrated over the contact plane, with the two signs kept apart**, because they mean
  opposite things and a mean cancels them. `sc_gap_vol` is the void, `sc_interlock` the
  interdigitated volume, both in A^3; `sc_gap_index` is the void over the retained contact area,
  the intensive form. On 1AO7 interlock is 1,089 A^3 against a 315 A^3 void.
- **Four more reading the same field by sign**: `sc_interlock_frac` (share of cells with a negative
  gap — 0.771 on 1AO7, against the 0.71 this module's own 60-crystal calibration reports),
  `sc_gap_depth` (how far in, over those cells alone), `sc_gap_height` (how far off, over the void
  cells alone) and `sc_gap_asym` (the balance of the two volumes, in [-1, 1]). `sc_gap_mean` is
  `interlock_frac` weighting the two, so these are the three numbers it collapses into one.
- **`co_pep` / `co_mhc`** — contact order (Plaxco 1998), the package's first sequence-separation
  descriptor. New on the MHC helices (R^2 0.275); **not** new on the peptide (R^2 0.519), because
  our peptides are overwhelmingly 9-mers so the length normalisation has almost no range: 619
  distinct values over 12,662 structures.
- **`partcoef_tcr` / `partcoef_pmhc`** — Di Paola's participation coefficient
  `P_i = 1 - sum_s (k_si/k_i)^2`, their most discriminative descriptor. Adopted with its verdict
  recorded rather than hidden: on this corpus it is **redundant**, R^2 0.730 and 0.720 on the
  incumbents, the pMHC side's nearest neighbour being `g_loop_overlap` from the same bipartite
  object. It costs nothing over the contact map already built.

### Changed
- `surface_complementarity` returns seven more keys (`gap_vol`, `interlock`, `gap_index`,
  `interlock_frac`, `gap_depth`, `gap_height`, `gap_asym`). Existing keys are unchanged.
- `descriptors("topology")` is 70 columns, 34 compositional / 22 geometric / 14 topological. The
  family is no longer identical to `footprint_topology_features()`; `test_features_families` now
  asserts the union rather than an equality that only forbade the family from growing.
- Units gain `A^3`. The vocabulary is closed by test, so a genuinely new physical unit is added
  there deliberately.

### Notes
- **Two corrections the audit proposed are falsified on this corpus and were NOT applied.** "Any
  sum-over-contacts potential is extensive in interface size" is false for ours — every `Phi_*` has
  R^2 <= 0.044 on interface size, `Phi_tcr_pep` 0.006 — so per-BSA normalisation solves a problem
  this package does not have. And "`S_cell` is interface size wearing a hat" does not hold: R^2
  0.303, *lower* than `D2_pep24`'s 0.441.
- The gap is the **channel** Jones & Thornton 1996 named, not their gap volume index, which divides
  a Voronoi gap volume by interface ASA. Ours is a raster height-field gap and the catalogue says
  so rather than citing their formula.

## [2.29.0] — 2026-09-01

**Eighteen descriptors, 123 -> 141: the footprint's evenness re-derived from the contact graph
itself, length-agnostic readings of the Calpha/Cbeta maps, and five quantities the package already
computed but never catalogued. Plus the package is laid out in named layers.**

### Added
- **Eight contact-graph descriptors** (`tcren.topology.graph`). The existing coverage measures bin
  contacts into a fixed 12/24-cell partition and the existing topology measures build a flag complex
  at an arbitrary 7/8 A radius; these read the 5 A bipartite incidence directly, so neither a
  partition nor a radius is chosen. `g_comp_frac` is the parameter-free counterpart of
  `fp_b0_frac_r7`, and `g_cyclo_frac` the size-free cyclomatic number the footprint docstring has
  said since 2.12 that the raw form is not. `g_loop_overlap` is new in kind: the mean Jaccard overlap
  of the six loops' partner sets, which reads a footprint collapsed onto a few residues.
- **Four spectral and two Cbeta descriptors of the interface distance maps.** `m_erank_*`/`m_gap_*`
  are functionals of the singular spectrum of `exp(-D^2 / 2 sigma^2)` at sigma = 5 A; `m_face_*` is
  the mean `d_Calpha - d_Cbeta` over contacting pairs, positive when side chains lean towards each
  other and negative when the backbones are close and the side chains are not.
- **Five descriptors promoted out of `tcren.topology.pose` into the catalogue**: `m_face_tp`
  (previously `sidechain_toward`), `ca_cb_agreement_tp`, `ca_cb_agreement_tm`, `degree_evenness_tp`
  and `frac_well_coordinated_tp`. All five were computed and returned since the pose layer was
  written and named in no registry, so `tcren features` could not emit them.
- **`scripts/audit_architecture.py`** — the import graph against the intended layering. `--check`
  exits non-zero on any module-level cycle, which is what a CI gate would use.

### Changed
- **`recognition.py` is three modules.** It had grown to 1,151 lines holding the catalogue,
  computing the interface block and dispatching batches, so a caller that only wanted to know what a
  column means imported all of it. Now `descriptors.catalogue` (data and selection, no arithmetic),
  `descriptors.compute` and `descriptors.table`. `tcren.recognition` re-exports all three and no
  existing import changes.
- **The flat middle of the package is four named sub-packages**: `docking` (renamed from `orient`,
  which said only that something gets moved -- it holds a native reimplementation of TCRdock's
  rigid-body parameterisation), `topology`, `energetics` and `mechanics`. Every old top-level name
  is a transparent re-export.
- `_iter_typed` moved out of `cli` and the batched-annotation helpers out of `paper.helpers`, both
  into `annotation.batch`. Six modules had been importing upwards into the app layer for them.
- The interface energy sum (`_interface_energy`, `_contact_weights`, `_phi_scale`) moved from
  `pipeline` to `energetics.scoring`, where four modules had been reaching up for it.
- `energetics.ddg` is `energetics.mutation`: as a submodule it shadowed the `ddg` function of the
  same name on the package.
- `pose.sidechain_toward` -> **`m_face_tp`**. One name for one number; it is identical to the
  catalogued column to 4.4e-16 over 196 crystals, and shipping both would recreate the duplicate the
  2026-07-28 audit removed.

### Note on length
Nine of the 141 carry a new `STATUS` entry recording a measured coupling to peptide or CDR3 length,
on 143 class I Native2026 crystals over both axes. This is not a deprecation: every one keeps most
of its variance after both lengths are regressed out (58 per cent in the worst case), and a
length-coupled column is useful beside a length-free one because a model can form the contrast that
cancels the shared part. What the entries buy is knowing which is which. The mechanism on the
peptide axis is the class I groove: closed at both ends, so a longer peptide must bulge, and a
bulged peptide sits closer to the receptor.

### Unchanged
`Q`, `T`, `Pi` and `S`, and every coefficient in `data/reliability_moments.json`. The blocks are
`Q_FEATURES_GEOM` and `T_FEATURES_TOPO`, neither of which gains a term.

### Note for anyone holding a feature table
`provenance.registry_digest()` changes, so `tcren recognize --features <table>` raises on any table
written by 2.28.0 or earlier. Recompute with `tcren features`.

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
