# tcren — roadmap

Where the project is going, and what each direction is waiting on. The single place for forward
plans: [STATUS.md](STATUS.md) is where the project *is*, [CHANGELOG.md](CHANGELOG.md) is what has
landed, and this file is what has not.

Last reviewed against `master` at v3.6.0 (2026-09-06); the previous review was 2026-09-03 at v3.0.0.

## Landed since the last review

- **v2.11.0** — TCRen2 became the **default** TCR:peptide potential, and the shipped matrix was
  re-derived on the 362 fully annotated αβ complexes of Native2026 (down from all 374).
  `derive-potential` now derives from αβ TCR:pMHC only, unconditionally.
- **v2.12.0** — **`P_native` shipped and is the recommended score**: three channels (geometry,
  footprint topology, contact energetics), each a conditional-linear-Gaussian Bayes network fitted
  by expectation maximization with no binding label, combined by adding log-odds; `tcren features`
  and `tcren recognize --features` are the two commands that reach it. The combiner zoo it replaces
  was deleted in the same release — see [OBSOLETE.md](OBSOLETE.md) for the list. **Reversed in
  v2.26.0**, which discarded `P_native` outright; the entry stays as the record of what v2.12.0
  shipped.
- **v2.12.1** — `paths.tcren_home()`. An installed wheel could not find its reference data at all
  before it, so tcren was usable only from a git checkout.

## Landed since that review — v2.13.0 – v3.0.0

Off this roadmap or off the manuscript's needs; [CHANGELOG.md](CHANGELOG.md) is the record and
[STATUS.md](STATUS.md) carries the condensed history.

- **v2.13.0** — `tcren.potts`, the contact map as a random variable rather than a property of one
  structure. `Π`, the interface energy referenced against the partition function, comes from here,
  and it is the third block of the composite below.
- **v2.15.0 – v2.16.0** — **a composite defined for one structure**: `S_free` =
  `Q/sd_Q + T/sd_T + (Π − μ)/sd_Π`, three fit-free directional blocks over the 374 Native2026
  crystals; `reliability.t_score`, the shape block, which loses 0.06 ROC-AUC on the balanced VDJdb
  panel when the epitope has no solved complex to template on, against `Q`'s 0.24; plus
  `inversion_flag` and `screening_yield`.
- **v2.22.0 – v2.24.0** — the presentation interface gets its own Hamiltonian, a second
  presentation potential, and `centred_potential`, the gauge that lets a pinned model reproduce a
  referenced score.
- **v2.26.0** — **`P_native` discarded**, with the whole v1 score block and every frozen-coefficient
  composite; `F_*` → `Phi_*`; `tcren.provenance` stamps every feature table with a catalogue digest
  and a table written under a different one is refused rather than silently scored.
- **v2.27.0** — the composite is `S` (`reliability.s_free` → `reliability.s_score`, column `S_free`
  → `S`), with no alias.
- **v2.28.0** — every out-of-fold-fitted read-out removed: `p_binder`, `available_links`,
  `correct_confidence`, `available_corrections` and `tcren diagnose`.
- **v2.29.0 – v2.30.0** — the catalogue 123 → 141 → **164** descriptors, and the flat middle of the
  package split into `docking`, `topology`, `energetics` and `mechanics`. This is where roadmap item
  7's **Lawrence–Colman shape complementarity** landed: `sc_shape` is their Sc computed on a raster
  rather than a dot surface, catalogued in the `topology` family. Its R² when regressed on all 141
  incumbent descriptors is 0.445 — the more familiar and the *less* novel of what 2.30.0 adopted,
  against the gap descriptors' 0.131 and 0.255 — and it was adopted anyway.
- **v3.0.0** — **`tcren.score`**: one frozen object and five read-outs, each defined for a single
  structure, with `tcren assess` and `tcren fit-holdout`. This is what replaced `P_native`, and the
  difference is the manifest — the 8,292 hold-out structures the fit used ship in the wheel and
  `fit-holdout` regenerates the shipped arrays from them bit for bit.

## Landed since that review — v3.1.0 – v3.6.0

- **v3.1.0** — a structure set is annotated **once**, not once per structure. `tcren scoring` over
  the 374-crystal `Native2026` set went 1,032 s → 145 s with identical output.
- **v3.2.0** — `tcren explain` / `score.residue_deltas`: which residue carries a score.
- **v3.3.0** — the shipped TCRen2 matrix re-derived on **359** of 374 crystals, three withdrawn PDB
  entries dropped. A provenance correction, not a change of instrument (*r* = 0.9996).
- **v3.4.0** — **the self-consistent arm became addressable.** `tcren2_em` and
  `potts_tcr_peptide_sc` resolve by name; `potts selfcon`, `potts fit --em`, `potts ensemble` and
  `potts.geom` all landed here; `tcren score --delta` and the per-template `-c` manifest with them.
- **v3.5.0** — `tcren potential-table`, so a figure can be drawn from the matrix the scoring
  commands resolve rather than from a path into the wheel.
- **v3.6.0** — **the Ca geometry prior ships and a contact map can be scored with no structure of
  the complex**: `geom_tcr_peptide`, `tcren potts geom`, `tcren potts predict`. Held out by
  structure, the predicted 15 Å envelope holds 7,399 of 7,687 crystal contacts (96.25 %) and the
  coupled model reads the crystal's indicator inside it at median ROC-AUC 0.8310 over 361
  structures. This is rung 2 of item 1 below.

## Open

Ranked by what each unblocks. Every row names what is missing, not what would be nice.

## 1. Place the interface with no crystal anywhere in the chain

**The ladder is what Cα geometry a reading is allowed**, and rung 3 is the open one:

| rung | Cα geometry from | crystal of this complex? | state |
|---|---|---|---|
| 1 | the complex's own coordinates | **yes** | `available_pairs`, shipped |
| 2 | superimposed crystals, graphon-averaged | no | `geom_tcr_peptide`, shipped in 3.6.0 |
| 3 | a parametric loop, no crystal coordinates | no | **open** |

Rung 2 already delivers "score from sequence"; what rung 3 changes is the parameterisation, and
with it the two things rung 2 measurably gives up. It places a pair about **1.5 Å further out**
than the crystal does (mean predicted d_Cα 12.46 Å against 10.87 Å on four receptors), so the
distance term deflates the layer and the expected contact count on 2b4 falls 24.6 → 4.9: the
**ordering** survives and the **level** does not. And a 12 × 12 averaged surface cannot represent a
loop whose tip sits far from the mean, which is exactly the receptor a template-free method is
wanted for.

Three pieces, two of which exist:

- **The chain model.** `geometry.reach_max` is already the closed-form maximum reach of a length-`L`
  loop held at its neck, from the discrete-Frenet virtual-bond chain, and `geometry.cdr3_internal_coords`
  already extracts the bonds, angles and pseudo-torsions to fit it against. What is missing is the
  forward direction: building Cα coordinates *from* internal coordinates.
- **The angle and torsion distributions.** The benchmark repo has an unwired stack for exactly this
  — `bench/scripts/cdr3_torsions.py` → `torsion_fit.py` → `models/torsion_vonmises.py` — written
  2026-09-01, in no `recompute.sh` target, with no output on disk. Check it before writing anything.
- **The placement.** The hard part, and the one with no prior art here: putting the built loop in
  the groove frame. `docking/align.py` and `data/canonical_frame.json` are the entry points.

The gate is the one rung 2 already has a number for: held-out recall of the crystal's contacts
inside the predicted envelope, and the area of `p_model` inside it
(`$TCREN2_BENCH/bench/eda/geom_prior_recall.py`). Rung 3 has to beat 0.9625 and 0.8310 to be worth
shipping, and it has to close the calibration gap to be worth preferring.

## 2. Make peptide stability a score, not a diagnostic

`tcren.dynamics` separates best from worst binders in 4/4 clones where the additive contact energy
fails, and the intra-peptide term behaves as Sewell's hypothesis predicts. It is still a research
readout: no `rmsf` column is in `recognition.DESCRIPTORS`, so `tcren features` does not emit it and
neither `recognize` nor `assess` can reach it. Four questions gate shipping it.

- Does it survive on **crystal** structures? Every number so far is from modelled CPL complexes, so
  the effect could be reading model quality rather than peptide physics. The Canonical2026 set is the
  control that exists already.
- Does it hold at more steps and more seeds? `rmsf` is the least noisy readout (CV 0.115 over six
  seeds) but 4000 steps and one seed is what the CPL numbers used.
- Is **contact + stability** worth shipping as a combined score? Within-clone z-sum lifts mean AUC
  0.721 → 0.826 and improves 5/7, but n = 7 clones cannot carry that (Wilcoxon p = 0.22). More clones,
  or a per-structure test.
- Does the **per-position** stability profile localise to the P3/P6 pair Dolton et al. Fig. S4 names?
  That would turn a global `rmsf` into a mechanism, and the 4C6 pair split (P = 0.0012) says the
  signal is there to find.

## 3. Repack inside the MC loop

`relax_interface` samples the backbone with χ held fixed, so a move that would be favourable *after*
the side chains relax is rejected. `repack` is now 6 ms, which is what makes a repack-per-cycle
affordable — and it is the difference between a stability probe and a real FlexPepDock analogue.

## 4. Side-chain construction

`repack` rotates the side chains a model **has**. `substitute_peptide` strips past Cβ by design, so
that path still returns 44 of 77 heavy atoms and cannot be scored at atom level at all. Needs ideal
internal geometry per residue type. AlphaFold/TCRmodel output is full-atom, so `repack` already
covers the common case — this is what the *substitution* path is waiting on. Rung 3 of item 1 needs the same machinery on the receptor side.

## 5. De novo peptide placement — PART 2 of the review

A fast in-house kernel for placing a peptide into an empty groove, benchmarked against FlexPepDock:
anchor-constrained CCD closure + DOPE MC + rotamer repack. `_refine` (rigid-body MC) and `_fold`
(CCD) already cover the *template-based* case; the packer and the backbone sampler are two of the
three pieces. Side-chain construction (4) is the third. See `refine/CPP_REWRITE.md`.

## 6. Surface topology, from descriptor to feature

Narrowed by v2.30.0, not closed. Nineteen `sc_*` descriptors built on the same
`topology.surface` height fields are catalogued now and reach the `tcren features` table, so the
module as a whole is no longer a dead end. What is still unreachable is the **epitope-shape**
read-out: the scalars exist and separate literature-named featureless from bulged epitopes
completely, but none of `relief` / `peak_to_valley` / `frac_above_ridge` / `phobic_centre` is in
`recognition.DESCRIPTORS`, so no `tcren features` table carries them and neither `Q`, `T`, `S` nor
any channel of the score set can see them. They are fit-free and z-scoreable against
`native_reference()`, which is exactly `cohort`'s premise.

Two open questions behind it: whether the **charge and hydropathy** channels carry immunogenicity
signal (Chowell et al. 2015 says TCR-contact hydrophobicity should), and whether the map distance
predicts **cross-reactivity** between epitopes rather than merely clustering copies of the same one.

## 7. Full-scale fold benchmark on aldan3

`scripts/fold_benchmark.sbatch`, n ≈ 374 with every oracle. FlexPepDock is minutes per structure — it
burned 21 min of CPU on six structures locally without finishing, which is why this is a cluster job
and not a local one.

## 8. Smaller, well-scoped

- **`2wbj`** is the single class-II Canonical2026 structure whose β-sheet core still fails to map
  (93/94).
- **`surface_distance` is an O(n²) Python loop** — fine at 374 maps, not at thousands.
- **Mouse class-II reference is sparse**; extend if a mouse class-II cohort ever needs it.
- **2D map polish**: an optional "contacting residues only" mode for less cluttered overlays.

---

## What each manuscript figure still needs from this library

The manuscript is being taken figure by figure. A library change should be able to name the figure
it serves; this is that map. Figure state lives in `$TCREN2_MS/FIGURES.md`, not here.

| figure | what it claims | what it reaches for here | still missing |
|---|---|---|---|
| 1 `fig:pipeline` | the five descriptor families and which read-out answers which question | schematic only — no library call | nothing |
| 2 `fig:potential` | the potential's genealogy, ΔΦ separation, the contact map read three ways | `potential-table`, `score --delta`, `potts contacts`, `potts geom`, `potts predict` | nothing |
| 3 `fig:peptide` | peptide ranking on combinatorial libraries, residue by residue | `tcren.ddg`, `cpl.predict_response`, `pipeline --delta` | nothing |
| 4 `fig:receptor` | receptor ranking, functionally validated | `tcren features`, `assess`, `score.binder_score`, `cohort.q_score` | unreviewed |
| 5 `fig:template` | receptor ranking under template scarcity | the same, plus the template covariate reportable | unreviewed |
| 6 `fig:dynamics` | kinetics and the licence to score one static structure | the three per-interface energies, `potts` marginals, the occupancy law | `dynamics.rmsf` is not in `recognition.DESCRIPTORS` — item 2 |
| 7 `fig:rescale` | ipTM / pLDDT diagnostics | `assess`, `score.confidence_residual`, `pose_score`, `reliability.af_band` | unreviewed |
| 8 `fig:decision` | the end-user decision procedure | the whole score set | unreviewed |
