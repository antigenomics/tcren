# tcren — roadmap

Where the project is going, and what each direction is waiting on. The single place for forward
plans: [STATUS.md](STATUS.md) is where the project *is*, [CHANGELOG.md](CHANGELOG.md) is what has
landed, and this file is what has not.

Last reviewed against `master` at v3.0.0 (2026-09-03); the previous review was 2026-08-27 at v2.12.1.

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

## Open

Ranked by what each unblocks. Every row names what is missing, not what would be nice.

## 1. Make peptide stability a score, not a diagnostic

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

## 2. Repack inside the MC loop

`relax_interface` samples the backbone with χ held fixed, so a move that would be favourable *after*
the side chains relax is rejected. `repack` is now 6 ms, which is what makes a repack-per-cycle
affordable — and it is the difference between a stability probe and a real FlexPepDock analogue.

## 3. Side-chain construction

`repack` rotates the side chains a model **has**. `substitute_peptide` strips past Cβ by design, so
that path still returns 44 of 77 heavy atoms and cannot be scored at atom level at all. Needs ideal
internal geometry per residue type. AlphaFold/TCRmodel output is full-atom, so `repack` already
covers the common case — this is what the *substitution* path is waiting on.

## 4. De novo peptide placement — PART 2 of the review

A fast in-house kernel for placing a peptide into an empty groove, benchmarked against FlexPepDock:
anchor-constrained CCD closure + DOPE MC + rotamer repack. `_refine` (rigid-body MC) and `_fold`
(CCD) already cover the *template-based* case; the packer and the backbone sampler are two of the
three pieces. Side-chain construction (3) is the third. See `refine/CPP_REWRITE.md`.

## 5. Surface topology, from descriptor to feature

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

## 6. Full-scale fold benchmark on aldan3

`scripts/fold_benchmark.sbatch`, n ≈ 374 with every oracle. FlexPepDock is minutes per structure — it
burned 21 min of CPU on six structures locally without finishing, which is why this is a cluster job
and not a local one.

## 7. Smaller, well-scoped

- **`2wbj`** is the single class-II Canonical2026 structure whose β-sheet core still fails to map
  (93/94).
- **`surface_distance` is an O(n²) Python loop** — fine at 374 maps, not at thousands.
- **Mouse class-II reference is sparse**; extend if a mouse class-II cohort ever needs it.
- **2D map polish**: an optional "contacting residues only" mode for less cluttered overlays.
