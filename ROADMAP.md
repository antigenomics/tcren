# tcren — roadmap

Where the project is going, and what each direction is waiting on. The single place for forward
plans: [STATUS.md](STATUS.md) is where the project *is*, [CHANGELOG.md](CHANGELOG.md) is what has
landed, and this file is what has not.

Last reviewed 2026-08-18, against `master`.

Ranked by what each unblocks. Every row names what is missing, not what would be nice.

## 1. Make peptide stability a score, not a diagnostic

`tcren.dynamics` separates best from worst binders in 4/4 clones where the additive contact energy
fails, and the intra-peptide term behaves as Sewell's hypothesis predicts. It is still a research
readout: nothing in `score`, `recognize` or `cohort` can reach it. Four questions gate shipping it.

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

The scalars exist and separate literature-named featureless from bulged epitopes completely, but
`tcren surface` is a dead end — none of `relief` / `peak_to_valley` / `frac_above_ridge` /
`phobic_centre` reaches the 34-feature `recognize` table or `cohort`'s fit-free scores, so nothing
downstream can use them. They are fit-free and z-scoreable against `native_reference()`, which is
exactly `cohort`'s premise.

Two open questions behind it: whether the **charge and hydropathy** channels carry immunogenicity
signal (Chowell et al. 2015 says TCR-contact hydrophobicity should), and whether the map distance
predicts **cross-reactivity** between epitopes rather than merely clustering copies of the same one.

## 6. Full-scale fold benchmark on aldan3

`scripts/fold_benchmark.sbatch`, n ≈ 374 with every oracle. FlexPepDock is minutes per structure — it
burned 21 min of CPU on six structures locally without finishing, which is why this is a cluster job
and not a local one.

## 7. Smaller, well-scoped

- **Lawrence–Colman shape complementarity** (`src/_geom/geom.cpp:13`) — the one hard interface
  descriptor still missing, and the surface ray-casting makes the surface normals cheap.
- **`2wbj`** is the single class-II Canonical2026 structure whose β-sheet core still fails to map
  (93/94).
- **`surface_distance` is an O(n²) Python loop** — fine at 374 maps, not at thousands.
- **Mouse class-II reference is sparse**; extend if a mouse class-II cohort ever needs it.
- **2D map polish**: an optional "contacting residues only" mode for less cluttered overlays.
