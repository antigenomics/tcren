---
name: tcren
description: tcren — TCR-pMHC contact potential (TCRen) pipeline; conventions and public API
---

# tcren Skills Guide

`tcren` reproduces and extends the TCRen contact-energy potential (Nat Comput Sci 2022)
on a pure-Python pipeline (structure parsing → contacts → TCR/MHC annotation → potential
derivation → epitope-ranking benchmarks). Annotation uses the `arda` package
(mmseqs2-backed), a runtime dependency published to PyPI as `arda-mapper` (imports as `arda`,
`>=2.5.7`) — no separate checkout and **no `ARDA_HOME`** (arda auto-fetches its reference and a
static mmseqs binary on first use). Repo-local `.venv` via uv, no conda (`bash setup.sh`).

## Batch annotation — never loop (mmseqs2 is the parallel layer)

**All structure annotation (TCR chain typing AND MHC allele mapping) must gather every
sequence first, make ONE batched mmseqs2 call, then map the output back for downstream
per-structure analysis.** mmseqs2 parallelises internally across threads — that is the
parallel layer; Python orchestration is a single call.

- Each per-structure annotate call pays a fixed ~825ms mmseqs2 process+index-load cost;
  a batch of 300 sequences costs the same ~930ms total.
- A `ProcessPoolExecutor(fork)` over structures **deadlocks** (fork after mmseqs2/BLAS
  spawn threads). A `ThreadPoolExecutor` runs but still pays the fixed cost N times.
- `paper/helpers.py::_batch_annotate` does TCR annotation for a whole dataset in 2 arda
  calls (human + mouse). MHC annotation uses the same pattern: `mhc.annotate_mhc_batch(structures)`
  — ONE mmseqs search over every candidate MHC chain, sliced back per structure.

Reference: `arda.annotate_sequences([(id, seq), ...])` — one call, threads internally.

## Threading model — annotation batched, threads only for structural ops

- **Annotation (TCR + MHC) is never Python-threaded and never per-structure.** It is one
  batched mmseqs2 call; mmseqs2 is the parallel layer (do NOT pass it a thread count). No
  `ProcessPoolExecutor`/`workers`.
- **Use threads ONLY for the embarrassingly-parallel, mmseqs-free stages:** structural
  alignment (Kabsch/SVD superposition), peptide mutation, relaxation, and rendering — i.e.
  pymol / Rosetta / FlexPepDock and figure generation. `orient.run_folder(threads=…)` threads
  the parse and the align+write stages (default `os.cpu_count()`); annotation between them is
  the single batched pass. `tcren orient -t N`. **`superimpose` is the same**: `run_superimpose`
  batch-annotates all inputs, then threads the mmseqs-free ensemble alignment + write
  (`superimpose(..., annotate=False)` on the pre-annotated structures). `tcren superimpose -t N`.
  `-s` accepts file / dir / .tar.gz / glob; `-o` is a directory, or a single structure file
  (one input) whose extension must match `--mmCIF`/`--compress` (validated by `_output_target`).

## Two orientation commands — `superimpose` vs `orient`

- **`tcren superimpose` / `tcren.orient.superimpose(s, db_dir=…)`** — bring a NEW structure into
  the canonical frame against a canonical *database* (default `data/Canonical2026`). It detects
  the input's MHC class + species, selects every DB member of that class+species (from the DB's
  `orient_metadata.json` — bundled in the package for the shipped `Canonical2026`, which
  `fetch-data` populates with structures only), superposes the query groove Cα onto each, and
  **averages** the rigid transforms (chordal/SVD mean rotation + mean translation) into one
  consensus placement. The matching DB subset is batch-annotated once and cached per process.
- **`tcren orient` / `tcren.orient.run_folder(...)`** — BUILD a canonical DB from native
  complexes using the per-class derived frame (how `Canonical2026` is produced). Not for orienting
  a single new structure — use `superimpose` for that. Writes `<out>/orient_metadata.json` (the
  format `superimpose` reads); `--metadata foo.csv` writes CSV instead.
- **HF upload is NOT a user command.** `--push-to-hub` was removed; maintainers run
  `scripts/push_canonical_to_hub.py` instead.

## Structure scoring — `tcren.run_pipeline` / `tcren scoring`

- **It is called `tcren scoring`, not `pipeline`** (renamed 2026-07; the old name errors with a
  pointer). It only scores: canonicalisation, region mapping and the Cα/contact/atom-distance
  matrices are `tcren annotate` / `superimpose` / `contacts`.
- `run_pipeline(structure, superimpose=True, db_dir=…)` → `PipelineResult`: import → annotate
  (alleles + chains + MHC groove) → superimpose onto the canonical DB (canonical Cα) → resmarkup
  + 5 Å contacts → per-interface energies. Scores: **TCRen** for TCR↔peptide, **MJ** for TCR↔MHC
  and peptide↔MHC, plus `total` (sum of the residue-pair potential over each interface's contacts).
- **Φ and ΔΦ with the TP / TM / PM breakdown come from this one command.** Columns `F_tcr_pep`,
  `F_tcr_mhc`, `F_pep_mhc`, `F_total`; `run_pipeline(…, reference_aa="A")` / `tcren scoring --delta`
  adds `dF_*` and `dF_total` (each interface's `tcren.ddg.reference_delta`; `dF_tcr_mhc` ≡ 0, the
  peptide is not in that interface). Use ΔΦ, not Φ, when each candidate carries its **own**
  generated pose. **Column names are shared with `tcren recognize`** — one vocabulary — but the
  **keys differ**: `scoring` emits `pdb.id`, `recognize`/`features` emit `complex.id`, so rename
  one before joining (`scoring --geometry` does that rename internally).
  (Manuscript notation: Φ for the energy, φ for a potential matrix entry, F for the
  binder-direction channel `−Φ_TCR:pep` in `cohort.f_score`.)
- `tcren scoring --geometry` appends the interface descriptors + `Q` by calling
  `recognition_table` + `cohort.q_score` — it does **not** reimplement them.
- `-s` accepts a file, directory, `.tar.gz`, quoted glob, `.txt`/`.list` manifest (one path per
  line, `#` comments), comma-separated list, or a repeated flag. `structure_paths` handles globs
  and manifests; `resolve_sources` does the comma/repeat splitting.

## Compiled extensions — `_align`, `_refine`, `_relax`, `_fold`, `_geom` (pybind11 / scikit-build-core)

- FIVE C++ exts, all in `CMakeLists.txt` (`pybind11_add_module` each, `install(TARGETS _align
  _refine _fold _geom _relax …)`): `src/_align/align.cpp` (MHC pseudoseq fitting alignment),
  `src/_refine/refine.cpp` (potential-guided peptide refinement), `src/_relax/relax.cpp` (DOPE
  interface energy for `tcren energy`/ΔΔG), `src/_fold/fold.cpp` (CCD loop closure) and
  `src/_geom/geom.cpp` (interface geometry for `tcren binder`; also `interface_clashes` and
  `contact_stability`, the clash + fragility kernels). Adding a sixth = same pattern; a new *function*
  in an existing module (as clashes/stability were added to `_geom`) needs no `CMakeLists.txt` change.
- The MHC-pseudosequence fitting-alignment hot path is a C++ ext (`src/_align/align.cpp`,
  `CMakeLists.txt`). Build backend is `scikit-build-core` (not hatchling); `uv pip install -e .`
  builds it once at install. Funcs: `fitting_score`, `best_hit` (GIL released over
  candidates), `align` (traceback). Scoring matches Bio.Align's fitting config EXACTLY (BLOSUM62,
  placed-gap open -11/extend -1, free target + end gaps), so `tcren.mhc.pseudo` falls back to
  Biopython transparently when the ext is absent. ~40 ms vs Bio 59 ms vs pure-Python 15 s for 4k
  candidates (a modest 1.5x — Bio's aligner is already C). `editable.rebuild = false` in
  pyproject: the ext is built once at `pip install -e .` (do NOT rebuild on import — that needs
  cmake on PATH at import time and breaks pytest/CI). CI (uv, `astral-sh/setup-uv`) installs
  `uv pip install --system -e . --no-deps` + explicit runtime deps (so arda-backed tests skip) and
  `uv pip install --system cmake ninja` to build the ext.

## Peptide substitution + refinement — `tcren.refine` (`tcren refine`)

- `substitute_peptide(structure, new_peptide)` — backbone-preserving identity swap on the peptide
  chain (keep N/Cα/C/O+Cβ, drop side-chain beyond Cβ); pure data-model, no atoms moved. Region
  markup is carried over onto the new residues — without it the contact map has null `pos.from/to`
  and every scorer fails. `score_peptides`
  is *virtual* (matrix lookup over the fixed contact map) — substitution is only needed to then refine.
- `refine_peptide(structure, restraint_w=0.5, …)` → `(structure, energy)`: knowledge-based rigid-body
  **Metropolis MC** of the peptide via the `_refine` C++ kernel. Energy = **DOPE** atom-level
  distance-dependent statistical potential (Shen & Sali 2006) over all peptide↔partner heavy-atom
  pairs (chain-agnostic, atom-class × distance lookup with linear interp) + **harmonic restraint to
  the input pose**. DOPE is used ONLY for refinement — deliberately decoupled from the TCRen/MJ
  potentials tcren *scores* with (no circularity: don't refine and score against the same quantity).
  DOPE's short-range bins are repulsive, so it supplies its own clash term (no separate clash). The
  restraint keeps it local (without it the search can drift to another favourable pocket). Partners =
  all non-peptide chains within a 12 Å shell. NOT physics MD; Rosetta FlexPepDock is the deferred path.
- DOPE data: `src/tcren/data/dope_potential.npz` (158 atom classes × 29 bins 0.75–14.75 Å), built by
  `scripts/build_dope.py` from the pymod/altmod MODELLER libraries; `tcren.refine._dope()` loads it.
- `tcren refine -s … -o … [--substitute PEP] [--steps N] [--restraint W]`. Native pose ≈ stays
  (RMSD ~0.2 Å); a buried/clashed peptide relaxes locally. Deterministic given `seed`.

## Clash detection + contact stability + register fix — `tcren.clashes`, `tcren.stability`, `tcren.refine.register`

QC for **generated** (AlphaFold/TCRmodel) complexes: their peptide-swap poses are routinely non-physical
(forced poses), which corrupts the contact energy the score reads.
- `interface_clashes(structure, tolerance=0.4, severe=0.6) -> ClashReport` — **native `_geom` kernel**
  (`_geom.interface_clashes`; the numpy path `_clash_pairs_numpy` is kept as the reference/fallback).
  Heavy-atom vdW overlaps (Bondi radii) between the peptide chain and its partners, broken down by
  partner `chain_type` (`by_partner`), with `n_clashes`/`n_severe`/`max_overlap`/`clash_score` + the
  worst residue pairs. `has_clashes(structure)` is the bool convenience. Self-check `python -m tcren.clashes`.
- `contact_stability(structure, cutoff=5.0, delta=1.0) -> StabilityReport` (`tcren.stability`) — TCR:peptide
  contact fragility read straight off the contact map (native `_geom.contact_stability`, numpy reference
  behind it). Per contact `margin = cutoff − dmin`; report has `n_contacts`, `mean_margin`, `frac_robust`,
  `frac_marg_lt1`, `exp_lost` (expected contacts lost under a `delta`-Å shift). A coordinate-only interface
  positional-confidence readout. Self-check `python -m tcren.stability`.
- `tcren recognize` appends both as output columns: `n_clashes`, `clash_score`, `exp_lost`, `mean_margin`,
  `frac_robust` (extra columns, **not** part of the 34 `RECOGNITION_FEATURES` the models consume).
- `check_register(model, reference=None) -> RegisterReport` — always reports the clash burden; with a
  correctly-registered `reference` (crystal / trusted pose) it adds the **anchor-Cα RMSD** in the
  MHC-groove frame (`peptide_rmsd`) → `wrong_register` (True/False; `None` without a reference).
  **A heavy clash burden alone is NOT a register call** — AF swap models are routinely clashy; register
  needs the reference. (The ila1 CPL lesson: forced pose → raw TCRen ROC 0.35, recovered to ≈0.77 on the crystal register.)
- `fix_register(model, template, engine="ccd") -> ModelResult` — re-threads the model's peptide sequence
  onto `template`'s correctly-registered backbone and re-refines via `model_peptide` (the
  FlexPepDock-functional path; needs the `_refine`/`_fold` kernel). Equal peptide length required.
  Template-free re-docking (CCD to canonical anchor targets) is the future extension.

## Poly-alanine reference score — `tcren.ddg.reference_delta` (geometry-normalized TCRen)

- `reference_delta(cm, peptide, pot, interface="tcr_peptide", reference_aa="A") -> float` = ΔΦ = Φ(peptide)
  − Φ(all-`reference_aa` peptide). It is the **full-peptide alanine scan** (== sum of `alanine_scan().ddG`)
  and subtracts the pose's identity-independent geometry baseline Φ(polyAla).
- **`interface="complex"` sums both peptide-bearing interfaces** (TCR:peptide with `potential`,
  peptide:MHC with `mhc_potential=`, default MJ), matching `cpl.response_matrix`'s cell convention.
  `ddg`, `neoantigen_ddg` and `reference_delta` all take it; CLI is `tcren ddg --interface complex`.
  **Use it to rank whole library peptides.** `tcr_peptide` alone is blind to presentation — on 1ao7
  the C-terminal anchor swap `LLFGYPVYV -> LLFGYPVYA` reads ΔΔG **0.0000** over `tcr_peptide` and
  **-0.9740** over the complex. The two effects are NOT separable in a library varying every
  position; report the per-interface terms beside the complex.
- **Use for GENERATED poses only.** On a *fixed* contact map ΔΦ = Φ − const → ranking unchanged
  (no-op); it differs only across candidates with their *own* structure (AF swap models). There it
  normalizes out the per-pose interface geometry: **rescues forced/wrong-register poses** whose geometry
  corrupts raw Φ (CPL ila1 TCR-ranking ROC **0.35 → 0.83**), but costs a little where the AF geometry is
  itself informative (mel5/mel8) — so it is a mode, not the default. Sweep: `scripts/ala_reference_sweep.py`
  (manuscript repo). It is **NOT** an affinity ΔΔG (dimensionless contact-preference, no free energy).
- **Ranking, not affinity (ATLAS).** Both raw Φ and ΔΦ are *within-receptor* peptide rankings: they order
  a peptide panel against one fixed TCR:pMHC (Garcia B*27:05 EC50, Spearman ρ≈0.75–0.9) but predict
  equilibrium binding *across* complexes only weakly — on ATLAS SPR raw Φ and ΔΦ track dG/Kd/koff/kon at
  |ρ|≤0.3 (cache-free recompute, `scripts/atlas_tcren.py` + `atlas_within_series.py`, manuscript repo).
  Turning ATLAS Garcia-like (thread a within-TCR panel on one pose) partly recovers the signal — proof
  it's the task setup (within-receptor ranking vs cross-complex affinity), not the potential.

## Potential decomposition + the MJ reference tables — `tcren.potential`

- **A contact matrix is not purely an interaction.** `Potential.decompose()` splits it exactly by
  double-centring: `e(a,b) = mean + H(a) + H(b) + J(a,b)`, `J` with zero marginals. `mean` and both
  `H` depend on one residue each, so an additive per-position model already absorbs them — **`J` is
  the only part a sum over positions cannot express.** Report `J`, not the raw sum, when the claim is
  "this is pairwise chemistry".
- `Potential.hydrophobicity_fit()` → `C0 + C1(q_a+q_b) + C2 q_a q_b` (Li–Tang–Wingreen, PRL 79:765,
  1997: the MJ matrix is nearly rank one). R² = **0.98** on `mj1996`, **0.85** on the bundled `mj`.
  Consequence to state when using such a potential: the interaction term is only `C2·q_a·q_b`, so it
  **cannot prefer one pair of side chains over another of equal hydrophobicity**.
- **Both refuse a directed potential.** TCRen is TCR→peptide; decomposing it is meaningless.
- `Potential.components()` returns that same split as **three scorable potentials**, because an
  interface score is a sum over contacts and the split carries through to it: `size` (the grand mean
  everywhere → `mean × contact count`, an *area* term), `comp` (`H(a)+H(b)`, degree-weighted
  composition) and `pair` (`J`, the interaction). Score a structure with each in turn to say which of
  the three a result is actually reading. **A matrix with no positive entries has a large negative
  mean, so its interface sum is dominated by the contact count** — that is how an interface-area
  effect ends up wearing a chemical name.
- `mj1996()` = MJ 1996 Table 3 raw contact energies (AAindex MIYS960101, transcribed). **Five of its
  210 unique pairs disagree with the AAindex record** by 0.04–0.28 (M–V, D–M, E–M, H–R, A–P; r =
  0.99978). Left untouched, pinned by `tests/unit/test_aaindex.py`; `aaindex("MIYS960101")` is the
  curated alternative. `mj_partition_energy()` = MJ 1985 one-body scale (MIYS850101, larger = more
  hydrophobic; **opposite sign convention** to a contact energy). Cross-check that validates both:
  r = +0.98 between the partition scale and the `q` axis recovered from the 1996 matrix.
- **The bundled `mj()` is Miyazawa–Jernigan 1999, identified 2026-08-29** — AAindex3 `MIYS990106`,
  400 of 400 cells exactly, runner-up off by 0.65. It is **not** 1985 and **not** 1996 Table 3, which
  is what the old "upstream table unrecorded" warning here said. `keskin()` is `KESO980101`, the
  solvent-mediated interfacial form, likewise 400/400 with the runner-up off by 2.77. Both files are
  left byte-for-byte untouched; what changed is that they can be cited. See `SOURCES`.
- **Reference state decides what a comparison measures.** `mj` (−0.079) and `betancourt` (−0.057) are
  *pair-contact* matrices with the one-body term removed; `keskin` (−3.547) and `mj1996` (−3.166) are
  *raw contact energies* that keep it. The like-for-like pairs are `mj` ↔ `betancourt` and `keskin` ↔
  `mj1996`; comparing across the groups compares reference states, not derivations.

### Every published contact matrix — `tcren.potential.aaindex`

`src/tcren/data/aaindex3.txt` is the whole AAindex3 flat file, bundled verbatim (47 records, 80 kB),
so adding a matrix to a comparison costs a string rather than a transcription.

- `catalogue()` → one row per entry with `kind`, `symmetric`, `n_missing`, `mean`, `min`, `max` and
  the citation fields. **Read `mean` for reference state** (above).
- `aaindex(acc)` → one entry as a `Potential`. It **refuses** the 2 contact-*count* tables
  (`TANS760102`, `MIYS960103`) and the 3 side-chain-*distance* tables (`BONM030104`–`BONM030106`),
  because scoring a contact map with a count table is a silent category error; `entry(acc)` still
  returns them deliberately. 42 of the 47 are usable energies.
- `identify(pot)` → accessions ordered by max |Δ|. **An identification needs an exact match AND a
  distant runner-up**; that is how `mj()` and `keskin()` were pinned down.
- Three `ZHAC*` entries are asymmetric by construction (row secondary structure vs column) and
  `decompose()`/`components()` refuse them.
- **AAindex's PMID field sometimes cites the paper that tabulated a matrix, not the one that derived
  it** (`MIYS850102` carries Bastolla 2001). Check the entry's own author/title/journal before
  citing.

## Contact-map Potts model — `tcren.potts` / `tcren potts`

Every other scoring path reads the contact map a structure **has**. This one models the map itself:
the random variable is the whole configuration `sigma`, one binary per *available* residue pair.

```
E(sigma) = - sum_a eta_a sigma_a - 1/2 sum_ab A_ab sigma_a sigma_b ,   P(sigma) = exp(-E)/Z
```

A **site** `a = (i,j)` is a receptor residue and a partner residue whose **Calpha** atoms lie within
`radius` (15 A); `sigma_a = 1` iff a heavy-atom contact formed within `cutoff` (5 A, the unchanged
TCRen definition). The availability mask is **Calpha-only on purpose** — a side-chain- or
type-aware radius would make the single-body field circular. That reference state is the point:
a TCRen potential is a Boltzmann inversion *conditioned on a contact existing*, so a residue that
could have reached the peptide and declined contributes nothing to it; here the non-event is the
observable.

```python
from tcren.potts import PottsModel, available_pairs, fit_potts, score_sites, contact_probabilities

pairs = available_pairs(structure)                 # or partner="mhc" (needs annotate_mhc first)
model = PottsModel.bundled()                       # or fit_potts(pairs, weights=...)
score_sites(pairs, model)                          # energy, log Z, log-likelihood, psi
contact_probabilities(pairs, model)                # p_independent / p_model / p_conditional
```

```
tcren potts fit      -s structures/ -o potts.json --balance both
tcren potts score    -s structures/ -o scores.tsv          # bundled model by default
tcren potts contacts -s complex.pdb -o contacts.tsv
```

**Fitting needs no partition function.** `P(sigma_a | rest)` is logistic in
`eta_a + sum_k K_k n_k(a)` with `n_k(a)` the count of contacting neighbours in coupling class `k`,
so the coupled fit is a weighted-binomial GLM with a few extra integer covariates and stays convex
(Besag pseudolikelihood; consistency arXiv:1010.0311, same recipe as plmDCA arXiv:1211.1281).
**Scoring does**, and its reference is exact: `ais_log_z` anneals only the coupling term from
`beta = 0`, where the model IS the uncoupled one and `log Z_0 = sum_a log(1+exp(eta_a))` in closed
form. **Always read `ais_ess`** — effective sample size out of `--particles`; small means the
schedule was too short.

**Gauge: penalise then project.** The design is over-parametrised and identified by an L2 ridge,
*then* projected to zero-sum. An L2 penalty silently picks its own gauge (plmDCA flags this), so
the projection is a separate step, never a constraint. After it, `J` is directly comparable with a
double-centred `Potential`.

**Three coupling families**, all on sequence offsets so no extra coordinates are needed:
`K(di,dj)` within one loop (`|di|,|dj| <= 2`, 12 classes), `L(|dj|, same chain?)` across two
hypervariable loops (6), and `M` for the same receptor residue against both partners (1, joint
models only). Everything else is asserted uncoupled.

**The result to know.** On the 362 alpha-beta Native2026 crystals every **axial** class is positive
and every **off-axis** class negative — `K(+1,0) = +0.79`, `K(0,+1) = +0.66` against
`K(+1,+1) = -0.82`, `K(+1,-1) = -0.81`. A made contact recruits its own sequence neighbours onto
the *same* partner residue and suppresses the diagonal one. **This is invisible in raw data**: the
diagonal offsets look positive unconditionally (odds ratio +1.43) and only turn negative once the
axial terms and the Calpha distance profile are held fixed. The couplings buy +505.7 nats of
pseudo-log-likelihood for 18 parameters.

**`--coupling-matrix` is the fair way to compare potentials.** It fixes `J` to one scale on a
bundled matrix, so competitors carry identical parameter counts and identical designs and their
pseudo-log-likelihoods compare directly; and because the matrix is double-centred it contributes
nothing to the one-body marginals, so the 40 field parameters absorb all composition content in
every fit and the comparison is about **pair structure alone**. Measured on Native2026, the ranking
**inverts** between interfaces: TCRen2 beats MJ by 103.3 nats on TCR:peptide, MJ beats TCRen2 by
35.5 nats on the TCR:MHC groove, and TCRen2's fitted scale falls 5.4-fold across that move
(+1.131 -> +0.209) while MJ's barely moves (+0.803 -> +0.974). That is the measurement behind
scoring `F_tcr_mhc` with MJ and reserving TCRen for TCR:peptide — **do not reuse TCRen2 on the
groove.**

**`pin_centred=False` is the gauge that reproduces a referenced score.** Double-centring is right
for *comparing* matrices and wrong for *reproducing* one. `reference_delta` is a difference of
one-body sums, so a centred pin re-injects `n_i * c(a)` — the position's contact count times the
potential's partner-residue column mean (s.d. 0.0668 on TCRen2, `n_i` from 1 to 54) — and the
identity fails. `fit_potts(..., coupling_matrix="tcren2", pin_centred=False)` pins the raw matrix,
so the coupling **is** the potential and any linear read-out of the field reduces to the potential's
own score up to the fitted scale.

**The free energy has three limits, and the interface picks one.** `dPhi/d eta_a = p_a`, so
`Phi = sum softplus(eta)` is an interaction sum weighted by contact probability, and a fixed contact
map is the `p in {0,1}` case. Hard contact reproduces `reference_delta`; smoothed (`p` free) is for
a *plastic* interface; saturated (`p` driven to 0/1 by the interface itself) is where a fixed map is
already exact. The groove is measurably in the saturated limit — over four 100 ns trajectories,
36 of the 38 engaged peptide positions (94.7%) reach a maximum groove-pair contact frequency above
0.98, against 14 of 38 (36.8%) on the receptor side. See `docs/potts.rst`.

**Two bundled models**, `PottsModel.bundled("potts_tcr_peptide")` (the default, 64,622 sites /
7,865 contacts) and `"potts_tcr_mhc"` (239,093 / 15,451 — the TCR makes twice as many contacts with
the MHC as with the peptide). Both carry the alpha-beta HARD RULE, as `derive-potential` does.

**Bound versus unbound** (`bound_unbound`, `count_profile`, 2.14.0). `E(empty) = 0`, so the scores
`score_sites` already emits decompose exactly:

```
-E(sigma_obs)  =  log Z  +  L(sigma_obs)
[binding log-odds]  [capacity]  [typicality]
```

`bound_unbound` gives three readings of the whole-interface two-state contrast from ONE Gibbs pass:
`df_empty` = `log(Z - 1)`, exact; `df_threshold` = `log[P(N>=x)/P(N<x)]`, in which `Z` cancels so no
AIS is needed; and `mu_star`, the chemical potential at which the model's mean contact count matches
the observed one (`nan` outside the sampled support — that is not a failure, it is refusing to
extrapolate). They are not competing estimates: no sampler reaches `N = 0` for a docked pose, so
`df_empty` must come from `log Z`. `count_profile` gives the pooled `F(N) = -log p(N)` landscape, so
a threshold is read off it rather than assumed.

**Constraining a statistic of the whole configuration.** `gibbs(..., observer=fn)` calls `fn` on
each kept draw with the `(chains, n)` configuration matrix — that is the hook for accumulating the
kind of statistic a Lagrange multiplier couples to (Jaynes; Tkacik et al. 2014's K-pairwise model is
this with `O` = the total activity). Because such a statistic depends on sigma only through a
low-dimensional summary, one sampling pass serves every step of a moment-matching fit: reweight with
`tilt_mean` rather than resampling. A LINEAR tilt in `N` is exactly a constant added to every field,
`E - mu N = -(eta + mu).sigma`, which is what makes the reweighting checkable against direct
simulation.

**Gotchas.** `partner="mhc"` needs the groove regions, which chain typing alone does not assign —
run `annotate_mhc` / `annotate_mhc_batch` first or every MHC residue is silently dropped for want
of a within-region coordinate. Distance-bin edges are **global** (base 0), never derived from a
frame's own minimum, or a model fitted on one set indexes the wrong coefficients on another.
`psi` (log-likelihood per available pair), not `log_lik`, is the column to compare across
interfaces of different size.

## Ring stacking — `tcren.stacking.ring_stacking` (geometry, NOT an energy)

- A contact potential scores a pair by identity alone, so two rings face-to-face at 3.5 Å score the
  same as the same two residues brushing past edge-on. `ring_stacking(source, cutoff=7.5,
  min_seq_sep=1)` measures what identity cannot carry: `centroid_distance`, `interplanar_angle`
  (0 = face-to-face, 90 = edge-to-face), `vertical` (gap between planes) and `lateral` (sideways
  slide). Parallel-displaced stack = small vertical + a few Å lateral.
- **Pro is in `RING_ATOMS`** although not aromatic — its pyrrolidine ring packs face-on against
  aromatics via CH–π, which is the interaction the module exists to measure. Trp is represented by
  its six-membered ring (the face that stacks).
- Returns **no energy** by design. It says the rings are or are not arranged the way a stack is;
  valuing that is left to whoever has a potential that can.

## Intra-peptide term — `tcren.intra_peptide_energy` (the contacts a chain makes with ITSELF)

- Every Φ in the package sums over contacts between two **different** chains, so a peptide held in its
  bound conformation by its own side chains scores the same as one that is not. This term is that
  omission, and it is **off everywhere by default** — `scope="inter"`, `peptide_internal=False`,
  `intra_weight=0.0` all reproduce the previous output byte-for-byte.
- `all_atom_contacts(..., scope="inter"|"intra"|"all")`; `peptide_internal_contacts(structure,
  cutoff=5.0, min_seq_sep=3)` is the intra case plus the separation filter. Same 5 Å as everywhere
  else, so an internal and an interface contact mean the same thing; `|i−j| ≥ 3` is what does the
  filtering, because sequence neighbours touch by covalent geometry — over `tests/assets/pdb` the
  totals are 18 contacts at `|i−j| ≥ 3` vs 134 at `|i−j| ≥ 2`, the jump being `i`/`i+2` pairs.
- `intra_peptide_energy(cm, pot, peptide=None, contact_weight="residue") -> float`; needs
  `ContactMap.from_structure(..., peptide_internal=True)` (which stores the pairs *beside* `contacts`,
  never in it). `peptide=` threads a candidate onto the structure's peptide positions.
- **The potential is symmetrised**, `(F + Fᵀ)/2`. An intra-chain pair has no from/to orientation, and
  which residue lands on which side is the contact table's canonical `(chain.id, residue.index)`
  ordering, not chemistry. Matters for TCRen (directed), no-op for MJ. **Default MJ**: TCRen is derived
  from TCR↔peptide contacts and says nothing about a chain's contacts with itself.
- Wired in at three layers: `score_peptides(..., intra_weight=w, intra_potential=)` /
  `tcren score --intra-weight` (score = Φ + w·E_intra); `run_pipeline(..., intra_weight=w)` /
  `tcren scoring --intra-weight` (reports `F_pep_int` raw, folds w·it into `F_total`; potential
  overridable via `potentials={"peptide_internal": …}`); `tcren recognize --full` (`F_pep_int`,
  `n_pep_int`, catalogued `involves_tcr=False` — the peptide's own conformation is cohort identity).
- **Expect a sparse term.** A canonical extended class-I 9-mer makes 0–2 internal contacts, so it only
  separates candidates where the peptide is genuinely bulged or self-packed. Untested as a ranking
  signal — it is exposed so the assumption can be measured rather than inherited.

## Interface mechanics — `tcren.mechanics` (koff/kinetics, NOT ΔG)

- The TCR↔pMHC contact map as a network of breakable Cα-anchored Hookean springs (per-contact
  stiffness from heavy-atom-pair multiplicity; default `weight="invdist2"` = multiplicity/dist²).
  Pure-numpy, single structure, no MD. Public API:
  - `interface_springs(structure, cutoff=8.0, weight='invdist2') -> InterfaceSprings(a, b, k, rest, axis)`
    — the TCR-side/pMHC-side Cα anchors, spring stiffnesses, rest lengths, and the unit docking axis
    (TCR→pMHC). Raises if no peptide chain.
  - `stiffness_tensor(structure, cutoff, weight) -> dict` — linear-response `K = Σ kᵢ ûᵢ⊗ûᵢ`:
    `S_tot` (trace), `K_tens` (along docking axis), `K_shear` (`S_tot − K_tens`), `aniso`
    (`K_shear/K_tens`), `lam_max`/`lam_min`, `n_spring`. All `nan` if < 3 springs.
  - `rupture(structure, direction='tensile'|'shear'|'auto', cutoff, weight, break_strain=0.5, steps=80) -> dict`
    — steered-unbinding cartoon (rigidly pull pMHC off, break springs past strain): `rupture_force`
    (peak resisting force), `rupture_work` (∫ force·displacement), `n_spring`, `break_strain`.
    `"auto"` = min-force of tensile/shear.
  - `coupling_residues(structure, cutoff=5.0) -> dict` — residues in both an intra-body scaffold
    contact and the interface: `couple_pep`/`couple_mhc`/`couple_tcr`, `couple_total`, `n_interface`.
  - `interface_mechanics(structure, ...) -> dict` — **the one definition of "the mechanics row"**:
    the union of the three above under their shipped defaults. `tcren mechanics` and
    `tcren recognize --mechanics` both go through it, so the two agree by construction.
- **For a cohort, prefer `tcren recognize --mechanics -t 0`** over running `recognize` and
  `mechanics` as two commands. Both need the same annotated structure, so the second command repeats
  the parse and both mmseqs searches and returns a second table (CSV, keyed `pdb.id` rather than
  `complex.id`) that then has to be joined; the flag reuses the annotation and costs only the
  mechanics arithmetic. The structure must be MHC-annotated either way — an unannotated MHC chain
  silently empties the TCR:MHC half of the spring network (this is what the 2026-07-28 sharding bug
  did: `couple_mhc` = 0 in 523/523 rows, 67 springs against 279).
- **Caveat: these track the dissociation off-rate koff / kinetic stability (Bell–Evans rupture
  resistance ~ r0.5 on ATLAS), NOT the equilibrium ΔG/Kd** — rupture reflects the dissociation
  barrier, not the well depth (physically apt for the TCR mechanosensor / catch bonds). Use them
  **between structures** (one value per complex) to rank/compare; do not pool many per-residue or
  per-spring rows from one structure as independent samples — that is pseudo-replication.
- Self-check (no PDB): `python -m tcren.mechanics`.

## Docking geometry — `tcren.orient.docking` + `tcren.orient.tcrdock_geometry`

- **Two interpretable angles** (existing): `docking_angles(structure) -> DockingAngles(crossing_angle,
  crossing_angle_signed, incident_angle, ...)`. `crossing_angle` = the groove-plane "scanning" angle,
  `incident_angle` = the tilt. Computed from the Vα→Vβ axis in the groove frame; no reference DB.
- **Full rigid-body pose** (new, `tcrdock_geometry.py`): `docking_geometry(structure) -> DockingGeometry(d,
  torsion, tcr_unit_y, tcr_unit_z, mhc_unit_y, mhc_unit_z)` — native reimplementation of **TCRdock**
  (phbradley/TCRdock, MIT, commit `c5a7af4`; see `THIRD_PARTY_NOTICES.md`). MHC + TCR symmetry stubs (β-sheet
  floor / Vα-Vβ two-fold), MHC-I core by BLOSUM-align to TCRdock's template, TCR core by conserved IMGT
  framework positions from arda region markup. Needs a chain-typed + MHC-annotated structure; **class-I only**
  (class-II raises).
- **What we use / provenance finding (2026-07-05, validated on 618 TCRvdb models):** the upstream AF/TCRmodel2
  annotation table's `scanning_angle` **is** reproducible (= `crossing_angle`, r≈0.88), but its
  **`pitch_angle` is NOT** any clean geometric angle (best correlate `d`, r≈0.42) and out-discriminates every
  clean docking feature on TCRvdb (macro-PR≈0.72 vs `d`≈0.64 / `torsion`≈0.62 / tilt≈0.58) → its extra signal
  is **AlphaFold-confidence contamination, not geometry**. Prefer the documented `d`/`torsion`/`crossing`/
  `incident` over the opaque upstream `pitch_angle`.

## Publication figures — `tcren.viz.pymol` (headless PyMOL)

- `render(scene, png, size=(1200,1200), dpi=300, corner='bottom-left', gizmo=True)` ray-traces a
  PyMOL script body. Scene presets: `overlay_scene(ids, canon_dir, limit=8)` (ensemble, side-on),
  `groove_scene(pid, canon_dir, surface=False)` (peptide in the cleft, top-down; `surface=True` is
  the histo.fyi look), `interface_scene(pid, canon_dir, cdr_resi)` (peptide + CDR loops).
- **`CANONICAL_AXES` is the naming of `tcren.orient.frame`, and the gizmo prints it**: x=`width`
  (groove width, α1↔α2, PC3), y=`N→C` (groove axis toward the peptide C-terminus, PC2), z=`TCR`
  (docking normal, MHC floor→TCR, PC1). Same three directions as SwiftTCR / TCR3d; the PC *ranking*
  differs because `orient.frame` fits the whole complex and they fit the MHC groove alone.
- The gizmo is a **separate render pass composited at pixel coordinates**, not a CGO projected into
  the corner. PyMOL's orthoscopic viewport does not span `dist * tan(fov/2)` — measured, it is out
  by ~25% — so the projection route puts it off-frame. Do not "fix" this by reintroducing the math.
- Two PyMOL gotchas the module already handles: `auto_zoom` is **on** by default, so loading a CGO
  or a label pseudoatom re-frames the camera onto it (the module sets it off first); and
  `cmd.get_view()[0:9]` is **column-major**, so the transpose that inverts it is `rot[i*3+j]`.
- **`residue_importance(structure, interface='tcr_peptide') -> pl.DataFrame`** — the exact
  decomposition of Φ: `chain.id, residue.index, residue.aa, region.type, n_contacts, phi`, sorted
  most-favourable first. `importance_scene(pid, canon_dir, imp, by='phi'|'n_contacts')` paints it
  (sticks on a ramp, everything else pale). **The phi ramp is centred on zero**, so blue/red mean
  favourable/unfavourable rather than less/more — a range-fitted ramp reddens the least-favourable
  residue even in an all-stabilising interface. **Each contact is attributed to BOTH residues it
  joins, so the per-residue phi sums to 2x the interface total** — an attribution, not a partition;
  a test pins that factor of two.
- Chain roles after `tcren orient`: A=Vα, B=Vβ, C=peptide, D=MHCα, E=MHCβ/β2m — `CHAIN_COLOURS`.
- Notebooks: `notebooks/pymol_canonical_figures.ipynb` (static gallery) and
  `notebooks/pymol_interactive.py` (**marimo**, `pip install 'tcren[marimo]'`). marimo is reactive,
  so a name may be defined by exactly one cell — import shared symbols once in the setup cell and
  thread them through the signatures, or the export fails with `MultipleDefinitionError`.
- Docs gallery page with rendered examples: `docs/gallery.rst` + `docs/_static/gallery/*.png`.
- Everything except `render`/`probe_rotation` is pure Python and unit-tested without PyMOL
  (`tests/unit/test_pymol_viz.py`); the two that shell out are `slow` + skip when pymol is absent.

## Contact typing — `tcren.contact_types` (two schemes; v1 is frozen)

- **`scheme="v2"` (default)**: `salt_bridge, hydrogen_bond, cation_pi, stacking, aromatic, hydrophobic,
  polar, vdw, other`. `other` now means only "too far to be anything", never "unrecognised" — it fell
  from **72.3% to 13.9%** of TCR:peptide contacts on five crystals. A contact carries *several* types
  (`is_<type>` booleans are independent); `contact.type` is only the top-priority label.
- **`scheme="v1"`** is the old five-type residue-level scheme, kept byte-for-byte. `recognition.py`
  **pins it** — the frozen classifiers were fitted on its `ct_*` counts. Do not change v1.
- What v2 fixed: apolarity is decided per *atom* (a carbon with no bonded N/O), not per residue —
  v1's residue set excluded Tyr, the commonest TCR interface residue. H-bonds reach 3.9 Å with
  donor/acceptor typing (two carbonyl oxygens are no longer an H-bond). `stacking.ring_stacking` is
  finally joined via `stacked_pairs(structure)`.
- `residue_pair_types(structure, interface)` types from **every** atom pair, not the closest one —
  a salt bridge whose nearest contact is two carbons is invisible otherwise. `type_weights(typed)`
  gives 0/1 weights that drop pure-proximity contacts (`tcren score --drop-untyped`).
- **Hydrogens are now filtered** in `all_atom_contacts`. This changes contacts and energies for
  H-bearing depositions (5jhd: +7 of 28 contacts, −58.5% F_tcr_pep) and breaks legacy-oracle parity
  on 5jhd/7qpj, recorded as a subset relation in the regression test.

## Feature families — `tcren features` (descriptors) vs `tcren recognize` (scores)

- **Two commands, two jobs.** `tcren features` reads structures and writes descriptors;
  `tcren recognize` turns a descriptor table into scores. The feature pass is the expensive half,
  so run it once and re-score for free:

  ```bash
  tcren features  -s <in> -i placement,interface,topology,energetics -o feats.tsv
  tcren recognize --features feats.tsv -o scores.tsv
  ```

  `scores.tsv` is `complex.id` + **`Q`** (fit-free interface quality), the three channel posteriors
  **`G`** / **`T`** / **`E`** (geometry / topology / energetics), **`P_native`**, and **`S_free`**
  with its calibrated **`p_binder`**, and nothing
  else. `tcren recognize -s <in>` is the *other* mode: it reads structures and writes the
  descriptor table with `p_real`, not `Q`/`P_native`.
- **Five families, split by invariance** (`tcren.recognition.DESCRIPTORS`, `FAMILIES`); a
  *family* is a slice of the descriptor table, a *channel* is one of the three networks `P_native`
  sums, and `P_NATIVE_POOL` is the map between them:
  `placement` (groove-frame pose — angles, TCRdock params, ride height/shift/offset, CDR3 frames;
  frame-**dependent**), `interface` (contact size + chemistry), `topology` (the *shape* of the
  contact set, size-free), `energetics` (Φ and ΔΦ), `kinetics` (spring network; off by default).
- `placement` + `interface` were one `geometry` family and `energetics` was `physics` until
  2026-08-24. Both retired names still resolve in `descriptors()`; the split is what lets the
  independence claim be stated — measured on VDJdb, topology ⟂ interface at |ρ| = 0.023 while
  topology–placement is 0.177 (0.448 on TCRvdb), because uniform coverage *is* ride height.
- **Only the requested families are computed.** `-i topology` never builds the energies.
- **Contact counts are `interface`, not `topology`** (`FOOTPRINT_SIZE_FEATURES`). A shape channel
  carrying the interface's size would correlate with the interface channel by construction.
- `P_native` (`tcren.cohort.p_native`) combines **three** channels — `geometry`, `topology`,
  `energetics` — each fitted as its own latent-class Bayes network by EM
  (`GaussianBNClassifier.fit_em`), their log-odds added. **No binder label enters.** EM learns each
  channel's sign, which is what makes the measured coupling `C*` unnecessary: on a cohort whose
  contact energy runs backwards the energetics coefficient simply comes out negative.
  - **`geometry` pools the `placement` and `interface` FAMILIES into one network** (`P_NATIVE_POOL`).
    Adding log-odds is the exact posterior only across channels that are conditionally independent
    given the class, and those two are the most dependent pair measured (|ρ| = 0.244). Pooling them
    is what the assumption requires; summing them as two terms counts the dependence twice, and
    measurably: the four-channel sum reads 0.817/0.668 (TCRvdb/VDJdb ROC) against the pooled
    three-channel 0.832/0.718.
  - `rule="flat"` instead pools every channel's features into one network. It holds the top of a
    ranking better where cohorts are small (VDJdb P@10 0.872 vs 0.812) and ranks worse overall
    (0.689 vs 0.718). Both are reported in the paper; neither dominates.
  - **`T` is just `p_native(channels=("topology",))`** — the shape channel read on its own. It is
    what replaced the hand-written `fp_score` z-sum.
  - EM is monotone **only with the DAG fixed**, which is the default; `relearn_structure=True`
    changes the model family between rounds and the likelihood can fall.
  - A mixture is identified only up to permutation — `orient_by` is what stops the two components
    swapping between runs, and `P_NATIVE_ORIENT` gives the per-channel default. **A leading `-`
    means lower-is-native**; no shipped channel needs it, because each orients on a column that
    already runs binder-upwards (`burial`, `D2_pep24`, `neg_energy` — the last is `-E`, so higher
    is more favourable). Orienting the energetics channel on a raw Φ column *would* need the `-`.
  - **Anchors are optional and semi-supervised** — `{row_index: 0|1}` over the caller's own rows,
    pinned at every E-step; those rows stay in the design matrix and are still scored. The default
    is `anchors=None`, a fully unsupervised fit. Anchoring a row you then score reads the label
    back out: an early draft did that and reported 0.83 where the honest number is 0.69.
  - Keep the feature count small: the BIC hill climb is quadratic (0.01 s at 18 features on 618
    rows, 1.7 s at 40, **45 s at 89**). `P_NATIVE_FEATURES` is the compact default, keyed by
    FAMILY; `_channel_columns(channel)` resolves a channel through the pool.

## Single-structure reliability — `tcren.reliability` / `tcren assess`

**`P_native` is not the score to ship.** `cohort.p_native` refits a latent-class model on every
call and **raises when a cohort has fewer rows than features**, so it is undefined for one structure
and its value depends on what else was scored alongside it. Use `S_free`.

```python
from tcren.potts import score_sites, available_pairs, PottsModel
from tcren.reliability import s_free, p_binder, af_band

potts_scores = score_sites(available_pairs(structure), PottsModel.bundled())  # writes neg_energy
v = s_free(feature_table, energy=potts_scores["neg_energy"])   # energy optional
p = p_binder(v, link="binder_bm|S_nat")                        # frozen out-of-fold Platt
b = af_band(meta["iptm"], reference="binder_bm|ipTM")          # the generator diagnostic
```

`S_free = Q/sd_Q + T/sd_T + (Pi - mu)/sd_Pi`. Three fit-free directional blocks `z(x)' C^-1 s` over
the Native2026 crystals, each divided by its own native spread.

**The outer transform is a DIVIDE, not a z.** A block score's native mean is 0 by construction, so
re-centring is a no-op; its variance is `s' C^-1 s`, which is not 1 — measured 1.43 (`Q`), 1.61
(`T`), 14.13 (`Pi`). Without the division the energy would carry ten times the weight of the
geometry. Equal weight *in native-sd units* is the claim.

`Pi` is `neg_energy` from `tcren.potts`, the interface energy read against the partition function
rather than a poly-alanine reference. It is the frozen choice because it is the **least redundant
with `Q`** of the five ways of spending `-E = log Z + L`: native Pearson +0.33, against +0.75 for
the contact count and +0.51 for `log Z`. Pass `energy=None` and the two-block form is returned —
still defined, and reported as such rather than imputed.

`t_score` is new and is the block that **survives without a template**: on the balanced VDJdb panel
`T` loses 0.06 ROC-AUC where the epitope has no solved complex, against `Q`'s 0.24. It needs
`q_score(..., signs=T_SIGNS)`, because the footprint's connected-component fraction at 7 A runs the
other way from the rest.

**Calibration is per-link and the name is a contract.** `p_binder`'s links were fitted out of fold
and each expects the score its name carries; handing a raw `S_free` to a `min rank%(...)` link is a
category error, not a rescaling. `available_links()` / `available_bands()` list what ships.

```
tcren features -s models/ -i placement,interface,topology,energetics -o feats.tsv
tcren potts   score -s models/ -o potts.tsv        # writes neg_energy, keyed on pdb.id
# join potts.tsv's neg_energy onto feats.tsv on the structure id, then:
tcren assess --features joined.tsv -o assessed.tsv
```

**Gotcha.** `tcren features` does not compute `neg_energy` — it comes from `tcren potts score` and
has to be joined. Skip the join and `assess` emits the two-block `Q + T` form; it says so in its
report, so read the first line rather than assuming three blocks.

**The forced-pose flag — why a confident model can be wrong in a readable way.** A generator pushed
into a confident but wrong pose does not build a random interface. To seat the chains it picks
residue pairs it believes are favourable, so the recognition energy comes out *good* — often better
than a genuine complex of the same epitope. The energy **inverts** under forcing instead of
degrading: on a 24-cohort forced-pose panel (1,707 structures) it reads macro ROC-AUC **0.4952** and
is below 0.5 in **15 of 24 cohorts**, where ipTM reads 0.6093.

```python
from tcren.reliability import inversion_flag, screening_yield
f = inversion_flag(feature_table, energy=potts_scores["neg_energy"])  # NaN without the energy
y = screening_yield(v, budget=0.10, prevalence=0.48)   # what testing the top 10% implies
```

`inversion_flag` is the energy block minus the mean of the two shape blocks, in native-sd units.
Large positive = the energy is vouching for a structure the footprint does not, which is the pattern
to distrust; a generator fakes favourable contacts far more easily than a well-formed footprint.
It ranks and triages — it is not calibrated, `p_binder` is.

`screening_yield` returns the cut only: how many structures, at what threshold and percentile, plus
`expected_hits` under a stated prevalence. **Enrichment over random is not returned** — it needs
labels the function does not have, and a NaN there would read like a measurement.

`assess` emits three blocks: reliability (`S_free`, `p_binder`), ranking within the set (rank,
percentile, expected precision at a recall budget), and — when the table carries ipTM — the
generator diagnostic (`af_band`, `p_nonbinder_af`, `s_free_roc_in_band`). The last column is the
actionable one: on the balanced VDJdb panel the **top ipTM decile is 26.2 %** [18.7, 35.5]
**non-binders**, and is also the band where `S_free` reads highest.

## Footprint shape — `tcren.footprint` / `tcren footprint`

- **Reach for it when the energy is at chance but the pose still looks wrong.** It reads the same
  contact map as a *shape*, not a sum: how evenly the six CDR loops spread their contacts, whether
  the germline/CDR3 division of labour holds, and whether the footprint is one connected patch.
  No potential, no reference structure, no fitted parameter.
  CLI: `tcren features -s <in> -i topology -o <out>`, then `tcren recognize --features <out>`
  for the cohort-standardised shape posterior `T`. (`tcren footprint` is the same code path,
  now hidden and superseded; its `--score` no longer emits the removed `fp_score` z-sum but the
  same `T` posterior, so there is no reason to prefer it.)
- **The MHC pass must run AFTER chain typing, and it is not optional.** `classify_chains` leaves an
  MHC chain typed generically `"MHC"`; `interface("tcr_mhc")` matches the supertype `annotate_mhc`
  assigns. Skip it and six of the twelve cells are structurally unreachable with no error —
  `p_germ_mhc` reads 0.06 instead of 0.78 and `H_cell` is normalised by ln 12 over a partition half
  of which can never be occupied. `cell_counts` now warns; that bug shipped once.
- **NaN must become null before any polars aggregation.** polars propagates NaN through
  `mean`/`std`, so one contact-free structure turned a whole z-scored channel to NaN, `fill_nan(0)`
  flattened it, and `fp_score` silently became its other channel alone (0.815 → 0.691 on TCRvdb).
- **No canonical orientation is needed** — every feature is invariant under rigid motion, which
  `test_footprint.py::test_every_feature_is_invariant_under_rigid_motion` pins by rotating and
  translating a complex and demanding an identical row. Only chain typing + CDR markup. MHC *region*
  markup is **not** used, so the "MHC needs a second pass or it silently empties" trap does not apply.
- **Annotation goes through `_iter_typed`/`iter_annotated_set`** — one mmseqs call per organism for
  the whole set. Do not call `classify_chains` per structure here (that is what `tcren surface`
  does, and it is an order of magnitude slower over a cohort).
- `footprint_features(s) -> dict` (29 features at the default two radii), `footprint_batch(paths_or_structures) -> pl.DataFrame`.
  For the cohort-standardised shape score use `cohort.p_native(table, channels=("topology",))` (`T`);
  the old `footprint_score` z-sum was removed at 2.12.
- **Coverage**: cells are the 6 CDR loops × {peptide, MHC} (12) or with the peptide split into
  thirds (24). `H_cell` is the normalised Shannon entropy, `D1`/`D2` the Hill numbers — `D2`, the
  *effective number of engaged cells*, separates better because it discounts the rare cells a decoy
  populates weakly. Refining the **peptide** side helps; refining the MHC side into helices/floor
  does not, which is why that partition is not offered.
- **Topology**: `_flag_betti` builds the flag complex on the contacted pMHC Cα at a radius —
  `fp_b0_*` patches, `fp_b1_*` holes (via a GF(2) rank of the triangle boundary), `fp_chi_*`.
  `h0_pers_ent` is the H₀ persistence entropy, whose barcode **is** the MST edge lengths, so no
  filtration library is needed. `b0` is most informative at 7 Å and `b1` at 8 Å; both ship.
- **The bipartite contact graph's cyclomatic number `E − V + C` is deliberately not offered.** With
  ~30 contacts among ~30 residues it is dominated by `E` and just tracks interface size; the patch
  count is scale-free instead. If someone asks for "the Betti number of the interface", this is the
  distinction to make.
- `n_loop_contacts` counts only what the partition sees — the six CDR loops. Framework contacts are
  excluded by construction, so it is smaller than the full interface count. The topology features
  are **not** restricted this way: they use every contacted pMHC residue.
- **It is not `n_contacts`, and the distinction is load-bearing.** `n_contacts` belongs to the
  `potts` family — the available pairs that engaged, 29 on 1ao7 against this module's 66 — and it
  is the column `tcren diagnose` standardizes against the frozen Potts moments. Through 2.19.0 both
  passes wrote the one name and the later pass won; a stale table carrying the loop tally as
  `n_contacts` is now refused by `reliability.correct_confidence` rather than corrected.

## Surface topology — `tcren.surface` / `tcren surface`

- **Reach for it when the question is about the pMHC alone**, with no TCR in the structure or before
  one is docked: how featureless is this epitope, which epitopes present a similar face, is the
  peptide bulged. A contact potential cannot answer any of these — it needs an interface that exists.
  CLI: `tcren surface -s <in> -o <out> [--compare] [--svg <dir>] [--cells <csv>] [--channel h]`.
- `surface_map(structure) -> SurfaceMap` — height + hydropathy + charge on a 64×32 raster over the
  groove; `surface_stats` gives `relief`, `peak_to_valley`, `frac_above_ridge`, `phobic_centre`;
  `surface_distance`/`surface_tree` compare epitopes (SURFMAP Manhattan distance + linkage).
- **The frame is refit per structure** — z from the groove-floor plane normal, **y from the peptide**
  (the floor's own principal axis is NOT the groove axis), origin on the peptide centroid. So maps
  compare without prealigning inputs, and `force_pca`/canonicalisation are not prerequisites.
- Heights come from **ray casting in the groove frame**, not Shrake-Rupley points: sphere sampling is
  fixed in global axes, so the map wobbled 1.35 Å (median cell) under a rigid rotation of the input.
- Validated on all 374 Canonical2026: `frac_above_ridge` 0.054 (8-mers) → 0.569 (13-mers); same-epitope
  maps closer than cross-epitope at P = 0.917.
- **Both faces, one frame.** `surface_map(structure, side="tcr")` rasters the receptor's V-domain
  *underside* instead of the groove — same frame, same grid, lowest point per cell instead of
  highest — and `surface_complementarity(pmhc_map, tcr_map)` compares them cell for cell:
  `shape_r` (positive = complementary), `charge_r`/`charge_product` (**negative** = complementary,
  plus meeting minus), `phobic_r`/`phobic_product` (positive), `gap_mean`/`gap_sd`, and the
  per-channel map distances `d_h`/`d_charge`/`d_phobic`. CLI: `--side tcr`, `--complementarity <csv>`.
- **The two limits on that comparison are calibrated, not assumed** (60 Native2026 crystals).
  `DEFAULT_EXTENT` is sized for a class-II 15-mer, far wider than any receptor footprint: over the
  full extent a TCR projection reaches only **0.741** of occupied pMHC cells *at any Z cutoff*, and
  the shortfall is nearly all at the far groove end (0.348 beyond y = +15 Å against 0.987 near
  y = 0). So `COMPARE_WINDOW = (12, 12)` Å crops first, and `MAX_GAP = 10` Å is where the Z curve
  has gone flat (0.895 at 4 Å → 0.951 at 10 Å → 0.962 uncapped). Peptide-owned cells are covered at
  **0.917** even over the full extent, so cropping discards no epitope surface.
- **The gap is mostly negative and that is correct**: median `h_tcr − h_pmhc` = −1.7 Å with 71% of
  cells interdigitated, because the two faces interlock rather than stack. The Z cutoff is one-sided.
- Figures: `viz.surface2d.render_surface_map(smap, channel)` → SVG string (hand-built, zero deps).

## Parsing speed — the PDB fast path

- `parse_structure` takes a **vectorised ATOM-record path** for `.pdb` at `model=0`: the records are
  sliced as one uint8 array and residue breaks found with `np.flatnonzero`. **3.3x** end to end
  (22.1 -> 6.6 ms per crystal), and it is what makes a dataset-scale pass tractable — Biopython was
  86% of the wall clock (19.0 s of 30.3 s over Native2026, building 2.1 M Atom objects to discard).
- **Exact or it does not run.** `_parse_pdb_fast` returns `None` — falling back to Biopython — on a
  blank element column, a short ATOM line, or an unparseable coordinate. Coordinates are read
  **through float32**, matching Biopython bit for bit rather than being more accurate: a PDB's three
  decimals come back as `59.42599869`, and every previously computed number depends on that.
  Verified by re-deriving all three potentials of `pot_realmock.py` from structures: max |delta| =
  0.0000000000 over 1,140 cells.
- The remaining cost is ~655 k `Atom` constructions per 120 structures. That is the data model, not
  the language — **a C++ extension would not help without making residues array-backed**, which is a
  redesign, not an optimisation.

## Rotamer-averaged contacts — `tcren.rotamers`

- `contact_probabilities(structure, interface) -> df[..., p]` and `soft_energy(structure, potential)`.
  Replaces the 0/1 contact indicator with `p_ij` averaged over χ rotamers, Boltzmann-weighted under
  **DOPE** (never under the potential being scored — that would be circular).
- Why: under a deliberately wrong χ1 the hard contact set keeps Jaccard 0.66 and the energy moves by
  |ΔΦ| = 0.524; the averaged map keeps 0.95 and moves 0.054, against energies of magnitude 0.4–2.2.
- Rotating everything past Cβ about Cα–Cβ **is** χ1 exactly. `max_chi=2` default (3^n rotamers).
  ~0.24 s/structure — see `refine/CPP_REWRITE.md` for when this needs to be C++ (MC loops, not one-shot).

## Backbone dynamics / peptide conformational stability — `tcren.dynamics` (the Sewell hypothesis)

- **Library only — there is no `tcren stability` command.** Call `peptide_stability` /
  `stability_table`, both exported from the top-level `tcren` namespace. Batch driver, with the CPL
  design and the AUC bookkeeping: `scripts/sewell_stability.py`.
- `peptide_stability(structure, intra_weight=1.0) -> Stability`: flexible-backbone Metropolis MC of
  the peptide's φ/ψ against DOPE. The readout is **not a better pose** — it is `rmsf` (ensemble
  spread), `drift` (how far the mean moves from the input) and `energy_gap`.
- **Why it exists**: a contact potential scores whichever conformation it is handed and cannot tell a
  peptide held in the TCR-facing conformation by its own side chains from one merely drawn there.
  Sewell (`suggestions/sewell.txt`): intra-peptide interactions stabilise the productive bulge, and
  "poor binders could perhaps still make many contacts but fail to stabilise" it.
- **`intra_weight` is the switch.** `1.0` includes the peptide's DOPE contacts with itself, `0.0`
  removes them; `stability_table` runs both at the same seed and reports `delta_rmsf`.
- Moves are torsional and exact — perturb one φ/ψ, rotate everything downstream. Both peptide termini
  are free so no loop closure is needed. `backbone_torsions` is the tree; **φ moves the side chain**
  (Cβ hangs off Cα, off the N–Cα axis), **ψ does not** (Cβ is on the Cα side of Cα–C).
- **Noise**: `rmsf` is the least noisy readout (CV 0.115 over seeds), `drift` 0.25, `energy_gap` 1.56
  — the last is a max over the trajectory, so do not use it for per-structure comparisons. Aggregate
  over structures; ~0.9 s per run at 4000 steps.
- **Not MD**: no solvent, no force field, no time. DOPE is knowledge-based, so `rmsf` compares
  *between* structures run with the same settings, not against an MD RMSF in Å.

## Side-chain repack — `tcren.repack` (native `_relax.repack`)

- `repack(structure, chains=("PEPTIDE",), max_chi=2) -> (structure, report)`. Places each side chain
  in the χ conformer DOPE prefers; `report` carries `n_conformers`, `energy`, `p_best` per residue.
- **Same input, same atoms, measured**: side-chain RMSD 4.131 → **2.364 Å in 6 ms**; OpenMM returns
  4.133 Å (unchanged) in 3103 ms, because a local minimiser cannot cross a torsional barrier.
  8/8 structures improved. `tcren refine --repack` runs it after the rigid-body MC.
- **It rotates side chains a model has; it does not build missing ones.** `substitute_peptide` strips
  past Cβ, so that path still returns 44 of 77 heavy atoms — side-chain *construction* is a separate,
  open roadmap row (`refine/CPP_REWRITE.md`).
- Conventions the kernel relies on: the input conformer is index 0 of every enumeration (so a repack
  can never be worse than its input), and weighting is **mean field** (each residue against its
  neighbours at their input conformation — coupled side chains are not resolved).

## Peptide position — `tcren.scoring`

- `peptide_positions(cm, structure)` adds `peptide.pos` (1-based), `peptide.aa`, `peptide.role`
  (anchor/tcr_facing). **Pass the structure**: the sequence comes from it, not from the contacts —
  the class-II 9-mer register heuristic needs the whole peptide (4ozg's gliadin core resolves to
  P1/P4/P6/P9 = positions 2/5/7/10, the published register).
- `position_weights(ann, "uniform"|"central"|"tcr_facing")` → `score_peptides(..., weights=)`.
  `position_profile` decomposes Φ along the peptide (sums exactly to the total); `central_strain`
  isolates the middle band.
- All scoring reweighting goes through one hook: `pipeline._contact_weights`.

## Wrong-TCR decoys — `tcren.shuffle` (`tcren shuffle`)

- `make_decoys(structures, n_per=10, within_class=True, seed=0)` / `graft_tcr(pmhc_source, tcr_source)` /
  `run_shuffle(dir, out, n=10, ...)`. Keep each **oriented** complex's pMHC intact, graft on a **different**
  complex's TCR (within-MHC-class derangement, no real pairing) → wrong-TCR-on-real-pMHC negatives. Real
  (label 1) vs decoy (label 0) trains a **label-free** TCR-recognition classifier; peptide:MHC energy is
  invariant under the graft (built-in control), TCR:peptide/TCR:MHC contacts are new.
- **Direct chain replacement, NOT `orient.substitute_tcr`.** Inputs must be canonically oriented (one common
  frame); the graft copies chains with no per-pair alignment, so each grafted TCR keeps its native MHC–TCR
  docking angle → the decoy set spans the real docking-angle variance (substitute_tcr would collapse every
  donor onto the host's MHC pose). CLI: `tcren orient -s natives/ -o oriented/ && tcren shuffle -s oriented/ -o shuffled/ --n 10`.
- **Finding (2026-07-05):** real-vs-shuffled is learnable at **AUC 0.876** (RF; F_tcr_pep the top feature).
  But a shuffled-trained (crystal) model does **not** transfer to AF-modeled TCRvdb (0.55–0.62 vs AF 0.79) —
  crystal→AF distribution shift. Use it as a label-free recognition prior / supplementary benchmark, not as a
  drop-in TCRvdb scorer.
- **BN classifier** (`tcren.recognition.GaussianBNClassifier`): pure-numpy conditional-linear-Gaussian Bayes
  net — DAG over standardized features (BIC hill-climb on within-class-centred data) + class `y` and MHC class
  as discrete parent nodes; classifies by the Gaussian log-likelihood ratio. `fit/predict_proba`, gzip-JSON
  `save/load`, `to_dot` (graphviz). Trained model shipped at `src/tcren/data/shuffle_bn.json.gz`; the
  reproducible appendix (train+eval, gnuplot ROC/PR + balanced metrics, graphviz BN, marginals) is
  `shuffle_bn/` in the technical appendix (`make`) — moved to the manuscript repo 2026-07-28,
  `2026-tcren/archive/tcren-appendix/`; still in this repo's git history. Decoys are regenerable (`tcren shuffle --seed 0 --n 10`); manifest committed,
  full PDBs belong on HF (351 MB).
- **Distribution-aware logistic** (`tcren.recognition.BayesianLogisticRecognizer` + `encode_features`): a
  *discriminative* alternative to the BN. `encode_features` maps each feature by its natural family — circular
  `dock_torsion` → (cos, sin) von-Mises stats, `chain_balance` → logit, counts/continuous linear, drops the
  duplicate `n_hbond` — then a Bayesian logistic (PyMC NUTS, weakly-informative or horseshoe prior) is fit and
  its posterior mean frozen into `src/tcren/data/shuffle_logistic.json.gz` (dep-light numpy `predict_proba`).
  Real-vs-shuffled 5-fold CV **ROC-AUC 0.885** (matches RF, > BN 0.865 / raw-logistic 0.870). On TCRvdb a
  *supervised* refit with the encoding gives **0.860** pooled (> AF 0.794, raw-feature logistic 0.855); the
  frozen real-vs-shuffled transfer does NOT carry (0.53, crystal→AF shift, same as the BN). Appendix
  `logistic_stan/` (`make PY=<pymc-venv>`; ROC/PR + posterior-forest gnuplot, encoding table).
- **`tcren recognize` / `recognition_features` (2026-07-06):** `recognition.recognition_features(struct)`
  ports the manuscript's 34-descriptor extractor into tcren (docking geometry + TCRen/MJ F & poly-Ala dF +
  contact tallies + biopython ΔSASA `burial` + `mhc_class_bin`) — verified **byte-exact** vs
  `canonical2026_features.csv` (burial max diff 4e-11). Uses `import_structure` (C-gene trimmed) to match
  training; **no `_geom` C-ext needed** (only arda for annotation). `frozen_recognizers()` loads both
  shipped models; `real_probability(rows)` → `{"logistic","bn"}` P(real). CLI `tcren recognize -s pdbs/ -o
  out.tsv` writes one TSV row/PDB = 34 descriptors + `p_real` + `p_real_bn` (`--features-only` skips models).
  The user-facing "one TSV for a/b/d" answer; koff joins it under `--mechanics` (2026-07-28), and only
  ddF (ala), which is per-residue rather than per-structure, stays its own command `tcren ddg`.
- **`--full` feature table (2026-07-13, audited 2026-07-28):** `recognition_features(struct, full=True)` /
  `tcren recognize --full` append the **18 CDR3-frame** descriptors
  (`cdr3{a,b}_{reach,ou,ow,on,au,aw,an,topep,ext}`, FramePose groove-frame projection — the `cdr3b_*`
  strain signal) → 52 features total. Tuples: `RECOGNITION_FEATURES` (34), `CDR3_FRAME_FEATURES` (18),
  `FULL_FEATURES` (52).
- **Descriptor audit (2026-07-28):** every energy column is `F_*` (`e_cdr12`/`e_cdr3a`/`e_cdr3b` →
  `F_cdr12`/`F_cdr3a`/`F_cdr3b`); the duplicate `e_tcr_mhc` and `ct_tp_hydrogen_bond` columns are gone
  (they equalled `F_tcr_mhc` and `n_hbond`); the **12 matrix-swap** columns
  (`{tcren,mj,d}_{tp,cdr12,cdr3a,cdr3b}`) were **removed** — `tcren_*` duplicated `F_*`, and MJ is not the
  potential used on TCR:peptide. New: `crossing_signed` (signed scanning angle, carries docking polarity)
  and `DESCRIPTORS` / `descriptors(family, tcr_only=)` — the catalogue giving each column's family
  (`placement`/`interface`/`topology`/`energetics`/`kinetics`/`score`, with `geometry` and
  `physics` surviving as aliases) and whether the receptor enters it. Five columns do not —
  `F_pep_mhc`, `dF_pep_mhc`, `mhc_class_bin`, `F_pep_int` and `n_pep_int`; they carry cohort identity, so receptor questions must use
  `tcr_only=True`. Frozen recognizers verified **bit-identical** through `_FROZEN_ALIASES`.
- **`--scores` — LEGACY, v1 reproduction only.** Emits the frozen `p_bind` (`binder.binder_score`)
  and `p_forced` (`recognition.forced_pose_score`). Both are fitted, neither is used anywhere in the
  TCRen2 manuscript, and `p_forced`'s coefficients are not re-derivable. Use `reliability.s_free`
  for binder ranking and `cohort.strain_z` for forced-pose grading; both are fit-free.
- **`-t/--threads` on `tcren scoring` and `tcren recognize` (2026-07-26):** both accept a file, a
  directory, a `.tar.gz`, a quoted glob or a `.txt` manifest; `-t N` runs N concurrent workers (`-t 0`
  = all cores). Cohort-relative scores (`Q`, `P_native`, and the legacy `q_bind`/`s_strain` under
  `--scores`) are still computed over the **whole** set, never per batch. `scoring` gains ~7.6x on 8 threads; `recognize` less (its cost is Python
  featurisation, not mmseqs), so batch its annotation rather than expecting linear scaling.

- **`cohort.q_coupled` / `cohort.coupling` — DEPRECATED at 2.12, superseded by `p_native` (2026-07-26):**
  `q_coupled(q, energy)` = `¼[1+erf(z(Q)/√2)]·[1+erf(r·z(ΔΦ)/√2)]` with `r = coupling(q, energy)`, the
  cohort correlation between the geometry and energy channels. Two Gaussian tail probabilities multiplied
  — binding needs both an interface and favourable residues in it — with the energy admitted in
  proportion to `r`, because `E[z(Q)|z(ΔΦ)] = r·z(ΔΦ)` is exactly the part of it that is evidence about
  interface nativeness. **No threshold, no softness constant, no fitted coefficient.** A forced cohort
  disarms itself: `r < 0` flips the energy's sign, `r ≈ 0` collapses its factor to the constant ½ and
  leaves geometry alone. `coupling()` alone is a label-free forced-pose diagnostic (TCRvdb: −0.25 on the
  template-forced GLCTLVAML poses, +0.48 on the clean YLQPRTFLL ones). On TCRvdb it reaches macro
  ROC 0.799 / PR 0.817 / P@10%recall 0.949, ahead of every TCRmodel2 confidence with no generative term.
  For receptor ranking pass the **TCR**-referenced ΔΦ (the peptide is fixed, so the peptide reference
  carries nothing); for peptide ranking pass `reference_delta`'s peptide-referenced ΔΦ.

- **`cohort.f_score` — the contact-energy channel (2026-07-24):** `f_score(table)` =
  `z(-(F_tcr_pep+F_tcr_mhc))`, binder-oriented (`cohort.F_TERMS`). F reads contact chemistry but is
  **pose-conditional**: it inverts on forced poses, which is the whole reason the shape channel
  exists. It **no longer feeds `P_native`**: since 2.17.0 the `energetics` channel draws on the
  `potts` family (`neg_energy`, `log_z`, `log_lik`), not on `F_TERMS`. Either way EM fits that
  channel's sign per cohort rather than being told it. Do not hand-combine `f_score` with Q — that
  is what `p_native` is for.

## MHC allele reference — `tcren build-mhc-ref`, built on demand

- **Not bundled in the wheel.** The curated allele FASTA is built from IMGT/HLA + UniProt mouse by
  `tcren build-mhc-ref`, which must be run **once after a `pip install`**; every command that
  annotates a structure needs it. It lands in `database/mhc/` under `paths.tcren_home()`.
- `paths.tcren_home()` is the single root for on-disk reference data: `$TCREN_HOME` when set,
  else the source checkout (recognised by its `pyproject.toml`), else `$XDG_CACHE_HOME/tcren`
  (`~/.cache/tcren`). Added at 2.12.1, because deriving the roots from `Path(__file__).parents[2]`
  resolved to `site-packages`' parent in an installed wheel and every annotate failed.
  `paths.data_dir()` is `data/` under it, overridable on its own with `$TCREN_DATA_DIR`.

## MHC mapping speed — `mhc.reference.reference_db()`

- `easy_search(query, reference_fasta())` rebuilt the 28k-allele target DB + k-mer index on EVERY
  call (~4.5 s). `reference_db()` runs `createdb` + **`createindex`** once into gitignored
  `data/mhc_cache` (the index is the real cost); `map_mhc` / `annotate_mhc_batch` pass that DB →
  repeated single-structure searches drop to ~0.9 s (5×). Rebuilds if the FASTA is newer.

## Annotation CLI — one `annotate`, no separate `mhc`

- `tcren annotate -s … [--regions all|tcr|mhc|peptide] [--pseudo]` emits ONE per-residue markup
  covering TCR (CDR/FR), MHC groove (HELIX/FLOOR) and peptide. `--regions` filters by chain class;
  `--pseudo` adds `MPS` rows. The old `tcren mhc` command was removed — its allele/class info is
  available in the library via `mhc.map_mhc` / `mhc.annotate_mhc`.

## MHC pseudosequence (MPS) — `tcren.mhc.annotate_pseudo`

- Marks the NetMHCpan 34-residue groove pseudosequence on an annotated structure (region `MPS`).
  Committed FASTAs `src/tcren/data/{mhci,mhcii}_pseudo.fa` (built by `scripts/build_pseudo_fasta.py`
  from NetMHCpan tables; unique seqs, header `<allele>|n=<count>`).
- The 34 positions are **scattered**, so mmseqs/local search can't find them (no shared k-mer).
  Instead each candidate 34-mer is threaded through the chain with a **fitting alignment** (free
  chain gaps; positions are N→C ordered) — ~0.1 s over all ~4k seqs, no prebuilt index. One best
  hit is chosen; class I marks MHCa only (never β2m), class II splits across MHCa+MHCb.
- Validation notebook `notebooks/mhc_pseudosequence_mps.ipynb`: MPS residues vs. 5 Å peptide
  contacts (~half are direct contacts; the rest line the groove toward the TCR).

## Structure output format — `--mmCIF` / `--compress`

- Every command that writes a structure (`orient`, `superimpose`) outputs plain `.pdb` by
  default; `--mmCIF` switches to `.cif`, `--compress` adds a trailing `.gz`. In the library:
  `structure_output_path(dir, id, mmcif=…, compress=…)` + `write_structure(s, path)` (dispatches
  PDB/mmCIF by suffix; a minimal `_atom_site` mmCIF loop that round-trips through Biopython).

## Fetching recent structures — `tcren fetch-recent` / `tcren.recent`

- `tcren fetch-recent [--discover --after YYYY-MM-DD]` → `data/pdb_recent/` (gitignored):
  downloads PDB ids (Native2026 seed; `--discover` adds an RCSB full-text TCR:pMHC search) as
  **mmCIF `.cif.gz`** (the PDB deprecates split `.pdb`; handles **extended >4-char ids**), then
  keeps only complexes with all **5 required chains** (MHCa + b2m/MHCb + peptide + TCR pair),
  validated by one batched annotation pass. `tcren.recent.{fetch_ids,discover_similar,native2026_ids}`.

## Paper-reproduction module (`tcren.paper`)

```python
from tcren.paper import (
    bootstrap, fetch_hf_structures, fetch_vdjdb, fetch_pdb_dates,
    copy_external_inputs, copy_legacy_results,
    contact_table,            # mir extract_contact_map replacement (per structure)
    annotate_structure_set,   # batched TCR annotation over a folder -> (contacts, markup)
    mhc_annotation,           # per-structure MHC allele + class over a folder
    compare,
)
```

- Notebooks live in `notebooks/natcompsci2022/`. HF structure sets are fetched (gitignored) into
  per-set folders directly under `notebooks/data/` — **`notebooks/data/Native2022`,
  `notebooks/data/Native2026`, `notebooks/data/PolyV2022`, `notebooks/data/Bobisse`,
  `notebooks/data/Bigot`** (no `structures/` wrapper). `Canonical2026` is Native2026 after
  `tcren orient`. All structures are gzipped (`*.pdb.gz`).
- Non-structure inputs + 2022 comparison baselines are **committed** under
  `notebooks/natcompsci2022/data_legacy/` (vdjdb, Birnbaum, MJ/Keskin, IEDB, epitope lists,
  `TCRpMHCmodels.tar.gz`, PDB dates, mir/R oracle) — never a pipeline input. `results_new/` is computed.
- Root `data/` — resolved through `paths.data_dir()`, i.e. `data/` under `paths.tcren_home()`,
  which is the checkout in a dev install and `~/.cache/tcren` from a wheel — holds the library
  dataset (gitignored structures): `Native2026` (orientation
  references), `Canonical2026` (the default `superimpose` database), `PDB_date.tsv`,
  `TCRen_potential.csv`. `Canonical2026`'s `orient_metadata.json` moved into
  `src/tcren/data/` — it must ship in the wheel, since `fetch-data` downloads only structures and
  an installed user has no repo `data/`. `setup.sh` runs `tcren fetch-data` at install to
  populate `Native2026` + `Canonical2026` from HF (or lazily on first `superimpose`/`orient`).
  Orientation references load 1ao7/1fyt from `data/Native2026` via `tcren.paths`. The numerical
  regression oracle (legacy mir/R outputs: `contact_maps_PDB.csv`, `tcren_am/tcren.txt`, the
  `example/` set) lives under `tests/assets/oracle/`; the legacy R/Java pipeline was deleted.

## Symmetric TCRen — `derive_tcren(..., symmetric=True)`

TCRen is **directed**: `from` is a TCR residue, `to` a peptide residue, so `F[a,b] != F[b,a]` and
the shipped matrix is 19x20 (classic drops the `from == "C"` row). MJ — used for the presentation
interfaces — is symmetric 20x20. `symmetric=True` folds the **raw counts** onto their transpose
(`potential.symmetrize_counts`, `N + N.T`) *before* the log-odds, so the marginals and hence the
expected term are rebuilt from the folded counts.

```python
from tcren.potential import derive_tcren, symmetrize_counts
sym = derive_tcren(contacts, variant="classic", symmetric=True)   # exactly symmetric, 20x20
```

- **Not a post-hoc average.** Averaging the finished energies leaves the asymmetric background in
  place; on Native2026 the two differ by 0.29 mean / 0.82 max. Symmetrise counts, not energies.
- **Cys is grafted, not dropped.** Free Cys is essentially absent from CDR loops — 4 of 8062
  Native2026 contacts (0.05%) are TCR-side Cys vs 32 (0.40%) peptide-side. The fold grafts the
  peptide-side observations onto the Cys row, so the symmetric matrix keeps a full 20x20 alphabet
  instead of dropping a column for lack of data. `drop_cys` is forced `False` when `symmetric`.
- **Why it is interesting:** it puts TCRen in MJ's functional form, making them directly
  comparable for the first time. They turn out to be nearly **uncorrelated** (Pearson +0.07,
  Spearman +0.10 over the 210 unordered pairs) — TCR:peptide recognition chemistry is not generic
  protein-folding contact chemistry, which is the argument for TCRen existing separately at all.
  Symmetric vs directed TCRen: r = +0.67, mean|diff| 0.29.
- Default is `symmetric=False`; the shipped `TCRen_potential.csv` is unchanged.

## Geometry: contacts, region pairs, docking angle

```python
from tcren.project2d import region_pair_contacts, region_pair_summary  # needs chain-typed + MHC-annotated structure
from tcren.orient import docking_angles
```

- `region_pair_summary(s, kind="closest"|"cb"|"ca")` — inter-chain contact counts for **every**
  region pair (CDR↔peptide, peptide↔MHC, TCR↔MHC, intra-TCR …), not just one interface. Three
  contact definitions: `closest` (5 Å closest heavy-atom pair — the original TCRen definition, the
  only kind that carries a `contact_type` bond classification), `cb` (8 Å Cβ, Cα fallback for Gly),
  `ca` (12 Å Cα). Region-pair labels are ordered canonically (direction-independent).
- Bond types come from the heavy-atom heuristic `project2d.classify_contact` (salt_bridge /
  hydrogen_bond / aromatic / hydrophobic / polar). The external `biotite.structure.hbond`
  (Baker-Hubbard) needs **explicit hydrogens** — it returns 0 on X-ray crystals (no H), so it is
  only useful on protonated / NMR / MD structures. Use the heuristic for crystal structures.
- `docking_angles(s)` — TCR crossing + incident angle from a groove frame built from the peptide
  principal axis + peptide→TCR normal (NOT the whole-complex PCA basis, which the Vα–Vβ spread
  contaminates). ~20–70° crossing for αβ; requires a peptide chain (γδ without peptide raises).

## Gotchas

- nbconvert: pass `--ExecutePreprocessor.kernel_name=python3` (or `=tcren-nb`) or cells silently don't run.
- MHC allele strings from the mapper carry full resolution (e.g. `HLA-A*02:608N`); for
  IEDB-style matching, truncate to 2-field group (`HLA-A*02`).
