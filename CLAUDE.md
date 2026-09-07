# CLAUDE.md — tcren

`tcren` is the library: structure -> features -> scores for alphabeta TCR:peptide-MHC complexes,
released to PyPI. See `README.md` for what it exposes, `skills/tcren/SKILL.md` for the public API,
`docs/` for the reference, `ROADMAP.md` and `STATUS.md` for what is in flight.

## Three repositories, one boundary (do not blur it)

| repo | role |
|---|---|
| this repo | the library. Every descriptor, score and fit lives here |
| `~/vcs/projects/2026-tcren2-code` | the benchmark: datasets, task definitions, evaluation (ROC/PR/CI), the ledger |
| `~/vcs/manuscripts/2026-tcren2-ms` | LaTeX plus a generated publication layer |

If you find yourself computing an AUC here, stop — it belongs in the benchmark. If you find yourself
re-implementing featurisation in the benchmark, stop — it belongs here.

## The score set — the 3.0.0 public surface

**One frozen object, five read-outs, every one defined for a single structure.** The transform, the
class means and the covariance are all frozen on a hold-out that ships in the wheel, so nothing is
estimated from the rows being scored.

| name | tier | what is estimated |
|---|---|---|
| `peptide_score` | 0 | nothing; the direction is fixed by the potential |
| `pose_score` | 1 | a covariance over hold-out binders — no negative, no label |
| `confidence_residual` | 1 | the same covariance, read as a conditional mean |
| `binder_score` | 2 | class means and covariances, from hold-out binder labels |
| `channel_scores` | 2 | the same object, marginalized to one descriptor family |

Channels: `placement`, `interface`, `shape`, `energetics`, `mechanics`. Two commands —
`tcren features` then `tcren assess`. `cohort.q_score`, `reliability.t_score` and
`reliability.s_score` stay as the **fit-free predecessor tier**; `S` composes with `binder_score`
rather than being replaced by it, and is reported beside it.

## The five benchmark blocks — what this library is built to serve (author, 2026-08-30)

**The manuscript's backbone is five benchmark blocks, chosen to cover the most common problems a
user will want to work with.** They are why the public API looks the way it does; a new capability
should name the block it serves.

| # | block | task | what the library must expose |
|---|---|---|---|
| 1 | **CPL** | peptide ranking for a fixed receptor | `tcren.ddg` (`ddg`, `neoantigen_ddg`, `reference_delta`), `tcren pipeline --delta` |
| 2 | **TCRvdb** | receptor ranking for a fixed epitope, functionally validated | `tcren features`, `tcren assess`, `tcren.score.binder_score`, and the fit-free tier `cohort.q_score` / `reliability.t_score` / `reliability.s_score` |
| 3 | **VDJdb panel** | receptor ranking under template scarcity | the same, plus the template covariate being reportable rather than inferred |
| 4 | **Kinetics / ergodicity** | the licence to score one static structure, and magnitude | the three per-interface contact energies, `tcren.potts` contact marginals, the occupancy law, and — since 3.6.0 — the same marginals with the Cα geometry *modelled*: `tcren potts geom`, `tcren potts predict --ensemble`, `geom_tcr_peptide_ens` |
| 5 | **ipTM / pLDDT diagnostics** | which confident complexes are not real, and whether the confidence is warranted | `tcren assess`, `tcren.score.confidence_residual`, `tcren.score.pose_score`, `reliability.af_band` |

Blocks 1-2 are the two ranking tasks; block 3 is the same task in the regime the applied question
lands in; block 4 licenses the whole approach; block 5 is the one a user reaches for first, holding
an AlphaFold model and no way to judge it.

**Nothing is fitted at score time, and the frozen fit is reproducible.** The 2.28.0 removal took out
every out-of-fold-fitted read-out, `correct_confidence` and the `tcren diagnose` command with it. Its
replacement is `score.confidence_residual` at a strictly better standing: tier 1, the reported ipTM
minus `E[logit ipTM | x]` under the binder Gaussian, with **no binder label anywhere**. The tier-2
read-outs (`binder_score`, `channel_scores`) do read hold-out labels, but the hold-out is named,
deposited and shipped — `score.holdout_manifest()` returns its 8,292 structures and `tcren
fit-holdout` regenerates the shipped arrays from them, which is exactly the contract `P_native` could
not offer. **Report the template split, never the pool.** On the 22-cohort VDJdb panel every arm
divides under template coverage, and a macro number over the two strata is a Simpson artefact.

**A contact map no longer needs a structure of the complex (3.6.0).** `available_pairs` reads
`d_Cα` off the coordinates, so every contact read-out was conditional on holding a solved complex.
`geom_tcr_peptide` models the distance matrix instead — a graphon average over the same 359
crystals `tcren2` is derived from — and `tcren potts predict --sequences` scores a map with **no
structure read at all**. Held out by structure the predicted 15 Å envelope holds **7,399 of 7,687
crystal contacts (96.25 %)** against the crystal's own envelope at 99.78 %, and the coupled model
reads the crystal's indicator inside it at **median ROC-AUC 0.8310** over 361 structures.

**Average over geometries, not over surfaces (3.7.0).** Scoring one map read off the *mean* surface
recovered less than half the contacts, and the fix is an order of operations rather than a
parameter. A contacting pair sits at the near end of the distribution being averaged — measured
**7.94 Å**, placed at **10.90 Å** over 7,651 contacting pairs in 364 crystals — and per-pair
averages need not describe a chain any loop could adopt, consecutive Cα atoms being a fixed virtual
bond apart. `predict_contact_map` / `tcren potts predict --ensemble` score under **each** reference
geometry and average the probabilities, which makes every sample a solved chain: **24.62 expected
contacts against 21.02 observed** (ratio 1.17) where the mean surface reads 9.80 (ratio 0.47), and
the contact ROC-AUC rises **0.9242 → 0.9325**. Both the level and the ordering improve, so this
supersedes the mean-surface reading rather than qualifying it. `geom_tcr_peptide_ens` ships the
2,154 kept surfaces; `geom_tcr_peptide` is untouched and still the default of `predicted_sites`.

**What it does not buy is a count specific to one receptor** — Pearson +0.23 against the observed
count, MAE 5.95 contacts against a population s.d. of 6.08 — because a linear average over
reference geometries returns the level of a typical interface. That is step 3 of the ladder in
`ROADMAP.md` §1, a parametric loop with no crystal coordinates anywhere, and it is open. Record:
`$TCREN2_BENCH/bench/eda/out/geom_ensemble.md`.

**Every TCR against every pMHC, and a geometry you can draw from (3.9.0).** `Canonical2026` puts
every complex in one frame, so a TCR from one crystal against a pMHC from another needs no
transform — `tcren potts cross` enumerates all **79,061** within-class pseudo-complexes against the
359 solved and writes **36,269,320** parquet rows in **31 s**, one numpy broadcast per
`(class, chain, region, L, M)` bucket and no loop over pairs. It buys **length coverage** (131
`(L, M)` combinations for TRA:CDR3 against 62 solved) and cell-level sample size; it does **not**
buy receptor:peptide dependence, the cross product being the product of the two marginals. The
clash filter is measured: 3.35 Å, just under the **3.366 Å** closest approach any solved complex
makes over 175,358 pairs.

**A pairing is `shuffle.graft_tcr`, and the frame is the whole contract.** Keep the host pMHC, drop
its TCR, put the donor's in — the graft applies *no* transform because the set is co-framed, and the
table matches a built chimera to 4.5e-5 Å over 40 of them. **On an unoriented set every pairing is
meaningless and the table looks perfectly normal**, which is why `residue_tables` measures the
peptide-centroid spread and refuses rather than scoring; `--superimpose auto` orients at ~3 ms a
structure. `Canonical2026` reads 3.61 Å untouched. **`tests/assets/pdb` reads 141.7 Å as deposited**
— the fixtures were never co-framed, and every cross test written before the guard existed was
exercising an invalid regime without failing. Distances are frame-invariant, so those tests were
valid arithmetic; they just were not testing the product.

`tcren potts chain` then fits a joint Gaussian over Cα **coordinates** rather than distance
surfaces. A mean of distance matrices describes no chain — consecutive Cα are 3.7985 Å apart,
s.d. 0.0668 over 642 CDR gaps — and recovering coordinates from distances is a non-convex EDM
completion, so working in coordinates makes the constraint exact (`geometry.enforce_bond`) instead
of relaxed. **Length stops being a joint problem**: `L` conditions the loop block and `M` the
peptide block separately. **The draw set brackets the truth** — over 540 held-out contacts the best
of 512 draws places each to **0.04 Å** — while the mean over draws overshoots by +3.25 Å, because a
contacting pair sits at the near end of the distribution being averaged. **So the generator is not
the bottleneck; the selection rule is**, and the one rule tried is worse (below).

**`g_orient` — does the residue turn toward the peptide (3.9.0).** The Cα geometry says a pair is
close enough; whether it *faces* the partner was not representable, the site table carrying
`d_heavy` and `d_ca` and no Cβ. Held at a fixed Cα distance the Cα→Cβ cosine separates contact from
non-contact at ROC-AUC 0.69–0.84 over 4–12 Å, and fitted it is worth **+787.0 nats for 8
parameters** (369 crystals, `J` pinned to `tcren2`) against the 3.8.0 exclusion class's +81.8 for
one. **Where it pays is a crystal**: held out over 24, the map on the complex's own Cα goes js
0.7052 → **0.7290** and ROC-AUC 0.9599 → **0.9691**. On a *drawn* geometry it is neutral, because
`cos_cb` is a deterministic function of the Cα trace and can only be as receptor-specific as the
trace. The Cβ is **reconstructed on both sides** (`geometry.virtual_cb`, 0.547 Å residual) so
fitting and scoring describe one quantity. Opt-in; every shipped model has `g_orient == []`.

**Two things measured and not to be re-tried (3.9.0).** *Annealing over drawn geometries costs*:
weighting draws by their own `log Z` on a **union** site set — so the site-count confound behind the
earlier failure is ruled out — still runs js 0.4342 → 0.3538 and E[contacts] 22.0 → 47.1 against
22.6 observed from β = 0 to β = 1, because `log Z` is a free energy and rises with occupancy. β
defaults to 0. *And the graphon registration mismatch is deliberate*: `_resample` and
`predicted_distances` are not adjoint, but making them so buys 0.19 Å on CDR2 and costs 0.49 Å on
CDR3 (+3.26 against the shipped +2.91 over 7,732 contacts, 370 leave-one-out crystals). `endpoint`
registration aligns loop **termini** across lengths, which is the operation an average performs.
**So the CDR2 deficit is not a registration artefact** — it survives every pairing at 5.6–5.8 Å
against CDR3's 2.3–2.8 Å.

**The average is for the ranking read-outs; the contact map itself wants a constrained geometry
(author, 2026-09-07).** *"Average is good for some metrics and readouts, but it is unacceptable for
actual contact map modelling."* The criterion for a map is that the expected count reproduce the
crystal's within **±3** and that a contact probability be Beta distributed — occupied nearly always
or nearly never, never p ≈ 0.5. Both hold on the crystal's own Cα and neither survives the flat
average, and the reason is **not** a calibration that could be tuned out: a logistic-power
recalibration fits γ = 1.009, so the map is already calibrated. It is **placement**, and
`potts.map_agreement` says so in one number: over 366 held-out crystals the Jensen–Shannon
similarity to the crystal's own map, level divided out, reads **0.4200** for the flat ensemble
against **0.7052** on the crystal's own Cα. Mechanically it is receptor-residue occupancy — at a
half-contact threshold the flat ensemble engages 15.25 receptor residues at E[degree] 1.35 contacts
against step 1's 10.53 at 1.77, where the crystal makes 11.28 at 1.86. **Compare arms, never an arm
against the crystal's hard count**: a probabilistic map spreads mass by construction and step 1 does
it too. Records: `$TCREN2_BENCH/bench/eda/out/geom_sharpen.md`, `geom_select.md` and
`contact_fidelity.md`.

**Bigot is block-3 tier; Bobisse is not.** Both score a poly-V placeholder peptide, but the
templates are opposite. Bigot's 14 receptors have **no deposited ternary complex** — 0 of 14 CDR3α
and 0 of 14 CDR3β occur in any of the 362 Native2026 crystals — so a pure-potential arm there is
outside the regime a contact potential is measured in (energetics reads 0.692 templated against
0.580 template-free, FINDINGS §2); report it case by case beside a geometry channel, never as a
peptide-ranking floor. **Bobisse's receptor is 6UK4** (CDR3α `AFMDSNYQLI` / CDR3β `ASSRTSPTDTQY`,
cognate KQWLVWLFL), so it is template-rich like CPL; its limit is being one ranking unit.

## Hard conventions

- **Scope is alphabeta TCR:pMHC and nothing else.** A TRA and a TRB chain, a standard-amino-acid
  peptide, a class I or class II groove. Out-of-scope structures are dropped, not widened into.
- **CDR3 is not junction.** `junction_aa` carries the Cys104/Phe118 anchors and is two residues
  longer than IMGT CDR3. Confirm which a caller means before any Pgen, clustering or matching work.
  **Adding the anchors is a lookup, not a concatenation**: Cys104 is the last residue of FR3, and
  the 118 anchor is the first residue of FR4, which opens with the J-region F/W-G-x-G motif.
  `"C" + cdr3 + "F"` fails both ways — the anchor is Trp on some J genes (302TIL alpha carries FR4
  `WGAGTK`, junction `CAFMDSNYQLIW`), and a CDR3 may itself end in F so the junction carries two
  (`ASEGVHEQF` → `CASEGVHEQFF`). Asking whether a sequence "ends in F" is not a test; the FGXG
  motif in FR4 is.
- **Every on-disk data root goes through `paths.tcren_home()`** — `$TCREN_HOME`, then the source
  checkout, then `$XDG_CACHE_HOME/tcren`. Deriving a path from the source layout breaks the
  installed wheel, silently, and exits 0 with an empty table.
- **A command whose every row carries an error exits non-zero.** That guard exists because the
  silent-success failure mode above cost a debugging session; keep it.
- **`pitch_angle` is generator-confidence leakage, never interface geometry.** It out-discriminated
  every clean docking angle for exactly that reason. `recognition.STATUS` flags it, and it is not in
  the coordinate set `tcren.score.transform.working_set` hands the frozen model.

## The descriptor whitepaper — `appendix/` (2026-09-02)

`appendix/descriptors.pdf` is the reference for what the catalogue *is*: the reduction chain, the
nine operators, every formula and its derivation, MHC class I/II, and all 164 rows with units and
`STATUS`. `make` in that directory rebuilds it. Three things not to relearn:

- **The hierarchy is a chain of quotients**, $S \to C \to B \to T$, and $B \twoheadrightarrow T$ is
  the step that loses information. Every functional of the cell tally is a functional of the
  biadjacency matrix and not the reverse, so **redundancy runs down the chain and never up it**.
  That orders a redundancy screen's output before it is run, and it is why the parameter-free
  `g_*` forms are the ones to prefer over the radius-tagged `fp_*` ones on parsimony grounds.
- **Nineteen descriptors are one formula.** The Hill number $D_q$ generates the whole `*_cell` /
  `*_loop` / `pep_cov_*` / `g_even_*` / `m_erank_*` / `partcoef_*` set. Two identities not in
  `STATUS`: the participation coefficient is Gini--Simpson $1 - 1/D_2$, and the participation ratio
  is $D_2/k$. **`m_erank_*` is order 2, not Roy \& Vetterli's order 1** — the code is correct and
  its docstring says so, but the literature name points at $D_1$, so any comparison against a
  published effective rank must use $q = 2$.
- **Class II lives only in the crystals** — 94 of 374, and every modelled set is entirely class I.
  `mhc_class_bin` is therefore constant on both receptor benchmarks, and no class II claim is
  checkable against a model. The gap channel is class-invariant (`sc_interlock_frac` 0.769 against
  0.762), so interdigitation is not a class I artefact.

`scripts/gen_appendix.py` emits everything factual about the catalogue; `--check` fails when it is
stale, as for `gen_descriptor_table` and `gen_family_graph`. Never edit `appendix/generated/*`.

## Working here

`./setup.sh` creates the venv and installs editable; `./setup.sh --test` runs the suite. Docs build
with `-W`, so a warning is a failure.

**Two scripts sit above `setup.sh` and the split between them is the point.** `./bootstrap.sh` is
the environment and the data — it wraps `setup.sh`, fetches every set of the HuggingFace deposit
(about 620 MB) and builds the MHC allele reference — and is **deliberately outside any time
budget**. `./reproduce.sh` assumes both are on disk and runs everything the manuscript prints:
four library stages here (106 s measured, 2026-09-05, macOS arm64), then a handoff to
`$TCREN2_BENCH/recompute.sh` and `$TCREN2_MS/make`. **Neither script has a `--skip-slow` flag and
neither should get one** — the benchmark repository's equivalent claims it drops only stages no
reported table reads, which is false for four of the eight it names. Both fail loudly, by name, on
a missing input. `docs/reproduction.rst` is the prose version. Update `skills/tcren/SKILL.md` and `docs/` in the same change
as any new public capability, and `CHANGELOG.md` before any release. A new descriptor also needs
`python scripts/gen_appendix.py` — the whitepaper's tables come from the catalogue, not from prose.
