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
