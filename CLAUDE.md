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

## The five benchmark blocks — what this library is built to serve (author, 2026-08-30)

**The manuscript's backbone is five benchmark blocks, chosen to cover the most common problems a
user will want to work with.** They are why the public API looks the way it does; a new capability
should name the block it serves.

| # | block | task | what the library must expose |
|---|---|---|---|
| 1 | **CPL** | peptide ranking for a fixed receptor | `tcren.ddg` (`ddg`, `neoantigen_ddg`, `reference_delta`), `tcren pipeline --delta` |
| 2 | **TCRvdb** | receptor ranking for a fixed epitope, functionally validated | `tcren features`, `tcren recognize`, `cohort.q_score`, `reliability.t_score`, `reliability.s_free` |
| 3 | **VDJdb panel** | receptor ranking under template scarcity | the same, plus the template covariate being reportable rather than inferred |
| 4 | **Kinetics / ergodicity** | the licence to score one static structure, and magnitude | the three per-interface contact energies, `tcren.potts` contact marginals, the occupancy law |
| 5 | **ipTM / pLDDT diagnostics** | which confident complexes are not real, and correcting the confidence | `tcren diagnose`, `reliability.correct_confidence`, `available_corrections` |

Blocks 1-2 are the two ranking tasks; block 3 is the same task in the regime the applied question
lands in; block 4 licenses the whole approach; block 5 is the one a user reaches for first, holding
an AlphaFold model and no way to judge it.

**Block 5 is the only fitted read-out.** `correct_confidence` takes the generator's confidence as a
prior and adds the structure as log-odds, learning four coefficients and freezing them — the same
standing as the Platt links. Every other score takes no label anywhere. Say so wherever it is
reported. Report the template split, never the pool: pooled over the panel the correction reads
-0.015, which averages +0.051 over the 6 template-covered cohorts against -0.039 over the 16 free
ones, and the pooled number is a Simpson artefact.

## Hard conventions

- **Scope is alphabeta TCR:pMHC and nothing else.** A TRA and a TRB chain, a standard-amino-acid
  peptide, a class I or class II groove. Out-of-scope structures are dropped, not widened into.
- **CDR3 is not junction.** `junction_aa` carries the Cys104/Phe118 anchors and is two residues
  longer than IMGT CDR3. Confirm which a caller means before any Pgen, clustering or matching work.
- **Every on-disk data root goes through `paths.tcren_home()`** — `$TCREN_HOME`, then the source
  checkout, then `$XDG_CACHE_HOME/tcren`. Deriving a path from the source layout breaks the
  installed wheel, silently, and exits 0 with an empty table.
- **A command whose every row carries an error exits non-zero.** That guard exists because the
  silent-success failure mode above cost a debugging session; keep it.
- **`pitch_angle` is generator-confidence leakage, never interface geometry.** `cohort.P_NATIVE_BANNED`
  makes naming it in a channel raise rather than silently fit.

## The descriptor whitepaper — `appendix/` (2026-09-02)

`appendix/descriptors.pdf` is the reference for what the catalogue *is*: the reduction chain, the
nine operators, every formula and its derivation, MHC class I/II, and all 164 rows with units and
`STATUS`. `make` in that directory rebuilds it. Three things not to relearn:

- **The hierarchy is a chain of quotients**, $S \to C \to B \to T$, and $B \twoheadrightarrow T$ is
  the load-bearing step. Every functional of the cell tally is a functional of the biadjacency
  matrix and not the reverse, so **redundancy runs down the chain and never up it**. That orders a
  redundancy screen's output before it is run, and it is why the parameter-free `g_*` forms are the
  ones to prefer over the radius-tagged `fp_*` ones on parsimony grounds.
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
with `-W`, so a warning is a failure. Update `skills/tcren/SKILL.md` and `docs/` in the same change
as any new public capability, and `CHANGELOG.md` before any release. A new descriptor also needs
`python scripts/gen_appendix.py` — the whitepaper's tables come from the catalogue, not from prose.
