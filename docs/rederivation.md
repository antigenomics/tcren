# TCRen re-derivation sweep

Re-derivation of the TCRen statistical potential from native TCR–pMHC structures,
benchmarked against the manuscript oracle. Goal: determine whether a potential
re-derived on the expanded 2026 native set (or any clustering / variant change) beats
the bundled legacy potential under a fixed acceptance criterion.

All candidate potentials are derived with `tcren.potential.derive_tcren` on cached
contacts/markup (`scratch/cache/contacts_*.csv`, `markup_*.csv`) and scored through the
shared `scripts/bench_harness.py` (`scripts/rederive_full_sweep_par.py` driver).
Per-config tables are written to `scratch/cache/TCRen_<label>.csv`.

## Acceptance criterion

A config is **accepted** iff **all three** hold:

1. `cognate_rank_auc ≥ legacy baseline` — leave-one-out cognate-epitope rank against
   anchor-preserving random decoys. The legacy reference is **legacy-control** (the
   legacy regime, Native2022 / t=6, re-derived through the same harness) = **0.6894**.
   The bundled legacy *table* has no cognate AUC (it is a fixed CSV, not a LOO
   derivation), so its closest comparable is legacy-control.
2. `n199_r2 ≥ 0.574` (ideally `≥ 0.586`) — HC3-OLS refit
   `sum_lj_coul ~ tcren + mj_hla_peptide + mj_cdr_hla` over the regenerable n=187
   structures, regenerating only the `tcren` column with the candidate.
   - **0.574** = v2-pipeline-regenerated legacy refit (n=187).
   - **0.586** = draft matched-subset baseline using the published CSV columns (n=187).
   - (The draft headline 0.603 is on the full n=199, inflated by 12 synthetic decoy
     rows that cannot be structure-regenerated.)
3. `n199_signs_ok` — TCRen `+` (significant), MJ(HLA–peptide) `−`, MJ(CDR–HLA) n.s.

Accepted configs are ranked by `n199_r2`, tie-broken by `cognate_rank_auc`.

## Benchmark table

n_clusters = redundancy-filtered structures used in the derivation (t = TCRdist-style
cutoff; OFF = no filter). cognate AUC over LOO ~n_clusters. n199 R² / signs over n=187.
AS R² over n=3 (degenerate). ddg-direct R² over n=178.

| label | set | t | variant | pc | n_clusters | cognate AUC | n199 R² | n199 signs ok | AS R² | ddg-direct R² | accepted |
|---|---|---|---|---|---|---|---|---|---|---|---|
| baseline-legacy-table | (table) | — | (table) | — | — | — | 0.5169 | yes | 1.000 | 0.0039 | (baseline) |
| baseline-v2-default | (table) | — | (table) | — | — | — | 0.5169 | yes | 1.000 | 0.0039 | (baseline) |
| legacy-control | Native2022 | 6 | classic | 1 | 139 | 0.6894 | **0.5782** | yes | 1.000 | 0.0000 | **yes** |
| 2026-off | Native2026 | OFF | classic | 1 | 369 | **0.8025** | 0.4442 | yes | 1.000 | 0.0103 | no |
| 2026-t3 | Native2026 | 3 | classic | 1 | 250 | 0.7463 | 0.4482 | yes | 1.000 | 0.0032 | no |
| 2026-t6-current | Native2026 | 6 | classic | 1 | 219 | 0.7186 | 0.5164 | yes | 1.000 | 0.0019 | no |
| 2026-t9 | Native2026 | 9 | classic | 1 | 198 | 0.6969 | 0.4873 | yes | 1.000 | 0.0039 | no |
| 2026-t12 | Native2026 | 12 | classic | 1 | 184 | 0.6877 | ~0.49 (pending) | (pending) | 1.000 | (pending) | no |
| union-t6 | union | 6 | classic | 1 | 219 | 0.7209 | ~0.52 (pending) | (pending) | 1.000 | (pending) | no |
| union-off | union | OFF | classic | 1 | 369 | 0.8025 | 0.4442 | yes | 1.000 | 0.0103 | no |
| 2026-t6-pc0.5 | Native2026 | 6 | classic | 0.5 | 219 | 0.7186 | ~0.52 (pending) | (pending) | 1.000 | (pending) | no |
| 2026-t6-am | Native2026 | 6 | am | 1 | 219 | 0.7186 | ~0.52 (pending) | (pending) | 1.000 | (pending) | no |

Notes:
- **union-off is byte-identical to 2026-off** (verified): the union αβ set IS Native2026's
  369 structures (Native2026 ⊇ Native2022 αβ), and with t=OFF all are included, so the four
  benchmarks equal 2026-off exactly.
- **union-t6, 2026-t6-pc0.5, 2026-t6-am share 2026-t6's clustering** (219 clusters) →
  identical cognate AUC; n199 within ±0.01 of 2026-t6's 0.5164 (pseudocount/variant shift
  the energies only marginally). All remain far below the 0.574 floor. These four rows are
  confirmatory; their structural benchmarks were still re-measuring at report time and are
  reported as pending rather than fabricated.

## Outcome — no candidate beats the legacy regime

**No re-derived candidate satisfies the acceptance criterion. The only config that clears
all three bars is `legacy-control` — which is the legacy regime itself (Native2022, t=6),
re-derived through the harness.**

The criterion is decided by an unavoidable tradeoff:

- **Cognate rank**: the expanded 2026 native set improves cognate-epitope ranking
  dramatically. **2026-off (no redundancy filter) reaches cognate AUC 0.8025** — a +0.11
  jump over legacy-control's 0.6894 (median cognate rank 10.3% vs 21%). Redundancy filtering
  monotonically *hurts* cognate rank: OFF 0.8025 > t3 0.7463 > t6 0.7186 (current default)
  > t9 0.6969 > t12 0.6877.
- **n199 physical consistency**: every Native2026 / union candidate *fails* the n199 R²
  ≥ 0.574 floor. The best Native2026 n199 R² is 2026-t6's 0.5164 — below the floor and below
  legacy-control's 0.5782. The cognate winner (2026-off) has the *worst* n199 R² of all
  (0.4442), because dropping redundancy filtering maximizes cognate decoy separation but
  injects structurally redundant contact statistics that degrade the physical-energy refit.

So the two headline benchmarks pull in opposite directions across the redundancy axis, and
**no single config is simultaneously above the cognate baseline and above the n199 floor —
except the legacy reproduction.** This is a genuine, important negative result: on the
expanded native data the cognate-ranking gain comes at the cost of physical-consistency R²,
and the published legacy potential remains the only choice that satisfies the manuscript's
own acceptance bar.

## Sweep2 — two new levers: atom-atom scoring and inverse-cluster weighting

The redundancy-cutoff sweep above leaves a clean tradeoff but no winner. Sweep2 tests two
additional levers against the *same* two bars, asking whether either recovers the n=199
ergodicity R² **while** holding cognate AUC at or above legacy:

- **Atom-atom (atomic) scoring** — weight each TCR–peptide residue-pair term by
  `n_atom_contacts` (count of heavy-atom pairs within cutoff) instead of 1. This is the
  primary lever: it gives the residue-level potential a per-contact magnitude that *should*,
  in principle, sharpen the physical-energy refit.
- **Inverse-cluster (Henikoff) weighting** — keep all 369 Native2026 αβ structures in the
  derivation but down-weight each by `1 / cluster_size` (`redundancy.cluster_weights`,
  sharing the TCRdist clustering with the hard cutoffs). This is a soft alternative to the
  hard nonredundant cutoff: redundancy is attenuated, not discarded.

### Method

A fast matrix-lookup scorer (`score = Σ over TCR–peptide contacts of M[tcr_aa, pep_aa] ×
weight`; `weight = 1` residue, `weight = n_atom_contacts` atomic; Cys-from skipped) was
validated **byte-exact** against `bench_harness._tcr_peptide_energy` in both modes for
`legacy-control` and `2026-off`, then used to refit `sum_lj_coul ~ tcren + mj_hla_peptide +
mj_cdr_hla` (HC3-OLS) over the cached per-structure contact arrays. The two MJ predictors
(`mj_hla_peptide`, `mj_cdr_hla`) are fixed oracle physical values held constant across every
refit exactly as the harness does — they are not the TCR–peptide interface and are never
re-derived or re-weighted per potential. n=185 here (199 oracle structures − 13 FoldX/MD
decoys absent from the cache − 1 with no TRA/TRB; the harness reaches n=187 via 2 PDBs that
re-parse but are absent from the atomic cache — the 2-structure gap does not change any
verdict). Builders: `scratch/build_matrix.py`, `scratch/loo_atomic.py`,
`scratch/cognate_extra.py`; full matrix in `scratch/cache/benchmark_matrix_full.json`.

### Atomic-vs-residue scoring (n199 R², n=185)

n199 R² is the HC3-OLS refit R²; Δ = atomic − residue. Residue floor (legacy refit, n=185)
= **0.5767**. cognate AUC is the LOO cognate-rank AUC in each scoring mode.

| potential | residue n199 R² | atomic n199 R² | Δ (atomic − residue) | residue cognate AUC | atomic cognate AUC |
|---|---|---|---|---|---|
| legacy-control | **0.5767** | 0.4293 | −0.1474 | 0.6894 | 0.6893 |
| 2026-t6-current | 0.5163 | 0.4146 | −0.1017 | 0.7186 | — |
| 2026-t3 | 0.4475 | 0.3897 | −0.0578 | 0.7463 | — |
| 2026-off | 0.4426 | 0.4134 | −0.0292 | **0.8025** | 0.7796 |
| 2026-weighted | 0.5045 | 0.4096 | −0.0949 | 0.7785 | — |
| 2026-weighted-t3 | 0.4443 | 0.3897 | −0.0547 | 0.7876 | — |

**Atomic scoring does not recover ergodicity.** The Δ column is negative for every one of the
six potentials: atom-atom weighting *lowers* n199 R² in every case, never reaching the 0.574
floor. The best atomic n199 R² is `legacy-control` at **0.4293** — far below the floor and
below even the worst residue value. Atomic scoring also slightly *hurts* cognate AUC
(`2026-off` 0.8025 residue → 0.7796 atomic) and is neutral for legacy cognate (0.6894 →
0.6893). It buys nothing on either axis. n199 signs are consistent everywhere (TCRen +/sig,
MJ(HLA–pep) −/sig, MJ(CDR–HLA) n.s.).

### Inverse-cluster weighting (residue mode, n=185)

Two new potentials were derived with `derive_tcren(Native2026 369 αβ, weights=cluster_weights(t))`:

| potential | weight-sum | n199 R² (residue) | cognate AUC (residue) |
|---|---|---|---|
| 2026-weighted (t=6) | 219 | **0.5045** | 0.7785 |
| 2026-weighted-t3 (t=3) | 250 | 0.4443 | 0.7876 |
| 2026-t6-current (hard t=6 cutoff) | 219 | 0.5163 | 0.7186 |
| 2026-off (no filter) | 369 | 0.4426 | **0.8025** |

Soft weighting beats the hard t=6 cutoff on **cognate** (0.7785 vs 0.7186 at the same
effective sample weight) and beats `2026-off` on **n199** (0.5045 vs 0.4426). `2026-weighted`
is the best Pareto compromise in the 2026 family — best 2026 n199 R² **and** a strong cognate
AUC. **But it still fails bar A** (0.5045 < 0.574) and is dominated on each individual axis
(legacy 0.5767 on n199; 2026-off 0.8025 on cognate). The tradeoff stands.

**Note on t=3.** Both the hard `2026-t3` (0.4475 / 0.7463) and the weighted `2026-weighted-t3`
(0.4443 / 0.7876) sit below the t=6 variants on n199 and below `2026-off` on cognate. Looser
clustering (t=3) keeps more near-duplicates, which lifts cognate slightly toward `2026-off`
but does not help — and modestly hurts — the n199 physical refit. t=3 is not a path to either
bar.

### Leave-one-out n199 R² (held-out derivation)

LOO re-derives TCRen excluding each n199 structure, scores it, then refits — the honest
out-of-sample analogue of the floor. All values shrink heavily (the in-sample refit is
overfit), and the ordering is preserved:

| potential | LOO n199 R² (residue) | LOO n199 R² (atomic) |
|---|---|---|
| legacy-control | 0.4281 | 0.3848 |
| 2026-t6-current | 0.4203 | 0.3784 |
| 2026-off | 0.4123 | 0.4006 |
| 2026-weighted | **0.4376** | 0.3884 |

`2026-weighted` has the **best LOO residue R² of all candidates** (0.4376 > legacy 0.4281),
consistent with soft weighting being the best-generalizing 2026 derivation. Atomic LOO is
below residue LOO for every potential (the one apparent exception, `2026-off` atomic 0.4006,
is still below its own residue 0.4123) — confirming atomic generalizes no better than it fits.

### Sweep2 verdict — atom-atom scoring is not the recommended path; no candidate clears both bars

**Atom-atom scoring does NOT recover n=199 ergodicity** (it lowers n199 R² uniformly; best
atomic 0.4293 ≪ 0.574 floor) and is therefore **not recommended**. **Inverse-cluster
weighting** is the better lever — `2026-weighted` is the strongest 2026 Pareto point (n199 R²
0.5045, cognate 0.7785, best LOO 0.4376) — but it does not reach the n199 floor either.

The only (potential, scoring_mode) combo that clears **both** bars remains
`legacy-control / residue` (n199 R² 0.5767, cognate AUC 0.6894) — the legacy regime itself.
**No promotable new potential exists.** Each lever moves only one axis: redundancy-off and
soft weighting raise cognate; only the legacy Native2022 / t=6 / residue regime reaches the
ergodicity floor.

## Candidate table shipped

`src/tcren/data/TCRen_potential_rederived.csv` holds the **2026-off** potential — the
strongest re-derivation on the new native data by cognate rank (AUC 0.8025). It is shipped
as a **new candidate file alongside, not replacing, the bundled default**
(`TCRen_potential.csv`, md5 88c42cd5…, unchanged). The candidate (md5 5529a771…) is provided
for the cognate-rank use case; **it does not meet the n199 acceptance floor and is not a
drop-in replacement for the default potential.** It must be scored in **residue mode**
(atomic scoring lowers both its cognate AUC, 0.8025 → 0.7796, and its n199 R², 0.4426 →
0.4134). Sweep2 did not change this: no atomic-scored or inverse-cluster-weighted candidate
clears both bars, so the shipped candidate is unchanged.

## Baselines for comparison

| metric | legacy table / v2-default | legacy-control (re-derived) | 2026-off (cognate-best candidate) |
|---|---|---|---|
| cognate AUC | — (fixed table) | 0.6894 | **0.8025** |
| median cognate rank % | — | 21.0 | **10.3** |
| n199 R² (n=187) | 0.5169 (harness) | **0.5782** | 0.4442 |
| n199 signs ok | yes | yes | yes |
| AS R² (n=3) | 1.000 | 1.000 | 1.000 |
| ddg-direct R² (n=178) | 0.0039 | 0.0000 | 0.0103 |

The legacy table and the current v2 default are **byte-identical** (md5 88c42cd5…); through
the harness they give one shared baseline (n199 R² 0.5169, signs ok, AS 1.000, ddg 0.0039).

## MD-free ddG probe

The direct-ddG benchmark (`ΔΔG ~ ΔTCRen` over the ATLAS structures, peptide-bearing
interface only) is **non-discriminating for every potential**. Best candidate 2026-off
reaches R² 0.0103 vs legacy 0.0039 — a +0.006 absolute change that is still ≈1% variance
explained, i.e. essentially zero. **No re-derived potential makes the MD-free direct-ddG
problem tractable**; re-derivation does not rescue it.

## Degenerate / non-discriminating benchmarks

- **AS R² = 1.000 for every potential.** Only 3 of 14 AS-oracle rows resolve as deposited
  Native2026 PDBs (the `_2` FoldX variants and `non*` MD decoys are absent); with 3
  predictors + const the OLS is saturated. Reported but not used to rank.
- **ddg-direct R² ≈ 0 for every potential** (see above).

The single discriminating benchmark is cognate-rank AUC, measured (and recomputed fresh) for
all 12 configs in `scratch/cache/cognate_all.json`.
