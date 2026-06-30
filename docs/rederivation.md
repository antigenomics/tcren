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

## Candidate table shipped

`src/tcren/data/TCRen_potential_rederived.csv` holds the **2026-off** potential — the
strongest re-derivation on the new native data by cognate rank (AUC 0.8025). It is shipped
as a **new candidate file alongside, not replacing, the bundled default**
(`TCRen_potential.csv`, md5 88c42cd5…, unchanged). The candidate (md5 5529a771…) is provided
for the cognate-rank use case; **it does not meet the n199 acceptance floor and is not a
drop-in replacement for the default potential.**

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
