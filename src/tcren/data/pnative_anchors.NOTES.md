# `pnative_anchors.csv` — the reference set `p_native` anchors on

**273 real TCR:pMHC complexes, 161 epitopes, 27 HLA allele groups, 201 unique CDR3 pairs**
(206 crystallographic, 67 AlphaFold models of VDJdb-curated pairs). One descriptor row each, over
the four families `placement, interface, topology, energetics`.

## Why it exists

A finite mixture is identified only up to permutation of its components, so `p_native` needs
something that fixes which component means *native*. Anchors are extra rows of the design matrix
whose responsibilities are pinned; they are never scored. Before 2.13.0 the library shipped no
anchors and fitted unsupervised, while the published benchmark numbers came from a protocol that
lived outside the package — so the two disagreed. This table makes them the same object.

## How it was selected

`bench/scripts/build_pnative_anchors.py` in the TCRen2 benchmark repository, from
`huggingface.co/datasets/isalgo/vdjdb_structure_models`. Deterministic; **no random seed exists**.

1. MHC class I, and the epitope must not appear in any benchmark cohort (23 epitopes, read from the
   cohort tables rather than typed, so adding a cohort cannot silently create leakage).
2. One row per (epitope, CDR3α, CDR3β) — clonal expansion is heavy in this source, one pair recurring
   34 times. Ties are broken on `tcr_pmhc_hash`, which is content-derived and therefore stable.
3. Capped per epitope so no epitope dominates: the largest is 2.9% of the set. Crystals supply the
   epitope breadth, the models the depth.
4. Rows with any missing descriptor in the compact feature set are dropped (10 of 283) — an imputed
   anchor is a fabricated constraint.
5. Asserted disjoint from every scored structure by hash, not only by epitope.

`epitope` and `source` are carried for provenance; neither is read when fitting.

## What this set does and does not reproduce

Measured 2026-08-27 on the two TCRen2 receptor-ranking benchmarks, macro ROC-AUC over cohorts:

| anchors | TCRvdb | VDJdb balanced |
|---|---|---|
| other cohorts of the same benchmark, **both classes** (the published protocol) | 0.832 | 0.718 |
| other cohorts of the same benchmark, **positives only** | 0.725 | 0.574 |
| **this table** (273 external positives) | 0.660 | 0.494 |
| none (unsupervised) | 0.607 | 0.497 |

So the published numbers rest on two things this table does not supply: **negative** anchors, and
anchors drawn from the **same generator and pipeline** as the structures being scored. Pinning which
mixture component means *native* needs only one class; the negatives in the published protocol are
doing discriminative work, which makes it semi-supervised transfer rather than anchoring.

This table is therefore offered as `anchors="auto"`, not as the default, and it is not the protocol
behind any published TCRen2 number.
