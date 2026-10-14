# Score mathematics and topology audit — 2026-10-14

Scope: code review and tiny synthetic checks; no scientific analyses rerun. Inspected score/{__init__,model,transform,fit,explain}.py, cohort.py, reliability.py, topology/graph.py, energetics/scoring.py, potential/model.py, mechanics/dynamics.py and their relevant tests/callers. Existing dirty work preserved. This is not complete coverage of every topology, mechanics, potential derivation or energetics module.

## Implemented confirmed fixes

1. **P1: wholly missing Q/T blocks became finite zero.** `src/tcren/cohort.py` q_score used `np.nansum`, giving 0 for every-NaN inputs; `reliability.s_score` consequently returned 0 for an entirely missing structure. Smallest repair preserves partial-row behavior but returns NaN for no finite coordinate. The shared zscore constant-reference guard also now preserves missing input instead of turning it into 0.
2. **P2: signs ignored in legacy equal-weight branch.** q_score returned before checking/applying signs. Moved sign validation before branch and applies signs in the mean as in the covariance path.
3. **P2: irrelevant missing descriptors suppressed marginal scores.** `ScoreModel.coordinates` formerly filtered all columns before binder_score removed non-receptor coordinates or channel_scores selected a family. Synthetic table with all descriptors=0.5 scored −14375.25720105; setting excluded Phi_pep_mhc=null yielded NaN, while dropping that column restored exactly the finite score. Added optional coordinate selection before row mask; receptor/channel readouts now use their own required columns. Actual required missingness remains NaN. No fitted constants changed.

Verification: regression tests failed before edits (three cohort tests; one score test). After edits: `.venv/bin/python -m pytest tests/unit/test_score.py -k 'not refitting and not reproduc' tests/unit/test_cohort.py tests/unit/test_reliability.py -q` => 43 passed, 2 deselected, 1.86 s. Each run bounded by subprocess timeout=120 seconds. Two refitting tests intentionally not run. Changes only src/tcren/cohort.py, src/tcren/score/__init__.py, tests/unit/test_cohort.py, tests/unit/test_score.py.

## Open findings requiring explicit disposition

### P1 — Tier-1 no-negative claim contradicts the fitted preprocessing

Evidence: `src/tcren/score/fit.py:63-78` fits `Transformer(names=cols).fit(X)` on ALL joined rows at line 69, before separating class 0 and 1 at line 74. `score/transform.py:124-146` fits nonlinear Yeo–Johnson lambdas and location/scales to that input. `score/__init__.py:166-172` (pre-edit numbers; pose_score docstring) explicitly says no negative enters. `confidence_residual` makes the same tier claim. Even though binder moments use only positive rows, their nonlinear coordinate system is learned from negative rows too; shrinkage toward identity and clipping make even scaling material. Binder membership itself also requires the labelled class split. This is not inference-time fitting, but it does defeat the stated negative-independent training contract.

Smallest safe disposition: state the actual training dependency, or separately design a binder-only transformation and regenerate/revalidate after author approval. Do not silently refit the released model or change standing claims/results. No new benchmark performed.

### P1 — No confidence training data silently creates a fake fitted residual

Evidence: `score/fit.py:82-92` initializes conf_mu/conf_cov/conf_var to zero and only estimates them when iptm exists and >200 finite positive confidence values remain. The fitter explicitly supports missing iptm. Those zero defaults are serialized normally; `score/__init__.py:confidence_residual` always computes a prediction and reports logit(iptm)-prediction. Hence a valid model without confidence training returns raw logit(iptm) as a purported geometry-conditional residual, instead of signalling unavailable. Also `_logit` clips ±inf to finite endpoints before `good = isfinite(ip)`, so infinite confidence can enter this fit as valid.

Smallest repair: persist/validate confidence-availability explicitly, return unavailable or raise for that readout, and test non-finite raw ipTM before clipping. Keep existing trained model unchanged. Model metadata currently has conf_var>0, so the shipped fit is not the missing-confidence example.

### P2 — Model catalogue digest is loaded but never checked

Evidence: `score/__init__.py:holdout_model` copies meta['catalogue_digest'] into ScoreModel without comparison. `cli.py:assess` runs `provenance.check(features_table)` against the installed catalogue, then loads arbitrary --model; it does not compare that model's digest. Thus a correctly stamped new feature table can be scored by a stale model sharing names but different meanings. The module docstring promises this raises. The current shipped digest equals registry_digest(), verified by inspection, so this is a guard failure, not a claim the shipped model is stale.

Smallest repair: compare model catalogue digest with current catalogue at loading/scoring boundary (including custom ScoreModel path if supported); one tampered-metadata test.

### P2 — Holdout joins can silently duplicate or omit training structures

Evidence: `score/fit.py:54-65` chooses first shared key then inner-joins with no uniqueness/cardinality/missing-key validation. Duplicate feature IDs or manifest IDs replicate rows; absent keys disappear. Only aggregate minimum size is checked. Invalid y values are cast to int and selected later without binary validation. This undermines a frozen manifest's reproduction contract and can bias class priors/epitope weighting while fitting succeeds.

Smallest repair: validate unique keys and labels, require manifest IDs be present (or explicit documented omission reporting), check join cardinality. No change to valid-input fitted arrays needed. Not exercised with actual holdout refit.

### P2 — Attribution residue filtering implicitly changes the score being explained

Evidence: `score/explain.py:137-138` passes `receptor=chain_types != ('PEPTIDE',)` to `_read`. `chain_types` is documented purely as a residue-reporting restriction. Asking to colour peptide residues therefore changes binder_score from the default receptor marginal to peptide-mode full coordinates, so score.full differs depending on which residues were requested. It is not simply a restriction of the same attribution.

Smallest repair: expose receptor mode separately and keep residue selection independent; test score.full invariance under target filtering. Existing whole-structure scientific outputs were not recomputed.

## Mathematical/physical qualifications, not asserted runtime bugs

- `topology/graph.py:79-81` claims exp(-D²/(2σ²)) tends to the binary 5 Å contact map as σ→0. For every nonzero distance it tends to zero, not an indicator D≤5. Correct the comment; no descriptor changes needed.
- `score/__init__.py:score_table` calls binder_score+logit(ipTM) a naive-Bayes combination. ipTM is not established as calibrated binder posterior or likelihood-ratio evidence, and independence is not established. The defined heuristic can be retained, but its probabilistic derivation is not justified by merely putting both terms on a logit scale. No request to replace its recorded empirical result.
- Confidence augmentation uses raw cross-covariance/variance with a shrunk descriptor block (`fit.py:78,91-92`). That does not guarantee the enlarged covariance is PSD for every fit. The CURRENT shipped conditional variance is positive: conf_var 0.436427574026391 minus explained 0.38871833530633054 = 0.04770923872006044. Therefore no actual shipped covariance failure is reported. Validate Schur complement for future fits rather than assuming joint consistency.
- `score/explain.py` asserts noncontact residues move no descriptor. That is not true for whole-chain geometry and Cα-distance-derived descriptors. Current output only reports contacted residues, which is defensible as a restricted attribution; claiming omitted residues have exactly zero contribution needs removal or proof.
- `cohort.q_score` docstring says it reduces to equal-weight mean for C=I, but the displayed/implemented formula is the sum; ranking is equivalent, numeric scale is k times the mean. Existing test checks correlation, not equality. Clarify prose, preserve frozen block scaling.
- mechanics/dynamics accurately disclaims MD/physical-time meaning. Its fixed partner shell is selected at initial coordinates (`_pack`), so it is an approximation when peptide moves outside the initial shell; no error bound is supplied. No sampler runs or model changes performed.

## Verified / not verified

Verified by code: Gaussian marginal means/covariance subblocks, conditional-mean algebra, orthogonal full-rank projection contract; graph cyclomatic count E−V+C and normalized Laplacian construction; directional contact-potential indexing in score_peptides and symmetric intra-peptide energy; transform missingness and clipping flow. Focused existing score tests passed.

Not certified: full topology/footprint/surface operators, potential derivation and redundancy estimator, rotamer mutation machinery, dirty new background scan/sidechain implementation (owned by another reviewer), native MC kernel, complete manuscript claims, calibration/generalization, scientific benchmark statistics, installation/release behavior. No new controls, nulls, sweeps, fitting or scientific reruns were launched.

## Follow-up implementation — supersedes open-findings statuses above

Parent authorized safe guards after reviewing the findings. The catalogue-digest, missing-confidence, and ambiguous-key/binary-label findings are now fixed in code; no scientific model was refitted or rewritten.

- `score/__init__.py`: `_logit` preserves NaN and ±inf as NaN; `ScoreModel.confidence_available` explicitly controls unavailable residuals; model loading checks registry digest and raises StaleTableError. Existing manually constructed models retain default availability=True. Legacy bundles without an availability field infer it from positive conf_var. Constant-confidence fits made by the updated fitter carry an explicit True flag and remain valid.
- `score/fit.py`: validates unique non-null keys on both inputs and binary non-null y before fitting; join declares 1:1 cardinality. Existing subset fitting is retained (the committed slice intentionally joins against the wider manifest), with unmatched feature/manifest counts and incomplete-row counts added to metadata. Confidence availability is recorded when >200 finite binder ipTM observations are available. All existing valid-input arithmetic is unchanged; only metadata and invalid-input handling differ.
- `topology/graph.py`: corrected the false sigma→0 hard-contact claim, without altering the Gaussian kernel.
- `cohort.py` docstring: corrected C=I to equal-weight sum, without altering numeric scale.

Tests: six new guard cases were observed failing before implementation. Combined score/cohort/reliability suite subsequently returned 49 passed, 2 model-refitting tests deselected, in 0.98 s. Three additional actual-fitter synthetic cases returned 3 passed in 2.58 s: 260 artificial observations, 210 positive/50 negative, two real catalogue descriptors, a 261-row manifest; cases have absent, infinite, and constant finite confidence. They verify availability and unmatched-row metadata. All subprocess test budgets were 120 s. These are new synthetic unit checks, not scientific reruns.

Exact files changed across both stages: src/tcren/cohort.py; src/tcren/score/__init__.py; src/tcren/score/fit.py; src/tcren/topology/graph.py; tests/unit/test_cohort.py; tests/unit/test_score.py. Parent owns documentation/release/archive integration. No attribution API change was made: the chain_types/receptor-mode coupling remains unresolved. Tier-1 pooled-transform contradiction likewise remains for author disposition; no fitted arrays or standing empirical conclusions were changed.

## MHC reference readiness repair — full-suite follow-up

Root cause: `data/mhc_cache/alleles_db.idx.dbtype` and `.idx.index` survived while the actual `.idx` data file was absent. `reference_db` previously checked only fresh `.dbtype` markers and therefore returned a broken database to both mapper and region-projection callers. Installed arda  mmseqs wrappers provide createdb/run/easy_search, not an index-ignore fallback; preserving the existing preindexed contract avoids rebuilding the prefilter on each annotation.

Changed only `src/tcren/mhc/reference.py` and `tests/unit/test_reference_db_concurrency.py` for this follow-up. `_db_complete` verifies sequence, header, and prefilter-index components, marker sizes, and data sizes against maximum offset+length in their index tables. Split numeric data parts are recognized. Existing build_lock retained. Build occurs in a temporary directory on the cache filesystem, createindex uses one thread, and incomplete output is rejected before publication. Previous derived DB files are moved into a local backup; readiness marker is removed first and new readiness marker published last. Existing concurrency test fakes now produce actual minimal database components instead of markers alone; the local FASTA fixture is synthetic.

Verification: three regression cases failed before implementation (missing index, truncated index, incomplete build accepted). Five readiness/concurrency tests passed after implementation. Actual local index rebuilt successfully in 3.848 s, `_db_complete=True`; old files backed up under `data/mhc_cache/.alleles_db-backup-2up2f072`. New index is about 915 MiB. Raw `database/mhc/alleles.aa.fasta` was neither changed nor downloaded. Actual `map_mhc` smoke on tests/assets/pdb/1ao7.pdb returned A=MHCa/MHCI and B=B2M/MHCI. Final readiness+region check: 8 passed in 0.90 s, 120 s budget. Actual index repair had a 180 s subprocess budget. Scoped git diff --check passed.

Limits: tests establish local thread contention behavior via the existing cross-process lock implementation, not independently verified multi-node/shared-filesystem semantics. Component checks catch missing/truncated files, not arbitrary same-size corruption or MMseqs-version incompatibility. No scientific benchmark rerun and no native-extension rebuild was performed by this agent.
