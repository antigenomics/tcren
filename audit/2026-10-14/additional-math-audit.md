# Additional mathematical audit — 2026-10-14

Read in full: mechanics/springs.py (352 lines at audit start), topology/literature.py (262), topology/pose.py (336), score/explain.py (149), score/model.py (183), score/transform.py (215). Traced catalogue units, score/CLI attribution callers, tests/unit/{test_mechanics,test_literature,test_residue_deltas,test_score}.py, and contact-layer construction. No claim that all equations in the repository or manuscript are certified. No scientific benchmarks or model refits were performed.

## Confirmed defects repaired after parent authorization

### P1/P2: class-II target chain identities collapsed

Original `topology/literature.py:138-149` discarded chain.id.to while retaining only target category ('mhc') and residue index. Deduplication merged separate chain-A/chain-B residues with the same numbering; `participation_coefficient` at original lines 176-180 then grouped those distinct residues together. `contact_order` at original lines 204-214 also treated the two chain numbering axes as a single sequence.

Independent synthetic evidence: two distinct MHC residues, both index 5 but each contacted by a different single CDR loop, incorrectly produced partcoef_pmhc=0.5 rather than 0. Both edges from one TCR residue to the two same-index targets collapsed to one at the shared frame source. Both regression tests failed before the fix.

Repair: retain `ct` in the shared contact frame and deduplication; group pMHC participation by target+chain+residue. Contact order uses each chain's contacted span and averages normalized separations over supported loop–chain pairs. Legacy public-helper input frames lacking ct retain one-chain-per-target interpretation. This preserves single-chain semantics, apart from possible floating arithmetic roundoff due to operation order, while correcting multi-chain semantics.

**Release gate:** these corrections can change class-II co_mhc and participation descriptors; no catalogue digest, frozen fit, reference moments, cached scientific feature tables or published results were changed. Author-approved reproduction/model compatibility disposition is required before release. The existing class-I modelled-set observations were not rerun or weakened.

### P2: shear direction was not the constrained maximum

Original `mechanics/springs.py:144-153` projected the largest full 3D eigenvector onto the plane, then fell back to the *smallest* eigenvector if that projection vanished. Projection of an eigenvector does not solve the constrained Rayleigh quotient. With K=diag(1,2,3), docking axis z, it selected x (stiffness=1) instead of y (stiffness=2).

Repair: construct an orthonormal basis U for the axis-orthogonal plane, solve the 2x2 eigenproblem U'KU, and return U times its leading eigenvector. The regression failed before repair and now passes. Tensile path unchanged; shear and auto rupture outputs can change. No measured results were recomputed. Same release/reproduction gate applies to those changed readouts.

### P2: attribution selection changed the score mode

Original `score/explain.py:137-138` derived `receptor` from whether chain_types exactly equalled ('PEPTIDE',). Thus a reporting filter switched binder/channel scoring from receptor marginal to full peptide-presentation mode. Synthetic descriptors=0.5 gave score.full −14375.257201051412 with all residues and −15880.371011596068 with PEPTIDE filtering.

Caller trace: `cli.py:explain` documents --chain-types as a filter and forwards it, with no score-mode option; `tests/unit/test_residue_deltas.py` existing call uses score='peptide', which ignores receptor mode. No residue_deltas calls were found in benchmark/manuscript Python/Markdown/RST sources searched. Notebook or external downstream users were not exhaustively searched.

Repair: explicit keyword `receptor: bool=True` on residue_deltas, forwarded to the readout independently of chain_types. Explicit False recovers peptide-presentation scoring. Synthetic regression uses real frozen-model readouts with fake descriptor generation and verifies selection invariance plus explicit False. No corpus structure attribution was recomputed. Parent owns skill/docs/CLI integration; CLI remains consistently receptor mode by default.

## Remaining findings / qualifications

### P1/P2 scientific units: mechanics outputs are not calibrated SI

`mechanics/springs.py:94-99` defines k as 1, atomic multiplicity, or multiplicity/d_heavy^2 with distances in Angstrom. There is no calibrated energy factor or conversion to SI. The tensor at 111-114 therefore inherits k's chosen arbitrary scale; with default invdist2, its raw numeric dimensions are Angstrom^-2, force from lines 203-204 has Angstrom^-1, and work at 210 is dimensionless up to the unspecified energy scale. `descriptors/catalogue.py:569-576` labels these N/m, N, and J. Those labels cannot be justified by the implementation alone. Smallest safe repair is accurate reduced/arbitrary-unit documentation, not inventing a conversion constant. Changing the catalogue also changes its digest; parent must reconcile the frozen-model provenance deliberately, not silently. No catalogue change made here.

Also catalogue rupture_force says 'weaker axis' but shipped default is tensile (`springs.py:157`, `interface_mechanics`), and aniso's prose at catalogue line 572 describes along/across although the formula is across/along. These are definition/prose mismatches, not evidence against recorded empirical rankings.

### P2: rupture horizon does not ensure full rupture

`mechanics/springs.py:193-210` fixes max displacement at 2*max(rest), independently of break_strain and orientation, although the docstring promises integration until all springs break. For an aligned unit-rest spring and break_strain=3, final extension is only 2 and it survives; synthetic rupture still returns finite peak/work without indicating incompletion. Even default strain=.5 is not guaranteed for an initially antiparallel spring: at displacement 2L its length returns to L. A sufficient orientation-independent upper bound is t>(2+break_strain)*max(rest), by |v+t d| >= t-|v|. Exact individual break locations could provide a tighter bound, but changing horizon also changes grid spacing and historical outputs. Record and obtain model/reproduction disposition before altering shipped numerics.

`rupture` accepts arbitrary direction strings as shear, does not validate positive integer steps or nonnegative finite break_strain, and has no boundary check for zero-distance heavy-atom contacts in invdist2. These are input-validation gaps. Peak/work are finite-grid approximations with springs removed before force evaluation at crossed thresholds, not exact continuous rupture integrals. Existing self-check/test helper duplicates the rupture loop, so its success alone cannot independently validate that integration.

### P2: Transformer nonfinite contract fails for direct use

`score/transform.py:104-120` clips +/-inf to finite endpoints for logit/Fisher/atanh and clamps -inf for count/log transforms, before transform's line 181 nonfinite check. Synthetic Transformer(chain_balance), fitted to [.2,.8], maps input +inf/-inf to approximately ±3.882488, not NaN; NaN remains NaN. Main ScoreModel.coordinates rejects raw nonfinite values before calling Transformer, so ordinary score_table is protected. Public exported Transformer is not. Smallest repair: apply the input finite mask to transformed outputs. No change made here.

### Attribution omits noncontact residues but does not prove their deltas zero

`score/explain.py:19-21` claims any noncontacting residue moves no descriptor. The actual target restriction is implemented, but whole-chain docking/C-alpha geometry descriptors can depend on noncontacting residues. Existing `test_noncontacting_residue_moves_nothing` checks only peptide_score, a contact sum; it does not establish the claim for binder/pose or every descriptor. Describe attribution as restricted to interface-contacting residues, or explicitly define/test a wider target policy; do not present the limited test as a universal invariant. `_drop_residue` retains annotations/sequence labels intentionally while deleting residue atoms; downstream length/geometry conventions need separate validation for each descriptor, not an assumed 'physical deletion' interpretation.

### Literature attribution must distinguish adapted descriptors

`contact_order` is a mean pairwise sequence separation among target residues reached by a loop, normalized by *contacted span*, not the full-chain-normalized intrachain-contact order introduced by Plaxco. It is a valid explicit adapted statistic, but calling it precisely the published formula would overstate equivalence. Surface code appropriately distinguishes its height-field gap from a Voronoi/ASA gap-volume index. Citation metadata was not independently retrieved in this pass; no new bibliographic claims are supplied.

### Joint/model boundary validation gaps

`score/model.py:37-40` is the usual unbiased reliability-weight covariance only when weights are nonnegative, sum to one, and sum(w²)<1. Internal fit_holdout supplies normalized valid weights, but Joint.fit accepts arbitrary arrays without checking these conditions. `anomaly`/`anomaly_on` use `k or len(lam)`, so k=0 means all rather than zero directions, and negative k slices unexpectedly. `project(m)` can create m names but only p projected columns when m>p. Unknown subset names silently disappear in `_idx`; an empty subset yields float indices. Positive definiteness/logdet sign is not validated before inversion. These are lower-priority direct-API input-domain gaps, not observed failures in the valid frozen scoring route.

The unweighted Ledoit–Wolf shrinkage coefficient combined with weighted moments is explicitly an approximation, not the exact weighted optimum. Preserve that distinction. The covariance/transform no-negative training contradiction and confidence-augmentation PSD qualification are already recorded in /tmp/tcren-score-audit.md; this pass does not resolve them.

## Positive algebra checks and their limits

- Springs tensor is sum k*u*u', symmetric PSD for k>=0 and nonzero rest lengths. Trace is sum k when u is unit. Axis quadratic form and shear trace complement are correct; shear is the sum over two in-plane stiffnesses, not one shear-axis eigenvalue. Coupling counts are literal intersections of correctly keyed residue sets.
- Participation formula 1−sum module-share² and pose degree evenness (sum degrees)²/(n*sum degree²) are correct for their selected contact sets. `_spearman` excludes nonfinite pairs and guards insufficient/constant data. Favorability=−J and margin=cutoff−distance give the intended positive tight/favorable sign. OLS slope signs follow that same convention. These sign checks do not establish empirical discrimination or causality.
- All branches of score/model.py and transform.py were read. Gaussian log odds use class covariance blocks and the common normalizing constant cancels. Orthogonal projection mu'=W'mu, S'=W'SW preserves the full-rank posterior and yields a marginal when truncated. Conditional mean algebra is correct for positive-definite blocks. `_wmean_cov` agreed with numpy.cov(aweights=w) on three synthetic weighted rows. A two-dimensional conditional mean example returned exactly 4.0 as derived.
- Yeo–Johnson branches match scipy on seven signed inputs and six lambdas −3,−1,0,1,2,3; maximum absolute difference across these cases was 1.82e-11. This is a small numerical check, not exhaustive overflow/extreme-domain certification. Circular cos/sin expansion, fixed transformations, fitted per-coordinate standardization, and final clipping are traced. 'Every transform monotone' and 'only yeo carries fitted parameters' are imprecise prose: circular embedding is not scalar monotonic, and all outputs receive fitted location/scale.

## Verification and changed scope

Four mechanics/literature regressions failed before the repairs. After changes: `.venv/bin/python -m pytest tests/unit/test_mechanics.py tests/unit/test_literature.py -m 'not slow' -q` => 18 passed, 4 slow deselected, 0.32 s (120 s bound). Attribution regression failed before repair and passed afterward: 1 passed, 5 deselected, 0.33 s (120 s bound). Scoped git diff --check passed. Numerical snippets above ran in the repository venv; no scientific dataset fitting or evaluation.

Changed this stage: src/tcren/mechanics/springs.py; src/tcren/topology/literature.py; src/tcren/score/explain.py; tests/unit/test_mechanics.py; tests/unit/test_literature.py; tests/unit/test_residue_deltas.py. Report: /tmp/tcren-additional-math-audit.md. No generated scientific outputs, frozen arrays, catalogue units/digest, or reported benchmark measurements changed.

Not reviewed exhaustively: topology/surface/footprint implementations, descriptor extraction for every masked residue, native dynamics/relax kernels, scientific benchmark metrics/protocols, manuscript derivations/citations, version migration and full release eligibility. Parent is consolidating the wider audit.

## Final narrow follow-up

Parent authorized the remaining direct-Transformer nonfinite repair and attribution wording correction. `score/transform.py:_apply_raw` now replaces nonfinite inputs with NaN before any clipping or trigonometric transform. Finite-input arithmetic is unchanged. A synthetic regression spans all nine transform families and checks +inf, -inf, NaN, and a finite value; it failed before repair (Anscombe transformed -inf to a finite value; circular branches also emitted warnings), then passed. `score/explain.py` now explicitly restricts attribution to interface-contacting residues and makes no universal zero claim for omitted residues.

Focused validation: `.venv/bin/python -m pytest tests/unit/test_score.py -k transform -q -W error` => 3 passed, 21 deselected, 0.82 s, 120 s subprocess budget. Scoped diff check passed. This closes those two findings above. No model-domain expansion, fitted-array changes, scientific refits, or benchmark runs.
