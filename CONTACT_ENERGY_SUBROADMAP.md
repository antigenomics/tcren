# Contact energy: sparse reference models and complete task audit

Author direction, 20 October 2026. Preserve original TCRen and all accepted
baselines. Improve discrimination of supported and mock modeled complexes with
contact information beyond a closest-atom 5 Å residue edge. Keep every measured
gain and loss, including poorly performing tasks. No assay-fitted combination
weights, score-direction selection, threshold scans or silent default promotion.

## Evidence and immediate decision

The frozen residue+atom layer was evaluated on 3808 models and both residue
potentials, with 7616 original scores exactly preserved. Its 248 case/potential/
metric comparisons include original and both completed surface references.
ILA1 and YLQ improve; SB27 loses. Both potentials gain on two of seven CPL
AUROCs and lose on five. Template-covered and template-scarce receptor tasks
have different outcomes and must remain separate. This is a candidate result,
not a replacement for original TCRen or a universal classification restoration.

The current reference has 20927 observed directed residue+atom cells from
48827 native and 211517 antibody peptide-target atom contacts: 3438 cells
are singletons, 1024 native-only, 10338 antibody-only and 9565 shared.
Physical contacts are correlated within complexes. These counts are not that
many independent structures or freely supported parameters. Adding distance,
angle, residue identity and context axes would multiply sparsity.

Completed follow-through: all nine exported energy columns now have the recorded
ATLAS/MD endpoint audit (72 Pearson/Spearman entries). Fresh contact-arm exports
retain 176 matched correlations on 179 successful native structures, with ten
explicit eligibility failures. On the common 178-structure MD population, the
atom layer lowers Pearson correlation from 0.453750 to 0.339511 for TCRen2.

The first legacy gate scored 382764 candidate/potential rows across 222 ranking
units using the bundled classic coefficients and current TCRen2. The reference
identity amendment below distinguishes that May 2026 table from the paper artifact. Median
within-unit AUROC changes from 0.927000 to 0.923963 for IEDB (204 units; 79 gains,
108 losses, 17 ties), and 0.738095 to 0.726190 for Bigot (14; 5 gains, 7 losses,
2 ties). Bobisse retains rank 1/79; all three corrected Birnbaum AUROCs improve.
The largest IEDB unit loss is 0.865 to 0.124. These results do not pass a blanket
legacy nonregression claim.

The notebook-protocol follow-up adds 630924 fresh scores across 315 ranking
units. Median AUROC changes from 0.9405 to 0.9375 for the exact shuffled-RNG
random pool (218 units; 96 gains, 111 losses, 11 ties), 0.793 to 0.766 for modeled
complexes (47; 21 gains, 24 losses, 2 ties), and 0.956 to 0.930 for their crystal
counterparts (47; 16 gains, 26 losses, 5 ties). All three original-protocol
Birnbaum AUROCs improve. Notebook model/crystal arms draw distinct historical
pools; each compared matrix uses identical candidates within its arm. Their
original 2022 population and LOO/LCO outputs remain separate evidence.
Shuffled-receptor scoring now has a paired-CDR3 CLI route. Its fresh 470-pair
comparison retains 130 native and 340 shuffled candidates over 130 templates;
all 260 donor CDR3 sequences match fresh annotation. Pooled summed AUROC changes
0.730430 to 0.724502; within-template AUROC has 16 gains, 14 losses and 100 ties.
The TRA and TRB pooled AUROCs improve separately; no adverse sum is concealed.

The original-paper gate is now source-linked and scored separately. Older Git
source constructs the literal ERGO export from the same candidate object as the
structural score; no execution record certifies it as the final published random
stream. On 137 recovered original-panel pools, median AUROC changes from 0.960
to 0.937 (38 gains, 84 losses, 15 ties). The 47 original model/crystal comparisons
now use identical recovered pools: models 0.783 to 0.769 (20/25/2), crystals 0.945
to 0.931 (13/30/4). Historical scope defects remain visible as diagnostics.

The original junction shuffle retains all 302 source pairs in its ledger; 282
are scoreable (132 native, 150 shuffled). Pooled summed AUROC changes 0.698965
to 0.680025; TRA improves 0.628662 to 0.637929, TRB falls 0.655884 to 0.643965.
Within-template AUROC has 4 gains, 9 losses and 62 ties over 75 evaluable pools;
57 positive-only pools remain explicit nulls. Actual adjacent Cys and F/W anchors
are recovered from annotation, never fabricated. Twenty unavailable pairs retain
source/topology/peptide reasons and are excluded identically from both matrices.

Of 136 original IEDB-eligible units, 28 complete standard-amino-acid pools can be
recovered without an unrecorded random sample: median AUROC 0.948745 to 0.922592
(6 gains, 17 losses, 5 ties; 18090 candidates per potential). Another 107 historical
samples have no recoverable sampling state; one complete pool contains 86
nonstandard controls and remains wholly unavailable. This subset is not the
full paper. Original Birnbaum SUM aggregation is now retained exactly, including
its duplicated native peptide: all three AUROCs improve, while precision and
median positive rank worsen in two units. No random seed is fitted to published
outputs.

Both contact export and candidate CLI scoring now use bounded coarse workers
with all allocated CPUs. On 128 raw models, serial/parallel scores and contact
exports agree exactly; cold runtime changes from 28.331 to 8.428 seconds
(3.361-fold). This does not establish full-cohort scaling. Virtual residue
substitutions have no mutant atom coordinates, so atom/surface corrections for
those candidates remain undefined rather than silently set to zero.

A completed support-only audit maps 184 observed residue/atom labels to 16
existing chemical signatures. Unordered pooling reduces 256 possible directed
cells to 136 unordered cells (135 observed). Native support covers 133 cells,
five seen in only one PDB; antibody support covers 135 cells, none seen in only
one PDB. This supports antibody augmentation and count-level symmetry, but does
not validate the existing heuristic atom flags as a force field. Carbonyl versus
aliphatic carbon, sulfur chemistry and histidine protonation remain mapping
questions. A subsequently specified 19-type topology partition separates carbonyl
sites and sulfur classes and declares histidine state unknown. Its support audit
finds 187 of 190 unordered cells observed: native 180 (six single-PDB cells),
antibody 187 (two single-PDB cells). The estimator, endpoint-composition reference,
unseen policy and aggregation were frozen before evaluation; endpoint composition
is explicitly not exposed opportunity.

The next numerical definition is frozen in [CHEMICAL_CORE_SPEC.md](CHEMICAL_CORE_SPEC.md).
It uses unordered actual-count enrichment, directed residue-conditioned endpoint
composition and an exact conditional partition. Unseen core W=1 gives logZ;
unmapped or zero-conditional-support queries abstain with zero. The optional
`chemical_core` CLI arm preserves the earlier `atom_count` definition and both
original residue potentials. Focused hand-counted, source-validation, CLI CPU
mode and actual serial/spawn checks precede any numerical resource estimation.
No assay result has been used to choose this estimator or its missingness rules.

## Completed symmetric-core evaluation

The frozen 19-type core now has 187 observed of 190 unordered cells from the
same 48827 native and 211517 antibody observations. Its 3808-model ranking run
preserves all 7616 original and directed atom-count scores exactly. Independent
completed-output checks cover all 248 metric entries, raw hashes, cohort
membership and physical edge sums. Original missing edge values are preserved:
42 of 176146 exported edge rows, representing eight models under each potential.
No score or model is promoted.

The 190-cell core is not the whole estimated reference: its frozen binder and
target conditional compositions each contain 128 populated cells in 20
normalized residue rows. All components share observations; these dimensions do
not establish independent fitted-parameter counts. The 400 residue-pair log
partitions are derived, not separately fitted. Their contributions vary with
contacted residue identities, so removing them after a loss changes the model.

For TCRen2 AUROC, the new core gives SB27 0.806968 (original 0.861510,
previous surface 0.927116, boundary surface 0.918672, directed atom 0.646254);
ILA1 0.717088 (0.591236, 0.741744, 0.680124, 0.964096);
YLQ 0.460255 (0.477427, 0.364654, 0.347459, 0.608226);
and GLC receptor ranking 0.370433 (0.379849, 0.431576, 0.436849, 0.422913).
The respective positive/negative model counts are 161/64, 161/164, 174/249 and
135/59. SB27 improves against the sparse directed atom layer but loses against
all three older references. ILA1 retains an original-score gain while giving up
most of the directed layer's gain. Both receptor-ranking epitopes lose against
the directed layer; the GLC loss also includes both surfaces and original.

Across seven CPL AUROCs, core versus original has four gains/three losses for
TCRen2 and five/two for the bundled classic TCRen. For six template-covered VDJdb cases,
both potentials have three gains/three losses versus original. For sixteen
scarce-template cases, both have eight/eight. Against the directed layer,
covered cases have two gains/four losses for TCRen2 and one/five for TCRen;
scarce cases have six/ten and five/ten/one tie respectively. These strata remain
separate. Full per-case AP and precision-at-recall gains and losses are retained
in private evaluation records alongside every AUROC.

The ranking run finished in 182.727 seconds on a 16-CPU allocation with sampled
peak process-tree RSS 6.736 GB and 81 live OS threads. Kernel runnable-thread
count reached 34; this does not verify a 16-thread numerical concurrency bound.
Native follow-through preserves 179 successes, ten failures and all original/atom
scores. Its 16 correlations use exact accepted populations: all six ATLAS Pearson
entries rise against original, but four of six Spearman entries fall (23 koff,
55 Kd/dG observations). All four MD correlations lose against original and both
surfaces (178 structures); TCRen2 Pearson changes 0.453750 to 0.378634,
Spearman 0.455755 to 0.372674. Core beats the directed atom layer in 12/16
native entries, losing all four koff entries. No direction is selected. Virtual legacy residue substitutions
still lack candidate atom coordinates; no atom corrections are imputed for them.

## Original-paper reference and remaining legacy panel

The spring R-only root table is byte-identical to the December 2022 artifact;
all 380 values match saved Figure 3. The bundled classic table is the May 2026
updated reference, with a different native population. Earlier numeric results
remain valid for that recorded table, but their `published2022` label did not
establish an original-paper comparison. Root, pipeline, full/LOO/LCO and
redundancy/debiasing variants retain separate identities; no default is replaced.

The recovered original CSV was scored explicitly, without fitting, on the same
191382 locked candidate instances in the first legacy batch. Current scores
were read unchanged. Across 204 IEDB units, median AUROC is 0.907500 for the
original table and 0.923963 for TCRen2 (101 gains, 92 losses, 11 ties). Bigot's
14 units change 0.750000 to 0.726190 (5/7/2); Bobisse remains rank 1/79.
Corrected Birnbaum has two AUROC gains and one loss across three units. These
are full frozen tables, not a reconstruction of paper LOO/LCO or sampled pools.
A separate saved-paper-score comparison confirms identical 602 Bigot candidates
and the same 5/7/2 AUROC outcome; Bobisse's 79 saved-paper LOO scores retain
perfect ranking. No current score, matrix, threshold or candidate was changed.

The remaining locked protocols now add 521641 original-reference scores and
7920 metric entries, with current and May-reference results preserved exactly.
For 137 recovered original random pools, original-paper versus current median
AUROC is 0.973 versus 0.937 (34 gains, 90 losses, 13 ties); for 47 modeled
complexes it is 0.848 versus 0.769 (17/29/1). The separate notebook crystal
protocol changes 0.970 to 0.930 over 47 units (9/33/5). The 28 recovered original
IEDB pools change 0.952149 to 0.922592 (6/19/3). Original Birnbaum SUM changes
0.922882 to 0.936294 in median AUROC, but has one gain and two losses across
three units. Its duplicated native occurrence remains summed twice. Pooled
CDR3-shuffle AUROC rises 0.684253 to 0.724502; original junction-shuffle AUROC
falls 0.718763 to 0.680025. The latter retains 20 unavailable pairs and 57
positive-only pools with undefined discrimination metrics. These losses prevent
a blanket legacy nonregression claim; no score direction or default is changed.

The paper also describes a large-scale VDJdb homology-model benchmark of 6150
paired TCR–epitope cases. It is distinct from the 47 modeled/crystal complexes.
The saved multi-method Figure 4c ranks lack TRA/HLA/model identity and include
method-dependent row counts and duplicate beta-CDR3/peptide keys. Its exact
candidate/method populations and source linkage remain open; do not certify the
whole legacy paper from the already completed smaller panels.

The current-task audit additionally closes 112 graded-response correlations for
all contact arms and 384 previously omitted archival-energy metric cells, using
saved scores and exact locked populations. Native-reference-weighted totals
remain excluded from final models even when they gain on CPL. Experimental
primary-screen code defines a GLC-over-YLQ log2FoldChange contrast. Exact export
lineage into the restricted benchmark table remains unverified; retain raw signed
associations and distinguish source-supported interpretation from a proved table
mapping. Never choose an epitope-specific sign from prediction accuracy.

## Deferred model-storage harmonization

Ongoing global VDJdb modeling and manuscript-specific selections have separate
owners. Keep a single canonical coordinate payload with versioned cohort/source
manifests; preserve the balanced positive/negative manuscript selection and its
labels. Audit exact byte identity before deduplication: matching model IDs can
refer to different coordinate payloads. Storage migration is deferred behind the
energy audit; no new complete mirror, file deletion, label revision or model
regeneration is part of this step. Private source paths, overlap manifests and
consumer migration gates remain in the research roadmap.

## Scientific priority — author clarification, 21 October2026

The unavailable 6150-case historical registry is a documented reproduction issue,
not a blocker for improving the scientific model. Use the available historical
coordinate sets and retain exact limits on what their comparisons establish;
do not wait for an absent registry or replace it with modern generated models.

Prioritize whether atom patterns add information conditional on an existing
residue contact. The completed contact-class statistics permit an exact
residue-contact-density versus within-residue atomic-occupancy attribution,
with unchanged populations and explicit opportunity/zero-contact states.
This is forensic evidence, not a trained contact filter or an accuracy claim.
Keep the original score and both frozen atom candidates intact.

The next reference question is whether reference counting and query aggregation
use the same physical observation unit. Existing fitting counts physical atom
contacts while querying averages atom terms per residue edge. Audit that declared
choice using native/antibody evidence before estimating another table; no
assay-selected weights, added distance regression or discarded losses follow.

## Work order and verification

Additional author-requested audit, 21 October: compare binder-only, shared and
nonbinder/mock-only residue and atom contacts on the locked 1706 VDJdb/TCRvdb
models. Keep exact epitope/MHC cases, template coverage, receptor regions and
both TCRvdb epitopes separate. Full observed atom inventories supply amino-acid/
atom availability denominators, including zero-contact sites. Equal-complex
densities and equal-case overall effects address class imbalance; paired-junction
clusters block repeated receptors. Retain exact support, missing opportunities,
all contrasts and non-estimable uncertainty. These are observed class sets, not
proof of chemically required/forbidden contacts. Do not score the discovery
labels with a filter learned from those same labels. The optional
`contact-class-audit` CLI exports these diagnostics without fitting a potential
or altering any accepted score; its plan is frozen before cohort extraction.

Completed audit: all 1706 models in 24 exact cases, with zero failures, 1676
paired-junction variance blocks and 93942 observed chemistry motifs. Final
evidence contains 2556705 feature/case rows and 375768 overall chemistry rows.
Among 468285 class-exclusive feature/case rows, 62241 lack counterclass endpoint
opportunities; 406044 retain eligible counterclass models with zero contacts.
These are descriptive motif/case counts, not independent contacts or predictive
gains. Fixed 128-feature blocks preserve the full cohort per motif and use the
CPU allocation in bounded worker windows; the completed 16-CPU run took 471.153
seconds with sampled peak tree RSS 9539665920 bytes. Private reports retain all
sets, effects, supports, uncertainty states and execution identities. The next
decision is whether a physically interpretable contact rule can be frozen for
independent evaluation; this audit does not authorize evaluating a learned
filter on its discovery labels or changing the accepted energy defaults.


1. **Publish safely.** Keep raw inputs, private provenance, host paths and private
   evaluation records out of the public payload and its ancestors. Preserve the
   private research history. Verify the public branch and both host tree hashes.
2. **Audit every existing task before another model.** Retain every energy score,
   unit, direction, potential identity, population, missing-value policy and metric.
   Include CPL, both TCRvdb epitopes, each VDJdb coverage/epitope stratum, ATLAS
   kinetics/affinity, MD interfaces and raw confidence diagnostics. Separate the
   earlier complete five-block run from the newer ranking-only contact arms.
   An untested score/task cell is a gap, never a successful audit or a zero effect.
3. **Restore the legacy 2022 gate.** Inventory cognate/random, cognate/IEDB,
   Birnbaum, Bobisse, shuffled receptors, modeled/crystal and Bigot tasks. Match
   the exact published potential, historical rederivation, full-set and LOO/LCO
   arms separately. Preserve the original per-unit ranking endpoints; pooled ROC,
   median within-unit ROC and rank-derived AUC are different quantities. Compare
   every current candidate with original legacy TCRen on identical candidate sets.
   Report material regression by task; do not lower floors or hide a bad task.
   No numerical meaning of “not much regression” is invented before author review.
4. **Resolve source/task discrepancies.** Existing Bobisse records contain one
   receptor and 79 candidate peptides over seven HLA/length templates; the stored
   shuffle task varies CDR3 donors. The primary Bobisse study also describes four
   Modeller receptor complexes on one fixed neoepitope/HLA, using alternative alpha
   chains. Recover their paired sequence definitions and historical coordinates
   before calling that separate receptor comparison evaluated. Older
   TCRpMHCmodels structures are a distinct generator regime.
5. **Reduce the atom representation before estimating.** Verify a published,
   chemically explicit atom mapping and its reuse terms. Audit support by unique
   complex and reference family, not only atom counts. Keep antibody and native
   source counts visible. Freeze the mapping, reference measure, weighting,
   unseen policy and aggregation before using assay labels. Do not average or
   symmetrize an existing energy matrix as a substitute for recomputing counts.
6. **Specify one minimal next candidate.** Chemical grouping may reduce the table
   below a 400-cell directed residue matrix: four directed chemical types give
   16 cells; 18 symmetric chemical types give 171, or 324 directed cells. These
   are table sizes, not independent degrees of freedom. Choose from actual support
   and chemistry; no menu is selected by assay accuracy. Preserve the original
   residue score and test whether the additional term supplies distinct evidence.
7. **Close execution gaps, then evaluate.** Reuse CLI-owned quantities and the
   existing bounded runtime. Extend missing candidate-aware legacy threading and
   all-task exports rather than execute old scripts that overwrite fixtures or
   reuse calculated rows. Fresh authorized scoring rebuilds from raw inputs, uses
   a frozen wheel/reference and all allocated CPUs with one bounded compute layer.
   Verify serial/parallel numerical equality, actual process-tree threads, scaling,
   full membership and all gains/losses. Do not rerun completed scores just for prose.

## Symmetry, transfer and possible score structure

A chemical pair can have exchange symmetry while its observed binder→target
statistics remain asymmetric. Pool unordered physical contacts once, with a
matching reference measure, if estimating a symmetric chemical core. Pooling
counts and then taking a logarithm differs from averaging two log-potentials.
Aliasing equivalent atom names within each amino acid is not cross-residue
chemical pooling and does not remove most of the current sparsity.

Antibodies can supply chemical support; they do not make pMHC pocket contexts,
paratope composition or template-generated poses exchangeable with TCRs. Their
larger contact total currently dominates the pooled reference. Any future
redundancy weighting or prior strength is an estimated/selected dependency to
expose and freeze, not an assay-tuned coefficient or invented effective sample size.

A useful *proposal*, not an implemented score, is

    E = E0 + Σ_edges mean_atom_pairs Δu_chem + Σ_interface_atoms ψ_environment.

E0 is the unchanged dimensionless residue contact score. Δu_chem is a dimensionless
chemical log-ratio residual, with a declared opportunity/reference state; it must
be conditional enough to avoid counting E0's composition signal twice. ψ is an
atom/site term for collective features such as burial with unsatisfied polarity,
coordination saturation or solvent response. Its units and coefficient must be
compatible with E0 before addition. Physical kcal/mol terms require an explicit
conversion and cannot be mixed numerically with log-odds. No ψ, new coefficient,
angular table or distance regression is fitted or evaluated by this document.

For joint observation z and context h, the exact identity is

    -log[P(z,h)/Q(z,h)] = -log[P(z)/Q(z)] - log[P(h|z)/Q(h|z)].

P is the observed reference distribution and Q the chosen opportunity/reference
measure, on identical support. A compact chemical-core plus directed-context
model approximates the second factor; an unrestricted conditional table restores
sparsity. A sum of independent marginal penalties is not automatically a joint
likelihood or a binding free energy. Explicit many-body site terms need evidence
and sufficient inputs; they are not justified simply by increasing model complexity.

A compact conditional version can preserve role-specific opportunity without
training another large residue+atom table. Let qB(c|a,h) and qT(d|b,h) be normalized
chemical-type opportunities for binder residue a and target residue b in context h.
With a shared W(c,d)=exp(-u(c,d)), define

    p(c,d|a,b,h) = qB(c|a,h) qT(d|b,h) W(c,d) / Z(a,b,h),
    Z(a,b,h) = Σ_c,d qB(c|a,h) qT(d|b,h) W(c,d).

The atom log-ratio correction is u(c,d)+log Z(a,b,h), with dimensionless u and
normalized reference opportunity. Symmetry of u need not symmetrize qB/qT. A
four-type symmetric core has 10 cells, versus 16 if directed; 16 types have
136 symmetric cells, and 18 types have 171. One global additive gauge cancels
after conditional normalization when the opportunity distributions are fixed.
This is an architecture proposal, not a fitted/frozen candidate or an exact
reformulation of the current table. Actual noncontact opportunities are needed;
the antibody contact table alone cannot supply them. Conditional normalization
and edge averaging still do not establish a physical or independent likelihood.

## Verified primary literature and implications

- [PIPER (arXiv:q-bio/0605018)](https://arxiv.org/abs/q-bio/0605018): chemical atom
  grouping and a symmetric pair matrix, but antibody transfer is imperfect. Its
  18-type mapping is delegated to ACP; the inspected preprint is not a complete
  mapping or an established redistribution license. Its 6.5 Å definition is not
  interchangeable with the present 5 Å counts.
- [Application of asymmetric statistical potentials to antibody–protein docking](https://doi.org/10.1093/bioinformatics/bts493):
  a directly relevant sparse-data example uses four pooled atom groups and preserves
  role asymmetry. It argues against symmetrizing all interface statistics merely
  because microscopic pair interactions can be symmetric.
- [Orientation-dependent hydrogen bonding](https://www.physics.rutgers.edu/~morozov/pdfs/Kortemme_Morozov_HB_JMB2003.pdf):
  donor–hydrogen–acceptor and acceptor-base geometry supply chemistry beyond a
  radial cutoff. Existing N/CA/CB orientation descriptors are different quantities;
  current contact-only tables may lack bonded-neighbor coordinates. No Rosetta
  weights, hydrogen reconstruction assumptions or physical energy units are imported.
- [Reference-ratio formalism (arXiv:1008.4006)](https://arxiv.org/abs/1008.4006):
  the denominator defines the statistical model's interpretation; contact endpoint
  frequencies do not measure exposed-atom opportunity.
- [Anisotropic potentials (arXiv:physics/0302009)](https://arxiv.org/abs/physics/0302009):
  distance and orientation can contain additional structural information, shown
  for fold decoys rather than these TCR assays. This is not evidence for adding
  poorly supported distance/angle regressions here.
- [Cooperative hydration (arXiv:2206.09070)](https://arxiv.org/abs/2206.09070) and
  [solvent correlations (arXiv:2306.09392)](https://arxiv.org/abs/2306.09392):
  solvent-integrated behavior need not be pair-additive. These motivate a separate
  environment audit; they do not provide calibrated TCR score coefficients.

## Completion gates

The goal remains active until the task/score coverage audit and legacy comparisons
are complete for available inputs, historical asset limitations are explicit,
and any chosen frozen candidate has
all requested matched comparisons. Source review or a complete ranking-only run
cannot certify this broader goal. Keep scientific records private where required;
no manuscript findings or package default changes without explicit author review.

## Within-residue atom-pattern candidate —21 October2026

The author authorized implementation and frozen evaluation of an edge-balanced
native/antibody atom-pattern potential. `ATOM_PATTERN_SPEC.md` owns its exact
observation measure, equations, support and comparison contract. Each reference
residue contact contributes one unit divided over its physical atom pairs;
directed conditional atom-mark log ratios are averaged once per query CDR edge
and added to unchanged original scores with coefficient one. No assay/source
weights, new distance regression, blacklist or default promotion.

The frozen resource precedes evaluation. Preserve every current ranking case,
potential, metric and template stratum; verify exact original-score identity and
membership before attaching both completed surface baselines. Native energy and
graded-response endpoints require their own accepted membership checks. A
ranking-only result cannot establish all-task completion or legacy threading
capability. Private evidence retains exact resource/input/output hashes and all
gains and losses; do not place private provenance or benchmark rows here.

## Reconciliation after the edge-balanced comparison —21 October2026

The author requests a completed-output mechanism audit and a correction-versus-
joint-score design discussion, without structure/scoring reruns or assay fits.
`ATOM_PATTERN_RECONCILIATION.md` is the public-safe equation/design record.
The private evidence ledger retains all resource-cell, source/marginal, residue-
edge, model and actual changed-pair attributions; no new candidate is evaluated.

The atom-count and edge-balanced arms have the same observations, equivalences
and support. Distinguish changed cell counts, pair normalizers and both directed
marginals. Native/antibody attribution is bookkeeping under its declared path
and count gauge, not an ablation or independently evaluated source-only model.
Existing full atom observations support exact reconstruction where present;
compressed CPL edge means leave some source terms unidentified. That data issue
does not block the rest of the scientific goal or authorize a structure rerun.

An important exact identity is E_joint^S=E_res^S+D_S when both factors use one
reference measure S. Keeping an older E0 instead gives an explicit residue-prior
difference E0-E_res^S. Also E_P[D|a,b]=-KL(P||qBqT), generally not zero: conditioning
alone does not remove residue-dependent mean preferences. Never interpret this
as evidence that assay-selected centering, another coefficient or dropping E0
would improve recognition. Such alternatives are not evaluated by this audit.

Next scientific order:

1. Retain original scores, both surfaces and all tested atom/core arms, with
   identical-population gains/losses. Do not promote a weaker arm by relabeling it.
2. Define the conditional state before estimating more cells: available atoms,
   count/multiplicity, whole contact pattern, receptor region and peptide site.
   State which quantities are actually observed in each reference/query source.
3. Specify the opportunity/reference denominator on the same support. Antibody
   transfer may support chemical interactions without supplying TCR/pMHC pocket
   marginals. Do not equate native-like packing with functional recognition in
   complexes generated by the same modeling pipeline.
4. Compare coherent joint rebuilding and a deliberately retained residue prior
   algebraically before any new evaluation. Atom summation versus edge averaging
   changes the observation unit; it is not a purely cosmetic implementation choice.
5. Freeze any later estimator/data extension separately before evaluation, retain
   all-task and legacy coverage, and fit no assay/source weights. The present
   authorization continues the analysis/design; it adds no new scientific run.
