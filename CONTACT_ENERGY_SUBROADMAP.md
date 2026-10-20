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
units using the actual published 2022 coefficients and current TCRen2. Median
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
Shuffled-receptor scoring now has a paired-CDR3 CLI route; fresh evaluation and
exact original-paper coordinate/protocol mapping remain open.

Both contact export and candidate CLI scoring now use bounded coarse workers
with all allocated CPUs. On 128 raw models, serial/parallel scores and contact
exports agree exactly; cold runtime changes from 28.331 to 8.428 seconds
(3.361-fold). This does not establish full-cohort scaling. Virtual residue
substitutions have no mutant atom coordinates, so atom/surface corrections for
those candidates remain undefined rather than silently set to zero.

## Work order and verification

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
   shuffle task varies CDR3 donors. The author recalls a Bobisse receptor-ranking
   panel: retain that as a source-identity question rather than silently substituting
   a different task. Older TCRpMHCmodels structures are a distinct generator regime.
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
four-type symmetric core has10cells, versus16 if directed;18symmetric types have
171. This is an architecture proposal, not a fitted/frozen candidate or an exact
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
are complete, missing producers are resolved, and any chosen frozen candidate has
all requested matched comparisons. Source review or a complete ranking-only run
cannot certify this broader goal. Keep scientific records private where required;
no manuscript findings or package default changes without explicit author review.
