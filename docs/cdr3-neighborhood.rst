CDR3 chemical-neighborhood diagnostics
======================================

These experimental rules preserve the prepared score as a comparator and do not
change ``screen`` defaults. Fresh ``export --quantity cdr3-neighborhood`` reads
all heavy-atom pairs within 5 Angstroms on CDR3:peptide and CDR3:MHC interfaces.
MHC includes the complete annotated MHC chains, rather than helices alone.
CDR3 excludes the conserved junction anchors.

Each physical residue edge contributes once. The existing v2 classifier supplies
salt-bridge, hydrogen-bond, cation-pi, stacking, aromatic and hydrophobic flags.
A supported atom witness uses a CDR3 side-chain atom and has neither a severe
Bondi-radius overlap (greater than 0.6 Angstrom) nor same-sign charged endpoints
within 4 Angstroms. These are distance/type heuristics, not angular hydrogen-bond
validation or binding labels. An edge can have both a supported witness and a
conflicting witness. No potential is estimated or missing energy entry filled.

For supported witnesses let ``d`` be atom distance and ``r`` the sum of Bondi
radii, both in Angstroms. Packing is ``min(d/r,r/d)`` in [0,1]. A residue edge
takes its maximum supported packing; unsupported edges contribute zero. This
avoids counting a residue interaction more often because it has more atoms.

Run ``neighborhood-rules --features features.tsv --models models.tsv -o scores.tsv``.
The model table supplies unique ``pdb.id``, ``iptm``, ``topology`` and
``confidence_weight`` columns. These reconstruct prepared as
``S = w*iptm + (1-w)*topology``. The command reads no outcome labels, preserves
all model rows, exports prepared and all ten fixed diagnostic rules, and records
input hashes and rule definitions. Finite inputs and bounded components are required.

Six soft rules retain ``w*iptm`` and replace the topology term by
``(1-w)*topology*Q``. Their respective Q values are: supported peptide-edge
fraction; mean best peptide-edge packing; one minus conflicted peptide-edge
fraction; one minus conflicted MHC-edge fraction; one minus conflicted MHC edges
at residues without peptide support divided by all MHC edges; and the product
of peptide/MHC conflict-free fractions. Empty conflict denominators imply Q=1;
empty peptide-support denominators imply Q=0. Local anchoring requires support
at the same physical CDR3 residue, not merely elsewhere in its chain.

Four hard gates set the *entire score* to zero on failure: neither CDR3 contacts
peptide; either CDR3 lacks any pMHC contact; neither CDR3 has chemical peptide
support; or either CDR3 lacks chemical pMHC support. Thus these hard-gate arms
also suppress confidence. No thresholds or mixture coefficients are selected
from query labels. Per-cohort nonregression must be evaluated explicitly;
no rule is a guaranteed-negative classifier.

.. automodule:: tcren.topology.neighborhood
   :members: neighborhood_counts, neighborhood_features, neighborhood_scores

Detailed contact evidence
-------------------------

``export --quantity cdr3-contacts --interface all`` retains every heavy-atom
CDR3 contact witness within 5 Angstrom to peptide or complete MHC, with separate
interface labels, independent v2 chemical flags and ``sidechain.from``.
Model-level contact enrichment must count each feature once per model.

Directed TCR–MHC potential
--------------------------

``derive-potential --structure-dir Canonical2026 --interface tcr_mhc --balance
both -t 0 -o TCRen_m.csv --contacts-out TCRen_m_contacts.csv`` derives the separate
asymmetric TCR–MHC matrix with the existing TCRen classic estimator. Rows are
TCR and columns MHC. The 5 Angstrom contact cutoff, pseudocount 1, separate
contact marginals, redundancy balancing, alpha-beta eligibility and withdrawn
entry exclusions match the current TCRen derivation. All TCR regions and complete
MHC chains use the existing contact-map definition; this is not CDR3-only.
As in classic TCRen, receptor cysteine is removed after normalization, leaving
19 TCR rows and 20 MHC columns. The contact output retains physical chain IDs
and exact structure weights. No score default is changed.

The companion peptide-to-MHC matrix uses ``--interface peptide_mhc --keep-cys``
and output ``TCRen-pm.csv``, with the same Canonical2026 source and balancing.
It retains all 20 peptide rows and 20 MHC columns. Its row/column marginals
are derived from peptide-to-MHC contacts, separately from TCRen_m.

Residue saturation surprisal
----------------------------

``saturation-export -s MANIFEST -o residues.parquet -t 0`` exports all CDR3
side chains, their canonical heavy-atom inventory sizes n, distinct contacted
atom counts k, completeness, junction-relative positions and target residues.
Contact is at most 5 Angstroms to peptide, complete MHC, or their union.
Junction-relative position 1 is the conserved Cys before CDR3. Gly has n=0;
missing side-chain atoms remain flagged, not counted as fully contacted.
``--native`` restricts to nonwithdrawn class-I alpha-beta complexes.
Annotation is batched once; coarse descriptor chunks use bounded workers.

``saturation-score --native native.parquet --features query.parquet -o RUN``
excludes exact anchor-excluded CDR3 pairs shared with queries, then estimates
p_full=(F+0.5)/(N+1) separately by amino acid and interface. N counts complete
native side chains with k>0 and F counts k=n. Complete engagement contributes
-log2(p_full) bits. The bounded full score divides the observed bit sum by the
maximum possible sum over those same supported contacting residues. It is an
information-weighted saturation fraction, not joint-event surprisal.

A second score uses one pooled atom probability per interface,
p=(sum(k)+0.5)/(sum(n)+1), over all complete native side chains including zeros.
For each complete query residue with k>0, its contribution is
-log2[Binomial(k;n,p)/(1-(1-p)**n)]. Its maximum over k=1,...,n normalizes the
sum to [0,1]. For n=1 the conditional distribution is deterministic and carries
zero information. Empty information produces zero, with explicit counts.
Unsupported full-contact reference cells abstain and are counted. Raw bits,
normalizers, native counts and residue-level contributions are retained.

The additional ``full_sum`` and ``binomial_sum`` scores map each raw bit sum S
to S/(1+S), using a fixed one-bit scale. These preserve raw-sum rankings and
retain absolute engagement information; the inventory-normalized fractions above
answer a different question about completeness relative to available information.

All quantities describe empirical contact rarity, not chemical satisfaction
or calibrated binding probabilities. Correlated atoms violate the binomial
independence approximation. High surprisal can reflect implausible contacts.
A bounded optional mixture is (ipTM+w1*prepared+w2*surprise)/(1+w1+w2), with
finite nonnegative w1,w2. Coefficients tuned on evaluated labels are in-sample
results. These commands change no screening defaults.

Mock-reference occupancy
------------------------

``saturation-mock --reference mocks.parquet --metadata mocks.tsv --features
queries.parquet --queries baseline.parquet -o NEW_DIR`` fits two experimental
references. Metadata must contain unique ``pdb.id``, ``epitope`` and zero ``y``
(computed mispairings). Benchmark epitopes and either exact anchor-excluded
CDR3 chain are excluded before fitting. Coordinate-duplicate exclusion belongs
to the benchmark's raw-source audit. Outputs retain exclusions and raw counts.

The balanced reference gives equal total weight to each epitope/chain type and
equal weight to each distinct CDR3 within that cell; repeated models of a CDR3
share its weight. Model/chain weights have mean one. Weighted frequencies retain
the half-count smoothing; they are not independent-observation counts.
``balance_mock_residues`` and ``fit_reference`` implement this route.

``tail_sum`` is S/(1+S), where S sums -log2 P(K>=k | K>0) for complete contacting
side chains, K~Binomial(n,q). The same pooled reference q includes zero-contact
residues. This upper tail rewards high occupancy rather than both improbable
extremes. Atom dependence remains an approximation; scores are not calibrated
binding probabilities. ``score_saturation`` exports the individual tail bits.
Arm/weight selection must exclude held-out evaluation labels. Earlier exploration
of those cohorts still prevents claiming untouched external validation.

``tail`` divides the same observed upper-tail bit sum by its maximum attainable
sum over **all** complete CDR3 side chains with more than one atom, including
noncontacting residues. For each such residue the denominator contributes
-log2 P(K=n | K>0). Ala has no conditional information; Gly has no side-chain
opportunity. An empty denominator scores zero. This fixed opportunity denominator
keeps scores in [0,1] without compressing large bit sums toward one. It does not
estimate solvent-accessible opportunities or add a fitted calibration constant.
