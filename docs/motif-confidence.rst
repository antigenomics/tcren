Confidence gates and local contact motifs
=========================================

These development commands extend :doc:`motif-mixture` without changing screening
defaults. All eligible source-labelled mocks remain, without caps or receptor
collapse. Accuracy is a full-fit development readout.

The smallest confidence gate is
``S=(1-q*(b+m))*I+q*b*B+q*m*M``, with ``b,m>=0`` and ``b+m<=1``.
I is generator ipTM, B is raw complex beta1/E and M is empirical motif enrichment.
All three weights remain nonnegative and sum to one. q is generator ipTM or
mean residue pLDDT/100, either per complex or averaged over the supplied unlabelled
(dataset,stratum,epitope) cohort. No peptide-specific coefficient is fitted.
Cohort means therefore depend on the supplied batch composition.

``motif-confidence-inputs`` reads supplied generator metadata and recovers missing
pLDDT from hash-verified model C-alpha values through ``mean_bfactor(atom="CA")``.
This interpretation requires generated structures: experimental B-factors are
not pLDDT. Original metadata values remain unchanged. Source flags, hashes and
recovery counts are exported.

``motif-confidence-fit`` fits only b and m, using equal-epitope macro AP on labelled
TCRvdb models. NLV, TTD and GIL AP floors come from the supplied previous-score
table, and the complete coefficient table marks feasible candidates. If none
passes, the best macro-AP candidate is exported with ``preserves_priority=false``;
it must not be described as satisfying preservation. Missing score components
use explicit TCRvdb-training means. VDJdb labels do not enter this fit.

Motif evidence
--------------

The original residue tail conditions on amino-acid type. Optional
``motif-fisher --joint-residue`` instead counts the joint event of that residue
type and at least the observed contact multiplicity within its chain/region.
The numerator remains the type-specific contact tail; the denominator covers
all complete residues in the region. The inclusive plus-one estimator handles
unseen types; a zero-contact statistic still gives p=1. Whole-model leave-out
control ranks use the same denominator, and every control graph remains.

``motif-pair-counts`` creates the smallest sequence-context subgraph: two adjacent
observed residues in one loop/peptide, connected by their peptide bond, as a
single focal patch. With m atom contacts to q distinct partner residues, added
beta1 is m-q. Shared partner residues are counted once; distinct physical atom
contacts remain distinct. Zero-contact patches and incomplete-atom flags remain.
Consecutive observed sequence indices within a loop are required; this assumes
the generated backbone is connected. Intrinsic amino-acid cycles cancel in the
increment. Pair histograms use the same reference machinery and full control
identity/exclusion audit as single-residue histograms.

``motif-residue-audit`` exports residue contributions and descriptive counts for
TTD arginines, GIL CDR3-beta RS/flanking residues, GLC V-encoded CDR1/2 loops and
peptide GIGI segments. These named motifs are interpretation tags, never scoring
bonuses or extra fitted predictors. CDR3 germline contributions are not inferred
from a fixed prefix length. Native crystal examples use the same mock reference
as a structural diagnostic, without interpreting crystal B-factors as confidence.

The implementation reuses Polars/NumPy kernels, existing empirical histograms,
``ap_columns`` and the shared bounded runtime. Coefficient chunks or raw-input
reading chunks form one coarse process layer; native kernels stay serial in
workers. No new dependencies or per-model annotation worker pool are introduced.

Two further simple options
--------------------------

``motif-confidence-fit --odds`` uses a standard log-odds offset gate,
``g=alpha*q/(alpha*q+(1-alpha)*(1-q))``, and
``S=(1-g)*I+g*((1-r)*B+r*M)``. Only alpha and r are fitted; the logit slope is one.
This is the algebra of intercept-only logistic recalibration, used here only as a
bounded weighting rule, not as calibrated binding probability. Its form follows
standard `logistic recalibration <https://onlinelibrary.wiley.com/doi/full/10.1002/sim.9921>`_.

``motif-fisher --positive-reference BOUND.json`` uses signed evidence
``F=2*sum(log(p_bound_tail/p_mock_tail))`` instead of negative-tail surprise.
The empirical whole-mock survival still maps F to bounded enrichment. These are
ratios of estimated tail frequencies, not exact graph likelihoods; no chi-square
p-value is assigned. A positive reference must declare source label 1 and pass
the same exact-pair/hash exclusion audit. Any epitope balancing retains every
eligible model and records weights/counts. The positive pool must exclude
benchmark pairings and conflicting negative-control pairings before assessment.
The reference cannot be derived from the benchmark labels being scored.

The motif audit uses biological examples for interpretation only: the GIL
`RS/notch interaction <https://www.nature.com/articles/ni942>`_, the GLC
`AS01 germline contribution <https://pubmed.ncbi.nlm.nih.gov/21124993/>`_, and
`shared Melan-A/BST2 contact motifs <https://doi.org/10.1016/j.cell.2023.06.020>`_.
All-residue and adjacent-residue score definitions stay epitope-independent.

Complementary confidence weights
--------------------------------

``--complementary`` removes the fixed Betti/enrichment ratio:
``wB=b*q``, ``wM=m*(1-q)``, ``wI=1-wB-wM``. With b,m,q in [0,1], the weights
are nonnegative and sum to one. This remains two global fitted coefficients.
``--rank-confidence`` optionally replaces observed confidence by its ordinary
inclusive empirical CDF in the unlabelled TCRvdb reference feature distribution.
For cohort-mean gates that reference contains one mean per TCRvdb epitope cohort,
including ambiguous-labelled feature rows. The reference values are measured
confidence statistics, not fitted peptide-specific coefficients. New cohorts use
the same frozen numeric CDF. Cohort means depend on batch composition.

``motif-fisher --condition-epitope`` retains the shared all-mock residue histograms
and every whole control graph, but calibrates a query against whole mocks with
the same peptide sequence. If no same-peptide mock exists, the full control pool
is used explicitly. Per-query null-model counts and the specific/fallback flag
are exported; sparse strata imply coarse empirical probabilities. No stratum is
capped or discarded. This option changes the background comparison, not the
number of fitted score coefficients.

Frozen score and interpretation
-------------------------------

Fit the complementary cohort-confidence route with::

   tcren motif-confidence-fit --features components.parquet --metadata confidence.tsv --gate cohort-plddt --complementary --rank-confidence -t 0 -o fit
   tcren motif-confidence-score --features fit/predictions.parquet --resource fit/weights.json -o scores.tsv

The score command needs no binding labels and fits nothing. Supply the full
unlabelled cohort whose mean defines the confidence context; an isolated complex
is a different context. Frozen CDF comparisons use twelve decimal places only
to stabilize numerical ties at reference means, not as a scientific cutoff.
``motif-confidence-report`` produces complete flat cohort/stratum comparisons
from a list of fit directories.

The residue audit optionally reads original ``--annotations`` JSON/JSON.gz and
marks the V-aligned part of CDR3 using recorded AIRR V-start/end boundaries.
Structure residue indices are zero-based; AIRR sequence boundaries are one-based
and inclusive. Unknown boundaries stay unknown. Definite V-encoded CDR1/2 loops
remain marked, without guessing a fixed CDR3 germline prefix.

Optional individual-complex context
-----------------------------------

``motif-confidence-proxy --features components.tsv --fields 2 -o scores.tsv``
applies the frozen allele confidence lookup to individual complexes or small
sets. Required columns are ``mhc``, ``iptm``, ``betti`` and ``enrichment``. One
field means the HLA allele group (for example A*02), two means A*02:01. Unknown
alleles use an explicitly flagged global fallback; support counts are exported.
The lookup predicts cohort-mean confidence, not binding labels. The packaged
lookup uses the complete unlabelled reference; leave-epitope-out evaluation is
reported separately. This command is opt-in; ordinary screening retains measured
cohort-mean pLDDT and its existing minimum cohort size.

In the recorded held-out exploration, allele context gives TCRvdb macro AP
0.759925 versus 0.768183 for the accepted cohort default, over five epitopes and
1,782 classified models. Both allele resolutions agree on this reference because
each represented first field has one second-field subtype. Individual-model
confidence gates performed worse. This is a practical substitute with documented
losses, not an equivalent reconstruction of an unavailable cohort.
