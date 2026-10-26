Three-proportion recognition score
==================================

The experimental final model is one line:

.. math::

   S=w_I I+w_B B+w_M M,\qquad w_I,w_B,w_M\geq0,\quad w_I+w_B+w_M=1.

I is generator ipTM, B is the existing whole-complex recognition graph beta1/E,
and M is bounded contact-motif enrichment. The graph is CDR:peptide +
peptide:MHC + CDR:MHC, without internal peptide edges. No prepared combination
is used as a predictor: ipTM and Betti enter separately, so confidence is not
counted twice. All components and their proportions are explicit. No polynomial
power, intercept, sequence predictor or extra fitted classifier enters.

Fisher motif evidence
---------------------

Use every eligible negative control outside the benchmark, with no cap,
downsampling or receptor collapse. The local residue graph and reference-tail
calculation are defined in :doc:`residue-surprise`. The default local statistic
is added molecular cycles, m-q, retaining repeated physical atom contacts.

For local mock tail probabilities p_r, Fisher's statistic is
``F=-2*sum(log(p_r))`` over all six CDR loops and the epitope. F is unbounded;
the independent-test combined p-value ``chi2.sf(F,2*k)`` is bounded, where k
counts supported local observations. This analytic p-value is exported only as
a diagnostic: contact subgraphs overlap and the empirical tails are discrete.
The independence/continuity assumptions of the usual Fisher test do not hold.
See the `SciPy combine_pvalues documentation <https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.combine_pvalues.html>`_.

The scoring channel instead uses the whole-control-graph empirical survival:
``p_F=(1+sum_j v_j*I(F_j>=F))/(1+sum_j v_j)`` and **M=1-p_F**. Here v_j is the
explicit positive model weight, normally one, and F_j is a whole control graph's
Fisher statistic. Control local ranks leave that entire model out through
vectorized subtraction of its sufficient counts. This preserves within-graph
dependence rather than multiplying independent evidence. Every eligible control
is retained; a control with no supported local tests has F=0. Unsupported query
motifs abstain, with full counts and support flags exported. This is estimated
reference enrichment, not a calibrated binding probability or a formal
exchangeability guarantee for heterogeneous model sources.

Weight fitting and assessment
-----------------------------

Fit only the three proportions, i.e. two free coefficients, by equal-cohort mean
average precision across the five labelled TCRvdb epitopes. A simplex grid at
one-percentage-point coefficient precision implements constrained fitting; it
is not a search over graph definitions, powers or extra model features. Exact
AP tie blocks are retained. Ties prefer fewer active channels, then greater
confidence weight. No classification threshold is fitted.

By default, fit one global triplet on all labelled TCRvdb rows, without folds.
Missing-component means use these training rows; missingness indicators remain
in output. Accuracy on TCRvdb is explicitly full-fit, in-sample performance.
VDJdb uses the frozen triplet and retains every cohort and template stratum.
All gains/losses and observation counts remain visible.

``motif-mixture-fit-cohorts`` instead fits one triplet independently for TCRvdb,
template-present VDJdb and template-absent VDJdb. Each fit maximizes equal-epitope
macro AP within its own declared cohort. It exports a single resource containing
the three triplets and cohort-specific missing-value means. All accuracy in
these cohorts is then in-sample. ``motif-mixture-score`` dispatches by explicit
dataset/template metadata and rejects unassigned or ambiguous cohort membership.
There is no learned gate or confidence-dependent weight function.

Optional ``--folds`` and ``--residues`` preserve the earlier diagnostic. Recorded
exact-chain-component folds are inputs, not resplit. Each fold uses training
rows only. ``motif-mixture-diagnose`` exports within-fold AP/AUROC and the score
shift caused by applying each weight vector to the same feature rows. Fold-mean
ranking metrics and pooled OOF metrics remain separate: pooling scores from
different fitted rankers can introduce cross-fold scale artefacts. This issue is
described for AUC in `Airola et al., 2009 <https://proceedings.mlr.press/v8/airola10a.html>`_.
No calibration transform or additional fitted coefficients enter the final score.

Existing batched annotation/search APIs own structural exports. ``recognition-betti``
reuses ``topology.graph.recognition_graph_features`` and exports raw
``complex_cyclo_frac=beta1/E`` separately from ``receptor_contacts_available``.
The raw fraction remains defined when only peptide:MHC contacts remain; the
production recognition fraction retains its existing availability semantics.
Existing vectorized Polars/NumPy kernels suffice. The global fit partitions the
coefficient simplex into coarse chunks in one bounded process layer. Optional
fold diagnostics use six coarse jobs. Native kernels stay serial inside workers.
No worker pool surrounds the sequence-search engine.

CLI
---

::

   tcren export -s structures.txt --quantity recognition-betti -t 0 -o topology.tsv
   tcren motif-fisher --features query-residues.parquet --reference mock.json --controls control-residues.parquet -o motifs
   tcren motif-mixture-fit --baseline labels.parquet --topology topology.tsv --motifs motifs/scores.tsv -t 0 -o fit
   tcren motif-mixture-score --features fit/predictions.parquet --resource fit/weights.json -o scores.tsv
   tcren motif-mixture-fit-cohorts --features fit/predictions.parquet -t 0 -o cohorts
   tcren motif-mixture-score --features cohorts/predictions.parquet --resource cohorts/weights.json -o cohort-scores.tsv

Fit outputs component/score tables, per-cohort AP/AUROC, full and fold-specific
proportions, imputation means, objective values and input hashes. A frozen
weights resource scores another table containing iptm,betti,enrichment without
refitting. Released screening defaults remain unchanged.

.. automodule:: tcren.contacts.motif_mixture
   :members:
