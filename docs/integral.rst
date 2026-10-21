Integrated receptor readout
===========================

``integral-reference-1`` is an experimental class-I receptor-ranking candidate.
It combines four visible contributions without fitting assay labels. Original
energies and structural descriptors are unchanged; this is not a new default.

Usage
-----

Generate descriptors from raw structures in one batched CLI call::

   tcren features -s native.txt -i interface,topology,energetics -t 0 -o native.tsv
   tcren integral-reference --features native.tsv -o reference.json
   tcren features -s models.txt -i interface,topology,energetics -t 0 -o models.tsv
   tcren integrate --features models.tsv --reference reference.json \
       --confidence iptm.tsv -o integral.tsv

The reference input must be verified native structures; class-II rows are counted
and excluded. The confidence TSV contains unique ``complex.id`` and ``iptm``
columns. Keys are joined explicitly; extra or duplicate IDs fail. Missing confidence
is retained. If the features already contain ipTM, omit ``--confidence``.
Input descriptors must use the bundled TCRen2 potential and 5-Angstrom contacts.
Retain the producing CLI provenance: a numeric table alone does not prove its
potential, structural origin, cutoff or coordinate identity.

Definition
----------

For each model, let :math:`\Phi` be ``Phi_tcr_pep`` and :math:`N` be
``n_contacts_tp`` from the same all-TCR residue contact map. Define
:math:`x_E=-\Phi/N` for :math:`N>0`. The other structural components are
:math:`x_G=\mathtt{sc\_shape}`, the height-field Pearson shape complementarity,
and :math:`x_T=\mathtt{g\_cyclo\_frac}`, the CDR:pMHC graph's cycle count
divided by its edge count. The graph includes peptide and MHC contacts.

For component :math:`k\in\{E,G,T\}`, let :math:`r_{kj}` be its finite native
class-I reference values and :math:`n_k` their count. The frozen mid-distribution is

.. math::

   R_k(x)=\frac{\sum_{j=1}^{n_k}\mathbf{1}[r_{kj}<x]
               +\tfrac12\sum_{j=1}^{n_k}\mathbf{1}[r_{kj}=x]}{n_k}.

The integrated score, with higher values in the declared favourable direction, is

.. math::

   I=\tfrac14\left(\mathrm{ipTM}+R_E(x_E)+R_G(x_G)+R_T(x_T)\right).

All four contributions lie in [0,1]. Equal weights are a convention and allow
compensation between components. The score is dimensionless, not a binding
probability or a free energy. Native empirical distributions are estimated
dependencies, explicitly recorded with support and hashes. They are never estimated
from the query batch, so adding unrelated queries cannot change a model's score.

Missingness and limitations
---------------------------

Each missing component contributes a fixed 0.5 and has a ``missing_*`` flag.
``n_components`` counts observed contributions; zero observed contributions yields
a null integral score. Zero contacts makes mean energy unavailable. Invalid domains,
infinities, noninteger or negative contact counts, nonzero energy with zero contacts,
and class-II query rows fail explicitly. Evaluate complete and partial scores with
their coverage stated; neutral imputation is not observed evidence.

Neither averaging energy per edge nor normalizing graph cycles guarantees
independence from interface size. Shape and graph terms can share information.
Reference percentile and ipTM scales are bounded but not identically calibrated.
Transfer accuracy and improvement over components require measured, matched-row
evaluation. The implementation does not infer an epitope-specific sign reversal.

API
---

.. automodule:: tcren.integral
   :members: build_reference, integral_scores

Separate-training logistic readout
----------------------------------

``integral-fit`` fits the explicitly separate ``integral-glm-1`` candidate.
It does not change the equal-mean reference or the default score. Its four inputs
are ipTM and the three raw components defined above, before native percentile
conversion. Every input, imputation mean and scale is fitted using training rows
only; query batches cannot change the resource.

For training epitope :math:`e`, class :math:`y` and cell count
:math:`n_{ey}`, each row has weight :math:`w_i=N/(2K n_{ey})`, where
:math:`N` is the training row count and :math:`K` the epitope count. Each epitope
must contain both real and constructed mock examples. This gives equal total
weight to epitope/class cells; it is not an assay-reliability weight.
For each component :math:`k`, the weighted mean :math:`\mu_k` and population
standard deviation :math:`s_k` use its finite training observations. A constant
component uses :math:`s_k=1`. Define :math:`z_{ik}=(x_{ik}-\mu_k)/s_k` when
observed and zero otherwise, with missingness indicator :math:`m_{ik}`.
The frozen logit and readout are

.. math::

   \eta_i=b+\sum_{k=1}^{4}\beta_k z_{ik}+\sum_{k=1}^{4}\gamma_k m_{ik},
   \qquad S_i=\frac{1}{1+\exp(-\eta_i)}.

The intercept :math:`b`, four component coefficients :math:`\beta_k` and four
missingness coefficients :math:`\gamma_k` are fitted with scikit-learn's L2
logistic regression, ``C=1``, ``solver=lbfgs``, ``tol=1e-8`` and
``max_iter=1000``. Convergence warnings fail fitting. These are learned
coefficients, not a fit-free physical energy. :math:`S_i` describes the balanced
real/mock training task; it is not a calibrated experimental binding probability.
No covariance whitening, interactions or epitope-specific terms are fitted.

The labels table requires unique ``complex.id``, ``epitope``, ``receptor_id``
and binary ``y``. Exclusion TSVs contain ``epitope`` and ``receptor_id`` columns
respectively. Both must be nonempty; any training overlap fails. The caller owns
biological identity, current positive-record verification, constructed mock
provenance, coordinate deduplication and validation-split integrity. The commands
cannot infer those from arbitrary caller-supplied identity strings.

.. code-block:: console

   tcren integral-fit --features train.tsv --labels labels.tsv \
       --excluded-epitopes heldout_epitopes.tsv \
       --excluded-receptors heldout_receptors.tsv -o model.json
   tcren integral-predict --features query.tsv --model model.json -o scores.tsv

Feature tables include ipTM when available; a missing ipTM column is treated as
missing, as in the equal-mean implementation. All four components must have some
observed training support, and every training row must have an observed component.
Prediction retains partial rows, reports component and missingness contributions,
and returns null when all components are missing. Class/domain checks and potential
identity checks match the equal-mean route. Fitting refuses existing output files.
The JSON contains coefficients, transformations, counts and input hashes, without
training labels or receptor identifiers.

.. automodule:: tcren.integral_glm
   :members: fit_integral_glm, integral_glm_scores
