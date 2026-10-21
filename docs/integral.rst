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
