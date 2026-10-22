Task-specific ranking scores
============================

Epitope ranking and TCR ranking are different tasks and do not share a score.
:mod:`tcren.ranking` owns the **epitope-ranking** score, where the peptide varies while the
receptor and MHC are held fixed. The primary :doc:`screening` route supports both
tasks with separate frozen references. For receptor ranking it keeps ipTM primary
and reports the structural components alongside the combination.

The score
---------

Two components, each standardized against a frozen reference and averaged with equal weights:

.. math::

   S \;=\; \tfrac{1}{2}\Bigl(
       \frac{-(\Delta\Phi_{TP} - \mu_{E})}{\sigma_{E}}
     + \frac{\mathrm{ipTM} - \mu_{C}}{\sigma_{C}}
   \Bigr)

where :math:`\Delta\Phi_{TP}` is the **structural** poly-alanine reference difference on the
TCR:peptide interface — the epitope positions truncated to alanine in 3D with the contact map
rebuilt, from ``scoring --delta --reference-mode structural`` — and ipTM is the generator's
interface pTM, supplied by the caller. :math:`\mu_E` is the mean of the raw
alanine difference, before changing its sign; :math:`\mu_C` is the raw ipTM
mean. Both positive scales are frozen standard deviations of the same raw columns.

The energy enters negated because lower :math:`\Phi` is the more favourable interface, so higher
:math:`S` is always the more binder-like model. :math:`S` is dimensionless and unbounded: a
ranking score, not a binding free energy and not a calibrated probability.

What the constants are, and are not
-----------------------------------

:math:`\mu` and :math:`\sigma` are the mean and standard deviation of each component over a
held-out training pool of 729 models across 14 off-panel epitopes. **No label was read to
produce them**, so they are a change of units, not fitted weights. They are shipped in
``tcren/data/epitope_rank_reference.json`` and the resource's own SHA-256 travels with every
result.

Nothing is fitted at scoring time, and nothing is normalized against the query cohort: a
model's score depends on its own descriptors and the frozen reference, never on the other rows
supplied beside it. Scoring one model alone and scoring it inside a batch give the same number.

Missing components
------------------

Each component carries a ``missing_*`` flag and the row carries ``n_components``. A component
that is absent is **excluded from the mean**, not imputed to a neutral midpoint. A row with no
observed component scores null rather than receiving a default, because an unmeasured row must
not land on the same scale as a measured one.

Using it
--------

.. code-block:: console

   $ tcren scoring -s models.txt --delta --reference-mode structural -o energies.csv
   $ tcren rank-epitopes --features descriptors.tsv -o scores.tsv

``--features`` needs ``complex.id``, ``dPhi_tcr_pep`` and ``iptm``.

.. autofunction:: tcren.ranking.epitope_rank
.. autofunction:: tcren.ranking.bundled_reference
