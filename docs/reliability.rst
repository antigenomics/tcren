Reliability: scoring one modelled structure
===========================================

.. currentmodule:: tcren.reliability

A generator's confidence is not a binding prediction. On the balanced VDJdb panel, models in
ipTM's **top decile** are still **26.2 %** [18.7, 35.5] non-binders — and that is exactly the band
where the coordinates carry the most information the confidence does not read. This module is the
read-out for one structure at a time.

``S_free``: three blocks, one divide
------------------------------------

.. math::  S_{\mathrm{free}} \;=\; \frac{Q}{\sigma_Q} \;+\; \frac{T}{\sigma_T}
           \;+\; \frac{\Pi - \mu_\Pi}{\sigma_\Pi}

Each block is a fit-free directional score :math:`z(x)^\top C^{-1} s` standardized against the
**Native2026 crystals**: :math:`Q` the interface geometry (:data:`tcren.cohort.Q_FEATURES_GEOM`),
:math:`T` the footprint *shape* free of its size (:data:`T_FEATURES_TOPO`), and :math:`\Pi` the
interface energy read against the partition function (``neg_energy`` from :mod:`tcren.potts`).

**The outer transform is a divide, not a z-score.** A block score's native mean is 0 by
construction, so re-centring does nothing; its *variance* is :math:`s^\top C^{-1} s`, which is not
1. Measured native spreads are 1.43 (:math:`Q`), 1.61 (:math:`T`) and 14.13 (:math:`\Pi`), so
without the division the energy would carry ten times the weight of the geometry. The claim the
formula makes is **equal weight in native-sd units**, and the division is what makes it true.

Why this and not ``P_native``
------------------------------

:func:`tcren.cohort.p_native` refits a latent-class model on every call and **raises when a cohort
has fewer rows than features** — so it is undefined for a single structure, and its value depends on
what else was scored alongside it. ``S_free`` fits nothing at call time. It is still emitted, and
still documented, as a cohort-refit score; it is no longer the recommended one.

Calibration, and what a probability costs
------------------------------------------

:func:`p_binder` maps a score through a **frozen out-of-fold Platt link** — leave-one-epitope-out on
the 22-cohort panel, within-epitope 5-fold on TCRvdb, coefficients the fold means. A probability is
a stronger claim than a rank, so read the expected calibration error beside it: the composed score
reaches ECE 0.020 on the panel where ipTM alone reads 0.065.

Each link's name is the score it expects. Passing a raw ``S_free`` to a ``min rank%(...)`` link is a
category error, not a rescaling; :func:`available_links` lists what is shipped.

The generator diagnostic
-------------------------

:func:`af_band` looks a confidence up in a frozen band table: how often a model *this* confident is
a non-binder, with a Wilson interval, and ``s_free_roc_in_band`` — what ``S_free`` still separates
inside that band. Bands are deciles of the benchmark's own confidence distribution, never scanned
for an effect. Values outside the range clamp to the end bands.

From the command line
----------------------

.. code-block:: console

   $ tcren features -s models/ -i placement,interface,topology,energetics -o feats.tsv
   $ tcren assess --features feats.tsv -o assessed.tsv \
       --link 'tcrvdb|S_nat' --band 'tcrvdb|ipTM'

   618 structures; S_free = Q + T in native-sd units, 618 finite
     p_binder via 'tcrvdb|S_nat'; mean 0.562
     top 50% of the set (309 structures): mean p_binder 0.737 against 0.562 overall
     generator diagnostic (tcrvdb|ipTM): 60 of 618 structures sit in the top confidence decile,
     where 15.3% [8.2%, 26.5%] of benchmark models are NON-binders and S_free still reads 0.773
     ROC-AUC

Add the energy term by joining ``tcren potts score``'s ``neg_energy``; without it ``assess`` emits
the two-block form and says so rather than imputing.

API
---

.. autofunction:: s_free
.. autofunction:: t_score
.. autofunction:: p_binder
.. autofunction:: af_band
.. autofunction:: reliability_reference
.. autofunction:: available_links
.. autofunction:: available_bands
.. autofunction:: moments
.. autofunction:: inversion_flag
.. autofunction:: screening_yield
.. autodata:: T_FEATURES_TOPO
.. autodata:: T_SIGNS
.. autodata:: PI_FROZEN
