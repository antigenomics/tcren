Reliability: scoring one modelled structure
===========================================

.. currentmodule:: tcren.reliability

A generator's confidence is not a binding prediction. On the balanced VDJdb panel, models in
ipTM's **top decile** are still **26.2 %** [18.7, 35.5] non-binders — and that is exactly the band
where the coordinates carry the most information the confidence does not read. This module is the
read-out for one structure at a time.

``S``: three blocks, one divide
------------------------------------

.. math::  S \;=\; \frac{Q}{\sigma_Q} \;+\; \frac{T}{\sigma_T}
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

Why this and not a cohort-refit posterior
------------------------------------------

The latent-class posterior this module replaced refitted on every call and **raised when a cohort
had fewer rows than features** — so it was undefined for a single structure, and its value depended
on what else was scored alongside it. Neither property survives contact with a user holding one
model. It was discarded in 2.26.0. ``S`` fits nothing at call time.

Calibration, and what a probability costs
------------------------------------------

:func:`p_binder` maps a score through a **frozen out-of-fold Platt link** — leave-one-epitope-out on
the 22-cohort panel, within-epitope 5-fold on TCRvdb, coefficients the fold means. A probability is
a stronger claim than a rank, so read the expected calibration error beside it: the composed score
reaches ECE 0.020 on the panel where ipTM alone reads 0.065.

Each link's name is the score it expects. Passing a raw ``S`` to a ``min rank%(...)`` link is a
category error, not a rescaling; :func:`available_links` lists what is shipped.

The generator diagnostic
-------------------------

:func:`af_band` looks a confidence up in a frozen band table: how often a model *this* confident is
a non-binder, with a Wilson interval, and ``s_roc_in_band`` — what ``S`` still separates
inside that band. Bands are deciles of the benchmark's own confidence distribution, never scanned
for an effect. Values outside the range clamp to the end bands.

From the command line
----------------------

.. code-block:: console

   $ tcren features -s models/ -i placement,interface,topology,energetics -o feats.tsv
   $ tcren assess --features feats.tsv -o assessed.tsv \
       --link 'tcrvdb|S' --band 'tcrvdb|ipTM'

   618 structures; S = Q + T in native-sd units, 618 finite
     p_binder via 'tcrvdb|S'; mean 0.562
     top 50% of the set (309 structures): mean p_binder 0.737 against 0.562 overall
     generator diagnostic (tcrvdb|ipTM): 60 of 618 structures sit in the top confidence decile,
     where 15.3% [8.2%, 26.5%] of benchmark models are NON-binders and S still reads 0.773
     ROC-AUC

Add the energy term by joining ``tcren potts score``'s ``neg_energy``; without it ``assess`` emits
the two-block form and says so rather than imputing.

Correcting the generator's confidence
-------------------------------------

``af_band`` says how often a confidence band is wrong. It does not say what to believe instead.
:func:`correct_confidence` does, by reading the confidence together with the coordinates:

.. math:: \mathrm{logit}\,P(\mathrm{binder}) = b_0 + b_c\,z(c)
          + b_S\,S + b_N\,N

with :math:`c` the generator's confidence, :math:`S` the single-structure binder
score and :math:`N` the observed contact count, both in native-sd units. It returns the corrected
probability **and its parts**, so a caller can see whether a number moved because of the generator
or because of the structure:

.. code-block:: text

   $ tcren diagnose --features feats.tsv --confidence iptm -o diagnosed.tsv
   618 structures corrected against 'tcrvdb|ipTM'
     the structure argues AGAINST 267 of 590 (45%); mean shift -0.048 nats, range [-4.31, +2.07]
     the five the generator is most confident about:
       67c6026f...  iptm 0.931  p_conf 0.792 -> p_corrected 0.935  (+1.33 nats)
       d3bcd432...  iptm 0.921  p_conf 0.777 -> p_corrected 0.765  (-0.07 nats)
       52f8a2cb...  iptm 0.919  p_conf 0.772 -> p_corrected 0.681  (-0.46 nats)

Two properties to state when reporting it.

**It is not fit-free.** Every other score in this module takes no label anywhere;
:func:`correct_confidence` learns four coefficients from labels and freezes them, exactly as
:func:`p_binder`'s Platt links do. The structural terms it reads are themselves fit-free, but the
weighting is not.

**It is validated where the epitope has structural precedent.** Leave-one-epitope-out on the
balanced VDJdb panel, the correction adds **+0.051** macro ROC-AUC to ipTM and **+0.068** to pLDDT
over the 6 cohorts whose epitope has a solved complex (n = 284), and *subtracts* about 0.04 over
the 16 that do not (n = 743). That is the template covariate everything in this framework divides
under: where no receptor has been co-crystallized with the peptide, nothing works, the generator's
own confidence included. Coefficients are rounded to one decimal, which costs under 0.003 macro
ROC-AUC.

API
---

.. autofunction:: s_score
.. autofunction:: t_score
.. autofunction:: p_binder
.. autofunction:: af_band
.. autofunction:: correct_confidence
.. autofunction:: reliability_reference
.. autofunction:: available_links
.. autofunction:: available_bands
.. autofunction:: available_corrections
.. autofunction:: moments
.. autofunction:: inversion_flag
.. autofunction:: screening_yield
.. autodata:: T_FEATURES_TOPO
.. autodata:: T_SIGNS
.. autodata:: PI_FROZEN
.. autodata:: CORRECTION_VALIDATED_ON
