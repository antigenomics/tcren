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
   $ tcren assess --features feats.tsv -o assessed.tsv --band 'tcrvdb|ipTM'

   618 structures; S = Q + T in native-sd units, 618 finite
     top 50% of the set (309 structures): mean S 0.737 against 0.562 overall
     generator diagnostic (tcrvdb|ipTM): 60 of 618 structures sit in the top confidence decile,
     where 15.3% [8.2%, 26.5%] of benchmark models are NON-binders and S still reads 0.773
     ROC-AUC

Add the energy term by joining ``tcren potts score``'s ``neg_energy``; without it ``assess`` emits
the two-block form and says so rather than imputing.

API
---

.. autofunction:: s_score
.. autofunction:: t_score
.. autofunction:: af_band
.. autofunction:: reliability_reference
.. autofunction:: available_bands
.. autofunction:: moments
.. autofunction:: inversion_flag
.. autofunction:: screening_yield
.. autodata:: T_FEATURES_TOPO
.. autodata:: T_SIGNS
.. autodata:: PI_FROZEN

Nothing here is fitted against a binding label
----------------------------------------------

Version 2.28.0 removed the last read-outs that were: the frozen Platt links behind ``p_binder`` and
the four-coefficient confidence correction behind ``tcren diagnose``, both of which were fitted
out of fold on the benchmarks. Every quantity this module now returns is a directional score or an
empirical band table, so a value computed today depends on the structure it was computed from and
on the 374 Native2026 crystals, and on nothing else. The removed read-outs and the numbers they
produced are recorded in the manuscript repository's ``LEGACY.md``.
