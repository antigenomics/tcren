Screening a cohort
==================

This is the entry point. Everything else in ``tcren`` is a toolkit — annotation and markup,
canonical orientation, contact listing, surface topography, descriptor extraction — and you reach
for it when you want a particular number. :command:`tcren screen` answers the question most users
arrive with: *rank these models for me*.

.. code-block:: console

   $ tcren screen -s models.tar.gz -o scores.tsv
   tcr cohort: 420 models, 1 peptide(s), 420 receptor(s); 420 scored

It is deliberately narrow about when it will answer.

It requires a cohort
--------------------

A **cohort** holds one side fixed across every structure. ``screen`` works out which side, and
refuses a set that is neither rather than scoring it as though it were one:

.. code-block:: text

   Do the structures form a cohort?
   │
   ├── same receptor, peptides vary ──► EPITOPE RANKING
   │       energy = -dPhi, the epitope substituted to alanine in 3D
   │       topology excluded
   │
   ├── same pMHC, receptors vary ─────► TCR RANKING
   │       energy = -dFT, the contacted CDRs substituted to alanine in 3D
   │       topology carried
   │
   └── neither ───────────────────────► refused, naming what differs
   
   Then: at least 30 models with a TCR:peptide interface?
   │
   ├── no  ──► ipTM and every component reported. NO combined score.
   └── yes ──► the above, plus the combined score.

The two are different tasks and do not share a score. **The alanine reference is taken on the side
that varies** — substituting the fixed side would reference away a constant. Topology is carried
only for TCR ranking because its direction is established there and inverted on peptide ranking.

The score
---------

Components are standardized against a frozen reference and combined by **weighted Stouffer**:

.. math::

   S \;=\; \frac{\sum_i w_i z_i}{\sqrt{\sum_i w_i^2}},
   \qquad
   z_i \;=\; d_i \,\frac{x_i - \mu_i}{\sigma_i}

where :math:`x_i` is component :math:`i`'s raw value, :math:`\mu_i` and :math:`\sigma_i` are its
frozen mean and scale, :math:`d_i \in \{+1, -1\}` is its fixed direction, and :math:`w_i` is its
weight. Higher :math:`S` is always more binder-like. :math:`S` is dimensionless and unbounded: a
ranking score, not a binding free energy and not a calibrated probability.

The confidence block carries the reference's ``default_weight`` and the structural components
share the remainder equally. That weight is **task-specific** — 0.8 for TCR ranking, where the
generator's ipTM is much the strongest single input, and 0.5 for epitope ranking, where the
energy is stronger.

Symmetric pooling does not work here, and was measured rather than assumed. Fisher's
:math:`-2\sum \ln p`, equal-weight Stouffer, logit-sum, the geometric mean and the minimum
percentile all score **below ipTM alone** on every receptor-ranking arm, because each treats
confidence and structure as exchangeable evidence when they are not.

Models the generator never placed on the groove
-----------------------------------------------

A model with **zero TCR:peptide contacts** is not a TCR:pMHC complex. The receptor was docked
somewhere else — on β2-microglobulin, or off the groove entirely — and its energy, shape and
topology are all computed over an empty interface. ``screen`` **discards** such models: they are
kept in the output, with ``no_interface`` true and the reason recorded, but they are not scored
and **they do not count toward the cohort floor**. A cohort that only reaches 30 by including
failed models has not reached it.

This matters more than it sounds. In the balanced VDJdb cohort RAQAPPPSW, 11 of 20
binder-labelled models have no interface at all, against 0 of 24 mocks — so over half that
cohort's positives are structures where the receptor never touched the peptide. That cohort's
ipTM reads 0.0479, which looks like a catastrophic inversion and is nothing of the kind: the
generator is correctly reporting that the binder-labelled models are bad, and the label disagrees
with the coordinates. Discarding is the honest response; reversing the score there would rank
"the generator failed to build this complex" as "this is a binder".

Flags: when not to trust a ranking
----------------------------------

A ranking can fail in a way an aggregate hides. On balanced VDJdb the cohort RAQAPPPSW has ipTM
**0.0479** — not a failure to rank but a near-perfect ranking *pointing backwards*. A score that
weights ipTM at 0.8 reproduces that inversion faithfully.

``screen`` cannot correct this. A cohort-level reweighting was measured and moved the aggregate by
at most +0.013 while making two of the four worst cohorts worse. What it *can* do is see it
coming, because the regime is detectable with no labels at all, and say so.

Every flag reads only the cohort's own descriptors.

``small_cohort``
    Fewer models than ``--min-cohort``. **Withholds the score** on both tasks.

``low_confidence``
    The cohort's mean ipTM is below the reference cut, 0.7237 — the 20th percentile of
    per-epitope mean ipTM over the 14 off-panel training epitopes. Cohort mean ipTM correlates
    **+0.7099** with ipTM's own accuracy across 29 evaluation cohorts, so this is the regime where
    the generator's ranking is least reliable. Advisory by default on both tasks.

``forced_pose``
    The cohort's mean interface strain, :func:`tcren.cohort.strain_z`, is above the reference cut
    of 0.5586. ``tcren`` has long recorded that the contact energy "reads real binding chemistry
    on well-modelled poses and *inverts* on forced ones". Advisory. It **abstains** — reporting
    neither true nor false — unless the strain was standardized against the same population the
    cut came from, because ``strain_z`` standardizes within whatever set it is given and an
    uncalibrated mean is not comparable to the cut.

``confidence_disagrees``
    ipTM and the structural block rank this cohort oppositely, so at most one of them is right
    about it. Advisory; on the evaluation cohorts it fires rarely and is not a reliable predictor
    of failure on its own.

Default flag behavior
^^^^^^^^^^^^^^^^^^^^^

Flags report by default on both tasks. The cohort floor and no-interface discard
still apply independently. Use ``--on-flag withhold`` to request withholding for a
low-confidence cohort; the default ``--on-flag report`` leaves its score unchanged.

Responding to a flag
^^^^^^^^^^^^^^^^^^^^

``--on-flag`` chooses what happens when a cohort is flagged ``low_confidence``. Measured over the
29 TCR-ranking evaluation cohorts before the no-interface discard, where the flag fires on 3:

====================  ==========  =========  =================  ============================
``--on-flag``         pooled ROC  pooled PR  cohorts below 0.5  worst cost on a good cohort
====================  ==========  =========  =================  ============================
``report`` (default)      0.7250     0.7049                  6  none, nothing changes
``structure``             0.7296     0.7063                  6  0.9735 → 0.9225
``geometry``              0.7273     0.7033                  6  0.9735 → 0.7316
``withhold``                   —          —                  —  a usable score is lost
``reverse``               0.7410     0.7162                  3  0.9735 → **0.0265**
====================  ==========  =========  =================  ============================

``reverse`` is the strongest rescue and the most dangerous. It takes the worst cohort,
RAQAPPPSW, from 0.0479 to 0.9521 and halves the number of cohorts ranking below chance — but if
the flag ever fires on a good cohort it inverts that cohort completely.

``structure`` is the middle ground: it drops ipTM and scores on the structural components alone,
which removes the confidence contribution but does not guarantee a correct ranking. In the
pre-discard evaluation it lifts RAQAPPPSW to 0.4667 — *to* chance rather than above it — and its worst misfire costs 0.05 rather than 0.95.

``report`` is the default. The only evidence that acting helps comes from the same panels this
score is reported against. Reversal is opt-in and is not recommended: RAQAPPPSW's apparent
inversion includes 11 of 20 binder-labelled models with no TCR:peptide contacts, against
0 of 24 mocks. Here confidence correctly identifies failed poses; reversing it would reward
modelling failure. These models are now discarded. The remaining cohort effects are open
research questions, not an established biological inversion mechanism.

A response rewrites only the combined score. Every component, its percentile and its missingness
flag stay as they were, so the evidence the response was applied to remains visible, and
``flag_response`` records which one ran.

When a score is withheld the components are still reported, and ``score_unavailable`` names
the reason. ``--on-flag report`` reports flags without changing the score. A flagged cohort
requires inspection of its models and provenance; a flag alone does not validate a corrective
response.

Why 30 models
-------------

Below a cohort of 30 the combined score is withheld and the components are reported alone. The
per-cohort noise in the margin over ipTM follows

.. math::

   |\Delta| \;\approx\; \frac{0.0765}{\sqrt{n}} + 0.0019

measured over 29 evaluation cohorts. At :math:`n = 20` that noise is 0.019 and at :math:`n = 30`
it is 0.016, against the largest effect the combination is known to produce, 0.035. Under about 30
models a score would be decoration: you could not tell it from ipTM even if it were better.

Raise or lower the floor with ``--min-cohort`` if you know what you are doing.

Where the combination helps, and where it does not
--------------------------------------------------

On balanced VDJdb, split by whether a native template exists for the complex:

================================  =======  =========  ========  ==========
stratum                            models      score      ipTM      margin
================================  =======  =========  ========  ==========
native template present               297     0.7299    0.6949     +0.0350
no native template                    792     0.4985    0.5513     −0.0528
================================  =======  =========  ========  ==========

Where a native template exists the coordinates are accurate enough for structural descriptors to
mean something, and the combination gains. Where none exists every component sits between 0.455
and 0.540 — chance — and the combination costs. pLDDT does **not** substitute as that gate, and
neither does surface topography.

So: trust the combined score where your complex has a template, and treat ipTM as the safer
primary where it does not.

What this does not claim
------------------------

- **A high score on a peptide-ranking panel is not by itself evidence of TCR specificity.** On
  CPL a combination built only from the pMHC face, with no TCR atom in it at all, reaches 0.917
  within clone — above ipTM alone. Much of that contrast is a property of the peptides rather
  than of recognition.
- **Neither evaluation panel tests epitope specificity directly.** Only 21 of 1,068 receptors in
  balanced VDJdb, and 4 of 1,648 TRB CDR3s in the five-epitope screen, appear against more than
  one epitope, so the panels almost never ask whether one receptor ranks differently for two
  epitopes.
- **The pooled margin over ipTM is small and within these benchmarks' resolution.** The ipTM
  baseline itself moves by 0.0014–0.0028 under an innocuous change of which rows are scored. The
  honest claim is that the combination *matches* ipTM overall and *adds* on the stratum above.

Everything else is a toolkit
----------------------------

:doc:`features` for the full descriptor table, :command:`tcren annotate` and
:command:`tcren contacts` for markup and contact listing, :command:`tcren superimpose` for
canonical orientation, :command:`tcren surface` for pMHC height and hydropathy maps,
:doc:`ranking` for the single-task epitope score, :doc:`kit` for the fit-free predecessor tier.

.. autofunction:: tcren.screen.screen_scores
.. autofunction:: tcren.screen.cohort_signature
.. autofunction:: tcren.screen.stouffer
.. autofunction:: tcren.screen.cohort_flags
.. autofunction:: tcren.screen.bundled_reference
