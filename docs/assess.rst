Assessing a modelled complex
============================

.. currentmodule:: tcren.score

A co-folding model builds a confident TCR:pMHC complex for any receptor--peptide pair, binding or
not. :mod:`tcren.score` reads the coordinates it produced and answers four separate questions about
them, and **every answer is defined for a single structure**: the transform, the class means and
the covariance are all frozen on a hold-out that ships with the package, so nothing is estimated
from the rows you are scoring and a score does not change depending on what was scored beside it.

Two commands
------------

.. code-block:: console

   $ tcren features -s models/ -i placement,interface,topology,energetics -o feats.tsv
   $ tcren assess --features feats.tsv -o scores.tsv

The first is the expensive pass and runs once. The second is arithmetic over that table.

Pass ``--peptide`` when the **peptide** is what varies across the structures being compared, as in a
combinatorial library or a mutational scan. Otherwise the five descriptors computed without the
receptor are marginalized out, because they are constant across every structure of one epitope on
one allele and a model reading them reaches the cohort's name without reading an interface.

What each column is for
-----------------------

.. list-table::
   :header-rows: 1
   :widths: 22 10 68

   * - column
     - tier
     - what it answers
   * - ``pose_score``
     - 1
     - Is this the kind of interface real complexes make? A one-class distance to the manifold
       hold-out binders occupy. **No negative and no binder label enter it** — the same standing as
       :func:`tcren.cohort.q_score`. This is the bad-pose channel.
   * - ``binder_score``
     - 2
     - Log-odds that the complex is a genuine recognition interface.
   * - ``channel_*``
     - 2
     - The same log-odds, marginalized to one descriptor family, so a number can be attributed.
   * - ``peptide_score``
     - 0
     - The poly-alanine-referenced recognition energy. Nothing is fitted in it. This ranks
       **peptides against a fixed receptor**, and it reads below chance on a receptor benchmark —
       a property of the reference frame, not a defect.
   * - ``confidence_residual``
     - 1
     - The reported ipTM minus what the coordinates say it should have been. A large positive
       residual is a model the generator is more certain of than its own geometry warrants.
   * - ``binder_iptm``
     - 2
     - ``binder_score`` plus ``logit(ipTM)``. Two log-odds added, no coefficient to fit, still
       defined for one structure. The recommended read when a confidence is available.

Higher is better throughout.

The five channels
-----------------

A marginal of a Gaussian is a sub-block of its covariance — exact, closed form, no re-fit — so
reading one family of descriptors costs an index and nothing else:

``placement``
   where the receptor sits in the groove frame.
``interface``
   how much interface it makes, and of what chemistry.
``shape``
   the footprint's shape, free of its size.
``energetics``
   the contact chemistry, in kT.
``mechanics``
   the interface read as a network of breakable springs.

They do not sum to ``binder_score`` and should not: the whole model also reads the correlations
*between* channels, which a per-channel view cannot show. What they give you is attribution — and,
sometimes, a better instrument. A channel beats the whole model where the whole model dilutes it:
on template-free cohorts ``channel_shape`` reads 0.637 against the full posterior's 0.615, and on a
combinatorial peptide library ``channel_energetics`` reads 0.700 against 0.542.

Which residue carries it
------------------------

The five channels say which *part of the structure* a score comes from. ``residue_deltas`` says
which **residue**::

    tcren explain -s model.pdb --score binder -o deltas.tsv

Each row is one interface residue, and ``delta`` is

.. math:: \Delta_i = L(x) - L(x_{\setminus i})

the score of the complex minus the score of the same complex with that residue's atoms removed.
Positive means the residue carries the score. It is defined for every read-out and every channel,
so one table colours a figure by any of them — pass it straight to
:func:`tcren.viz.pymol.importance_scene`.

``peptide_score`` is the exception and does not need it: that score is a sum over contacts, so
:func:`tcren.energetics.scoring.position_profile` already returns an **exact** per-position
decomposition that sums back to the score. Use the exact one where it exists. Leave-one-out is for
the four read-outs that have no exact split, because each is a function of ~149 whole-structure
scalars through a full covariance and no residue owns a share of one.

Cost is one descriptor pass per interface residue — a few minutes for a complex. Chain typing and
the MHC call are mmseqs searches and run **once** for the whole structure, never per mask. A masked
row whose descriptors go undefined is reported as a null delta rather than imputed.

Reproducing the frozen model
----------------------------

The coefficients are frozen, and the inputs they were frozen against are named:

.. code-block:: console

   $ tcren fetch-data                                   # the structure sets the manifest names
   $ tcren features -s <those structures> -o hold.tsv
   $ tcren fit-holdout --features hold.tsv -o refit.npz

``refit.npz`` matches the shipped model bit for bit. :func:`holdout_manifest` returns the 8,292
structures with their dataset, epitope, label and ipTM. Refit on your own hold-out by passing
``--manifest``, then read it back with ``tcren assess --model``.

The predecessor tier
--------------------

:func:`tcren.cohort.q_score`, :func:`tcren.reliability.t_score` and
:func:`tcren.reliability.s_score` are the fit-free directional scores this layer generalises, and
they are still supported. ``S`` leads the functionally validated receptor screen on its own, and it
*composes* with ``binder_score`` rather than being replaced by it.

API
---

.. autofunction:: pose_score
.. autofunction:: binder_score
.. autofunction:: channel_scores
.. autofunction:: confidence_residual
.. autofunction:: peptide_score
.. autofunction:: residue_deltas
.. autofunction:: score_table
.. autofunction:: holdout_model
.. autofunction:: holdout_manifest
.. autoclass:: ScoreModel
   :members:
.. autofunction:: tcren.reliability.artefact_directions
