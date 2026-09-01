A kit for AI-generated TCR–pMHC structures
==========================================

AlphaFold and TCRmodel2 will seat *any* TCR–peptide–MHC candidate in a plausible low-energy pose,
binding or not. The generator's own confidence reports whether a plausible interface could be
built; it does not report whether the chemistry and shape across that interface are those of a real
complex. TCRen2 reads the second question from the coordinates alone, with no binding label and
nothing the generator emits. This page is the decision procedure.

Two commands
------------

The expensive pass — parse, annotate, contact map, descriptors — runs once:

.. code-block:: console

   tcren features   -s models/ -i placement,interface,topology,energetics -o feats.tsv
   tcren recognize  --features feats.tsv -o scores.tsv

``scores.tsv`` carries ``Q`` (interface geometry), ``T`` (footprint shape), and ``S_free`` with its
calibrated ``p_binder``. Join your generator's ``iptm`` /
``plddt`` on the structure-file stem if you want to compose with them; they are not structural
quantities, so tcren does not compute them.

For a per-structure decision rather than a score table, :doc:`reliability` documents ``tcren
assess``, which adds the ranking and the generator diagnostic to the same input.

The three questions the kit answers
-----------------------------------

**1. Does this receptor bind this epitope?**
   ``S_free`` (:func:`tcren.reliability.s_free`) is the recommended answer: three fit-free
   directional blocks — geometry ``Q``, footprint topology ``T``, and the interface energy read
   against the partition function — each divided by its own native spread, so they carry equal
   weight in native-sd units. Nothing is fitted at score time, so it is **defined for a single
   structure** and its value does not depend on what else you scored alongside it.
   :func:`tcren.reliability.p_binder` turns it into a probability through a frozen out-of-fold
   Platt link.

   The cohort-refit posterior that used to sit beside it was discarded in 2.26.0. It refitted a
   latent class per call, so it was undefined for one structure and its numbers depended on which
   rows the fit was anchored on — neither property survives contact with a user who has one model.

**2. Did the generator have a template, and does that matter?**
   It matters enormously, and this is the result the method exists for. Split the VDJdb benchmark
   by whether *some* receptor has already been co-crystallized with that peptide, and every score
   that reads the interface collapses when the template goes:

   .. list-table::
      :header-rows: 1
      :widths: 40 20 20 20

      * - macro ROC-AUC
        - template-covered
        - template-free
        - lost
      * - generator ipTM
        - 0.692
        - 0.555
        - 0.136
      * - interface geometry ``Q``
        - 0.729
        - 0.497
        - 0.232
      * - shape channel ``T``
        - 0.756
        - 0.608
        - 0.148

   The shape of a footprint — how evenly six loops spread their contacts, whether the touched
   surface is one patch or several — is invariant under the rigid-body placement a co-folding model
   is optimizing, which is why it still says something once that model has produced a confident
   pose. Note also that the generator's confidence does *not* fall to warn you: a pose built
   without a template is scored confidently and wrongly.

**3. Is the recognition signal intrinsic, or an artefact of the generator?**
   Intrinsic. On experimental crystals, scoring the native epitope against wrong-epitope decoys of
   the same length on the *fixed* crystal contact map — no AlphaFold, no re-docking — ranks the
   native above the decoys.

Composing with the generator's confidence
-----------------------------------------

The generator's confidence and the structure are two read-outs of the same complex, and composing
them is worth measuring rather than assuming. Report the confidence on its own, the structural score
on its own, and the two together, on identical rows — a combination quoted without both of its parts
cannot be checked.

What the kit does *not* claim
-----------------------------

- **Not affinity.** TCRen2 ranks specificity. What a static interface reads of *dynamics* is
  specific and not uniform: the rupture work tracks the dissociation rate, ride height and coverage
  entropy track the equilibrium free energy where the rupture work does not, and alanine ΔΔG stays
  with molecular dynamics. See :doc:`features` and ``tcren mechanics``.
- **Not a substitute for reporting template coverage.** Nothing computable from the model announces
  which regime a cohort is in — the generator's confidence least of all — so template availability
  is a covariate to report, not one to infer. It needs only a PDB lookup on the peptide.
- **Not a number you can carry between versions.** Every table ``tcren features`` writes is stamped
  with the descriptor catalogue that produced it, and ``tcren recognize --features`` refuses a table
  written under a different one. A score is only meaningful against the catalogue it was computed
  from; recompute rather than re-read.
