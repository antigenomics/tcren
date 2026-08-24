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

``scores.tsv`` carries ``Q`` (interface geometry), the three channel posteriors ``G``, ``T``, ``E``,
and ``P_native``. Join your generator's ``iptm`` / ``plddt`` on the structure-file stem if you want
to compose with them; they are not structural quantities, so tcren does not compute them.

The three questions the kit answers
-----------------------------------

**1. Does this receptor bind this epitope?**
   ``P_native`` (:func:`tcren.cohort.p_native`) is the posterior of a latent class over three
   channels — geometry, footprint topology, and contact energetics — fitted by expectation
   maximization on the cohort you pass. **No binding label enters the fit.** On a two-epitope
   TCRvdb screen it reaches macro ROC-AUC 0.832 and PR-AUC 0.849; on a 22-cohort balanced
   real-versus-mock VDJdb benchmark, 0.718 and 0.685.

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
      * - ``P_native``
        - 0.721
        - **0.716**
        - **0.005**

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

Worth measuring, and worth measuring honestly. A plain logistic on ``P_native``, ipTM and pLDDT,
fitted and read in sample as a demonstration of complementarity rather than a ranking claim, adds:

- on TCRvdb, **nothing resolvable** — ΔROC +0.008 [−0.009, +0.027], ΔPR +0.012 [−0.012, +0.036];
- on VDJdb, a small PR gain only — ΔPR +0.034 [+0.005, +0.054], ΔROC's interval containing zero.

The generator's confidences rank well on their own; on these cohorts they carry little the
structure does not already say.

What the kit does *not* claim
-----------------------------

- **Not affinity.** TCRen2 ranks specificity. What a static interface reads of *dynamics* is
  specific and not uniform: the rupture work tracks the dissociation rate, ride height and coverage
  entropy track the equilibrium free energy where the rupture work does not, and alanine ΔΔG stays
  with molecular dynamics. See :doc:`features` and ``tcren mechanics``.
- **Not a substitute for reporting template coverage.** Nothing computable from the model announces
  which regime a cohort is in — the generator's confidence least of all — so template availability
  is a covariate to report, not one to infer. It needs only a PDB lookup on the peptide.
- ``P_native`` is **cohort-relative**: it standardizes and fits over the set you pass. Score a whole
  batch of candidates together, not one at a time. For a single structure, use ``Q`` against the
  shipped crystal reference (:func:`tcren.cohort.native_reference`).
