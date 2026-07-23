A kit for AI-generated TCR–pMHC structures
==========================================

AlphaFold / TCRmodel2 will seat *any* TCR–peptide–MHC candidate in a plausible low-energy pose — a static
"forced" pose that is not a Gibbs sample. The generator's own confidence (ipTM) reads largely as pose
canonicality / template coverage. TCRen adds an **intrinsic, generator-orthogonal** read of the same
structure: is the interface *natural*, and does it *engage* like a binder? The two are complementary, and
combining them is synergistic. This page is the decision procedure.

One command produces every score:

.. code-block:: console

   tcren recognize --full --scores -s models/ -o kit.tsv

giving, per structure: the 65 interface descriptors (:doc:`features`), the binder score ``p_bind``, the
forced-pose flag ``p_forced``, and the wrong-TCR flag ``p_real`` / ``p_real_bn``. Join your AlphaFold ``iptm``
(it is not a structural quantity, so tcren does not compute it) to the table on the structure-file stem.

The three questions the kit answers
-----------------------------------

**1. Is this AF model trustworthy, or "too good to be true"?**
   ``p_forced`` (:func:`tcren.recognition.forced_pose_score`) is a frozen strain classifier trained only on
   provenance (crystal-natural vs AF-forced, 5-fold AUC 0.762). High ``p_forced`` = the interface is
   stretched / thin / one-sided — residues placed to *mimic* a good contact energy rather than to bind.
   Among AF-*confident* poses (high ipTM), the high-strain ones are enriched for non-binders — the QC signal
   ipTM cannot give you.

**2. Does this TCR bind this epitope?**
   ``p_bind`` (:func:`tcren.binder.binder_score`) is the AF-orthogonal binder score (TCRvdb
   **raw-label** macro AUC 0.796, pooled 0.810; AF ipTM 0.794 / 0.793). Use it to screen many TCRs
   against one epitope. Label denoising is a separate algorithm and is not benchmarked here.

**3. Combined call — the synergy.**
   :func:`tcren.recognition.kit_score` = ``z(p_bind) + z(iptm)`` over the cohort — a fixed, no-fit
   combination. On TCRvdb raw labels it beats **either score alone** at precision:

   .. list-table::
      :header-rows: 1
      :widths: 34 22 22 22

      * - score
        - macro-PR
        - P@10% recall
        - P@20% recall
      * - AF ipTM
        - 0.782
        - 0.861
        - 0.816
      * - tcren ``p_bind``
        - 0.804
        - 0.912
        - 0.873
      * - ``kit_score`` (``p_bind`` + ipTM)
        - **0.847**
        - **0.969**
        - **0.939**

   Δ macro-PR vs ipTM = **+0.065** (95% CI [+0.022, +0.100], P(Δ>0)=1.00) for the no-fit
   z-sum shown above. A CV-honest leave-epitope-out logistic on the same two inputs confirms it
   more conservatively at +0.041 ([+0.005, +0.076], P=0.99) — that is a *different estimator*,
   not this row.

   .. note::
      ipTM is the **weakest** of AlphaFold's three confidences on this task. Against global pLDDT
      (macro-PR 0.808) the margin is ``+0.039``, not ``+0.065``. Quote the baseline you measured
      against.

   The combination also **corrects AF's errors**:
   strain flags AF false-positives among confident poses (AUROC 0.633), and ``p_bind`` rescues AF
   false-negatives among under-confident poses (0.732 vs ipTM 0.697).

``kit_score`` is cohort-relative (``z`` standardizes over the set you pass) — score a whole batch of AF
models together, not one at a time.

What the kit does *not* claim
-----------------------------

- **Not "beat AF" everywhere.** On the harder VDJdb-AF real-vs-mock task, combining does not beat ipTM
  (macro 0.639 vs 0.656); there TCRen's contribution is the interpretable forced-pose *gradient*
  (crystal < AF-real < AF-decoy), not a discrimination win.
- **Not affinity.** TCRen ranks specificity, not Kd/ΔG/koff (see the note on the landing page).
- ``kit_score`` needs the generator ipTM as input; the purely structural scores (``p_bind``, ``p_forced``,
  ``p_real``) do not, and also work on crystals with no generator at all.
