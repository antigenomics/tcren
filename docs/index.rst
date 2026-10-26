tcren: rank and inspect TCR–pMHC models
=======================================

Rank modelled TCR–pMHC complexes and inspect the confidence, contact enrichment
and topology behind each score. Start with one fixed pMHC and candidate receptors,
or one fixed receptor and candidate peptides. These are separate ranking tasks.

.. code-block:: console

   $ tcren screen -s models.tar.gz --confidence confidence.tsv -o scores.tsv

The default TCR score uses ipTM, recognition graph cycles and mock-calibrated
contact enrichment. Its weights depend on the unlabelled cohort mean pLDDT.
Scores are dimensionless ranking readouts, not affinity estimates or binding
probabilities. See :doc:`screening` for the complete contract.

.. image:: _static/screening-workflow.svg
   :alt: Choose the ranking task, calculate components, then check score availability.
   :width: 100%

Learn with a worked example
---------------------------

.. toctree::
   :maxdepth: 1

   notebooks/confidence_motif_score
   gallery

Complete a task
---------------

.. toctree::
   :maxdepth: 1

   getting-started
   reproduction
   ranking
   md-validation
   performance

Understand the result
---------------------

.. toctree::
   :maxdepth: 1

   kit
   fit-free
   motif-confidence
   motif-mixture
   residue-surprise
   beyond-the-contact-sum
   background

Look up definitions and APIs
----------------------------

.. toctree::
   :maxdepth: 1

   screening
   mhc-annotation
   features
   descriptor_table
   potentials
   contact-penalties
   prediction
   associations
   evaluation
   modules
   cdr3-neighborhood
   ATOM_PATTERN_SPEC
   INVENTORY_PATTERN_SPEC
   CHEMICAL_CORE_SPEC
   experimental-ranking
   oracle
   assess
   potts
   reliability

.. toctree::
   :hidden:

   notebooks/tcrvdb_frozen_score_reproduction

   notebooks/complementarity_map_2d
   notebooks/pocket_cdr_3d
   notebooks/canonical_frame_figures
   notebooks/pymol_canonical_figures
   notebooks/contact_thresholds_and_bondtypes
   notebooks/mhc_pseudosequence_mps
   notebooks/example_gil_a02_rs_motif
   notebooks/surface_topology

* :ref:`genindex`
* :ref:`modindex`
