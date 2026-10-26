tcren: rank and inspect TCR–pMHC models
=======================================

Start with a cohort of modelled αβ TCR:peptide–MHC complexes. ``tcren screen``
identifies the ranking task, computes a structural alanine reference on the side
that varies, and reports generator confidence alongside energy, geometry and topology.
The frozen task-matched combination is a dimensionless ranking readout.
It is not an affinity estimate or a calibrated binding probability.

.. image:: _static/screening-workflow.svg
   :alt: Cohort input branches into peptide or receptor ranking, followed by structural reference, component reporting and cohort checks.
   :width: 100%

Start with the ranking task
---------------------------

* **One receptor, candidate peptides:** use ``screen`` for separately modelled
  candidate complexes, or ``score`` to thread candidates on one fixed template.
  These are different structural protocols. ``rank-epitopes`` accepts an existing
  structural-energy/ipTM table for the single-task frozen readout.
* **One pMHC, candidate receptors:** use ``screen``. Keep ipTM primary and inspect
  the separate structural components, cohort flags and missingness.
* **A heterogeneous set or one model:** use ``features``, ``contacts`` or
  ``scoring`` to inspect its quantities. A set without a fixed side is not a ranking
  cohort; small cohorts retain components without a combined score.

.. code-block:: console

   $ tcren screen -s models.tar.gz -o scores.tsv

ipTM is supplied by the generator's metadata or ``--confidence``; tcren does not
compute it. Frozen references ship with the package. No parameter is fitted to a
query cohort. Models with no TCR:peptide contact remain in output with a discard
reason and do not count toward the default 30-model floor.

Read :doc:`getting-started` for installation and input formats, :doc:`screening`
for the component equations and output contract, and :doc:`kit` for inspection
and interpretation. :doc:`fit-free` documents estimated reference dependencies
and distinguishes current scoring from historical fitted readouts.

.. toctree::
   :maxdepth: 2
   :caption: Start and rank

   getting-started
   screening
   residue-surprise
   motif-mixture
   motif-confidence
   ranking
   kit

.. toctree::
   :maxdepth: 2
   :caption: Inspect structures and explain scores

   features
   descriptor_table
   potentials
   contact-penalties
   prediction
   associations
   md-validation
   beyond-the-contact-sum
   background
   gallery
   performance
   modules

.. toctree::
   :maxdepth: 2
   :caption: Reference dependencies and reproduction

   fit-free
   cdr3-neighborhood
   ATOM_PATTERN_SPEC
   INVENTORY_PATTERN_SPEC
   CHEMICAL_CORE_SPEC
   reproduction
   oracle
   assess
   potts
   reliability

.. toctree::
   :maxdepth: 1
   :caption: Tutorials

   notebooks/complementarity_map_2d
   notebooks/pocket_cdr_3d
   notebooks/canonical_frame_figures
   notebooks/pymol_canonical_figures
   notebooks/contact_thresholds_and_bondtypes
   notebooks/mhc_pseudosequence_mps
   notebooks/example_gil_a02_rs_motif
   notebooks/surface_topology
   notebooks/confidence_motif_score
   notebooks/tcrvdb_frozen_score_reproduction

Indices
-------

* :ref:`genindex`
* :ref:`modindex`
