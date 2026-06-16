tcren documentation
====================

``tcren`` is a Python re-implementation of the TCRen method for structure-based
prediction of T-cell-receptor recognition of epitopes. It parses TCR–peptide–MHC
structures, annotates TCR chains via `arda <https://github.com/antigenomics/arda>`_,
computes residue contacts, derives and applies residue-level statistical potentials,
and scores candidate epitopes.

.. toctree::
   :maxdepth: 2
   :caption: Contents

   getting-started
   modules

.. toctree::
   :maxdepth: 1
   :caption: Tutorials

   notebooks/complementarity_map_2d
   notebooks/pocket_cdr_3d
   notebooks/canonical_frame_figures
   notebooks/contact_thresholds_and_bondtypes
   notebooks/tcren_analysis

Indices
-------

* :ref:`genindex`
* :ref:`modindex`
