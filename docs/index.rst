tcren documentation
====================

``tcren`` is a Python re-implementation (and extension) of the TCRen method for structure-based
prediction of T-cell-receptor recognition of epitopes. From one TCR–peptide–MHC structure
(experimental or modelled) it parses and annotates the complex — TCR chains via
`arda <https://github.com/antigenomics/arda>`_, MHC chains mapped against a curated reference and
the groove partitioned — orients it into a canonical frame, computes residue contacts, and scores
every candidate epitope with a residue-level statistical potential derived from TCR:pMHC crystal
structures.

Where the original TCRen scored only TCR↔peptide contacts, this version scores all three interfaces
(TCR↔peptide with TCRen, TCR↔MHC and peptide↔MHC with Miyazawa–Jernigan) for the full binding
picture, and adds mutation ΔΔG, binder classification, pose refinement, and interface mechanics.

What tcren does
---------------

* **Score & rank epitopes** — ``score`` / ``rank`` / ``pipeline``: TCRen energy per candidate, a
  percentile rank against a random background, and the three-interface breakdown + total.
* **Mutation ΔΔG** — ``ddg``: alanine scans and neoantigen substitutions on the native contact map
  (virtual-matrix, no re-docking).
* **Binder classification** — ``binder``: binder vs non-binder for AlphaFold/TCRmodel2 models from
  AF-orthogonal interface geometry.
* **Interface feature table** — ``recognize``: one flat per-structure table — 35 core interface
  descriptors + ``p_real`` / ``p_real_bn`` by default (``--full`` expands to 65 features;
  ``--scores`` adds the frozen ``p_bind`` / ``p_forced``; see :doc:`features`).
* **Annotation & contacts** — ``annotate`` / ``contacts``: TCR CDR/FR, MHC groove helices/floor and
  peptide markup; multi-layer (5/8/12 Å) contact tables.
* **Canonical orientation** — ``orient`` / ``superimpose``: one common MHC frame, docking angles,
  reverse-dock detection.
* **Peptide substitution & refinement** — ``refine``: backbone-preserving substitution plus a
  DOPE-scored Monte-Carlo pose refinement (with CCD/OpenMM/ProMod3/FlexPepDock engines).
* **Potential derivation** — ``derive-potential``: re-derive the TCRen potential (classic/AM/LOO,
  with non-redundancy filtering) from a structure set.
* **QC, mechanics & maps** — steric-clash and register checks, an interface spring-network /
  rupture model, and 2D complementarity maps + 3D pocket/CDR views.

.. note::

   **Ranking, not affinity.** TCRen ranks peptide/TCR *specificity* for a given receptor; it is not a
   binding-affinity model. On the ATLAS SPR benchmark neither the raw contact energy nor its
   poly-alanine difference (:func:`tcren.ddg.reference_delta`) predicts Kd/ΔG/koff/kon (ρ ≤ 0.3 in
   magnitude). The
   one affinity-adjacent quantity a static structure predicts is the off-rate koff, via interface
   mechanics (:mod:`tcren.mechanics`) — not the contact sum.

.. toctree::
   :maxdepth: 2
   :caption: Contents

   getting-started
   oracle
   features
   gallery
   kit
   performance
   modules

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
   notebooks/tcren_analysis

Indices
-------

* :ref:`genindex`
* :ref:`modindex`
