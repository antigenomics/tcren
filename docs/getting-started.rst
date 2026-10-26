Install and rank models
=======================

Install the package and references
----------------------------------

Use Python 3.10 or later in an isolated environment:

.. code-block:: console

   $ python -m venv .venv
   $ source .venv/bin/activate
   $ python -m pip install tcren
   $ tcren build-mhc-ref
   $ tcren info

Binary wheels contain the native kernels. ARDA supplies TCR chain annotation
and obtains its reference/search backend on first use. The MHC sequence reference
is prepared by ``build-mhc-ref``; this is a one-time reference download, not
score fitting. Reference preparation needs network access. Prepare references
before moving to an offline compute environment.

For a source checkout, run ``bash setup.sh`` and activate ``.venv``. Re-run setup
after changing native code. Do not mix a development import with an installed-wheel
validation run.

.. _mhc-reference:

MHC annotation
--------------

MHC class, chain role and groove annotation support the structural contact graph.
Exact allele identity is not a default score input. A cropped structure sequence
may match several alleles equally well; its representative match does not replace
known experimental restriction. Pseudosequence extraction is a separate optional
annotation operation, not a prerequisite for the ranking formula.

.. _data-locations:

Reference locations
-------------------

``TCREN_HOME`` selects the reference root when set. Otherwise a source install
uses its checkout, while a wheel uses ``XDG_CACHE_HOME/tcren`` (normally
``~/.cache/tcren``). MHC sequences live in ``database/mhc/`` and search indexes in
``data/mhc_cache/``. ``TCREN_DATA_DIR`` overrides the structure-data directory.
These are reference assets; calculated results are written to the requested output.
Use ``tcren fetch-data`` for the native/canonical structure sets needed by
orientation and relevant examples. The default score's frozen mock calibration
is packaged and does not require those mock structures at scoring time.

Prepare a receptor-ranking cohort
---------------------------------

Supply separately modelled complexes containing paired alpha/beta TCR chains,
a peptide and MHC. Use one peptide–MHC background with candidate receptors.
Keep all unlabelled candidate models in the confidence context, including any
whose experimental outcomes are ambiguous.

A confidence TSV has one row per input structure stem:

.. code-block:: text

   complex.id    iptm    plddt
   model_001     0.82    93.4
   model_002     0.76    91.8

Use actual tab separators. ipTM is on [0,1] and pLDDT is on [0,100]. The example
values illustrate the schema, not a dataset. TCRen reads generator confidence;
it does not compute ipTM. Never use crystal B factors as pLDDT. The set's
``metadata.tsv`` can provide confidence instead of an explicit table.

.. code-block:: console

   $ tcren screen -s models.tar.gz --confidence confidence.tsv -t 0 -o scores.tsv

``-s`` accepts PDB/mmCIF, gzipped structures, directories, tar.gz archives,
quoted globs, a text manifest of paths or repeated arguments. Batch the inputs
in one invocation. ``-t 0`` uses the allocated CPU budget; a positive number
bounds it. Avoid a separate annotation process per model.

Check and interpret the output
------------------------------

#. Confirm all expected structure IDs are present and unique.
#. Read ``score_unavailable`` and ``no_interface`` before sorting scores.
#. Sort finite ``screen_score`` values from highest to lowest.
#. Keep raw ipTM, topology, enrichment, pLDDT context and missing-value flags.
#. Save the output provenance sidecar with the table.

The default requires 30 models with an interface; smaller cohorts retain
components without a combined score. For individual or small sets, use the
explicit confidence-context workflow in :doc:`motif-confidence`. A score is a
ranking readout, not a thresholded binding call. See :doc:`kit` for interpretation
and :doc:`screening` for exact equations, imputation and flags.

Rank peptides or inspect a structure
------------------------------------

For separately modelled peptides against one fixed receptor, use ``screen`` with
the epitope-ranking task. To thread peptides on a single fixed template instead:

.. code-block:: console

   $ tcren score -s template.pdb -c candidates.txt -o ranked.csv

These protocols have different contact maps. A frozen existing epitope-energy/ipTM
table can be scored with ``rank-epitopes``; see :doc:`ranking`.

For structural quantities and annotation:

.. code-block:: console

   $ tcren contacts -s models.tar.gz --interface tcr_peptide -o contacts.csv
   $ tcren scoring -s models.tar.gz --delta --reference-mode structural -o energy.csv
   $ tcren annotate -s model.cif.gz --regions all -o regions.csv

The structural alanine reference truncates side chains and rebuilds contacts.
It is distinct from virtual substitution on a fixed map. Use :doc:`features`
for descriptor exports, :doc:`potentials` for potential definitions and
:doc:`reproduction` for reproducible calculations.
