Reproduce structural calculations
=================================

Use a pinned installed TCRen wheel for calculations and retain its version and
artifact hash. Dataset manifests define the scientific population; TCRen owns
annotation, descriptors, scoring and evaluation. Plotting and manuscript builds
consume recorded tables and do not silently rerun calculations.

Freeze the inputs
-----------------

Record raw-coordinate hashes, generator confidence, exact paired receptor/pMHC
identities, supplied outcomes and the selected reference resource. Keep ambiguous
outcomes in unlabelled cohort context, excluding them only from binary evaluation.
Keep MHC class, epitope and template coverage explicit. Never infer a label from
an identifier, confidence value or published aggregate validation fraction.

Reference assets and final accepted results are distinct from intermediate
calculations. For a fresh reproduction, rebuild required descriptors from raw
inputs; do not silently resume from old feature rows. Use a new output directory
and retain command lines, installed version, failures and all expected input IDs.

Reproduce receptor scores from coordinates
------------------------------------------

Prepare a separate manifest and confidence table for each fixed-pMHC cohort:

.. code-block:: console

   $ tcren screen -s cohort.txt --task tcr --confidence confidence.tsv -t 0 -o screen.tsv

Confidence uses unique ``complex.id``, ``iptm`` and ``plddt`` columns. The default
uses the measured cohort-mean pLDDT and the frozen three-component rule described
in :doc:`screening`. Do not substitute the legacy two-weight reference or a
single-model proxy when reproducing the default. Keep withheld/no-contact rows
and audit row membership before evaluating finite scores.

For an exact replay from already exported components, use
``motif-confidence-score`` with the recorded frozen resource and complete
unlabelled cohort context. The notebook
:doc:`notebooks/tcrvdb_frozen_score_reproduction` demonstrates that route.
The retained input bundle is not distributed with the package; that notebook is
an optional supplied-data replay, not a self-contained public-data tutorial.
Replay verifies application of a rule; it does not establish fresh annotation
or contact extraction from raw coordinates.

Reproduce other structural quantities
-------------------------------------

.. code-block:: console

   $ tcren scoring -s models.txt --delta --reference-mode structural -o peptide-energy.csv
   $ tcren export -s models.txt --quantity receptor-reference --reference-mode structural -o receptor-energy.tsv
   $ tcren export -s models.txt --quantity recognition-betti -o graph.tsv
   $ tcren export -s models.txt --quantity contact-atoms --interface tcr_peptide -o contacts.parquet
   $ tcren export -s models.txt --quantity atom-inventory --interface all -o atoms.parquet

Use ``score`` for fixed-template candidate peptide threading, ``derive-potential``
for an explicitly specified potential derivation, and ``md``/``mechanics`` for
trajectory/mechanical readouts. Their command help defines required inputs.
Do not replace structural alanine references with fixed-map virtual substitution.
Do not treat original, updated, full-reference, leave-one-out and leave-cluster-out
potentials as the same baseline. :doc:`fit-free` and ``SOURCES`` describe the
parameter and annotation-reference dependencies.

Evaluate matched exported rows
------------------------------

Join outcomes and scores by exact IDs in the dataset-owning project, checking
one-to-one correspondence. The resulting table may contain several score columns.
A JSON file maps each score column to its favourable direction, ``1`` for higher
or ``-1`` for lower. Then calculate metrics through the CLI:

.. code-block:: console

   $ tcren evaluate --results rows.tsv --scores directions.json --groups dataset,stratum,epitope --endpoint y --kind binary --mask common --reference iptm --curves -o evaluation

``--mask common`` keeps the same finite population across compared scores;
``pair`` explicitly evaluates each score on its own available rows. Retain the
membership output, exact positive/negative counts and every comparator delta.
Use continuous evaluation for a measured continuous endpoint:

.. code-block:: console

   $ tcren evaluate --results rows.tsv --scores directions.json --groups dataset,epitope --endpoint response --kind continuous --mask common -o correlations

Request ``--pearson-test`` only when reproducing an analysis that includes that
test. A rank score is not a binding probability, and an adaptive development
comparison is not an untouched validation panel. No extra threshold search or
new fitting is implied by evaluation.

Run a dataset-specific campaign
-------------------------------

The source checkout's ``reproduce.sh`` delegates to an explicitly supplied driver:

.. code-block:: console

   $ TCREN_REPRO_PYTHON=/absolute/env/bin/python TCREN_REPRO_DRIVER=/absolute/project/driver.py ./reproduce.sh --help

The driver owns raw-source selection and the list of CLI stages; it must not
reimplement the scientific calculations. The wrapper contains no default personal
checkout paths or private data locations. Full campaigns must preserve each
block's input population, reference version and output contract. Completed tables
can subsequently be formatted into figures without rerunning the campaign.

Verify before claiming reproduction
-----------------------------------

Check exact IDs/counts first, then numerical outputs at their declared tolerance,
then aggregate metrics and all gains/losses. Report missing raw sources or failed
stages explicitly. A successful package test, frozen-score replay, documentation
build or figure render alone does not establish full scientific reproduction.
Review provenance sidecars before sharing: they can contain local paths and
private source identifiers even when the score table uses neutral IDs.
