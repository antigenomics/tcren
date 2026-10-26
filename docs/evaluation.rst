Evaluation and reproduction commands
====================================

These commands calculate explicitly requested summaries from supplied tables or
raw structures. They do not fit or replace the default recognition score.
Dataset selection, experimental identities and access permissions remain the
caller's responsibility. See :doc:`reproduction` for the end-to-end protocol.

Matched score evaluation
------------------------

``evaluate`` accepts a row table, a JSON mapping of score columns to directions
(``1`` higher-is-better; ``-1`` lower-is-better), identity/grouping columns and
an outcome column. Binary outcomes are 0/1 with missing outcomes retained in
membership; continuous outcomes use their supplied numerical scale.

.. code-block:: console

   $ tcren evaluate --results rows.tsv --scores directions.json --groups epitope --keys id --endpoint y --kind binary --mask common -o evaluation

``common`` uses one finite mask across scores and requested comparators;
``pair`` uses each score/outcome/comparator's finite mask. Outputs contain exact
used/excluded counts, metrics and membership. ``--reference`` compares every arm
with one named score; ``--baseline-scores`` instead names a prior column per arm.

Binary readouts are AUROC, average precision and maximum precision at recall
at least 0.1/0.5. Optional ``--curves`` exports tied-threshold PR/ROC points.
Continuous readouts are Pearson and Spearman correlations; ``--pearson-test``
explicitly adds the two-sided zero-correlation test. Unavailable/degenerate
comparisons remain unavailable. Output directories must be new.

Descriptive diagnostics
-----------------------

``evaluate-diagnostics --operation NAME`` exposes the existing descriptive
operations. Use its help for required column/group arguments:

* ``grouped-summary``: mean, median, extrema or within-group centering.
* ``distribution-summary``: finite-value counts and distribution summaries.
* ``context-correlations``: declared feature/context correlations.
* ``backbone-deformation``: matched-window differences with circular torsions.
* ``backbone-extrema``: declared backbone-window extrema.
* ``regression-bands``: an ordinary least-squares line with 95% Student-t
  mean-response and prediction bands; at least three finite pairs are required.
* ``wild-type-affinity``: the documented supplied ATLAS wild-type source transform.
* ``rank-selection``: explicit maximum-confidence selection against a supplied
  outcome, with first-tie handling.

These are observations and descriptive summaries. Request a fit, test or selection
only when it belongs to the analysis being reproduced. Grouping and missingness
must match that analysis; atom/window counts are not independent biological samples.

Exact protocol reproduction
---------------------------

``ranking-summary --config config.json -o output`` preserves the established
publication-table and tied-metric conventions. The config explicitly names input
panels, score directions, response columns/matrices and boolean inclusion fields;
there are no personal checkout or dataset-path defaults. This specialized route
requires ``Rscript`` with ``data.table``, ``jsonlite`` and ``Rcpp``. It is separate
from ordinary screening and the Python ``evaluate`` command.

``cutoff-calibration`` exposes the fixed native calibration protocol with its
existing cutoffs, decoys and paired bootstrap; it requires ``Rscript`` and
``data.table``. This is explicit calibration reproduction, not score-time work or
an automatic new control. Input source/chemistry and previous-comparison tables
are required. The packaged R resource preserves the established random-number
and interval conventions.

``potential-clusters`` exports explicit reference eligibility and cluster
membership. ``restore-modified-peptides`` restores a declared residue from source
deposits using observed backbone correspondence. ``substitute-peptides`` applies
an explicit template/id/peptide variant table through one batched annotation and
the existing substitution API; it performs no Monte-Carlo refinement.

.. code-block:: console

   $ tcren substitute-peptides --structures templates/ --variants variants.tsv --output variants/ --threads 0

``historical-candidates`` and ``historical-markup`` expose explicitly configured
candidate populations and residue membership for the existing historical protocols.
They do not generate new controls during ordinary scoring. Original/full/LOO/LCO
arms remain distinct and their selection is declared in configuration.

Read each command's ``--help`` before running a historical protocol. Preserve
source/resource hashes, accounting, exact membership and the original protocol
parameters. Optional R requirements do not apply to normal TCR ranking.

Python API
----------

.. automodule:: tcren.evaluation
   :members:

.. automodule:: tcren.evaluation_diagnostics
   :members:

.. automodule:: tcren.diagnostics
   :members:
