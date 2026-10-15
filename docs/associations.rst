pMHC feature–performance associations
====================================

Evaluate recorded descriptors/scores against ipTM, then relate their performance
to pMHC surface properties. This is descriptive evaluation, not a fitted score.
MHC classes, template coverage, datasets and tasks are never pooled.

Derive features
---------------

.. code-block:: console

   tcren export -s structures.txt --quantity pmhc -o pmhc.tsv -t 1

Files, directories, globs, archives and text manifests are supported. Full stems
identify rows; duplicates fail. Searches are batched. Accounting and annotation
sidecars retain input hashes, searched-index identities and failures. Empty or
all-failed exports exit nonzero. No potential is needed or used.

The receptor is removed before constructing the groove frame or calculating
exposure. Coordinates remain unchanged: a receptor-stripped bound model is not an
unliganded experimental pMHC. Peptides require standard amino acids; class II
requires both annotated groove chains. Core/register boundaries are not guessed.

.. list-table:: Feature definitions
   :header-rows: 1

   * - Column
     - Definition and units
   * - ``relief`` / ``peak_to_valley``
     - Standard deviation / range of peptide-owned surface heights, Å.
   * - ``frac_above_ridge``
     - Peptide-cell fraction above the 90th percentile of helix surface heights.
   * - ``area_frac_peptide``
     - Peptide-owned / all occupied map cells; projected area, not SASA.
   * - ``phobic_mean`` / ``phobic_centre``
     - Mean surface Kyte–Doolittle value, whole peptide / central spatial third.
   * - ``charge_mean``
     - Mean surface formal side-chain charge (histidine 0.1).
   * - ``peptide_sasa_A2`` / ``sidechain_sasa_A2``
     - Peptide heavy-atom / side-chain solvent-accessible area in pMHC, Å².
   * - ``sasa_hydropathy``
     - Side-chain-SASA-weighted Kyte–Doolittle value.
   * - ``peptide_length``
     - Number of resolved peptide residues.

SASA uses Biopython Shrake–Rupley with 100 sphere points, its element radii and a
1.4 Å solvent probe, after transforming coordinates into the pMHC groove frame.
Atoms are filtered only when they cannot occlude any peptide expanded sphere.
With exposed side-chain areas A_i (Å²) and Kyte–Doolittle values h_i:

.. math::

   H = \frac{\sum_i A_i h_i}{\sum_i A_i}.

Glycine has no heavy side-chain atoms. Zero total area gives an undefined H.
The existing ``pep_free_frac`` and whole-interface ``burial`` are different,
receptor-dependent quantities. Surface maps use the existing 64×32 raster over
40×50 Å. Raster and sphere discretization are approximations without a universal
error bound. Neither static quantity measures flexibility.

Compare scores
--------------

Results require unique ``id``, ``dataset``, ``task``, ``cohort``, ``mhc_class``
(``I``/``II``), ``template_coverage`` (explicit ``unknown`` if unavailable),
binary ``y``, ``iptm`` and score columns. Define receptor cohorts by peptide–MHC
pairing, not peptide alone; peptide panels are a separate task. Features require
unique ``id`` plus selected numeric columns. Rename exported ``pdb.id`` to ``id``
through an exact validated mapping. CSV, TSV and Parquet are supported.

The scores JSON supplies fixed directions, for example:

.. code-block:: json

   {"Phi_tcr_pep": -1, "m_erank_tp": 1, "plddt": 1}

Energy is lower-is-better. Geometry uses a documented fixed direction, never a
direction chosen from labels. ``pitch``/``pitch_angle`` are rejected as leakage.
Exclude historical fitted scorers from final-model manifests.

.. code-block:: console

   tcren associate --results recorded.csv --features pmhc_by_id.tsv \
     --scores scores.json --feature-columns relief,frac_above_ridge,sasa_hydropathy \
     -o associations

Output must be a new directory. ``metrics.tsv`` contains paired values/deltas,
directions, positive/negative and used/excluded counts. ``membership.tsv`` retains
every row's inclusion per arm. ``features.tsv`` contains label-blind medians,
quartiles and missing counts. ``associations.tsv`` contains descriptive Spearman
correlations with score performance, matched ipTM performance and their difference.
``provenance.json`` retains input hashes, settings and software identity.

Algorithm and limits
--------------------

For each cohort and score s, apply its fixed direction d_s in {-1, 1}. Keep rows
where s and ipTM are finite; all labels must be finite binary values. Calculate
both metrics on this exact set, then subtract ipTM's metric. Missing unrelated
scores or surface features never change this paired set.

ROC-AUC uses half-credit for ties. Average precision is the non-interpolated
precision–recall step integral. Precision at recall 0.1/0.5 uses the established
benchmark convention: maximum precision at thresholds attaining that recall.
These fixed reporting endpoints do not tune a scorer. Empty/single-class pairs
have undefined metrics.

Summarize each feature once across all supplied cohort rows, without labels or
confidence filtering. Correlate those medians with performance using Pearson
correlation of average ranks (Spearman), one equally weighted cohort per
observation, strictly within dataset/task/class/coverage strata. Fewer than three
finite cohorts or a constant variable produces an explicit unavailable status.
Matched ipTM associations can differ across arms because their support differs.

Repeated models are not independent epitopes. Label-blind model summaries still
describe receptor-conditioned bound conformations. A CPL median describes a
candidate-peptide panel rather than a fixed pMHC. Associations can reflect allele,
length, generator or cohort composition, not causation. No fitted coefficients,
calibration, significance tests or multiplicity-adjusted discovery claims are
produced. Class-II structures without labelled ipTM rows establish feature
coverage only.

For N rows, S scores, F features and G cohorts, per-cohort sorting is bounded by
O(S N log N) for metrics and O(F N log N) for summaries. Associations cost
O(S F G log G). Explicit membership output contains N S rows. NumPy/Polars,
sklearn metrics and SciPy ranks supply the computations, without new dependencies
or worker pools. Tables/output stay in memory; split large analyses by independent
dataset/task strata. Structural export cost is separate and depends on atom counts.

API
---

.. automodule:: tcren.topology.pmhc
   :members:

.. automodule:: tcren.associations
   :members: associate, read_table
