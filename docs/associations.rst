pMHC feature–performance associations
=====================================

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

Joint scalar features
---------------------

Use repeatable ``--feature-product LEFT:RIGHT`` options to append raw signed
per-model products before computing cohort summaries:

.. code-block:: console

   tcren associate --results recorded.tsv --features pmhc_by_id.tsv \
     --scores directions.json --feature-columns relief,frac_above_ridge,sasa_hydropathy \
     --feature-product relief:sasa_hydropathy \
     --feature-product frac_above_ridge:sasa_hydropathy -o joint_associations

The additional columns are ``relief__x__sasa_hydropathy`` (Å times the
Kyte–Doolittle index) and ``frac_above_ridge__x__sasa_hydropathy`` (dimensionless
fraction times that index). ``joint_features.tsv`` exports every supplied feature
row and product; provenance records the source pairs. Do not additionally list
generated product names in ``--feature-columns``; they are appended automatically.
Sources must be distinct numeric columns and names must be unique. Missing or
nonfinite sources, or an overflowing product, give an undefined product, never
zero imputation. Missing products never remove score/ipTM classification rows.

For model i with shape feature B_i and exposed hydropathy H_i, J_i = B_i H_i.
The cohort summary is median_i(J_i), not median_i(B_i) median_i(H_i). There is no
centering, cohort normalization, label-selected weight, threshold or fitting.
Declare the source features before inspecting joint results. A product association
describes joint scalar context; it is not evidence of spatial co-localization,
an interaction independent of the component associations, causation, or improved
classification by a new scorer. Compare product/component associations on identical
evaluable cohorts and preserve weaker as well as stronger associations. Product
construction uses one vectorized Polars expression per pair, O(P N) work for
P pairs and N models, without additional workers or dependencies.

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
   :members: associate, read_table, add_feature_products

Spatial co-localization from exact cells
----------------------------------------

``tcren export -s structures.txt --quantity pmhc-cells -t 4 -o cells.tsv``
exports occupied cells and ``cells.tsv.summary.tsv`` from one batched annotation
phase. No potential is needed. Raw-coordinate, residue, searched-index and grid
provenance sidecars apply. The default 64-by-32 grid, fixed extent, 1.4 Å probe,
KD field and padded eight-neighbour smoothing reuse ``surface_map``. The receptor
is removed before the frame is fitted. Removing it does not reverse induced fit.
The original topmost-atom source identifies cell ownership; smoothing can mix
neighbouring peptide/MHC fields, so these are map-level measurements.

``topology.pmhc.pmhc_cells`` owns the coordinate export;
``topology.pmhc.spatial_features`` owns the vectorized table reduction. For finite
peptide-owned cells, let h be height in Å, k the Kyte–Doolittle field and r the
existing 90th percentile of MHC-helix heights. Let e = max(h-r, 0). The summary
retains population covariances Cov(e,k) (``ridge_excess_hydropathy_cov``) and
Cov(h,k) (``height_hydropathy_cov``), in Å times KD index. A positive covariance
means relatively more hydrophobic map values occur at higher cells. Neither
covariance is a calibrated binding or pose-quality score. Their marginal means
and E(e*k) are retained: E(e*k) = E(e)*E(k) + Cov(e,k). This decomposition separates
spatial co-location from a product of marginal summaries, without fitted weights.

Counts distinguish finite peptide heights from common finite h/k cells. The
conditional above-ridge mean is missing when no such cells exist; covariance is
zero when excess height is identically zero. A missing helix ridge leaves
ridge-dependent quantities undefined, including the above-ridge count. Classes
I and II and derived allele calls are explicit; derived calls are not experimental
allele assignments. No class-II core/register is inferred. Duplicate cells, null
keys and inconsistent per-model metadata fail. All source cells stay inspectable.

The reduction takes O(C) expected time and O(C+S) memory for C occupied cells in
S structures, using grouped Polars expressions. Coordinates and annotation retain
the existing batch API's resource contract. Use the summary's unique ``pdb.id``
renamed to ``id`` as the input to ``associate``; retain classes/template strata and
original pair-specific score/ipTM rows. No new fitted score, cutoff or direction
selection is introduced. A single bound-conformation map cannot establish loop
formation, strain relative to an apo ensemble, or ground-truth forced docking.

Pose-reference validation must verify physical chain roles, peptide identity,
observed receptor/MHC sequence spans, construct substitutions and reference
coordinate identities before treating a docking metric as structural truth.
A total over six interfaces can include TCR-internal and peptide:MHC interfaces;
report the four recognition interfaces separately. Comparisons to noncognate or
substituted references remain reference-deviation diagnostics. Preserve all
exclusions, exact per-descriptor/ipTM rows and class/template strata.

The receptor-independent pMHC map is invariant to a rigid TCR displacement with
fixed peptide/MHC coordinates. Its spatial covariances therefore supply context
for receptor-relative coverage, gap, clash and placement measurements rather than
identifying that displacement alone. A singleton verified comparison does not
establish a detector's sensitivity or specificity. The CPL evidence audit is
recorded in ``audit/2026-10-16/CPL_POSE_TRUTH.md``.

Saved predictions and receptor-relative pose checks
---------------------------------------------------

The existing ``features -i placement,interface,topology --no-metadata`` CLI emits
placement, contact coverage, graph topology, gap and steric diagnostics without a
fitted score. Use one manifest and a fixed CPU budget, then compare each descriptor
with the same reference metric and both recorded confidence definitions on
identical finite rows. Preserve undefined descriptors and the full row ledger.

TCRmodel2 describes an interface-specific ipTM calculation that represents TCR and
pMHC as two partners. A run may record both overall and TCR-pMHC-specific ipTM;
retain both rather than inferring their equivalence. By contrast, a multimer DockQ
mean may include internal receptor and peptide:MHC interfaces. Keep the four
recognition interfaces and their contact/RMSD measurements explicit.

``clash_score`` is the sum of peptide-partner heavy-atom overlaps beyond the
existing tolerance, in Angstrom, including TCR and MHC partners. It is neither a
receptor-only displacement measure nor a calibrated energy. Structural agreement,
packing and binding labels answer different questions. Multiple saved predictions
of one complex permit a within-target comparison, not population detector
sensitivity or specificity. Keep clone, class and template strata separate and
record exact-reference template inclusion independently.

Before interpreting external DockQ, constrain physical receptor, peptide and MHC
roles. Permissive automatic sequence mapping can exchange peptide and MHC chains;
reproducing a historical total does not validate that assignment. Keep the raw
mapping, biological role agreement, copy membership and each interface value.
An unbound receptor supplies no receptor:pMHC reference placement: 5NMD contains
unbound 868, whereas 4MNQ contains an affinity-enhanced ILA1 receptor, not SB27.
Existing sequence/construct and class/template gates still apply.

Permissive historical mappings can also combine crystallographic copies and
change the number of nonempty interfaces in the total. Retain that count; do not
assume every historical multimer total averages the same six physical interfaces.
