Final fit-free calculations
===========================

The final scoring boundary permits fixed empirical statistical-potential tables.
Their observed contact counts, pseudocounts, weighting and reference convention
remain estimated dependencies and must accompany the table identity. No Gaussian,
GLM, fitted Potts or fitted transition readout belongs to the final specification.
Q/T/S also use estimated reference moments or correlation weights and remain
historical. These APIs and their measured results are preserved.

Searched annotation references
--------------------------------

ARDA 2.36.2 or later validates the source reference and actual MMseqs database
payloads by content, preserving legacy indexes as evidence. Source-root discovery
recognizes ARDA itself rather than an unrelated enclosing project.
``scoring``, ``features``, ``recognize`` and ``export`` retain ``searched_indexes``
in provenance sidecars: resolved database prefixes and SHA256 identities of the
searched sequence, header and index files, including split payloads. These are
captured around annotation, not inferred from package versions after scoring.
Mid-search changes and mixed generations across input chunks fail explicitly.
An empty identity record does not certify annotation provenance.

Contact sums and references
---------------------------

For a set of observed residue contacts C and a selected pair-potential table V,
the contact sum is Phi = sum over (i,j) in C of V(a_i,a_j). The identities a_i and
a_j are the observed amino acids. The virtual reference difference is
DeltaPhi = Phi(observed) - Phi(reference), with C unchanged. Lower values follow
the energy direction. These are dimensionless statistical-potential sums, not a
calibrated binding free energy. The manuscript reference uses structural alanine truncation: preserve backbone
and C-beta coordinates, remove the remaining side-chain atoms, construct C-beta
for glycine using the existing substitution geometry, and rebuild contacts.
This differs from the historical virtual reference and can change rankings.

``scoring --delta`` references peptide identities. ``dPhi_peptide_sum`` is the
unweighted TCR:peptide plus peptide:MHC reference sum. Historical ``Phi_total`` and
``dPhi_total`` use estimated native-reference scale weights and are excluded from
final scoring; they must not substitute for the unweighted sum. ``export`` references receptor
identities, or exports existing exact decompositions. No fitted replacement is
introduced:

.. code-block:: console

   tcren scoring -s structures/ --delta --reference-mode structural --no-superimpose --tcr-peptide-potential tcren2 -o peptide.csv
   tcren export -s structures/ --potential tcren2 --quantity receptor-reference --reference-mode structural -o receptor.tsv
   tcren export -s structures/ --potential tcren2 --quantity position -o positions.tsv
   tcren export -s structures/ --potential tcren2 --quantity residue -o residues.tsv
   tcren export -s structures/ --potential tcren2 --quantity contacts -o contacts.tsv
   tcren export -s structures/ --potential tcren2 --quantity placement -o placement.tsv

Use ``--reference-mode structural`` explicitly for the manuscript. The default
``virtual`` mode remains available for exact historical reproduction. Output rows
and provenance identify the mode; do not mix modes in a comparison. Structural
peptide references require alanine and use the selected cutoff and contact weighting.
The peptide reference contact map is built once per structure and reused across
both peptide interfaces. Potential sums then operate on the contact tables;
SASA and coordinate-based descriptors still require their structural calculations.
When running multiple outer workers, set ``POLARS_MAX_THREADS=1`` and the BLAS
thread limit before starting Python to avoid nested thread pools.

The receptor export substitutes alanine in all, CDR, CDR1/2 or CDR3 receptor
positions, separately for TCR:peptide and TCR:MHC and their explicit sum.
Position exports use ``position_profile``; residue exports use
``residue_importance``. Leave-one-out ``explain`` is not an equivalent decomposition.
Each contact contributes to both residue endpoints, so the sum of endpoint
contributions is twice the contact sum. Distances and placement are in angstroms.

Identity and accounting
-----------------------

Exports retain full source stems, physical chain IDs, sequential residue indices,
author residue numbers and insertion codes. CDR regions exclude conserved junction
anchors. One batched annotation is used per species plus one MHC annotation batch.
Files, directories, glob/manifests and supported structure archives are accepted.

Every input receives an accounting row, including parsing, annotation and export
errors. Empty/all-failed runs exit nonzero and keep diagnostics. Mixed runs retain
all failures for explicit acceptance by their consumer; a benchmark must verify
its full key set rather than silently discard failed rows.

Sidecars retain residue assignments, chain annotations, installed versions,
reference hashes, the exact potential tables and their hashes. A custom potential
has no inferred derivation or unit calibration. Source paths and annotations may
be private; review them before sharing generated outputs.

Historical assessment availability
----------------------------------

``score_table`` exposes a finite-value ``<readout>_available`` flag per readout,
``assessment_available`` for binder-score availability, and a JSON list of
``missing_descriptors``. Absent columns are marginalized; missing values in
present required columns remain missing. No imputation or refit is performed.
``inspect-models`` omits fitted assessment unless ``--historical-assessment`` is
explicitly requested.

Documentation builds require the ``docs`` extra and the external Pandoc executable
for nbsphinx. IPython supplies the notebook syntax lexer.

Molecular publication consumers
-------------------------------

Scene producers read exact ``position``, ``residue``, ``contacts``, ``placement``
and ``pmhc`` exports. They do not substitute leave-one-out ``explain`` for a
position decomposition or calculate new contact scores during drawing. Physical
chain, author residue number and insertion code identify selections; region
membership comes from annotation. Preserve raw-coordinate, searched-reference,
potential, camera and scene-geometry identities and check the displayed atom
colours against their authoritative rows. Inspect the final PNGs and composed
PDF pages. Keep class II, epitopes and template strata separate. Historical
images and manifests remain archived when labels or pixels change.
