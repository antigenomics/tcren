Generated-model inspection
==========================

``tcren inspect-models`` connects native AF3/TCRmodel2 confidence outputs to the same
annotation, geometry, descriptor and residue-contact APIs used on deposited structures.
It serves the ipTM/pLDDT diagnostics block. No separate analysis package is required.

Prepare a TSV manifest with unique ``complex.id`` and ``structure`` columns. Optional
columns are ``generator`` (``af3``, ``tcrmodel2``, ``boltz2``, ``tfold-tcr``,
``af2-multimer``, ``protenix`` or ``rf3``), ``confidence``,
``detail_confidence``, ``producer_structure`` and ``model_key``. Relative paths resolve against the manifest.
For AF3, point ``confidence`` at the summary-confidence JSON and
``detail_confidence`` at the full-confidence JSON. For TCRmodel2, point ``confidence``
at ``statistics.json`` and give the literal ``model_key`` (for example ``ranked_0``).
Native references omit confidence. Use a unique ID across jobs: two files named
``ranked_0.pdb`` are two different models.

.. code-block:: bash

   tcren inspect-models manifest.tsv analysis/ --threads 2

The output directory contains:

* ``models.tsv``: manifest identity, literal generator scalars, interface contact counts,
  descriptor availability and all four default descriptor families, including docking angles.
* ``residues.tsv``: TCR CDR/framework, peptide and MHC residue markup from one batched
  annotation pass. Residue indexes are zero based within each coordinate chain.
* ``chains.tsv``: chain roles, sequences, original author residue IDs and the library's
  allele/germline annotation. ``cdr3`` is the arda CDR3 region, not anchor-inclusive junction.
  Allele/germline values are sequence-reference annotations, not experimentally established
  genotypes; a truncated MHC groove does not necessarily distinguish full-resolution alleles.
* ``contacts.tsv``: directed, annotated residue pairs, closest heavy-atom distances,
  heavy-atom pair counts, and available global/chain-pair/token-pair confidence.
* ``assessment.tsv``: historical fitted readouts only with ``--historical-assessment``;
  disabled by default. Existing frozen tcren readouts on complete alphabeta TCR:pMHC
  structures; ``confidence_residual`` and ``binder_iptm`` are added where ipTM exists.
  The frozen model's calibration/hold-out provenance remains the one shipped by tcren;
  these readouts are not a new validation of transfer to AF3 or of binder discrimination.
  ``assessment_status`` and ``missing_descriptors`` distinguish unavailable complete-case
  readouts from numerical values without imputing missing interface evidence.
  Each emitted readout has a corresponding ``*_available`` boolean, including
  ``confidence_residual_available`` and ``binder_iptm_available`` where ipTM exists.
  Without the historical opt-in, the assessment table is empty and complete models
  report ``assessment_status=historical_not_requested``.
* ``angles.tsv``: existing ``docking_angles`` and ``tcr_placement`` API values, including
  signed crossing and incident angles and CDR height/shift/offset. These readouts retain their
  existing definitions; this export does not add a new descriptor to the catalogue.

Every table has a provenance sidecar. Parsing, annotation, confidence mismatch or
featurisation failures raise: an empty or failed batch must not look successful.
TCR-only and pMHC-only inputs receive annotation and contacts, with
``descriptor_status=incomplete_construct``; ternary descriptors and assessment stay
unavailable. No missing chain is invented.

Interpret the evidence separately
---------------------------------

Global ipTM is repeated per model/contact for joining; it is not a contact probability.
AF3 chain-pair ipTM and minimum PAE retain chain order and direction. Full AF3 confidence
JSON supplies separate token-pair PAE and contact probability. Protein token identities
are matched by coordinate chain ID and author residue number, never by a guessed sequential
index. Multiple tokens per residue (ligands/modifications) are outside this inspection surface.
AF3's native token contact probability is the distogram's <8 Å event; it is not
the probability of tcren's default inclusive 5 Å minimum-heavy-atom contact.
Keep the native probability, computed contact definition and distance cutoff separate.
The current TCRmodel2 statistics output supplies scalar confidence only; deleted upstream
PAE arrays are not reconstructed. Numeric strings are accepted while original metric names
are retained. A fused pMHC may report pTM and interface pLDDT without ipTM: absence stays null.
Do not compare these as interchangeable confidence estimates across construct types.

Contacts are symmetrised: each residue pair appears in both directions. An interface
selection in ``ContactMap`` chooses one biological direction. Empty interfaces remain
explicit in ``models.tsv`` even when there are no rows to write in the contact table.
The contact export defaults to 5 Å; ``--cutoff`` changes that export and its counts.
Descriptor definitions retain their catalogue thresholds (including the standard 5 Å
contact definition), regardless of the export cutoff. An angle describes rotation;
CDR placement and contact counts expose displacement that an angle alone cannot detect.

Library use
-----------

.. code-block:: python

   import polars as pl
   from tcren.prediction import inspect_predictions

   tables = inspect_predictions(pl.read_csv("manifest.tsv", separator="\t"), threads=1)
   models, contacts = tables["models"], tables["contacts"]
   # Filter contacts by complex.id, chain.type.from/to, region.type.from/to or residue.index.
   # assessment.tsv joins to models on complex.id; contacts retain the same model key.

``recognition_table(..., annotate=False)`` reuses already annotated ``Structure`` objects
for callers that also need contacts/markup. ``tcren.validation.structures`` provides
sequence-verified role mapping, coordinate diagnostics and chain-pair contact coverage.
These proximity diagnostics are not MolProbity clashscores. ``tcren.validation.md`` reads
recorded GROMACS XVG datasets and plots short-run readouts without running MD again.
Whole-system potential energies from different solvent boxes are not binding energies;
20 ps execution checks do not establish equilibration.

``compatible_reference_models(models, reference_roles, construct_roles)`` selects
queries covering every explicit reference role, including complete-complex queries
for a pMHC-only reference. The caller supplies the construct-to-role mapping and
uses ``map_sequence_roles`` to verify actual coordinate sequence spans. This
selection helper does not run a campaign comparison. In a partial-reference fit,
``all`` must mean all matched components; unmatched receptor coordinates must not
be scored against a pMHC-only reference (issue #22).

For multiple descriptor workers, call the API from an importable script inside an
``if __name__ == "__main__":`` guard, or use the CLI above. The current spawn backend
cannot import stdin/notebook entrypoints; use one worker there. Issue #21 tracks an
interactive-safe parallel backend without forking after annotation/BLAS threads.

Building this documentation locally
-----------------------------------

Install the ``docs`` extra and Pandoc (``brew install pandoc`` on macOS). The
notebook pages require this external document converter, as does the existing docs CI.
Then run ``.venv/bin/sphinx-build -W -b html docs docs/_build/html``.

Retaining TCRmodel2 residue confidence
--------------------------------------

TCRmodel2 normally removes its raw prediction results during finalisation. For a running
local producer, ``tcren.validation.tcrmodel2.export_generated_confidence`` exports
its PAE and per-residue pLDDT before that cleanup, retaining checksums of the original
result and unrelaxed PDB.
This reads trusted locally generated pickles, not arbitrary interchange files. Exported
JSON records exact producer-chain sequences and sequential residue positions.

Point ``detail_confidence`` at the retained JSON. The inspection API verifies its ranking
confidence against the requested rank and remaps tokens by unique exact sequence spans
when postprocessing renames chains and applies AHO author numbering. Contacts then carry
``token_pair_pae`` and ``residue_plddt.from/to``. No AF2 chain-pair ipTM or residue-contact
probability is invented. ``ranked_detail_file`` selects retained results by unique exact
producer ranking confidence, refusing ties instead of guessing. Earlier completed runs
whose raw results were removed remain scalar-only.

MD and controlled AF3 follow-up
-------------------------------

Validation helpers audit literal AF3 chains, MSA query sequences and template index bounds;
paired-MSA removal preserves other search outputs as an explicit inference-only experiment.
``tcren.validation.md.component_relative_rmsd`` fits a selected component and measures
its partner in that frame, using corresponding backbone coordinates in Å after PBC correction.
Atom selections must be one-dimensional integer indices; fractional, boolean and
out-of-range indices are rejected before coordinate selection.
``minimum_potential_frame`` requires exact saved-frame/energy timestamps in ps and kJ/mol.
It selects whole-system potential, never binding energy. ``close_atom_pairs`` exports atom
identities and source ``altloc.from/to`` labels for nonlocal heavy-atom proximity pairs,
with ``alternate_location_involved`` and possible-disulfide flags. Bond classification
is explicitly unknown; adjacent-residue exclusions are a sequence heuristic, not a
bond graph. Different altloc letters across residues do not establish conformer
incompatibility. All parsed conformers still contribute to contact geometry; writers
retain their existing one-conformer selection and preserve the selected label. This is
not a MolProbity clashscore. See issues #25 and #26. Missing original MD settings must remain explicit.

The native inspection also exports ``interfaces.tsv``: directional AF3 chain-pair confidence
with observed directed-contact counts, including zero-contact pairs. AF2 pair matrices remain
unavailable; global/named confidence is retained in models. The library helper
``confidence_chain_pairs`` exposes this table without treating confidence as contact probability.

Modern predictor controls
-------------------------

``tcren.validation.prediction_inputs.af3_to_boltz_input`` exports literal protein
sequences and retained unpaired MSAs to a Boltz YAML/CSV input. All CSV pairing keys
are -1; lowercase insertions are preserved. No templates, paired MSAs, steering or
affinity setting is silently transferred. The output directory must be new, and
the provenance records source and MSA checksums. This is an inference control,
not a fresh database search. Native confidence ingestion for additional providers
is tracked in issue #28; unsupported fields must remain unavailable.

``map_sequence_roles(..., allow_extended=True)`` explicitly permits an exact
requested domain inside a longer output chain. It returns coordinates only for
that span and records both reference/output offsets. The default remains strict.
``geometry_and_contacts`` records the actual ``cutoff_A`` and generic pair counts;
legacy ``*_le4_5A`` aliases appear only for a 4.5 Å calculation (issue #27).

``inspect_md_trajectory`` streams a GROMACS multi-model PDB after explicit PBC
correction. It requires TITLE timestamps, exact polymer sequences and matching
backbone atom identities. Supply named components as chain-ID lists and a
two-column potential-energy XVG. It returns partner-relative RMSD rows in Å/ps,
minimum saved-frame metadata and that frame's PDB text. The reference structure
must be specified explicitly; selection remains whole-system potential energy.

Native Boltz-2 confidence
-------------------------

Use ``generator=boltz2`` with the writer's ``confidence_MODEL.json`` and an
explicit ``producer_structure`` pointing to its processed structure NPZ.
``detail_confidence`` points to ``pae_MODEL.npz``; the corresponding
``plddt_MODEL.npz`` must exist beside it. Only standard all-protein inputs with
all chains valid are supported. The adapter validates chain order, unique IDs,
residue spans, full coordinate sequences and matrix dimensions. It preserves
directional chain-pair ipTM and the original ``confidence_score``; the latter
is never renamed ipTM. Boltz's 0–1 token pLDDT is explicitly converted to the
0–100 pLDDT units used in residue-contact exports. No contact probability or
AF3-specific minimum-PAE matrix is invented. All evidence files are checksummed.

``compare_same_frame_energy_forces`` compares one identical-coordinate rerun
against a numerical reference, requiring equal energy/force timestamps and three
force columns per atom. It reports total-potential difference (kJ/mol), vector
force RMS/relative error and percentile/max errors (kJ/mol/nm). The caller must
verify coordinate/topology identity and System-group atom order. These numerical
readouts do not establish affinity or structural quality.

For ``generator=tfold-tcr``, retain native scalar scores, ordered token/chain
sequences and per-residue pLDDT from the producer in one JSON bundle. The reader
requires explicit sequence indexing and ``probability_0_1`` pLDDT units, verifies
full coordinate sequences, and exports pLDDT in 0–100 units. Global ipTM remains
global; unsupported chain-pair confidence, PAE and contact probability remain
unavailable. Producer execution and generation of the retained bundle are separate
from library inspection.

``generator=af2-multimer`` accepts the analogous retained native AF2 bundle,
requiring ``score_0_100`` pLDDT units and preserving a supplied token PAE matrix.
No chain-pair ipTM is imputed. The generic Multimer baseline records the official
source and retained-MSA/no-template inputs separately from TCRmodel2.

MD preprocessing must make molecules whole, establish a compact first-frame
complex image, then unwrap continuously using that coherent first-frame
reference. A TPR containing wrapped starting coordinates can anchor ``nojump``
to broken chains. The reader rejects consecutive Cα gaps above 8 Å in this
complete-chain canonical MD scope (issue #29). If the first frame contains CRYST1,
``validate_initial_complex_image`` checks that every pair of chain Cα centroids is
in a nearest periodic image, using the full triclinic cell. A closer lattice image
raises an error before RMSD reporting. The returned ``initial_image_status`` is
``nearest_image_chain_centroids`` on success or ``unavailable_no_unit_cell`` when
CRYST1 is absent. This compact-centroid check neither unwraps coordinates nor
establishes docking correctness. It is applied only to the first frame, preserving
later displacement and unbinding. Independently fitting partners,
or clustering every later frame, must not conceal relative motion/unbinding.
``energy_drift`` reports NVE total-energy slope per ns and per atom; it does not
treat correlated frames as independent replicates or normalize by a potential zero.

``inspect_md_trajectory(..., contact_cutoff=5.0, docking=True)`` additionally
reports zero-inclusive chain-pair contact counts, observed residue-pair occupancy
over saved frames and native docking/placement descriptors for each frame. The
reader verifies heavy-atom identity before comparing contacts. Docking annotates
the reference once and transfers exact-sequence region markup onto each frame's
own residues. Occupancy is a fraction of saved frames, not an independently
sampled binding probability; never-observed residue pairs are omitted.
The caller owns the whole-molecule/compact-first-frame/continuous-unwrapping
preprocessing, must retain its commands and source hashes, and must provide the
explicit RMSD reference. The reader returns first and final saved coordinates,
exact frame times and the minimum-potential frame. Post-MD coordinates have no
newly inferred predictor confidence; keep the original prediction confidence separate.

``validation.prediction_comparison.compare_peptide_contact_sets`` compares
directed inspection contact tables against an explicit predicted comparator
over exact TRA/TRB/peptide sequences. It reports shared, reference-only and
candidate-only residue pairs and Jaccard agreement. A contact-free complete
matching complex is distinct from a missing or different sequence, which remains
unavailable. Both tables must use the same inclusive heavy-atom cutoff. Empty
unions remain undefined. Agreement is not native accuracy or binding evidence.
AF2 retained confidence tokens may use input chain IDs that differ from the PDB
writer output. The contact join maps exact full-chain sequences with a complete
one-to-one correspondence; it rejects ambiguous renamed identical chains and
sequence differences, including on contact-free structures. PAE direction and
token pLDDT are preserved without rewriting the producer bundle.

Protenix and RF3 summary confidence
-----------------------------------

``generator=protenix`` and ``generator=rf3`` accept native summary JSON scalar
values, retaining producer names such as ``iptm``, ``ptm`` and ``ranking_score``.
A ranking score is never substituted for ipTM. Native chain-pair arrays and
residue/token confidence remain unavailable until version-pinned producer fixtures
verify their correspondence to coordinates. These adapters reject
``detail_confidence`` and ``producer_structure`` instead of assuming AF3 semantics.
``models.tsv`` records
``confidence_token_mapping_status=unavailable_unverified_token_mapping``.

The scalar schemas were checked against the official
`Protenix confidence producer <https://github.com/bytedance/Protenix/blob/main/protenix/model/sample_confidence.py>`_
and `RF3 output writer <https://github.com/RosettaCommons/foundry/blob/production/models/rf3/src/rf3/inference_engines/rf3.py>`_.
Synthetic tests verify scalar preservation and rejection of unsupported detail inputs;
they do not establish a completed prediction campaign or validated token mapping.
