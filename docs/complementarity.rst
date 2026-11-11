Complementarity browser
=======================

The rewritten ``tcren.complementarity`` module owns one molecular graph and a
standalone local dashboard, independent of vdjdb-web. It changes no screening
scores. Atoms own coordinates; residue points reference their C-alpha atoms.

.. code-block:: console

   tcren complementarity -s models.list -o browser --task tcr -t 4
   python -m http.server --bind 127.0.0.1 --directory browser

Open ``http://127.0.0.1:8000``. Assets are local; compressed scenes load on
selection (at most 16 structures per overlay to bound browser memory). All input
structures remain individually available. Output must be a new directory. Inputs accept a structure, directory,
path manifest or tar archive. Archive identities retain member paths.
``--scratch`` selects an existing directory for one temporary archive member.
Repeat ``--exclude-id`` to omit explicitly named inputs; omissions and input
hashes are recorded in ``excluded.json``. Unknown exclusion IDs fail explicitly.
Parsing/annotation errors stop the run. Geometric errors are recorded in
``errors.json`` and cause a nonzero exit after remaining scenes export.

Ranking frames
--------------

``--task tcr`` groups identical peptide/MHC sequences and aligns annotated MHC
groove C-alpha anchors to the first encountered reference. ``--task epitope``
groups identical alpha/beta sequences and aligns framework C-alpha anchors.
Different constructs receive separate groups. Use ordered manifests for stable
reference selection. Coordinates and RMSD are in Angstroms.

The reference peptide's N-to-C terminal vector defines +x. The groove-floor plane
normal toward the peptide, orthogonalized against x, defines +z. The right-handed
y axis points toward the side camera. Side views show x/z, top views x/y.
Overlay members inherit one reference origin without independent recentering.

For TCR ranking, ``--origin center-of-mass`` uses peptide heavy-atom element masses, selecting one
alternate atom per name (blank/A first). Hydrogen/deuterium are excluded.
``--origin apex`` uses the reference peptide C-alpha with greatest height above
the groove plane, sequence order breaking ties. Missing, collinear and nonfinite
landmarks fail explicitly. For epitope ranking the origin instead uses the fixed
TCR variable-domain heavy-atom COM, or its groove-facing CDR C-alpha apex (lowest
reference z). Framework and CDR atoms define that variable domain; constant
domains do not determine its center. Scenes record the origin component explicitly.
Scenes record reference, transform, origin, anchor
count and RMSD. Alignment is geometric normalization, not evidence of binding.
Matched anchors also retain their chain role and resolved-sequence ordinal.
Execution retains compact reference landmarks rather than whole complexes for
every group; completed annotation batches are released before loading the next.

Nodes and contacts
------------------

Scenes retain all typed alpha/beta, peptide, MHC and beta-2-microglobulin polymer
atoms. CDR1–3, peptide and groove helices/floor are initially visible; chains,
frameworks and other regions remain independently selectable. Missing annotations
are reported. Residue positions resolve the owning C-alpha atom;
author numbers, insertion codes and sequence indices remain separate. Missing
C-alpha residues retain null coordinates and a warning. Heavy-atom nodes retain
alternate locations and distinguish repeated names.

Contacts are all inter-chain heavy-atom pairs within inclusive ``--cutoff``
(default 5 Angstroms). Different nonblank alternate conformers are excluded.
Directly covalent atom pairs are excluded. Residue multiplicity aggregates the
same atom table, counting witnesses rather than residue pairs or energy.
Default displayed groups are TCR–peptide, TCR–MHC and peptide–MHC; alpha–beta
contacts are optional. Unit edges, linear width (0.5 pixels per atom-pair contact)
and parallel witness curves sharing the same endpoints are available.
Atom mode displays individual witnesses.

Named amino-acid templates define covalent connectivity, with coordinate
validation: backbone distances in (0.4, 1.9] Angstroms, side chains in (0.4, 2.3].
Adjacent peptide C–N bonds require consecutive resolved sequence positions.
MSE replaces methionine SD with SE. Cysteine SG–SG at 1.7–2.3 Angstroms is an
explicit disulfide, including between chains. Unknown side-chain templates and
unavailable named bonds are flagged. Atoms outside the canonical template are
retained with an explicit untemplated-atom warning rather than invented bonds.
Modified-residue connectivity may therefore be incomplete. Bond orders are not assigned.
Missing backbone atoms are flagged; reported covalent degrees describe the
observed, modelled connectivity rather than unseen chemical neighbours. Residue
backbone traces derive from peptide bonds whose residues both have C-alpha atoms.
Alternate conformers remain distinguishable and are not simultaneous occupancy
evidence. Dashed edges are non-covalent inter-chain distance contacts.

Global atom counts include covalent degree, contact count/degree, partner residue
and partner chain counts. Residues include atom membership, contacting atoms,
atom-pair contact count, distinct residue-partner degree, covalent residue degree
and partner-chain count. Zeros remain explicit. The inspector separately labels
counts after chain, region and contact-group filters; zoom and pan do not change
these counts. Filters never overwrite global counts.
Hover/click exposes identities, source keys, ordinals and counts. Summary tables,
chain/region/edge selection, zoom, pan and label scope support inspection.
Side z and top y camera offsets are independent. ``Fit visible nodes`` adjusts
the shared zoom and camera offsets; ``Reset camera`` restores the reference-origin
view. These actions never recenter molecular coordinates or change counts.
Default peptide labels refer to the first selected structure, avoiding repeated
labels in overlays; select ``Peptide · every structure`` to label each member.
All points remain independently inspectable in either mode.
Each projection displays at most 250 labels, prioritizing peptide nodes and then
nodes with contacts after filters. A visible notice reports omitted labels.
The label limit never removes molecular points, edges, table rows or counts.

Compressed data
---------------

``tables/`` contains typed Zstandard-compressed Parquet and matching TSV.gz
partitions for structures, atoms, residues, bonds, atom contacts and aggregated
residue contacts. Partitions hold
at most eight structures. IDs are structure-local; join on ``structure_id`` plus
the local ID. Residue coordinates resolve ``ca_atom`` in the atom table.
``tables/manifest.json`` records schemas and row counts. JSON.gz is the browser
interchange serialization of the same model; the exporter performs no second
contact or degree calculation. Small JSON files retain input/reference manifests,
errors and run parameters. All runtime assets are local; no CDN is required.

``summary.parquet``, ``summary.tsv.gz`` and ``summary.json`` contain per-structure
counts and export totals. Run ``tcren complementarity-summary browser`` to check
table row counts and contact-incidence sums and obtain a compact JSON report
without annotation or contact recomputation. This verifies serialization and
count consistency; it does not establish biological activity or pose correctness.

Legacy and Python APIs
----------------------

``--legacy`` invokes unchanged ``project2d``/``viz.svg2d`` output, writing one SVG
per structure. New ranking frames and dashboard do not apply. Historical cutoff
and projection restrictions still apply. Existing APIs remain unchanged.

``canonical_frame`` constructs a reference; ``fit_frame`` fits matched landmarks.
``prepare.reference_frame`` and ``prepare.ranking_frame`` accept annotated
structures. ``scene.build_scene`` creates JSON-compatible scenes with schema
``tcren.complementarity/2``. ``dashboard.write_dashboard`` streams compressed
scenes and the viewer into a new directory.

Project-owned extensions
------------------------

Dataset-specific PWM thresholds, experimental matrices and motif/anchor joins
belong in consuming projects. They are not computed by this CLI. The browser
dispatches ``tcren:render`` with ``{scenes, atomic}`` after drawing; plotted points
carry ``data-structure-id`` and ``data-residue-id``. JavaScript consumers can reuse
``loadScene``, ``nodesFor``, ``choices``, ``updateStructureOptions``,
``fitCameraCoordinates`` and ``resetCameraCoordinates`` from
``viewer.js``. This permits local
overlays without duplicating molecular coordinate or contact calculations.
An HTML ``data-max-overlay`` attribute can lower the detailed selection limit;
the library always caps it at 16. Large overview rendering belongs in consumers.
Consumers can add ``display_label`` to index entries. Long structure identifiers
are abbreviated in the selection summary, with full identities in tooltips,
inspectors and all exports.

Metadata and table-owned surface/contact views
----------------------------------------------

Use an explicit metadata identity and path rather than inferring identity from a
filename. TSV, compressed TSV, CSV and Parquet inputs are supported; relative
paths resolve against the metadata file. IDs must be unique and nonempty.
An optional archive-member column points to individual structure members without
unpacking an entire archive. Outcome columns are never used to compute geometry.

.. code-block:: bash

   tcren complementarity --metadata structures.tsv --id-column id \
       --path-column structure_path -o browser --task tcr -t 4
   tcren complementarity-summary browser
   tcren complementarity-views --exports browser -o views --grid 64
   tcren complementarity-plot --tables views --structure-id native-1 \
       -o figures/native-1 --contact-mode pairs
   python -m http.server 8000 --directory views

``--structures`` and ``--metadata`` are mutually exclusive. For archived paths,
add ``--member-column member``. Use ``--task epitope`` for peptide topography in
the TCR-fixed frame. Existing verified browser exports can be used directly as
the input to ``complementarity-views``; no reannotation or contact detection is
performed by that command.

The views export has three typed long tables, partitioned by structure, in
matching Zstandard Parquet and TSV.gz. Composite identifiers always include
``structure_id``. The schemas and SHA256 checksums are in ``manifest.json``.

* ``topography_residues``: variable-partner residue identity, author numbering,
  region, C-alpha coordinates, atom-center z extrema, facing atom apex and KD
  hydropathy. C-alpha may be null; residue IDs join the original graph tables.
* ``topography_cells``: every fixed-grid x/y cell, occupied flag, exposed surface
  z, owner atom/residue and chemistry. Empty cells and unknown hydropathy stay null.
* ``contact_heatmap``: dense residue pairs for six TCR loop facets and three
  partner facets, C-alpha distance, resolution flag and original noncovalent
  atom-contact multiplicity. Missing C-alpha is null, never zero distance.

The surface is the unsmoothed sphere envelope with Bondi radii and a 1.4-Angstrom
probe by default. TCR ranking uses the lower face of alpha/beta variable domains;
epitope ranking uses the peptide upper face. Both retain the owning scene's
canonical frame. Heights are signed canonical z, not bulge above an independently
fitted groove. Grid ordering is ``(ny,nx)``; ``--extent x0 x1 y0 y1`` defaults to
an 80-by-80-Angstrom window. A finite window can omit peripheral atoms.
Summary fields ``occupied_boundary_cells`` and ``surface_may_be_clipped`` flag
occupied edge cells conservatively; inspect flagged extents before comparing
surface descriptors. The surface uses the union of exported heavy-atom spheres,
including alternative conformers; it is not an occupancy-weighted conformer model.

Facet columns are alpha CDR1/2/3 followed by beta CDR3/2/1; within every facet,
residues retain their natural N-to-C order. Partner rows are peptide and the
two groove helices: class-I alpha1/alpha2 or class-II alpha1/beta1. Empty facets
are recorded in the index. Contacting-pair masking uses existing heavy-atom
witnesses, not a C-alpha cutoff. ``--metric atom_contact_count`` instead shows
the witness count. Repeat ``--hide-row`` or ``--hide-column`` to hide facets.

Static PNGs read Parquet. The interactive table browser loads one selected
structure. PNG rendering requires the optional ``tcren[viz]`` dependencies.
The browser loads the selected
structure through JSON caches generated by reading those same written Parquet
partitions; cache metadata includes source paths and SHA256 hashes. Neither
consumer calculates heights, distances or contacts. ``summary.parquet``,
``summary.tsv.gz`` and ``summary.json`` contain descriptive counts and
surface ranges computed with Polars from the written tables. Mean hydropathy
is an occupied-cell mean over known KD values, not a residue mean. Comparisons
must retain canonical reference groups, endpoint definitions and uncertainty;
these tables alone do not establish binding or predictive performance.

Heatmap inspection
------------------

Contact-only filtering can hide every pair in a facet or structure. Empty
facets explicitly report no contacting pairs; unavailable values and absent
facets remain distinct. Switch to all-pair C-alpha distances to inspect the
resolved geometry. No contacts alone does not establish an experimental
nonbinder. Short display labels preserve full structure identities for joins.
Static figures and browser values remain consumers of exported Parquet tables.
Inspect occupied-grid-boundary warnings before interpreting whole surfaces.
