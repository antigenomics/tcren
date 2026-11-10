Complementarity browser
=======================

The rewritten ``tcren.complementarity`` module exports contact scenes and a local
dashboard without changing screening scores.

.. code-block:: console

   tcren complementarity -s models.list -o browser --task tcr -t 4
   python -m http.server --bind 127.0.0.1 --directory browser

Open ``http://127.0.0.1:8000``. Assets are local; compressed scenes load on
selection (at most 16 structures per overlay to bound browser memory). All input
structures remain individually available. Output must be a new directory. Inputs accept a structure, directory,
path manifest or tar archive. Archive identities retain member paths.
``--scratch`` selects an existing directory for one temporary archive member.
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
Overlay members inherit one reference origin without peptide recentering.

``--origin center-of-mass`` uses peptide heavy-atom element masses, selecting one
alternate atom per name (blank/A first). Hydrogen/deuterium are excluded.
``--origin apex`` uses the reference peptide C-alpha with greatest height above
the groove plane, sequence order breaking ties. Missing, collinear and nonfinite
landmarks fail explicitly. Scenes record reference, transform, origin, anchor
count and RMSD. Alignment is geometric normalization, not evidence of binding.

Nodes and contacts
------------------

Scenes include CDR1–3 on both receptor chains, peptide and annotated MHC groove
helices/floor. Missing annotations are reported. Residue positions are C-alpha;
author numbers, insertion codes and sequence indices remain separate. Missing
C-alpha residues retain null coordinates and a warning. Heavy-atom nodes retain
alternate locations and distinguish repeated names.

Contacts are all inter-chain heavy-atom pairs within inclusive ``--cutoff``
(default 5 Angstroms). Different nonblank alternate conformers are excluded.
Residue multiplicity counts these witnesses, not residue pairs or energy.
Default displayed groups are TCR–peptide, TCR–MHC and peptide–MHC; alpha–beta
contacts are optional. Unit edges, width proportional to square-root count
(capped at eight pixels), and parallel witness edges are available.
Atom mode displays individual witnesses.

Solid residue traces require consecutive resolved sequence indices and a C-alpha
distance <=4.5 Angstroms. Atom backbone traces use N–CA, CA–C, C–O/OXT and adjacent
C–N distances in (0.4, 1.9] Angstroms with compatible alternate locations.
Within-residue atom links use distances in (0.4, 1.2 times summed covalent radii]
for C/N/O/S/P/Se, with compatible alternate conformers; they are geometric bond
suggestions, not verified bond orders. Backbone links are highlighted separately.
These traces are not a complete chemical bond graph. Dashed edges are contacts.
Hover shows identities; region/edge selection, zoom and label scope control density.

Legacy and Python APIs
----------------------

``--legacy`` invokes unchanged ``project2d``/``viz.svg2d`` output, writing one SVG
per structure. New ranking frames and dashboard do not apply. Historical cutoff
and projection restrictions still apply. Existing APIs remain unchanged.

``canonical_frame`` constructs a reference; ``fit_frame`` fits matched landmarks.
``prepare.reference_frame`` and ``prepare.ranking_frame`` accept annotated
structures. ``scene.build_scene`` creates JSON-compatible scenes with schema
``tcren.complementarity/1``. ``dashboard.write_dashboard`` streams compressed
scenes and the viewer into a new directory.
