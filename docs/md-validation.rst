MD trajectory validation
========================

``tcren.validation.md`` reads recorded, PBC-corrected GROMACS output.
Coordinates use angstroms, timestamps use ps, whole-system potential uses
kJ/mol, and forces use kJ/mol/nm. These measurements do not estimate affinity.

Installed CLI
-------------

Use a new output directory for each export. For an already PBC-corrected,
continuous multi-model PDB, with explicit ``TITLE t=`` times in ps::

   tcren md inspect --trajectory protein-continuous.pdb --reference first.pdb \
     --potential potential.xvg --components components.json --mhc-class I \
     --pbc-whole --docking -t 1 -o md-inspection

``components.json`` maps nonempty, disjoint component names to lists of chain IDs,
for example ``{"TCR": ["A", "B"], "peptide": ["C"], "MHC": ["D"]}``.
For class II, pass ``--mhc-class II`` and include **both** groove chains in MHC
(for example ``["D", "E"]``). Class is declared by the caller and recorded;
it is not inferred from the number of chains. Exports never pool classes.
The name ``all`` is reserved for the whole-backbone fit. Each component fit
measures every component without a second fit of the measured partner.
The correspondence includes chain, sequence index, author residue number,
insertion code, atom name and alternate-location label; repeated identities fail.

Outputs are ``backbone_atoms.tsv`` (explicit atom/component correspondence),
``rmsd.tsv``, ``contact_counts.tsv``, ``contact_occupancy.tsv``,
``docking.tsv`` (header only when disabled), first/final/minimum-potential
PDBs, saved-frame times and minimum-frame JSON, and ``provenance.json``.
The latter records source/output hashes, units, class, component mapping,
CPU ceiling, searched annotation indexes when docking is enabled, and initial-image
check status. Failed inspection is marked
``failed`` rather than complete; existing directories are never overwritten.
No predictor confidence is inherited by evolved coordinates.

The ``--pbc-whole`` flag acknowledges input preparation; it does not perform
or prove unwrapping. Prepare whole molecules, a compact first-frame anchor and
continuous images with the GROMACS wrapper below. An absent unit cell is
explicitly reported as ``unavailable_no_unit_cell``. No repair is guessed.

The CPU ceiling accepts ``-t 1`` (serial), ``-t 0`` (all allocated) or an
explicit count within the allocation. Frame inspection and numerical kernels
are serial; only the optional single-reference annotation phase uses additional
cores, in one TCR query batch and one MHC batch, never a search per frame.
The public ``inspect_md_trajectory`` API enforces the same allocation and
native-kernel limits. No environment exports are required.

Atom-level proximity inspection uses the existing coordinate APIs::

   tcren md quality -s predicted.pdb --mhc-class II --cutoff 2.0 -o atom-quality

This writes ``coordinate_diagnostics.json`` and ``close_atom_pairs.tsv``,
including author residue/atom identities, alternate locations and possible
disulfides, with source/output hashes in ``provenance.json``. The strict cutoff
is in angstroms. Same-residue/adjacent same-chain exclusions are heuristic;
bond classification remains explicitly unknown. No clashscore is inferred.

For same-coordinate numerical reruns, including PME-order 4 versus 8::

   tcren md compare-forces --reference-energy order8-energy.xvg \
     --reference-forces order8-forces.xvg --energy order4-energy.xvg \
     --forces order4-forces.xvg --atoms 1000 --topology common.tpr \
     --coordinates same-frame.gro --same-system --mhc-class I -o force-comparison

Use the actual System-group atom count. XVGs contain one equal timestamp;
energies have two columns and forces have time plus three columns per atom.
``comparison.json`` preserves the existing API's energy difference and force
errors; provenance retains all four XVGs and declared common-input hashes.
``--same-system`` records a caller assertion, not proof of the XVGs' origin.
Preserve both original rerun commands, MDP settings and input checks separately.
These serial diagnostics do not run MD, establish PME convergence or estimate
binding energy. They do not change an active Stratos simulation.

Existing APIs and trajectory preparation
----------------------------------------

* ``component_relative_rmsd(reference, coordinates, components)`` fits each
  named component and measures every partner in that same frame.
* ``inspect_md_trajectory(trajectory_pdb, reference_pdb, components, energy_path)``
  streams multi-model PDB frames, verifies chain sequences and atom identities,
  and returns first, final and minimum-potential saved-frame coordinates.
  Optional ``contact_cutoff=5.0`` reports heavy-atom residue contacts, including
  zero-contact chain pairs and occupancy over saved frames. ``docking=True``
  annotates the reference once and reports frame-specific docking geometry.
* ``minimum_potential_frame`` requires matching saved-frame/energy timestamps;
  it neither interpolates coordinates nor selects an unsaved energy minimum.
* ``compare_same_frame_energy_forces`` compares identical-coordinate reruns.
  The caller must establish matching topology, coordinates and atom order.
* ``energy_drift`` reports NVE drift per ns and per atom without treating
  correlated frames as independent replicates.

Use ``python scripts/inspect_gromacs_trajectory.py --help`` for the GROMACS
export wrapper. It makes molecules whole, constructs a compact first-frame
anchor, then applies continuous unwrapping before analysis. This preserves
later relative motion. The wrapper requires an explicit GROMACS binary,
TPR, trajectory, two-column potential XVG, component-to-chain JSON mapping,
and a new output directory; command logs and source hashes accompany outputs.
The library rejects consecutive C-alpha gaps above 8 angstroms, but this guard
cannot detect every incorrect periodic image: input preparation remains required.

``tcren.validation.structures`` supplies sequence-verified component mapping,
coordinate diagnostics and explicit nonlocal close-atom pairs. Proximity pairs
are not a MolProbity clashscore. Original simulations and private campaign
orchestration are outside this library change.
