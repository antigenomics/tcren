MD trajectory validation
========================

``tcren.validation.md`` reads recorded, PBC-corrected GROMACS output.
Coordinates use angstroms, timestamps use ps, whole-system potential uses
kJ/mol, and forces use kJ/mol/nm. These measurements do not estimate affinity.

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
