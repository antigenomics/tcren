"""Read recorded GROMACS validation output; short trajectories do not establish affinity."""

import re
from pathlib import Path

import numpy as np


def read_xvg(path, *, columns=None):
    """Read one numeric XVG dataset, requiring finite values and strictly increasing time.

    Header/legend lines are skipped; multi-dataset files are rejected. Units remain the units
    recorded by GROMACS (ps for these campaign trajectories). No implicit time conversion.
    """
    rows = [
        line.strip()
        for line in Path(path).read_text().splitlines()
        if line.strip() and not line.lstrip().startswith(("#", "@"))
    ]
    values = np.loadtxt(rows, ndmin=2)
    if not len(values) or values.shape[1] < 2 or not np.isfinite(values).all():
        raise ValueError(f"Empty/nonfinite XVG dataset: {path}")
    if columns is not None and values.shape[1] != columns:
        raise ValueError(f"Expected {columns} XVG columns: {path}")
    if (np.diff(values[:, 0]) <= 0).any():
        raise ValueError(f"XVG time must increase strictly: {path}")
    return values


def plot_validation(directories, labels, output):
    """Plot whole-system potential, temperature, pressure, density and backbone RMSD."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    if len(directories) != len(labels) or not directories:
        raise ValueError("one label required per MD directory")
    plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(
        len(directories),
        5,
        figsize=(16, 2.7 * len(directories)),
        layout="constrained",
        squeeze=False,
    )
    columns = [
        (1, "Potential energy", "kJ/mol", None),
        (4, "Temperature", "K", 300),
        (5, "Pressure", "bar", 1),
        (6, "Density", "kg/m³", None),
    ]
    for row, (directory, label) in enumerate(zip(directories, labels)):
        energy = read_xvg(Path(directory) / "energy.xvg", columns=7)
        rmsd = read_xvg(Path(directory) / "backbone-rmsd.xvg", columns=2)
        color = ["#4477AA", "#AA3377", "#228833"][row % 3]
        for column, (index, title, unit, target) in enumerate(columns):
            ax = axes[row, column]
            ax.plot(energy[:, 0], energy[:, index], color=color, lw=1)
            if target is not None:
                ax.axhline(target, color="#555555", ls="--", lw=0.8)
            ax.set_ylabel(unit)
            ax.grid(alpha=0.2)
            if row == 0:
                ax.set_title(title)
        ax = axes[row, 4]
        ax.plot(rmsd[:, 0], rmsd[:, 1], color=color, lw=1)
        ax.set_ylabel("nm")
        ax.grid(alpha=0.2)
        if row == 0:
            ax.set_title("Backbone RMSD")
        axes[row, 0].set_title("")
        axes[row, 0].set_title(f"{label}\nPotential energy", loc="left", fontweight="bold")
        for ax in axes[row]:
            ax.set_xlabel("Production time (ps)")
    fig.suptitle(
        "Short MD execution validation: 20 ps production per system\n"
        "Energy is whole-system potential; RMSD uses each production-start reference; "
        "short runs do not establish equilibration",
        fontsize=12,
    )
    output = Path(output)
    fig.savefig(output.with_suffix(".png"), dpi=180)
    fig.savefig(output.with_suffix(".pdf"))
    plt.close(fig)


def component_relative_rmsd(reference, coordinates, components):
    """Backbone RMSD in Å after component fits; coordinates must be made PBC-whole.

    Atom correspondence and component indices are explicit, identical for both arrays.
    No independent fit of the measured component is applied after fitting its partner.
    """
    from ..docking._transform import kabsch

    reference = np.asarray(reference, dtype=float)
    coordinates = np.asarray(coordinates, dtype=float)
    if reference.shape != coordinates.shape or reference.ndim != 2 or reference.shape[1] != 3:
        raise ValueError("Expected corresponding N × 3 coordinate arrays")
    if not np.isfinite(reference).all() or not np.isfinite(coordinates).all():
        raise ValueError("Coordinates must be finite")
    indices = {name: np.asarray(values, dtype=int) for name, values in components.items()}
    for values in indices.values():
        if not len(values) or len(np.unique(values)) != len(values):
            raise ValueError("Empty/repeated component atom indices")
        if values.min() < 0 or values.max() >= len(reference):
            raise ValueError("Component atom index out of range")
    rows = []
    fits = {"all": np.arange(len(reference)), **indices}
    for fit, selection in fits.items():
        rotation, translation, _ = kabsch(coordinates[selection], reference[selection])
        aligned = coordinates @ rotation + translation
        for component, measured in indices.items():
            value = np.sqrt(np.mean(np.sum((aligned[measured] - reference[measured]) ** 2, axis=1)))
            rows.append(
                {
                    "fit_component": fit,
                    "measured_component": component,
                    "atoms": len(measured),
                    "rmsd_A": float(value),
                }
            )
    return rows


def minimum_potential_frame(energy_path, frame_times_ps, *, tolerance_ps=1e-4):
    """Select the saved frame with lowest matching whole-system potential (kJ/mol).

    The energy XVG must have time and potential columns only. No interpolation,
    implicit units, or binding-energy interpretation is permitted.
    """
    energy = read_xvg(energy_path, columns=2)
    times = np.asarray(frame_times_ps, dtype=float)
    if times.ndim != 1 or not len(times) or not np.isfinite(times).all():
        raise ValueError("Expected finite saved frame times in ps")
    if (np.diff(times) <= 0).any() or tolerance_ps < 0:
        raise ValueError("Frame times must increase and tolerance must be nonnegative")
    right = np.searchsorted(energy[:, 0], times)
    right = np.minimum(right, len(energy) - 1)
    left = np.maximum(right - 1, 0)
    selected = np.where(abs(energy[left, 0] - times) <= abs(energy[right, 0] - times), left, right)
    if (abs(energy[selected, 0] - times) > tolerance_ps).any():
        raise ValueError("One or more saved frames has no matching potential-energy timestamp")
    frame = int(np.argmin(energy[selected, 1]))
    return {
        "frame_index": frame,
        "time_ps": float(times[frame]),
        "potential_kJ_mol": float(energy[selected[frame], 1]),
        "selection_scope": "saved frames with exact energy-time correspondence",
    }


def _iter_pdb_frames(path):
    """Stream GROMACS multi-model PDB frames with explicit TITLE time in ps."""
    lines, time = [], None
    with Path(path).open() as handle:
        for line in handle:
            if line.startswith("TITLE"):
                match = re.search(r"\bt=\s*([-+\d.eE]+)", line)
                if match:
                    time = float(match.group(1))
            if line.startswith("ENDMDL"):
                if time is None or not lines:
                    raise ValueError("Trajectory frame lacks explicit TITLE time/coordinates")
                yield time, "".join(lines) + "END\n"
                lines, time = [], None
            elif line.startswith(("ATOM  ", "HETATM", "TER   ", "CRYST1")):
                lines.append(line)
    if lines:
        raise ValueError("Trajectory ends without ENDMDL")


def inspect_md_trajectory(
    trajectory_pdb,
    reference_pdb,
    components,
    energy_path,
    *,
    contact_cutoff=None,
    docking=False,
    organism="human",
    threads=1,
):
    """Stream PBC-corrected GROMACS frames, reporting partner-relative RMSD.

    Components explicitly map names to chain IDs. Exact backbone atom keys and
    polymer sequences must match the reference, preventing silent atom-order
    correspondence. GROMACS must make coordinates PBC-whole before this call;
    this reader neither unwraps coordinates nor claims a binding-energy estimate.
    Returns RMSD rows, whole-system minimum-potential saved-frame metadata and
    that frame's PDB text. Coordinates are Å and times are ps. Optional heavy-atom
    residue contacts report occupancy over saved frames, including zero-contact
    chain-pair counts. Optional docking readouts annotate the reference once and
    transfer sequence-validated region markup onto each frame's own residues.
    """
    from ..structure import parse_structure
    from .structures import parse_pdb_text

    reference = parse_structure(reference_pdb, keep_hydrogens=False)
    if threads < 1 or (
        contact_cutoff is not None and (not np.isfinite(contact_cutoff) or contact_cutoff <= 0)
    ):
        raise ValueError("Positive threads and finite positive contact cutoff required")
    if docking:
        from ..annotation import classify_chains
        from ..annotation.batch import annotate_batch
        from ..mhc import annotate_mhc_batch

        records = annotate_batch([reference], organisms=(organism,), threads=threads)
        classify_chains(
            reference, organism=organism, autodetect_species=False, precomputed_records=records[0]
        )
        annotate_mhc_batch([reference], threads=threads)

    def backbone(structure):
        for chain in structure.chains:
            ca = np.asarray([r.ca for r in chain.residues if r.ca is not None])
            if len(ca) > 1 and (np.linalg.norm(np.diff(ca, axis=0), axis=1) > 8).any():
                raise ValueError(
                    f"Broken backbone/PBC in chain {chain.chain_id}: consecutive CA gap >8 A"
                )
        return {
            (chain.chain_id, residue.seq_index, atom.name): atom.coord
            for chain in structure.chains
            for residue in chain.residues
            for atom in residue.atoms
            if atom.name in {"N", "CA", "C", "O"}
        }

    expected_sequences = {c.chain_id: c.sequence() for c in reference.chains}
    heavy_keys = {
        (c.chain_id, r.seq_index, a.name, a.element)
        for c in reference.chains
        for r in c.residues
        for a in r.atoms
    }
    atoms = backbone(reference)
    keys = sorted(atoms)
    xyz = np.asarray([atoms[k] for k in keys])
    if not keys:
        raise ValueError("Reference has no backbone atoms")
    selections = {}
    for name, chains in components.items():
        if not chains or not set(chains) <= expected_sequences.keys():
            raise ValueError(f"Unknown/empty component chains: {name}")
        selections[name] = [i for i, key in enumerate(keys) if key[0] in chains]
    times, rows, counts, docking_rows = [], [], [], []
    first_text, last_text = None, None
    occupancy = {}
    for frame_index, (time, text) in enumerate(_iter_pdb_frames(trajectory_pdb)):
        if first_text is None:
            first_text = text
        last_text = text
        structure = parse_pdb_text(text, source=f"frame-{frame_index}")
        if {c.chain_id: c.sequence() for c in structure.chains} != expected_sequences:
            raise ValueError(f"Frame {frame_index}: chain sequences differ from reference")
        coordinates = backbone(structure)
        if coordinates.keys() != atoms.keys():
            raise ValueError(f"Frame {frame_index}: backbone atom identities differ from reference")
        times.append(time)
        for row in component_relative_rmsd(xyz, [coordinates[k] for k in keys], selections):
            rows.append({"frame_index": frame_index, "time_ps": time, **row})
        if contact_cutoff is not None:
            from itertools import combinations

            from ..contacts.geometry import all_atom_contacts

            if {
                (c.chain_id, r.seq_index, a.name, a.element)
                for c in structure.chains
                for r in c.residues
                for a in r.atoms
            } != heavy_keys:
                raise ValueError(
                    f"Frame {frame_index}: heavy atom identities differ from reference"
                )
            contacts = all_atom_contacts(structure, cutoff=contact_cutoff)
            pair_counts = {}
            for row in contacts.iter_rows(named=True):
                pair = (row["chain.id.from"], row["chain.id.to"])
                pair_counts[pair] = pair_counts.get(pair, 0) + 1
                key = (pair[0], row["residue.index.from"], pair[1], row["residue.index.to"])
                entry = occupancy.setdefault(
                    key,
                    {
                        "frames_contacting": 0,
                        "minimum_distance_A": row["dist"],
                        "maximum_contact_distance_A": row["dist"],
                    },
                )
                entry["frames_contacting"] += 1
                entry["minimum_distance_A"] = min(entry["minimum_distance_A"], row["dist"])
                entry["maximum_contact_distance_A"] = max(
                    entry["maximum_contact_distance_A"], row["dist"]
                )
            counts.extend(
                {
                    "frame_index": frame_index,
                    "time_ps": time,
                    "chain.id.from": a,
                    "chain.id.to": b,
                    "cutoff_A": float(contact_cutoff),
                    "residue_pairs": pair_counts.get((a, b), 0),
                }
                for a, b in combinations(sorted(expected_sequences), 2)
            )
        if docking:
            from dataclasses import asdict, replace

            from ..docking import docking_angles, tcr_placement

            for chain in structure.chains:
                source = reference.chain(chain.chain_id)
                chain.chain_type = source.chain_type
                chain.chain_supertype = source.chain_supertype
                chain.allele_info = source.allele_info
                chain.regions = [
                    replace(
                        region, residues=[chain.by_seq_index(r.seq_index) for r in region.residues]
                    )
                    for region in source.regions
                ]
            structure.complex_species = reference.complex_species
            structure.cell_type = reference.cell_type
            docking_rows.append(
                {
                    "frame_index": frame_index,
                    "time_ps": time,
                    **asdict(docking_angles(structure)),
                    **asdict(tcr_placement(structure)),
                }
            )
    minimum = minimum_potential_frame(energy_path, times)
    selected_text = next(
        text
        for i, (_, text) in enumerate(_iter_pdb_frames(trajectory_pdb))
        if i == minimum["frame_index"]
    )
    contacts = [
        {
            "chain.id.from": a,
            "residue.index.from": i,
            "chain.id.to": b,
            "residue.index.to": j,
            "cutoff_A": float(contact_cutoff),
            **entry,
            "frames": len(times),
            "saved_frame_occupancy": entry["frames_contacting"] / len(times),
        }
        for (a, i, b, j), entry in sorted(occupancy.items())
    ]
    return {
        "rmsd": rows,
        "minimum_potential": minimum,
        "selected_frame_pdb": selected_text,
        "first_frame_pdb": first_text,
        "last_frame_pdb": last_text,
        "saved_frame_times_ps": times,
        "contact_counts": counts,
        "contact_occupancy": contacts,
        "docking": docking_rows,
    }


def compare_same_frame_energy_forces(reference_energy, reference_forces, energy, forces, *, atoms):
    """Numerical fidelity from identical-coordinate GROMACS reruns.

    Each XVG contains one explicit timestamp; forces are three columns per atom
    in the same System-group order (kJ/mol/nm). Energies are kJ/mol. The caller
    must verify coordinate/topology identity; no structural or affinity judgment
    is inferred from a force error or total-system energy difference.
    """
    if not isinstance(atoms, int) or atoms < 1:
        raise ValueError("Positive integer atom count required")
    reference_e = read_xvg(reference_energy, columns=2)
    candidate_e = read_xvg(energy, columns=2)
    reference_f = read_xvg(reference_forces, columns=1 + 3 * atoms)
    candidate_f = read_xvg(forces, columns=1 + 3 * atoms)
    arrays = [reference_e, candidate_e, reference_f, candidate_f]
    if any(len(a) != 1 for a in arrays) or any(a[0, 0] != arrays[0][0, 0] for a in arrays):
        raise ValueError("Require one identical explicit energy/force timestamp")
    baseline = reference_f[0, 1:].reshape(atoms, 3)
    delta = candidate_f[0, 1:].reshape(atoms, 3) - baseline
    distances = np.linalg.norm(delta, axis=1)
    rms_reference = float(np.sqrt(np.mean(np.sum(baseline**2, axis=1))))
    rms_error = float(np.sqrt(np.mean(distances**2)))
    return {
        "time_ps": float(reference_e[0, 0]),
        "atoms": atoms,
        "potential_difference_kJ_mol": float(candidate_e[0, 1] - reference_e[0, 1]),
        "reference_force_rms_kJ_mol_nm": rms_reference,
        "force_error_rms_kJ_mol_nm": rms_error,
        "force_error_relative_rms": rms_error / rms_reference if rms_reference else None,
        "force_error_p95_kJ_mol_nm": float(np.quantile(distances, 0.95)),
        "force_error_max_kJ_mol_nm": float(distances.max()),
    }


def energy_drift(energy_path, *, atoms, discard_ps=0.0):
    """Fit NVE total-energy drift in kJ/mol/ns, normalized by atom count.

    No normalization by an arbitrary potential zero, independent-replicate
    confidence interval, or automatic timestep acceptance is implied.
    """
    if not isinstance(atoms, int) or atoms < 1 or not np.isfinite(discard_ps) or discard_ps < 0:
        raise ValueError("Positive atom count and finite nonnegative discard time required")
    values = read_xvg(energy_path, columns=2)
    values = values[values[:, 0] >= discard_ps]
    if len(values) < 3 or values[-1, 0] == values[0, 0]:
        raise ValueError("Need at least three distinct energy timestamps")
    centered = values[:, 0] - values[:, 0].mean()
    slope = float(centered @ (values[:, 1] - values[:, 1].mean()) / (centered @ centered))
    return {
        "atoms": atoms,
        "samples": len(values),
        "start_ps": float(values[0, 0]),
        "end_ps": float(values[-1, 0]),
        "drift_kJ_mol_ns": 1000 * slope,
        "drift_per_atom_kJ_mol_ns": 1000 * slope / atoms,
        "energy_std_kJ_mol": float(values[:, 1].std()),
        "energy_range_kJ_mol": float(np.ptp(values[:, 1])),
        "end_minus_start_kJ_mol": float(values[-1, 1] - values[0, 1]),
    }
