import numpy as np
import pytest

from tcren.validation.md import component_relative_rmsd, minimum_potential_frame


def test_contact_report_cutoff_matches_actual_boundary(tmp_path):
    from tcren.validation.structures import geometry_and_contacts

    pdb = tmp_path / "boundary.pdb"
    pdb.write_text(
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C\n"
        "ATOM      2  CA  GLY B   1       4.800   0.000   0.000  1.00 20.00           C\nEND\n"
    )
    identity = [{"role": "left", "output_chain": "A"}, {"role": "right", "output_chain": "B"}]
    _, strict = geometry_and_contacts(pdb, identity)
    _, wider = geometry_and_contacts(pdb, identity, cutoff=5.0)
    assert strict[0]["residue_pairs_le4_5A"] == 0
    assert wider[0]["residue_pairs"] == 1
    assert wider[0]["cutoff_A"] == 5.0
    assert "residue_pairs_le4_5A" not in wider[0]


def test_extended_sequence_comparison_is_explicit_and_reports_span(monkeypatch):
    from tcren.validation import structures

    coords = np.arange(18).reshape(6, 3)
    monkeypatch.setattr(structures, "read_chains", lambda _: [("D", "AACDEF", coords)])
    with pytest.raises(ValueError, match="sequence match"):
        structures.map_sequence_roles("control", ["MHC"], {"MHC": "CDE"})
    mapped, identity = structures.map_sequence_roles(
        "control", ["MHC"], {"MHC": "CDE"}, allow_extended=True
    )
    np.testing.assert_array_equal(mapped["MHC"][0], coords[2:5])
    assert identity[0]["output_start"] == 2
    assert identity[0]["reference_start"] == 0


def test_streamed_md_keeps_atom_identity_and_exact_energy_times(tmp_path):
    from tcren.validation.md import inspect_md_trajectory

    atoms = []
    for i, (chain, xyz) in enumerate(
        [("A", (0, 0, 0)), ("A", (1, 0, 0)), ("B", (0, 1, 0)), ("B", (0, 1, 1))], 1
    ):
        atoms.append(
            f"ATOM  {i:5d}  CA  ALA {chain}{(i - 1) % 2 + 1:4d}    {xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00 20.00           C\n"
        )
    text = "".join(atoms)
    reference = tmp_path / "ref.pdb"
    reference.write_text(text + "END\n")
    trajectory = tmp_path / "trajectory.pdb"
    trajectory.write_text(
        "TITLE     t= 0.000\nMODEL        1\n"
        + text
        + "ENDMDL\n"
        + "TITLE     t= 10.000\nMODEL        2\n"
        + text
        + "ENDMDL\n"
    )
    energy = tmp_path / "potential.xvg"
    energy.write_text("0 -10\n10 -20\n")
    result = inspect_md_trajectory(trajectory, reference, {"TCR": ["A"], "MHC": ["B"]}, energy)
    assert result["minimum_potential"]["frame_index"] == 1
    assert result["saved_frame_times_ps"] == [0.0, 10.0]
    assert result["first_frame_pdb"] == result["last_frame_pdb"]
    assert len(result["rmsd"]) == 12
    assert max(row["rmsd_A"] for row in result["rmsd"]) < 1e-10
    displaced = "".join(
        line[:38] + f"{11.0:8.3f}" + line[46:] if line[21] == "B" else line for line in atoms
    )
    trajectory.write_text(
        "TITLE     t= 0.000\n"
        + text
        + "ENDMDL\n"
        + "TITLE     t= 10.000\n"
        + displaced
        + "ENDMDL\n"
    )
    result = inspect_md_trajectory(
        trajectory, reference, {"TCR": ["A"], "MHC": ["B"]}, energy, contact_cutoff=1.1
    )
    assert [r["residue_pairs"] for r in result["contact_counts"]] == [1, 0]
    assert len(result["contact_occupancy"]) == 1
    assert result["contact_occupancy"][0]["saved_frame_occupancy"] == 0.5
    trajectory.write_text(trajectory.read_text().replace(" CA  ALA B", " N   ALA B"))
    with pytest.raises(ValueError, match="atom identities"):
        inspect_md_trajectory(trajectory, reference, {"TCR": ["A"]}, energy)
    trajectory.write_text(
        "TITLE     t= 0.000\n"
        + text.replace("   1.000   0.000   0.000", " 100.000   0.000   0.000")
        + "ENDMDL\n"
    )
    with pytest.raises(ValueError, match="Broken backbone/PBC"):
        inspect_md_trajectory(trajectory, reference, {"TCR": ["A"]}, energy)


def test_same_frame_force_comparison_validates_units_dimensions_and_time(tmp_path):
    from tcren.validation.md import compare_same_frame_energy_forces

    paths = [
        tmp_path / name
        for name in ["reference-energy.xvg", "reference-force.xvg", "energy.xvg", "force.xvg"]
    ]
    for path, text in zip(paths, ["0 -100\n", "0 3 4 0\n", "0 -99\n", "0 3 4 1\n"]):
        path.write_text(text)
    result = compare_same_frame_energy_forces(*paths, atoms=1)
    assert result["potential_difference_kJ_mol"] == 1
    assert result["force_error_relative_rms"] == pytest.approx(0.2)
    paths[-1].write_text("1 3 4 1\n")
    with pytest.raises(ValueError, match="timestamp"):
        compare_same_frame_energy_forces(*paths, atoms=1)


def test_energy_drift_reports_ps_to_ns_and_per_atom_units(tmp_path):
    from tcren.validation.md import energy_drift

    p = tmp_path / "energy.xvg"
    p.write_text("0 10\n1 12\n2 14\n3 16\n")
    result = energy_drift(p, atoms=2)
    assert result["drift_kJ_mol_ns"] == pytest.approx(2000)
    assert result["drift_per_atom_kJ_mol_ns"] == pytest.approx(1000)
    with pytest.raises(ValueError, match="timestamps"):
        energy_drift(p, atoms=2, discard_ps=2)


def test_component_fit_retains_partner_displacement():
    xyz = np.random.default_rng(42).normal(size=(16, 3))
    moved = xyz.copy()
    moved[:8] += [5, 0, 0]
    rows = component_relative_rmsd(xyz, moved, {"TCR": range(8), "MHC": range(8, 16)})
    row = next(r for r in rows if r["fit_component"] == "MHC" and r["measured_component"] == "TCR")
    assert row["rmsd_A"] == pytest.approx(5)
    own = next(r for r in rows if r["fit_component"] == "TCR" and r["measured_component"] == "TCR")
    assert own["rmsd_A"] == pytest.approx(0, abs=1e-12)


def test_minimum_frame_uses_saved_frames_and_exact_times(tmp_path):
    path = tmp_path / "potential.xvg"
    path.write_text('@ title "Potential"\n0 -10\n1 -100\n2 -20\n')
    result = minimum_potential_frame(path, [0, 2])
    assert result["frame_index"] == 1
    assert result["potential_kJ_mol"] == -20
    with pytest.raises(ValueError, match="timestamp"):
        minimum_potential_frame(path, [0, 1.5])


def test_invalid_components_are_rejected():
    with pytest.raises(ValueError, match="out of range"):
        component_relative_rmsd(np.zeros((4, 3)), np.zeros((4, 3)), {"TCR": [4]})
