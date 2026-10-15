import json
import numpy as np
import pytest

from tcren.validation.md import component_relative_rmsd, minimum_potential_frame
from tcren.validation.prediction_inputs import audit_af3_input, disable_af3_pairing
from tcren.prediction import PredictionConfidence, confidence_chain_pairs


def test_boltz_control_retains_insertions_and_prevents_header_pairing(tmp_path):
    import csv
    from tcren.validation.prediction_inputs import af3_to_boltz_input

    source = tmp_path / "af3.json"
    source.write_text(
        json.dumps(
            {
                "dialect": "alphafold3",
                "modelSeeds": [42],
                "sequences": [
                    {
                        "protein": {
                            "id": "A",
                            "sequence": "ACD",
                            "unpairedMsa": ">query\nACD\n>UniRef100_shared\nAc-D\n",
                            "templates": [],
                        }
                    },
                    {
                        "protein": {
                            "id": "B",
                            "sequence": "EFG",
                            "unpairedMsa": ">query\nEFG\n>UniRef100_shared\nE-G\n",
                            "templates": [],
                        }
                    },
                ],
            }
        )
    )
    report = af3_to_boltz_input(source, tmp_path / "boltz")
    with (tmp_path / "boltz/A.csv").open() as handle:
        rows = list(csv.DictReader(handle))
    assert rows == [{"key": "-1", "sequence": "ACD"}, {"key": "-1", "sequence": "Ac-D"}]
    assert report["templates"] == "none transferred"
    with pytest.raises(FileExistsError):
        af3_to_boltz_input(source, tmp_path / "boltz")


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
    cell = f"CRYST1{20:9.3f}{20:9.3f}{20:9.3f}{90:7.2f}{90:7.2f}{90:7.2f} P 1\n"
    trajectory.write_text(trajectory.read_text().replace("TITLE", cell + "TITLE"))
    result = inspect_md_trajectory(
        trajectory, reference, {"TCR": ["A"], "MHC": ["B"]}, energy, contact_cutoff=1.1
    )
    assert result["initial_image_status"] == "nearest_image_chain_centroids"

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


def test_controlled_input_preserves_sequences_templates_unpaired(tmp_path):
    source = tmp_path / "original.json"
    protein = {
        "id": "A",
        "sequence": "ACD",
        "unpairedMsa": ">q\nACD\n>x\nA-D\n",
        "pairedMsa": ">q\nACD\n",
        "templates": [],
    }
    source.write_text(
        json.dumps(
            {"dialect": "alphafold3", "modelSeeds": [42], "sequences": [{"protein": protein}]}
        )
    )
    target = tmp_path / "controlled.json"
    disable_af3_pairing(source, target, {"A": "ACD"}, name="arm", seeds=[42])
    output = json.loads(target.read_text())["sequences"][0]["protein"]
    assert output == {**protein, "pairedMsa": ""}
    protein["unpairedMsa"] = ">q\nACE\n"
    source.write_text(
        json.dumps(
            {"dialect": "alphafold3", "modelSeeds": [42], "sequences": [{"protein": protein}]}
        )
    )
    with pytest.raises(ValueError, match="MSA query"):
        audit_af3_input(source, {"A": "ACD"})


def test_template_audit_checks_actual_cif_polymer_bounds(tmp_path):
    cif = "data_test\n_entry.id TEST\n_entity.pdbx_description 'test protein'\nloop_\n_entity_poly_seq.entity_id\n_entity_poly_seq.num\n1 1\n1 2\n1 3\n#\nloop_\n_atom_site.label_seq_id\n1\n2\n3\n#\n"
    p = tmp_path / "input.json"
    data = {
        "dialect": "alphafold3",
        "modelSeeds": [42],
        "sequences": [
            {
                "protein": {
                    "id": "A",
                    "sequence": "ACD",
                    "templates": [{"queryIndices": [0], "templateIndices": [3], "mmcif": cif}],
                }
            }
        ],
    }
    p.write_text(json.dumps(data))
    with pytest.raises(ValueError, match="polymer length"):
        audit_af3_input(p, {"A": "ACD"})
    data["sequences"][0]["protein"]["templates"][0]["templateIndices"] = [2]
    p.write_text(json.dumps(data))
    audit = audit_af3_input(p, {"A": "ACD"})
    assert audit["chains"][0]["template_details"][0]["polymer_residues"] == 3
    assert audit["chains"][0]["template_details"][0]["mapped_resolved_residues"] == 1


def test_invalid_components_are_rejected():
    with pytest.raises(ValueError, match="out of range"):
        component_relative_rmsd(np.zeros((4, 3)), np.zeros((4, 3)), {"TCR": [4]})


def test_chain_pair_confidence_does_not_require_contacts_or_invent_af2_pairs():
    c = PredictionConfidence(
        "summary",
        "af3",
        chain_pairs={
            ("A", "C"): {"chain_pair_iptm": 0.1},
            ("C", "A"): {"chain_pair_iptm": 0.2},
            ("A", "A"): {"chain_pair_iptm": 0.9},
        },
    )
    rows = confidence_chain_pairs(c).to_dicts()
    assert [r["chain_pair_iptm"] for r in rows] == [0.1, 0.2]
    assert confidence_chain_pairs(PredictionConfidence("stats", "tcrmodel2")).height == 0


def test_md_plot_separates_label_from_observable_and_exports_both_formats(tmp_path, monkeypatch):
    import matplotlib.figure
    from tcren.validation.md import plot_validation

    directory = tmp_path / "md"
    directory.mkdir()
    (directory / "energy.xvg").write_text("0 -10 0 0 300 1 1000\n1 -11 0 0 301 1 1001\n")
    (directory / "backbone-rmsd.xvg").write_text("0 0\n1 0.1\n")
    saved = []
    original = matplotlib.figure.Figure.savefig

    def check(figure, path, **kwargs):
        assert figure.axes[0].get_title() == ""
        assert figure.axes[0].get_title(loc="left") == "Synthetic complete complex\nPotential energy"
        saved.append(path.suffix)
        return original(figure, path, **kwargs)

    monkeypatch.setattr(matplotlib.figure.Figure, "savefig", check)
    plot_validation([directory], ["Synthetic complete complex"], tmp_path / "validation")
    assert saved == [".png", ".pdf"]
    assert (tmp_path / "validation.png").stat().st_size > 0
    assert (tmp_path / "validation.pdf").stat().st_size > 0


@pytest.mark.parametrize("gamma", [90.0, 60.0])
def test_initial_periodic_image_check_handles_triclinic_cells(gamma):
    from tcren.validation.md import validate_initial_complex_image
    from tcren.structure.model import Atom, Chain, Residue, Structure

    cell = f"CRYST1{20:9.3f}{20:9.3f}{20:9.3f}{90:7.2f}{90:7.2f}{gamma:7.2f} P 1\n"
    def structure(offset):
        return Structure("synthetic", [Chain(cid, [Residue(0, 1, "", "A", "ALA", (Atom("CA", "C", np.asarray(xyz, dtype=float)),))])
                                      for cid, xyz in [("A", [0, 0, 0]), ("B", [offset, 1, 0])]])
    assert validate_initial_complex_image(cell, structure(1)) == "nearest_image_chain_centroids"
    with pytest.raises(ValueError, match="periodic-image mismatch"):
        validate_initial_complex_image(cell, structure(21))
    assert validate_initial_complex_image("", structure(21)) == "unavailable_no_unit_cell"


def test_partial_reference_selects_complete_queries_without_scoring_unmatched_roles():
    from tcren.validation.structures import compatible_reference_models

    models = [{"construct": kind} for kind in ["tcr", "pmhc", "complete"]]
    roles = {"tcr": ["TRA", "TRB"], "pmhc": ["peptide", "MHC"],
             "complete": ["TRA", "TRB", "peptide", "MHC"]}
    selected = compatible_reference_models(models, ["peptide", "MHC"], roles)
    assert selected == models[1:]
    with pytest.raises(ValueError, match="at least one"):
        compatible_reference_models(models, [], roles)


def test_inventory_visits_nested_collections_and_archive_members(tmp_path):
    import io
    import tarfile
    from tcren.validation.inventory import coordinate_stream

    for name in ["first", "second/nested"]:
        directory = tmp_path / name
        directory.mkdir(parents=True)
        (directory / "coordinate.pdb").write_text("END\n")
    with tarfile.open(tmp_path / "third.tar.gz", "w:gz") as archive:
        member = tarfile.TarInfo("nested/coordinate.pdb")
        member.size = 4
        archive.addfile(member, io.BytesIO(b"END\n"))
    rows = list(coordinate_stream(tmp_path))
    assert len(rows) == 3
    assert len({name for name, _ in rows}) == 3
    assert all(text == "END\n" for _, text in rows)


@pytest.mark.parametrize("chain_id", ["../outside", "/outside", "..\\outside", "", ".", "..", "C:outside"])
def test_input_chain_ids_cannot_escape_export_directory(tmp_path, chain_id):
    from tcren.validation.prediction_inputs import af3_to_boltz_input

    source = tmp_path / "input.json"
    source.write_text(json.dumps({"dialect": "alphafold3", "modelSeeds": [1], "sequences": [
        {"protein": {"id": chain_id, "sequence": "ACD", "templates": []}}]}))
    output = tmp_path / "output"
    with pytest.raises(ValueError, match="safe single-component"):
        af3_to_boltz_input(source, output)
    assert not output.exists()


@pytest.mark.parametrize("selection", [[0.9, 1.1], [True, 1], [False, True], [[0, 1]], [0.0, 1.0]])
def test_component_indices_reject_silent_numeric_coercion(selection):
    with pytest.raises(ValueError, match="one-dimensional integers"):
        component_relative_rmsd(np.zeros((4, 3)), np.zeros((4, 3)), {"TCR": selection})
