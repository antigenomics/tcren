import numpy as np
import pytest

from tcren.validation.replay import compare_replay
from tcren.validation.structures import map_sequence_roles


SEQUENCES = dict(A="ACDEFG", B="HIKLMN", C="PQRSTV", D="WYADFH")
NAMES = dict(
    A="ALA",
    C="CYS",
    D="ASP",
    E="GLU",
    F="PHE",
    G="GLY",
    H="HIS",
    I="ILE",
    K="LYS",
    L="LEU",
    M="MET",
    N="ASN",
    P="PRO",
    Q="GLN",
    R="ARG",
    S="SER",
    T="THR",
    V="VAL",
    W="TRP",
    Y="TYR",
)


def pdb(path, *, rotation=None, translation=None, tcr_shift=0):
    lines = []
    serial = 1
    for chain, (role, seq) in enumerate(SEQUENCES.items()):
        for i, aa in enumerate(seq):
            point = np.array([i, chain * 10 + (i % 2), i % 3], dtype=float)
            if role in "AB":
                point[0] += tcr_shift
            if rotation is not None:
                point = point @ rotation
            if translation is not None:
                point += translation
            lines.append(
                f"ATOM  {serial:5d}  CA  {NAMES[aa]} {role}{i + 1:4d}    "
                f"{point[0]:8.3f}{point[1]:8.3f}{point[2]:8.3f}  1.00 50.00           C"
            )
            serial += 1
    path.write_text("\n".join(lines) + "\nEND\n")
    return path


def test_global_rigid_transform_does_not_change_replay_rmsd(tmp_path):
    ref = pdb(tmp_path / "ref.pdb")
    model = pdb(
        tmp_path / "model.pdb",
        rotation=np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]]),
        translation=[7, 9, 11],
    )
    result = compare_replay(ref, model, SEQUENCES)
    assert len(result["rows"]) == 18
    assert max(r["rmsd_A"] for r in result["rows"]) < 1e-8
    assert all(r["matched_ca"] == 6 for r in result["spans"])


def test_partner_fit_does_not_erase_tcr_displacement(tmp_path):
    result = compare_replay(
        pdb(tmp_path / "ref.pdb"), pdb(tmp_path / "model.pdb", tcr_shift=8), SEQUENCES
    )
    row = next(
        r
        for r in result["rows"]
        if r["fit_component"] == "MHC" and r["measured_component"] == "TCR"
    )
    assert row["rmsd_A"] == pytest.approx(8)
    row = next(
        r
        for r in result["rows"]
        if r["fit_component"] == "TCR" and r["measured_component"] == "TCR"
    )
    assert row["rmsd_A"] < 1e-8


def test_extended_native_terminal_loss_is_explicit_and_opt_in(tmp_path):
    path = tmp_path / "native.pdb"
    seq = "CDEFGHIK"
    lines = [
        f"ATOM  {i + 1:5d}  CA  {NAMES[aa]} X{i + 1:4d}    "
        f"{float(i):8.3f}{float(i % 2):8.3f}{float(i % 3):8.3f}  1.00 50.00           C"
        for i, aa in enumerate(seq)
    ]
    path.write_text("\n".join(lines) + "\nEND\n")
    query = {"A": "ACDEFG"}
    with pytest.raises(ValueError, match="Absent/ambiguous"):
        map_sequence_roles(path, "A", query)
    _, identity = map_sequence_roles(path, "A", query, allow_extended=True)
    assert identity[0]["reference_start"] == 1
    assert identity[0]["output_start"] == 0
    assert identity[0]["matched_residues"] == 5


def test_explicit_spans_exclude_declared_substitution(tmp_path):
    ref = pdb(tmp_path / "ref.pdb")
    model = pdb(tmp_path / "model.pdb")
    ref.write_text(ref.read_text().replace("GLU A   4", "LYS A   4"))
    spans = {
        "A": [
            dict(output_chain="A", input_start=0, output_start=0, length=3),
            dict(output_chain="A", input_start=4, output_start=4, length=2),
        ]
    }
    result = compare_replay(ref, model, SEQUENCES, reference_spans=spans)
    assert result["spans"][0]["input_indices"] == [0, 1, 2, 4, 5]
    assert max(r["rmsd_A"] for r in result["rows"]) < 1e-8
    spans["A"][0]["length"] = 4
    with pytest.raises(ValueError, match="sequence differs"):
        compare_replay(ref, model, SEQUENCES, reference_spans=spans)


def test_self_check_records_alert_without_repairing_model(tmp_path):
    from tcren.validation.replay import inspect_replay

    ref = pdb(tmp_path / "ref.pdb")
    model = pdb(tmp_path / "model.pdb", tcr_shift=8)
    original = model.read_bytes()
    result = inspect_replay(ref, model, SEQUENCES, alerts={"MHC/TCR": 5})
    assert result["status"] == "requires_review"
    finding = next(r for r in result["findings"] if r["kind"] == "rmsd_alert")
    assert finding["observed_A"] == pytest.approx(8)
    assert model.read_bytes() == original
    with pytest.raises(ValueError, match="threshold"):
        inspect_replay(ref, model, SEQUENCES, alerts={"MHC/TCR": float("nan")})


def test_contacts_use_shared_query_indices_not_author_or_output_offsets():
    import polars as pl
    from tcren.validation.replay import compare_replay_contacts

    chains = []
    mappings = {}
    for name in ["old", "new"]:
        mappings[name] = []
        for role, kind in [("A", "TRA"), ("B", "TRB"), ("C", "PEPTIDE")]:
            start = 1 if name == "old" and role == "A" else 0
            seq = SEQUENCES[role][start:]
            chains.append(
                {"complex.id": name, "chain.id": role, "chain.type": kind, "sequence": seq}
            )
            mappings[name].append(
                dict(
                    role=role,
                    output_chain=role,
                    reference_start=start,
                    output_start=0,
                    matched_residues=len(seq),
                )
            )
    contacts = []
    for name, index, aa in [("old", 1, "D"), ("new", 2, "D"), ("new", 0, "A")]:
        contacts.append(
            {
                "complex.id": name,
                "chain.id.from": "A",
                "chain.id.to": "C",
                "chain.type.from": "TRA",
                "chain.type.to": "PEPTIDE",
                "residue.index.from": index,
                "residue.index.to": 1,
                "residue.aa.from": aa,
                "residue.aa.to": "Q",
            }
        )
    result, coverage = compare_replay_contacts(
        pl.DataFrame(contacts), pl.DataFrame(chains), SEQUENCES, mappings, reference_id="old"
    )
    assert result["contact_jaccard"].to_list() == [1, 1]
    assert coverage["common_query_indices"]["A"] == [1, 2, 3, 4, 5]
    assert coverage["excluded_residue_pairs"]["new"] == 1
    contacts[0]["chain.id.from"] = "X"
    with pytest.raises(ValueError, match="chain differs"):
        compare_replay_contacts(
            pl.DataFrame(contacts), pl.DataFrame(chains), SEQUENCES, mappings, reference_id="old"
        )

    contacts[0]["chain.id.from"] = "A"
    contacts[2]["residue.index.from"] = -1
    with pytest.raises(ValueError, match="source sequence/index"):
        compare_replay_contacts(
            pl.DataFrame(contacts), pl.DataFrame(chains), SEQUENCES, mappings, reference_id="old"
        )


def test_ambiguous_truncated_extended_chain_cannot_be_hidden_by_unique_chain(tmp_path):
    from tcren.validation.inventory import compatible_sequence_span

    assert compatible_sequence_span("CDEFGCDEFG", "ACDEFG", terminal_budget=1) is None
    with pytest.raises(ValueError, match="Ambiguous"):
        compatible_sequence_span("CDEFGCDEFG", "ACDEFG", terminal_budget=1, reject_ambiguous=True)

    lines = []
    for chain, sequence in [("X", "CDEFGHIK"), ("Y", "CDEFGCDEFG")]:
        for i, aa in enumerate(sequence):
            lines.append(
                f"ATOM  {len(lines) + 1:5d}  CA  {NAMES[aa]} {chain}{i + 1:4d}    "
                f"{float(i):8.3f}{float(i % 2):8.3f}{float(i % 3):8.3f}  1.00 50.00           C"
            )
    path = tmp_path / "ambiguous.pdb"
    path.write_text("\n".join(lines) + "\nEND\n")
    with pytest.raises(ValueError, match="Ambiguous"):
        map_sequence_roles(path, "A", {"A": "ACDEFG"}, allow_extended=True)


def test_confidence_selected_pose_outlier_is_flagged_without_native_selection(tmp_path):
    from tcren.validation.replay import replay_ensemble_consistency

    models = {str(i): pdb(tmp_path / f"{i}.pdb") for i in range(3)}
    models["selected"] = pdb(tmp_path / "selected.pdb", tcr_shift=25)
    result = replay_ensemble_consistency(
        models, SEQUENCES, threshold_A=10, confidence_selected_id="selected"
    )
    assert result["dominant_fraction"] == pytest.approx(0.75)
    assert result["confidence_selected_in_dominant_cluster"] is False
    assert result["review_required"]
    assert result["clusters"][0]["members"] == ["0", "1", "2"]
    tied = replay_ensemble_consistency(
        {"0": models["0"], "selected": models["selected"]},
        SEQUENCES,
        threshold_A=10,
        confidence_selected_id="selected",
    )
    assert tied["dominant_cluster"] is None
    assert tied["confidence_selected_in_dominant_cluster"] is None
    assert tied["review_required"]


def test_extreme_clash_is_flagged_even_for_identical_reference(tmp_path):
    from tcren.validation.replay import inspect_replay

    model = pdb(tmp_path / "overlap.pdb")
    lines = model.read_text().splitlines()
    # Move one TRB CA onto a TRA CA, retaining exact chain/sequence identity.
    lines[6] = lines[6][:30] + lines[0][30:54] + lines[6][54:]
    model.write_text("\n".join(lines) + "\n")
    result = inspect_replay(model, model, SEQUENCES, alerts={"all/all": 5})
    assert not any(r["kind"] == "rmsd_alert" for r in result["findings"])
    assert any(r["kind"] == "extreme_interchain_clash" for r in result["findings"])
    assert result["minimum_interchain_distance_A"] == 1
    with pytest.raises(ValueError, match="Extreme-clash"):
        inspect_replay(
            model, model, SEQUENCES, alerts={"all/all": 5}, minimum_interchain_distance_A=3
        )


def test_producer_coordinate_identity_ignores_metadata_but_not_atom_frame(tmp_path):
    from tcren.validation.structures import same_atom_coordinates

    reference = pdb(tmp_path / "ref.pdb")
    metadata = pdb(tmp_path / "metadata.pdb")
    metadata.write_text("REMARK producer timestamp differs\n" + metadata.read_text())
    assert same_atom_coordinates(reference, metadata)
    moved = pdb(tmp_path / "moved.pdb", translation=[1, 2, 3])
    assert not same_atom_coordinates(reference, moved)
    metadata.write_text(metadata.read_text().replace(" CA ", " CB ", 1))
    assert not same_atom_coordinates(reference, metadata)
