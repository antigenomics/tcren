"""Retain alternate-location evidence without changing coordinate/contact inputs."""
from dataclasses import replace

import numpy as np
from polars.testing import assert_frame_equal

from tcren.contacts.geometry import all_atom_contacts
from tcren.structure import io
from tcren.validation.structures import close_atom_pairs


def test_altloc_metadata_keeps_exact_atoms_contacts_and_writer_selection(tmp_path, monkeypatch):
    records = [("A", "CA", "", -4., 0.), ("A", "SG", "A", 0., 0.),
               ("A", "SG", "B", 0., .3), ("B", "CA", "", 6., 0.),
               ("B", "SG", "A", 1.8, 0.), ("B", "SG", "B", 1.9, .2)]
    lines = []
    for serial, (chain, name, altloc, x, y) in enumerate(records, 1):
        lines.append(f"ATOM  {serial:5d} {io._atom_name_field(name)}{altloc or ' '}CYS {chain}   1    "
                     f"{x:8.3f}{y:8.3f}{0.:8.3f}  0.50  0.00          {name[0]:>2}")
    path = tmp_path / "alternate.pdb"
    path.write_text("\n".join(lines) + "\nEND\n")
    fast = io.parse_structure(path)
    assert io._parse_pdb_fast(path.read_bytes(), "alternate", True) is not None
    atoms = [a for c in fast.chains for r in c.residues for a in r.atoms]
    assert [(a.name, a.altloc) for a in atoms] == [(r[1], r[2]) for r in records]
    expected = np.array([[r[3], r[4], 0.] for r in records], dtype=np.float32).astype(float)
    np.testing.assert_array_equal([a.coord for a in atoms], expected)
    with monkeypatch.context() as patched:
        patched.setattr(io, "_parse_pdb_fast", lambda *a, **k: None)
        reference = io.parse_structure(path)
    reference_atoms = [a for c in reference.chains for r in c.residues for a in r.atoms]
    assert [(a.name, a.altloc) for a in reference_atoms] == [(a.name, a.altloc) for a in atoms]
    np.testing.assert_array_equal([a.coord for a in reference_atoms], expected)

    # The old representation differed only by missing labels: identical order,
    # coordinates, atom multiplicity and exact computed contact table must survive.
    without_labels = replace(fast, chains=[replace(c, residues=[replace(r, atoms=tuple(
        replace(a, altloc="") for a in r.atoms)) for r in c.residues]) for c in fast.chains])
    assert_frame_equal(all_atom_contacts(fast), all_atom_contacts(without_labels))
    rows = close_atom_pairs(path)
    assert len(rows) == 4
    assert {(r["altloc.from"], r["altloc.to"]) for r in rows} == {
        ("A", "A"), ("A", "B"), ("B", "A"), ("B", "B")}
    assert all(r["alternate_location_involved"] and r["possible_disulfide"] for r in rows)
    assert {r["bond_classification"] for r in rows} == {"unknown"}

    for suffix, writer in (("pdb", io.pdb_lines), ("cif", io.cif_lines)):
        output = tmp_path / f"chosen.{suffix}"
        output.write_text("\n".join(writer(fast)) + "\n")
        chosen = [a for c in io.parse_structure(output).chains for r in c.residues for a in r.atoms]
        assert [(a.name, a.altloc) for a in chosen] == [("CA", ""), ("SG", "A")] * 2
        np.testing.assert_array_equal([a.coord for a in chosen], expected[[0, 1, 3, 4]])
