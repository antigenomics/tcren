"""Physical atom identities and symmetry for the opt-in contact observations."""
from dataclasses import replace

import numpy as np
import polars as pl
import pytest
from polars.testing import assert_frame_equal

from tcren.contacts.geometry import all_atom_contacts
from tcren.structure.model import Atom, Chain, RegionMarkup, Residue, Structure


def complex_fixture():
    def residue(aa, atoms, index=0):
        return Residue(index, 10 + index, "A", aa, "LYS" if aa == "K" else "ASP",
                       tuple(Atom(name, name[0], np.array(xyz, dtype=float), alt)
                             for name, xyz, alt in atoms))
    lys = residue("K", [("N", (-1, 1, 0), ""), ("CA", (0, 0, 0), ""),
                        ("CB", (1, 0, 0), "A"), ("CB", (9, 0, 0), "B"),
                        ("NZ", (0, 1, 0), "")])
    asp = residue("D", [("N", (3, 1, 1), ""), ("CA", (3, 0, 0), ""),
                        ("CB", (3, 0, 1), ""), ("OD1", (5, 1, 0), "")])
    helix = residue("D", [("CA", (0, 3, 0), "")])
    floor = residue("D", [("CA", (0, 2, 0), "")], 1)
    chains = [Chain("T", [lys], "TRA", "TRAB", regions=[RegionMarkup("CDR3", 0, 0, "K", [lys])]),
              Chain("B", [], "TRB", "TRAB"), Chain("P", [asp], "PEPTIDE", "PEPTIDE"),
              Chain("M", [helix, floor], "MHCa", "MHCI", regions=[
                  RegionMarkup("HELIX_A1", 0, 0, "D", [helix]),
                  RegionMarkup("GROOVE_FLOOR", 1, 1, "D", [floor])])]
    return Structure("synthetic", chains)


def test_contact_observations_identity_geometry_and_unchanged_legacy():
    from tcren.contacts.observations import contact_observations, summarize_contact_observations
    s = complex_fixture()
    before = all_atom_contacts(s)
    atoms = contact_observations(s)
    from tcren.exports import export_structure
    assert_frame_equal(export_structure(s, "contact-atoms", None, None, interface="all"), atoms)
    assert set(atoms["interface"]) == {"tcr_peptide", "tcr_mhc", "peptide_mhc"}
    assert set(atoms["mhc_class"]) == {"I"}
    assert "B" not in set(atoms["altloc.from"]) | set(atoms["altloc.to"])
    assert atoms["dist"].max() <= 5
    assert atoms.filter(pl.col("dist") == 5).height > 0
    assert atoms.filter((pl.col("chain.id.from") == "M") &
                        (pl.col("residue.index.from") == 1)).is_empty()
    assert atoms["residue.index.pdb.from"].null_count() == 0
    assert atoms["insertion_code.from"].unique().to_list() == ["A"]
    assert atoms.select("atom.index.from", "atom.index.to").is_duplicated().sum() == 0
    assert atoms["chem.a"].to_list() == [min(a, b) for a, b in zip(
        atoms["chem.from"], atoms["chem.to"])]
    tp = atoms.filter(pl.col("interface") == "tcr_peptide")
    assert tp["phi.from_rad"].drop_nulls().n_unique() == 1
    assert tp["phi.to_rad"].drop_nulls().n_unique() == 1
    assert_frame_equal(all_atom_contacts(s), before)
    reversed_chains = replace(s, chains=list(reversed(s.chains)))
    other = contact_observations(reversed_chains)
    keys = ["chain.id.from", "residue.index.from", "atom.from", "chain.id.to", "atom.to"]
    cols = keys + ["chem.a", "chem.b", "dist", "phi.from_rad", "phi.to_rad", "omega_rad"]
    assert_frame_equal(atoms.select(cols).sort(keys), other.select(cols).sort(keys))
    # Renaming the physical chains reverses endpoint order, not just traversal order.
    renamed = replace(s, chains=[replace(c, chain_id={"T": "A", "P": "Z"}.get(c.chain_id, c.chain_id))
                                 for c in s.chains])
    swapped = contact_observations(renamed, interface="tcr_peptide").sort("chem.a", "chem.b")
    original = tp.sort("chem.a", "chem.b")
    for name in ("theta", "phi"):
        np.testing.assert_allclose(original[f"{name}.from_rad"].to_numpy(),
                                   swapped[f"{name}.to_rad"].to_numpy(), equal_nan=True)
    np.testing.assert_allclose(original["omega_rad"].to_numpy(), swapped["omega_rad"].to_numpy(), equal_nan=True)
    summary = summarize_contact_observations(atoms)
    assert summary.filter(pl.col("level") == "atom")["n_contacts"].sum() == atoms.height
    assert summary.filter(pl.col("level") == "residue")["n_contacts"].sum() == 3


def test_contact_observations_empty_glycine_and_missing_geometry():
    from tcren.contacts.observations import contact_observations, summarize_contact_observations
    s = complex_fixture()
    gly = replace(s.chain("T").residues[0], aa="G", resname="GLY")
    s.chain("T").residues = [gly]
    s.chain("T").regions = [RegionMarkup("CDR3", 0, 0, "G", [gly])]
    atoms = contact_observations(s, interface="tcr_peptide")
    for side in ("from", "to"):
        g = atoms.filter(pl.col(f"residue.aa.{side}") == "G")
        assert set(g[f"atom.{side}"]) <= {"N", "CA", "C", "O", "OXT"}
    assert atoms["omega_rad"].null_count() == atoms.height
    s.chain("T").residues = [replace(gly, atoms=tuple(replace(a, coord=np.zeros(3)) for a in gly.atoms))]
    assert contact_observations(s)["omega_rad"].null_count() == contact_observations(s).height
    empty = contact_observations(s, cutoff=0.00001)
    assert empty.is_empty() and "omega_rad" in empty.columns
    assert summarize_contact_observations(empty).is_empty()
    with pytest.raises(ValueError, match="interface"):
        contact_observations(s, interface="invented")
    s.chain("M").regions = []
    with pytest.raises(ValueError, match="helix region annotation"):
        contact_observations(s)


def test_opt_in_indices_match_flattened_heavy_atom_coordinates():
    s = complex_fixture()
    old = all_atom_contacts(s, atom_pairs=True)
    actual = all_atom_contacts(s, atom_pairs=True, atom_indices=True)
    assert_frame_equal(actual.drop("atom.index.from", "atom.index.to"), old)
    xyz = np.array([a.coord for c in s.chains for r in c.residues for a in r.atoms
                    if a.element != "H"])
    distances = np.linalg.norm(xyz[actual["atom.index.from"]] -
                               xyz[actual["atom.index.to"]], axis=1)
    np.testing.assert_allclose(distances, actual["dist"])


def test_geometry_matches_existing_scalar_convention_and_rigid_transform():
    from tcren.contacts.observations import _angle, _dihedral, contact_observations
    from tcren.refine.sidechain import bond_angle, dihedral
    points = np.array([[1., 1, 0], [0, 0, 0], [1, 0, 0], [1, 1, 1]])
    assert _angle(*[p[None, :] for p in points[:3]])[0] == pytest.approx(bond_angle(*points[:3]))
    assert _dihedral(*[p[None, :] for p in points])[0] == pytest.approx(dihedral(*points))
    s = complex_fixture()
    original = contact_observations(s)
    for c in s.chains:
        c.residues = [replace(r, atoms=tuple(replace(a, coord=a.coord[[1, 2, 0]] + 7) for a in r.atoms))
                      for r in c.residues]
    moved = contact_observations(s)
    cols = ["dist", "omega_rad", "theta.from_rad", "theta.to_rad", "phi.from_rad", "phi.to_rad"]
    np.testing.assert_allclose(original.select(cols).to_numpy(), moved.select(cols).to_numpy(), equal_nan=True)


def test_residue_summary_keeps_binder_target_roles_when_physical_order_changes():
    from tcren.contacts.observations import contact_observations, summarize_contact_observations
    structure = complex_fixture()
    # Physical TCR K→peptide D is stored peptide-first for these chain IDs.
    original = summarize_contact_observations(contact_observations(structure))
    renamed = replace(structure, chains=[
        replace(c, chain_id={"T": "A", "P": "Z"}.get(c.chain_id, c.chain_id))
        for c in structure.chains])
    assert_frame_equal(original, summarize_contact_observations(contact_observations(renamed)))
    residue = original.filter(pl.col("level") == "residue")
    assert residue.filter(pl.col("interface") == "tcr_peptide").select("chem.a", "chem.b").rows() == [("K", "D")]
    assert residue.filter(pl.col("interface") == "tcr_mhc").select("chem.a", "chem.b").rows() == [("K", "D")]
    assert residue.filter(pl.col("interface") == "peptide_mhc").select("chem.a", "chem.b").rows() == [("D", "D")]


def test_directed_views_preserve_physical_edges_for_role_selection():
    from tcren.contacts.table import directed_views
    from tcren.contactmap import ContactMap
    structure = complex_fixture()
    physical = all_atom_contacts(structure, atom_indices=True)
    views = directed_views(physical)
    assert views.height == 2*physical.height
    assert_frame_equal(views.head(physical.height), physical)
    reverse = views.tail(physical.height)
    assert reverse["atom.index.from"].to_list() == physical["atom.index.to"].to_list()
    assert reverse["atom.index.to"].to_list() == physical["atom.index.from"].to_list()
    assert reverse["residue.aa.from"].to_list() == physical["residue.aa.to"].to_list()
    assert reverse["residue.aa.to"].to_list() == physical["residue.aa.from"].to_list()
    interface = ContactMap.from_structure(structure).interface("tcr_peptide")
    assert interface.height == 1
    assert interface.select("residue.aa.from", "residue.aa.to").rows() == [("K", "D")]
