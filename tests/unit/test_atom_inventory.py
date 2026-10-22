"""Observed full-context atoms share contact selection and stable source identities."""
from dataclasses import replace

import numpy as np
import polars as pl
import pytest

from tcren.contacts import observations
from tcren.structure.model import Atom, Chain, RegionMarkup, Residue, Structure


def complex_fixture():
    def residue(index, aa="A", atoms=None):
        return Residue(index, 10 + index, "", aa, "ALA" if aa == "A" else "MSE",
                       tuple(atoms or [Atom("CA", "C", np.array([index, 0., 0.]))]))
    cdr = residue(0, atoms=[Atom("CA", "C", np.array([0., 0., 0.])),
                           Atom("CB", "C", np.array([1., 0., 0.])),
                           Atom("CB", "C", np.array([2., 0., 0.]), "B"),
                           Atom("CB", "C", np.array([3., 0., 0.]), "A"),
                           Atom("H", " h ", np.zeros(3)), Atom("D", "D", np.zeros(3))])
    framework = residue(1, "X")
    helix, floor = residue(0), residue(1)
    return Structure("synthetic", [
        Chain("ignored", [residue(0)], None),
        Chain("T", [cdr, framework], "TRA", "TRAB", regions=[
            RegionMarkup("CDR3", 0, 0, "A", [cdr]),
            RegionMarkup("FR4", 1, 1, "X", [framework])]),
        Chain("B", [], "TRB", "TRAB"),
        Chain("P", [residue(7), residue(11)], "PEPTIDE", "PEPTIDE"),
        Chain("M", [helix, floor], "MHCa", "MHCI", regions=[
            RegionMarkup("HELIX_A1", 0, 0, "A", [helix])]),
        Chain("b2m", [residue(0)], "B2M", "MHCI")])


def test_inventory_full_context_and_shared_conformer_selection():
    assert hasattr(observations, "atom_inventory")
    structure = complex_fixture()
    table = observations.atom_inventory(structure)
    assert table.height == 8
    assert table["atom.source_index"].to_list() == [1, 4, 5, 6, 7, 8, 9, 10]
    assert table.filter(pl.col("atom.name") == "CB").select("x", "altloc", "altloc.selected").row(0) == (3., "A", "A")
    assert table.filter(pl.col("chain.id") == "T")["n_altloc_omitted"].to_list() == [2, 2, 0]
    assert table.filter(pl.col("chain.type") == "PEPTIDE")["peptide.pos"].to_list() == [1, 2]
    assert table.filter(pl.col("residue.aa") == "X")["residue.resname"].to_list() == ["MSE"]
    assert table.filter(pl.col("role") == "b2m")["contact_eligible"].to_list() == [False]
    assert table.filter(pl.col("region.type") == "FR4")["contact_role"].to_list() == [None]
    _, metadata = observations._selected_atoms(structure)
    eligible = table.filter(pl.col("contact_eligible"))
    assert eligible["atom.source_index"].to_list() == metadata["atom.source_index"].to_list()
    assert eligible["altloc.selected"].to_list() == metadata["altloc.selected"].to_list()
    assert table.select("chain.id", "residue.index.pdb", "insertion_code", "atom.name").is_duplicated().sum() == 0


def test_inventory_keeps_glycine_extras_but_excludes_from_contacts():
    structure = complex_fixture()
    chain = structure.chain("T")
    gly = replace(chain.residues[0], aa="G", resname="GLY")
    chain.residues[0] = gly
    table = observations.atom_inventory(structure)
    assert table.filter(pl.col("atom.name") == "CB")["contact_eligible"].to_list() == [False]
    assert table.filter(pl.col("atom.name") == "CB")["contact_role"].to_list() == ["tcr"]


@pytest.mark.parametrize("coord", [np.array([np.nan, 0., 0.]), np.array([1., 2.]), np.array([None, 0, 0])])
def test_inventory_rejects_malformed_full_context_coordinates(coord):
    structure = complex_fixture()
    framework = structure.chain("T").residues[1]
    structure.chain("T").residues[1] = replace(framework, atoms=(replace(framework.atoms[0], coord=coord),))
    with pytest.raises(ValueError, match="finite three-dimensional"):
        observations.atom_inventory(structure)


@pytest.mark.parametrize("duplicate", ["chain", "residue", "author"])
def test_inventory_rejects_duplicate_physical_keys(duplicate):
    structure = complex_fixture()
    chain = structure.chain("T")
    if duplicate == "chain":
        structure.chains.append(replace(chain))
    elif duplicate == "residue":
        chain.residues.append(replace(chain.residues[0], pdb_index=99))
    else:
        chain.residues.append(replace(chain.residues[0], seq_index=99))
    with pytest.raises(ValueError, match="duplicate"):
        observations.atom_inventory(structure)


def test_inventory_scope_validation_and_typed_empty():
    structure = complex_fixture()
    for chain in structure.chains:
        chain.residues = [replace(r, atoms=()) for r in chain.residues]
    empty = observations.atom_inventory(structure)
    assert empty.is_empty()
    assert empty.schema["x"] == pl.Float64
    assert empty.schema["contact_eligible"] == pl.Boolean
    structure.chain("P").residues[0] = replace(structure.chain("P").residues[0], aa="X")
    with pytest.raises(ValueError, match="standard-amino-acid peptide"):
        observations.atom_inventory(structure)


@pytest.mark.parametrize("change, message", [
    ("receptor", "alpha-beta"), ("extra_peptide", "exactly one"),
    ("mhc_class", "consistent annotated"), ("missing_mhcb", "two-chain"),
])
def test_inventory_rejects_out_of_scope_complex(change, message):
    structure = complex_fixture()
    if change == "receptor":
        structure.chain("T").chain_type = "TRG"
    elif change == "extra_peptide":
        structure.chains.append(replace(structure.chain("P"), chain_id="P2"))
    elif change == "mhc_class":
        structure.chain("M").chain_supertype = None
    else:
        structure.chain("M").chain_supertype = "MHCII"
    with pytest.raises(ValueError, match=message):
        observations.atom_inventory(structure)


def test_inventory_retains_mhcii_beta_context_without_helix_annotation():
    structure = complex_fixture()
    structure.chain("M").chain_supertype = "MHCII"
    structure.chains.append(replace(structure.chain("M"), chain_id="N", chain_type="MHCb", regions=[]))
    table = observations.atom_inventory(structure)
    beta = table.filter(pl.col("chain.type") == "MHCb")
    assert beta.height == 2
    assert beta["role"].to_list() == ["mhc", "mhc"]
    assert beta["contact_role"].null_count() == 2
    with pytest.raises(ValueError, match="helix region annotation"):
        observations.contact_observations(structure)
