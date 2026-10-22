"""Hand-counted contact motifs and full-availability opportunities."""
from dataclasses import replace

import numpy as np
import polars as pl
import pytest
from polars.testing import assert_frame_equal

from tcren.contacts import class_features
from tcren.structure.model import Atom, Chain, RegionMarkup, Residue, Structure


def fixture(model="one", shift=0.):
    def residue(index, aa, xs):
        return Residue(index, 10 + index, "", aa, "ALA", tuple(
            Atom(name, "C", np.array([x, 0., 0.])) for name, x in xs))
    cdr = residue(0, "A", [("CA", 0.), ("CB", 1.)])
    fr = residue(1, "A", [("CA", 100.)])
    pep1 = residue(4, "D", [("CA", 2. + shift), ("CB", 3. + shift)])
    pep2 = residue(8, "D", [("CA", 200. + shift)])
    helix = residue(0, "D", [("CA", 5. + shift)])
    floor = residue(1, "D", [("CA", 0.)])
    return Structure(model, [
        Chain("T", [cdr, fr], "TRA", "TRAB", regions=[RegionMarkup("CDR3", 0, 0, "A", [cdr])]),
        Chain("B", [], "TRB", "TRAB"),
        Chain("P", [pep1, pep2], "PEPTIDE", "PEPTIDE"),
        Chain("M", [helix, floor], "MHCa", "MHCI", regions=[RegionMarkup("HELIX_A1", 0, 0, "D", [helix])]),
        Chain("b", [floor], "B2M", "MHCI")])


def assemble(*outputs):
    return [pl.concat([output[i] for output in outputs]).unique() for i in (2, 3, 4, 5)]


def test_hand_counts_availability_and_residue_collapse():
    output = class_features.features_for_structure(fixture())
    inventory, pairs, counts, binders, targets, dictionary = output
    assert inventory.height == 9
    assert pairs.height == 6  # Four CDR:peptide atom pairs, two CDR:MHC-helix pairs.
    assert pairs["dist"].max() == 5.  # Inclusive cutoff.
    assert pairs.select("binder.atom.source_index", "target.atom.source_index").is_duplicated().sum() == 0
    rows = class_features.opportunity_rows(counts, binders, targets, dictionary, ["one"]).join(dictionary, on="feature")
    chem_atom = rows.filter((pl.col("view") == "chemistry") & (pl.col("level") == "atom") &
                            (pl.col("target.chain.type") == "PEPTIDE"))
    assert chem_atom["C"].sum() == 4
    assert chem_atom["O"].sum() == 6  # Both selected CDR atoms × all three peptide atoms.
    chem_res = rows.filter((pl.col("view") == "chemistry") & (pl.col("level") == "residue") &
                           (pl.col("target.chain.type") == "PEPTIDE"))
    assert chem_res.select("C", "O").rows() == [(1, 2)]
    site_res = rows.filter((pl.col("view") == "site") & (pl.col("level") == "residue") &
                           (pl.col("target.chain.type") == "PEPTIDE"))
    assert site_res.select("C", "O", "target.peptide.pos").rows() == [(1, 1, 1)]
    assert binders["N_binder"].sum() == 10  # Per view: three atoms and two residues, including FR.
    assert "UNASSIGNED" not in set(dictionary["binder.region.type"])  # No FR contact yet.
    assert rows.filter(pl.col("C") > pl.col("O")).is_empty()


def test_zero_contact_model_has_opportunities_and_absent_amino_acid_has_none():
    original = class_features.features_for_structure(fixture())
    far = class_features.features_for_structure(fixture("far", 1000.))
    missing_structure = fixture("missing", 1000.)
    missing_structure.chain("P").residues = [replace(r, aa="E") for r in missing_structure.chain("P").residues]
    missing = class_features.features_for_structure(missing_structure)
    assert far[1].is_empty() and far[2].is_empty() and far[5].is_empty()
    counts, binders, targets, dictionary = assemble(original, far, missing)
    rows = class_features.opportunity_rows(counts, binders, targets, dictionary, ["one", "far", "missing"])
    far_rows = rows.filter(pl.col("pdb.id") == "far")
    assert far_rows.height == rows.filter(pl.col("pdb.id") == "one").height
    assert far_rows["C"].sum() == 0 and far_rows["O"].min() > 0
    missing_peptide = rows.join(dictionary, on="feature").filter(
        (pl.col("pdb.id") == "missing") & (pl.col("target.chain.type") == "PEPTIDE"))
    assert missing_peptide.is_empty()


def test_framework_contacts_included_and_floor_b2m_excluded():
    structure = fixture()
    chain = structure.chain("T")
    framework = chain.residues[1]
    chain.residues[1] = replace(framework, atoms=(replace(framework.atoms[0], coord=np.array([2., 0., 0.])),))
    inventory, pairs, _, _, _, dictionary = class_features.features_for_structure(structure)
    assert pairs.height == 9
    assert "UNASSIGNED" in set(dictionary["binder.region.type"])
    targets = inventory.filter(pl.col("atom.source_index").is_in(pairs["target.atom.source_index"].implode()))
    assert set(targets["chain.type"]) == {"PEPTIDE", "MHCa"}
    assert targets.filter(pl.col("role") == "b2m").is_empty()
    assert targets.filter((pl.col("chain.type") == "MHCa") & (pl.col("residue.index") == 1)).is_empty()


def test_altloc_and_rigid_translation_preserve_features():
    structure = fixture()
    chain = structure.chain("T")
    residue = chain.residues[0]
    chain.residues[0] = replace(residue, atoms=residue.atoms + (
        Atom("CB", "C", np.array([100., 0., 0.]), "B"),
        Atom("CB", "C", np.array([1., 0., 0.]), "A")))
    original = class_features.features_for_structure(structure)
    for c in structure.chains:
        c.residues = [replace(r, atoms=tuple(replace(a, coord=a.coord + np.array([3., 5., 8.])) for a in r.atoms)) for r in c.residues]
    moved = class_features.features_for_structure(structure)
    for i in range(1, 6):
        assert_frame_equal(original[i], moved[i])
    assert original[1].height == 6


@pytest.mark.parametrize("cutoff", [0., -1., np.nan, np.inf])
def test_invalid_cutoff_fails(cutoff):
    with pytest.raises(ValueError, match="cutoff"):
        class_features.features_for_structure(fixture(), cutoff=cutoff)


def test_duplicate_counts_and_invalid_C_fail():
    _, _, counts, binders, targets, dictionary = class_features.features_for_structure(fixture())
    with pytest.raises(ValueError, match="duplicate"):
        class_features.opportunity_rows(pl.concat([counts, counts]), binders, targets, dictionary, ["one"])
    with pytest.raises(ValueError, match="C.*O"):
        class_features.opportunity_rows(counts.with_columns(pl.lit(100, dtype=pl.Int64).alias("C")), binders, targets, dictionary, ["one"])


def test_empty_inventory_returns_typed_empty_tables():
    structure = fixture()
    for chain in structure.chains:
        chain.residues = [replace(r, atoms=()) for r in chain.residues]
    with pytest.raises(ValueError, match="MHC helix heavy atoms"):
        class_features.features_for_structure(structure)
    outputs = class_features.features_for_structure(fixture())
    rows = class_features.opportunity_rows(*[table.head(0) for table in outputs[2:]], ["one"])
    assert rows.schema == {"pdb.id": pl.String, "feature": pl.String, "C": pl.Int64, "O": pl.Int64}


@pytest.mark.parametrize("table_index,count", [(2, -1), (3, 0), (4, -1)])
def test_nonpositive_input_counts_fail(table_index, count):
    output = list(class_features.features_for_structure(fixture()))
    column = {2: "C", 3: "N_binder", 4: "N_target"}[table_index]
    output[table_index] = output[table_index].with_columns(pl.lit(count, dtype=pl.Int64).alias(column))
    with pytest.raises(ValueError, match="positive"):
        class_features.opportunity_rows(*output[2:], ["one"])


def test_opportunity_overflow_fails_before_integer_product():
    _, _, counts, binders, targets, dictionary = class_features.features_for_structure(fixture())
    binders = binders.with_columns(pl.lit(2**62, dtype=pl.Int64).alias("N_binder"))
    targets = targets.with_columns(pl.lit(5, dtype=pl.Int64).alias("N_target"))
    with pytest.raises(ValueError, match="overflow"):
        class_features.opportunity_rows(counts, binders, targets, dictionary, ["one"])


def test_missing_mhc_helices_fail_instead_of_becoming_no_opportunities():
    structure = fixture()
    for chain in structure.by_type('MHCa', 'MHCb'):
        chain.regions = []
    with pytest.raises(ValueError, match='MHC helix'):
        class_features.features_for_structure(structure)
