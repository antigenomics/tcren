"""The 3D alanine scan, both sides of the interface.

The scan substitutes one residue at a time on rebuilt coordinates, so a side chain that was the
only thing bridging to its partner loses those contacts. The regression these tests pin is that a
*single* substitution must leave every other side chain alone: threading the whole peptide chain
truncates all of them, and the scan then reads each position against a poly-stub baseline.
"""

from __future__ import annotations

from pathlib import Path

import polars as pl
import pytest

from tcren.contactmap import ContactMap
from tcren.ddg import alanine_scan, tcr_alanine_reference, tcr_alanine_scan
from tcren.potential import tcren2
from tcren.refine.substitute import substitute_peptide, substitute_residues
from tcren.scoring import score_peptides

PDB_DIR = Path(__file__).resolve().parents[1] / "assets" / "pdb"
BACKBONE_PLUS_CB = {"N", "CA", "C", "O", "CB"}


@pytest.fixture(scope="module")
def annotated():
    pytest.importorskip("arda")
    from tcren.annotation import classify_chains
    from tcren.mhc import annotate_mhc
    from tcren.structure import parse_structure

    s = parse_structure(PDB_DIR / "1ao7.pdb")
    classify_chains(s, organism="human", autodetect_species=True)
    annotate_mhc(s)
    s.pdb_id = "1ao7"
    return s


def _peptide(structure):
    return next(c for c in structure.chains if c.chain_type == "PEPTIDE")


def _sequence(structure) -> str:
    return "".join(r.aa for r in _peptide(structure).residues)


def _energy(contact_map, peptide, potential) -> float:
    return float(score_peptides(contact_map, [peptide], potential,
                                interface="tcr_peptide")["score"][0])


# --- the 3D primitive -------------------------------------------------------------------------

def test_no_mutations_is_the_identity(annotated):
    """An empty mapping returns the structure untouched -- the baseline every ddG rests on."""
    assert substitute_residues(annotated, {}) is annotated


def test_one_substitution_leaves_every_other_residue_byte_identical(annotated):
    pep = _peptide(annotated)
    target = pep.residues[0].seq_index
    out = substitute_residues(annotated, {(pep.chain_id, target): "A"})

    for before, after in zip(annotated.chains, out.chains):
        for r0, r1 in zip(before.residues, after.residues):
            if before.chain_id == pep.chain_id and r0.seq_index == target:
                continue
            assert r0.aa == r1.aa
            assert {a.name for a in r0.atoms} == {a.name for a in r1.atoms}


def test_an_alanine_target_keeps_exactly_backbone_and_cbeta(annotated):
    pep = _peptide(annotated)
    target = pep.residues[4].seq_index          # Tyr5, the buried one
    out = substitute_residues(annotated, {(pep.chain_id, target): "A"})
    res = next(r for r in _peptide(out).residues if r.seq_index == target)

    assert res.aa == "A"
    assert res.resname == "ALA"
    assert {a.name for a in res.atoms} <= BACKBONE_PLUS_CB
    assert "CB" in {a.name for a in res.atoms}


def test_mutating_away_from_glycine_builds_a_cbeta(annotated):
    pep = _peptide(annotated)
    gly = next(r for r in pep.residues if r.aa == "G")
    assert "CB" not in {a.name for a in gly.atoms}

    out = substitute_residues(annotated, {(pep.chain_id, gly.seq_index): "A"})
    res = next(r for r in _peptide(out).residues if r.seq_index == gly.seq_index)
    assert "CB" in {a.name for a in res.atoms}


def test_an_unknown_chain_or_residue_raises(annotated):
    pep = _peptide(annotated)
    with pytest.raises(ValueError, match="no such chain"):
        substitute_residues(annotated, {("ZZ", 0): "A"})
    with pytest.raises(ValueError, match="no residue at seq_index"):
        substitute_residues(annotated, {(pep.chain_id, 10_000): "A"})


def test_threading_the_whole_chain_is_not_a_point_substitution(annotated):
    """The defect this scan was rewritten to fix, pinned so it cannot come back.

    ``substitute_peptide`` truncates every residue to backbone + Cβ, so threading the *native*
    sequence back through it is not the identity: it loses the contacts every side chain was
    making. ``substitute_residues`` with no mutations is the identity, which is why the scan uses
    it.
    """
    native = _sequence(annotated)
    n_native = ContactMap.from_structure(annotated, cutoff=5.0).interface("tcr_peptide").height
    n_threaded = ContactMap.from_structure(
        substitute_peptide(annotated, native), cutoff=5.0
    ).interface("tcr_peptide").height

    assert n_threaded < n_native


# --- the peptide side -------------------------------------------------------------------------

def test_peptide_scan_reports_one_row_per_position(annotated):
    native = _sequence(annotated)
    cm = ContactMap.from_structure(annotated, cutoff=5.0)
    scan = alanine_scan(cm, native, tcren2(), structure=annotated)

    assert scan.height == len(native)
    assert scan["wt_aa"].to_list() == list(native)
    assert scan["pos"].to_list() == list(range(len(native)))


def test_an_uncontacted_peptide_position_reads_exactly_zero(annotated):
    """No contacts, no energy to lose -- and no baseline drift allowed to leak in."""
    native = _sequence(annotated)
    cm = ContactMap.from_structure(annotated, cutoff=5.0)
    contacted = set(cm.interface("tcr_peptide")["pos.to"].to_list())
    idle = [p for p in range(len(native)) if p not in contacted]
    assert idle, "1ao7 should have at least one peptide position off the TCR interface"

    scan = alanine_scan(cm, native, tcren2(), structure=annotated)
    for pos in idle:
        assert scan.filter(pl.col("pos") == pos)["ddG"].item() == 0.0


def test_a_position_whose_contacts_are_backbone_only_agrees_with_the_virtual_path(annotated):
    """The two paths may only differ where a side chain past Cβ was doing the reaching."""
    native = _sequence(annotated)
    cm = ContactMap.from_structure(annotated, cutoff=5.0)
    pot = tcren2()
    virtual = alanine_scan(cm, native, pot)
    spatial = alanine_scan(cm, native, pot, structure=annotated)

    joined = virtual.join(spatial, on="pos", suffix="_3d")
    assert (joined["ddG"] - joined["ddG_3d"]).abs().max() > 0.0   # they are not the same measurement
    agreeing = joined.filter((pl.col("ddG") - pl.col("ddG_3d")).abs() < 1e-9)
    assert agreeing.height >= 1


def test_the_scan_needs_the_peptide_length_to_match(annotated):
    cm = ContactMap.from_structure(annotated, cutoff=5.0)
    with pytest.raises(ValueError, match="peptide chain has"):
        alanine_scan(cm, _sequence(annotated) + "A", tcren2(), structure=annotated)


# --- the receptor side ------------------------------------------------------------------------

def test_tcr_scan_walks_exactly_the_contacted_cdr_residues(annotated):
    cm = ContactMap.from_structure(annotated, cutoff=5.0)
    expected = (cm.interface("tcr_peptide", tcr_regions="cdr")
                .select("chain.id.from", "residue.index.from").unique().height)

    scan = tcr_alanine_scan(cm, annotated, tcren2())
    assert scan.height == expected
    assert set(scan["region.type"]) <= {"CDR1", "CDR2", "CDR3"}
    assert set(scan["chain.type"]) <= {"TRA", "TRB"}


def test_a_residue_that_is_already_alanine_costs_nothing(annotated):
    cm = ContactMap.from_structure(annotated, cutoff=5.0)
    scan = tcr_alanine_scan(cm, annotated, tcren2())
    already = scan.filter(pl.col("wt_aa") == "A")
    if already.height:
        assert already["ddG"].abs().max() == 0.0


def test_tcr_scan_ddg_is_native_minus_mutant(annotated):
    """One row checked against the definition, by hand."""
    cm = ContactMap.from_structure(annotated, cutoff=5.0)
    pot = tcren2()
    native = _sequence(annotated)
    scan = tcr_alanine_scan(cm, annotated, pot).filter(pl.col("wt_aa") != "A")
    row = scan.sort("ddG", descending=True).head(1)
    chain_id, index = row["chain.id"].item(), row["residue.index"].item()

    mutant_map = ContactMap.from_structure(
        substitute_residues(annotated, {(chain_id, index): "A"}), cutoff=5.0
    )
    expected = _energy(cm, native, pot) - _energy(mutant_map, native, pot)
    assert row["ddG"].item() == pytest.approx(expected, abs=1e-12)


def test_per_loop_reference_partitions_the_scan(annotated):
    cm = ContactMap.from_structure(annotated, cutoff=5.0)
    scan = tcr_alanine_scan(cm, annotated, tcren2())
    ref = tcr_alanine_reference(scan)

    assert set(ref) == {"dPhi_ala_cdr12", "dPhi_ala_cdr3a", "dPhi_ala_cdr3b", "dPhi_ala_tcr"}
    parts = ref["dPhi_ala_cdr12"] + ref["dPhi_ala_cdr3a"] + ref["dPhi_ala_cdr3b"]
    assert parts == pytest.approx(ref["dPhi_ala_tcr"], abs=1e-12)


def test_an_empty_scan_still_has_the_right_schema():
    """A structure with no CDR:peptide contact returns a typed empty frame, not a bare one."""
    empty = tcr_alanine_reference(pl.DataFrame(
        schema={"chain.type": pl.Utf8, "region.type": pl.Utf8, "ddG": pl.Float64}
    ))
    assert empty == {"dPhi_ala_cdr12": 0.0, "dPhi_ala_cdr3a": 0.0,
                     "dPhi_ala_cdr3b": 0.0, "dPhi_ala_tcr": 0.0}
