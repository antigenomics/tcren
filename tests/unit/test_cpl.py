"""Unit tests for CPL response-matrix prediction (:mod:`tcren.cpl`).

Two layers. Most checks run on a tiny hand-built contact map whose energies are analytically
checkable, so the referencing algebra is pinned exactly rather than to a tolerance. A second layer
runs the whole path on a real deposited complex, because the toy map cannot catch the failure this
module is most exposed to: a peptide:MHC interface that comes out empty because the structure was
never MHC-annotated, which would silently zero every anchor cell.
"""

from __future__ import annotations

import numpy as np
import polars as pl
import pytest

from tcren.contactmap import ContactMap
from tcren.cpl import (AA20, equimolar_effect, mutation_effect, position_scan, response_matrix)
from tcren.potential import Potential
from tcren.scoring import score_peptides

ASSETS = __import__("pathlib").Path(__file__).resolve().parents[1] / "assets" / "pdb"


def _toy_potential() -> Potential:
    """Directed toy potential over a four-letter alphabet, distinct values so sums are traceable."""
    vals = {("A", "A"): 1.0, ("A", "K"): -2.0, ("A", "G"): 0.1, ("A", "L"): 0.7,
            ("L", "A"): 0.5, ("L", "K"): 3.0, ("L", "G"): 0.2, ("L", "L"): -1.3,
            ("G", "A"): 0.3, ("G", "K"): 0.9, ("G", "G"): -0.4, ("G", "L"): 2.2,
            ("K", "A"): -0.6, ("K", "K"): 1.1, ("K", "G"): 0.8, ("K", "L"): -1.9}
    rows = [{"residue.aa.from": f, "residue.aa.to": t, "value": v} for (f, t), v in vals.items()]
    return Potential(name="toy", matrix=pl.DataFrame(rows), alphabet=("A", "L", "K", "G"))


def _toy_contact_map() -> ContactMap:
    """TCR touches peptide positions 0 and 2; the MHC touches position 1 only.

    Position 1 is therefore a pure ANCHOR -- the case where the TCR term is identically zero and the
    summed score has to degrade to presentation alone without producing a NaN.
    """
    contacts = pl.DataFrame({
        "chain.type.from": ["TRA", "TRB", "PEPTIDE"],
        "chain.type.to": ["PEPTIDE", "PEPTIDE", "MHCa"],
        "residue.aa.from": ["A", "L", "G"],
        "residue.aa.to": ["G", "K", "L"],
        "region.type.from": ["CDR3", "CDR3", "PEPTIDE"],
        "residue.index.from": [10, 20, 1],
        "residue.index.to": [0, 2, 55],
        "region.start.from": [8, 18, 0],
        "region.start.to": [0, 0, 0],
        "pdb.id": ["toy", "toy", "toy"],
    })
    return ContactMap(pdb_id="toy", contacts=contacts, peptide_length=3)


@pytest.fixture
def toy():
    cm, pot = _toy_contact_map(), _toy_potential()
    return cm, pot, response_matrix(cm, "GGK", tcr_potential=pot, mhc_potential=pot)


# ---------------------------------------------------------------- shape and membership
def test_matrix_covers_every_contacting_position_and_all_twenty_residues(toy):
    _cm, _pot, rm = toy
    assert rm.positions == (1, 2, 3)          # 1-based; all three positions touch something
    assert rm.aa == AA20
    assert np.asarray(rm.phi).shape == (3, 20)
    assert rm.peptide == "GGK"


def test_interface_class_names_the_contact_not_the_potential(toy):
    """Position 2 touches only the groove, so it is an anchor; the others face the receptor."""
    _cm, _pot, rm = toy
    assert rm.interface_class == ("receptor", "anchor", "receptor")


def test_anchor_row_carries_presentation_only(toy):
    """At an anchor the summed rule must degrade to presentation alone -- not to NaN, not to zero.

    ``phi`` is a WHOLE-PEPTIDE energy, so the receptor term at an anchor row is not zero: it is the
    constant the other positions contribute. Constant is the invariant that matters, because it is
    exactly what both references subtract away, leaving the groove as the only thing that varies.
    """
    _cm, _pot, rm = toy
    anchor = rm.row_of(2)
    tcr_row = np.asarray(rm.phi_tcr)[anchor]
    assert np.ptp(tcr_row) == pytest.approx(0.0)           # constant: contributes to no difference
    assert np.ptp(np.asarray(rm.phi_mhc)[anchor]) > 0      # the groove does vary across residues
    for ref in ("equimolar", "wild_type"):
        mhc_only = np.asarray(rm.phi_mhc)[anchor]
        base = mhc_only.mean() if ref == "equimolar" else mhc_only[rm.column_of(rm.wild_type_at(2))]
        assert np.allclose(np.asarray(rm.referenced(ref))[anchor], base - mhc_only)


def test_phi_is_the_sum_of_the_two_interfaces(toy):
    _cm, _pot, rm = toy
    assert np.allclose(np.asarray(rm.phi),
                       np.asarray(rm.phi_tcr) + np.asarray(rm.phi_mhc))


def test_phi_cell_equals_scoring_the_threaded_peptide(toy):
    """Every cell must BE the shipped score of the corresponding threaded peptide, not a proxy."""
    cm, pot, rm = toy
    for pos in rm.positions:
        for aa in ("A", "L", "K", "G"):
            threaded = rm.peptide[:pos - 1] + aa + rm.peptide[pos:]
            expect = (float(score_peptides(cm, [threaded], pot, interface="tcr_peptide")["score"][0])
                      + float(score_peptides(cm, [threaded], pot, interface="peptide_mhc")["score"][0]))
            got = np.asarray(rm.phi)[rm.row_of(pos), rm.column_of(aa)]
            assert got == pytest.approx(expect), f"cell ({pos}, {aa})"


# ---------------------------------------------------------------- the two reference states
def test_wild_type_reference_is_zero_at_the_template_residue(toy):
    _cm, _pot, rm = toy
    for pos in rm.positions:
        assert mutation_effect(rm, pos, rm.wild_type_at(pos), reference="wild_type") == 0.0


def test_equimolar_reference_makes_every_row_sum_to_zero(toy):
    """That IS the equimolar reference: the row is centred on its own twenty-residue mean."""
    _cm, _pot, rm = toy
    assert np.allclose(np.asarray(rm.referenced("equimolar")).mean(axis=1), 0.0, atol=1e-12)


def test_the_two_references_differ_by_a_per_position_constant(toy):
    """The whole point of the distinction: within a row they differ by one number, not by shape."""
    _cm, _pot, rm = toy
    diff = np.asarray(rm.referenced("wild_type")) - np.asarray(rm.referenced("equimolar"))
    assert np.allclose(diff, diff[:, [0]])                      # constant along each row
    assert not np.allclose(diff, diff[[0], :])                  # but NOT the same constant per row


def test_equimolar_reference_keeps_the_template_residue_informative(toy):
    """Under the wild-type reference the template cell is a forced zero; here it is a measurement."""
    _cm, _pot, rm = toy
    vals = [equimolar_effect(rm, p) for p in rm.positions]
    assert any(abs(v) > 1e-9 for v in vals)


def test_positive_means_favourable_on_both_references(toy):
    """Lower energy is a better binder, so `reference - candidate` must put the best residue on top."""
    _cm, _pot, rm = toy
    phi = np.asarray(rm.phi)
    for ref in ("equimolar", "wild_type"):
        eff = np.asarray(rm.referenced(ref))
        for row in range(phi.shape[0]):
            assert np.nanargmin(phi[row]) == np.nanargmax(eff[row])


# ---------------------------------------------------------------- the three query forms
def test_mutation_effect_reads_one_cell_of_the_matrix(toy):
    _cm, _pot, rm = toy
    for ref in ("equimolar", "wild_type"):
        got = mutation_effect(rm, 3, "L", reference=ref)
        assert got == pytest.approx(np.asarray(rm.referenced(ref))[rm.row_of(3), rm.column_of("L")])


def test_position_scan_returns_all_twenty_best_first(toy):
    _cm, _pot, rm = toy
    df = position_scan(rm, 1)
    assert df.height == 20
    assert set(df["aa"].to_list()) == set(AA20)
    assert df["effect"].to_list() == sorted(df["effect"].to_list(), reverse=True)
    assert df["is_wt"].sum() == 1


def test_equimolar_effect_is_the_cost_of_giving_the_position_to_the_mixture(toy):
    """It must equal the equimolar cell of the residue being given up, by definition."""
    _cm, _pot, rm = toy
    assert equimolar_effect(rm, 1) == pytest.approx(mutation_effect(rm, 1, rm.wild_type_at(1)))
    assert equimolar_effect(rm, 1, "K") == pytest.approx(mutation_effect(rm, 1, "K"))


def test_equimolar_effect_is_the_mean_over_twenty_minus_the_residue(toy):
    """Spelled out against the definition rather than against another accessor."""
    _cm, _pot, rm = toy
    row = np.asarray(rm.phi)[rm.row_of(3)]
    assert equimolar_effect(rm, 3, "A") == pytest.approx(row.mean() - row[rm.column_of("A")])


def test_to_frame_emits_one_row_per_cell_with_both_references(toy):
    _cm, _pot, rm = toy
    df = rm.to_frame()
    assert df.height == 3 * 20
    assert {"effect_equimolar", "effect_wild_type", "phi_tcr", "phi_mhc"} <= set(df.columns)
    assert df.filter(pl.col("is_wt")).height == 3
    one = rm.to_frame("equimolar")
    assert "effect_wild_type" not in one.columns


# ---------------------------------------------------------------- failure modes
def test_a_position_that_contacts_nothing_raises_rather_than_reading_flat(toy):
    """Silent tolerance is the failure mode: no data must not be reported as no preference."""
    _cm, _pot, rm = toy
    with pytest.raises(KeyError, match="makes no contact"):
        rm.row_of(9)


def test_unknown_residue_and_unknown_reference_raise(toy):
    _cm, _pot, rm = toy
    with pytest.raises(KeyError):
        rm.column_of("B")
    with pytest.raises(ValueError, match="reference must be"):
        rm.referenced("polyalanine")


def test_wrong_peptide_length_raises(toy):
    cm, pot, _rm = toy
    with pytest.raises(ValueError, match="length"):
        response_matrix(cm, "GGKA", tcr_potential=pot, mhc_potential=pot)


def test_a_complex_with_no_peptide_contacts_raises(toy):
    _cm, pot, _rm = toy
    empty = ContactMap(pdb_id="bare", peptide_length=3, contacts=pl.DataFrame({
        "chain.type.from": ["TRA"], "chain.type.to": ["MHCa"],
        "residue.aa.from": ["A"], "residue.aa.to": ["L"], "region.type.from": ["CDR3"],
        "residue.index.from": [10], "residue.index.to": [55],
        "region.start.from": [8], "region.start.to": [0], "pdb.id": ["bare"]}))
    with pytest.raises(ValueError, match="no peptide contacts"):
        response_matrix(empty, "GGK", tcr_potential=pot, mhc_potential=pot)


# ---------------------------------------------------------------- a real complex
@pytest.mark.skipif(not (ASSETS / "1ao7.pdb").exists(), reason="1ao7 asset missing")
def test_real_complex_end_to_end():
    """The whole path on a deposited TCR:pMHC, including the MHC pass the anchors depend on."""
    from tcren import parse_structure
    from tcren.annotation import classify_chains
    from tcren.mhc import annotate_mhc

    s = parse_structure(str(ASSETS / "1ao7.pdb"), pdb_id="1ao7")
    classify_chains(s, organism="human")
    annotate_mhc(s)
    rm = response_matrix(ContactMap.from_structure(s, cutoff=5.0))

    assert len(rm.peptide) == 9                       # 1ao7 carries the HTLV-1 Tax nonamer
    assert rm.peptide == "LLFGYPVYV"
    assert 0 < len(rm.positions) <= 9
    assert np.isfinite(np.asarray(rm.phi)).all()
    # the groove must actually have been annotated -- an unannotated structure yields an empty
    # peptide:MHC interface and a silently all-zero presentation term
    assert np.abs(np.asarray(rm.phi_mhc)).sum() > 0
    # and the referencing identities must survive real data
    assert np.allclose(np.asarray(rm.referenced("equimolar")).mean(axis=1), 0.0, atol=1e-9)
    for pos in rm.positions:
        assert mutation_effect(rm, pos, rm.wild_type_at(pos), reference="wild_type") == 0.0
