"""Unit tests for the fast ΔΔG engine (S4).

Uses the same tiny hand-built contact map / potential as the scoring tests so the
energies are analytically checkable. No external oracle CSV (new method); the
checks are the ΔΔG identities and per-position structure.
"""

from __future__ import annotations

import polars as pl
import pytest

from tcren.contactmap import ContactMap
from tcren.ddg import alanine_scan, ddg, neoantigen_ddg, reference_delta
from tcren.potential import Potential
from tcren.scoring import score_peptides


def _toy_potential() -> Potential:
    vals = {("A", "A"): 1.0, ("A", "K"): -2.0, ("L", "A"): 0.5, ("L", "K"): 3.0,
            ("A", "G"): 0.1, ("L", "G"): 0.2}
    rows = [{"residue.aa.from": fr, "residue.aa.to": to, "value": v}
            for (fr, to), v in vals.items()]
    return Potential(name="toy", matrix=pl.DataFrame(rows), alphabet=("A", "L", "K", "G"))


def _toy_contact_map() -> ContactMap:
    # TCR 'A' contacts peptide pos 0; TCR 'L' contacts peptide pos 2.
    contacts = pl.DataFrame(
        {
            "chain.type.from": ["TRA", "TRB"],
            "chain.type.to": ["PEPTIDE", "PEPTIDE"],
            "residue.aa.from": ["A", "L"],
            "residue.aa.to": ["G", "G"],
            "region.type.from": ["CDR3", "CDR3"],
            "residue.index.from": [10, 20],
            "residue.index.to": [0, 2],
            "region.start.from": [8, 18],
            "region.start.to": [0, 0],
            "pdb.id": ["toy", "toy"],
        }
    )
    return ContactMap(pdb_id="toy", contacts=contacts, peptide_length=3)


def _toy_complex_potential() -> Potential:
    """A full 4x4 potential, so a peptide residue can sit on either side of a contact.

    `_toy_potential` only defines rows for the TCR residues A and L. The peptide:MHC interface puts
    the PEPTIDE on the ``from`` side, so scoring it needs rows for the peptide's own residues too.
    """
    aa = ("A", "L", "K", "G")
    rows = [{"residue.aa.from": f, "residue.aa.to": s, "value": v}
            for f, r in zip(aa, ([1.0, -2.0, 0.5, 3.0], [0.1, 0.2, -1.5, 2.5],
                                 [-0.7, 1.2, 0.3, -2.1], [2.2, -0.4, 1.7, 0.6]))
            for s, v in zip(aa, r)]
    return Potential(name="toy4", matrix=pl.DataFrame(rows), alphabet=aa)


def _toy_complex_map() -> ContactMap:
    """The toy map plus a peptide:MHC contact, so both peptide-bearing interfaces score.

    `_toy_contact_map` has an EMPTY peptide:MHC interface, which is fine for the receptor-only
    identities above but makes any presentation assertion vacuously true. Here MHC groove residue
    'K' contacts peptide position 1 -- the position the TCR does not touch -- so a mutation there
    moves the presentation term and nothing else. Note the peptide is on the ``from`` side of a
    peptide:MHC row, which is how `ContactMap.interface` selects it.
    """
    c = _toy_contact_map().contacts
    mhc = pl.DataFrame({
        "chain.type.from": ["PEPTIDE"], "chain.type.to": ["MHCa"],
        "residue.aa.from": ["G"], "residue.aa.to": ["K"],
        "region.type.from": ["PEPTIDE"], "residue.index.from": [1],
        "residue.index.to": [45], "region.start.from": [0], "region.start.to": [40],
        "pdb.id": ["toy"],
    })
    return ContactMap(pdb_id="toy", contacts=pl.concat([c, mhc.select(c.columns)]),
                      peptide_length=3)

def test_ddg_native_vs_native_is_zero():
    cm, pot = _toy_contact_map(), _toy_potential()
    assert ddg(cm, "AGK", "AGK", pot) == 0.0


def test_ddg_matches_independent_two_calls():
    cm, pot = _toy_contact_map(), _toy_potential()
    native, mutant = "AGK", "AGA"
    e_native = float(score_peptides(cm, [native], pot)["score"][0])
    e_mutant = float(score_peptides(cm, [mutant], pot)["score"][0])
    assert ddg(cm, native, mutant, pot) == pytest.approx(e_native - e_mutant)


def test_reference_delta_is_ddg_to_polyalanine():
    cm, pot = _toy_contact_map(), _toy_potential()
    assert reference_delta(cm, "AGK", pot) == pytest.approx(ddg(cm, "AGK", "AAA", pot))


def test_reference_delta_equals_alanine_scan_sum():
    # ΔΦ = Φ(real) − Φ(polyAla) equals the sum of the per-position native→Ala ΔΔGs (Φ is a contact sum).
    cm, pot = _toy_contact_map(), _toy_potential()
    scan_sum = float(alanine_scan(cm, "AGK", pot)["ddG"].sum())
    assert reference_delta(cm, "AGK", pot) == pytest.approx(scan_sum)


def test_reference_delta_is_constant_offset_on_fixed_map():
    # On ONE contact map, ΔΦ(p) = Φ(p) − const, so ΔΦ(p1) − ΔΦ(p2) == Φ(p1) − Φ(p2) (ranking-invariant).
    cm, pot = _toy_contact_map(), _toy_potential()
    d1, d2 = reference_delta(cm, "AGK", pot), reference_delta(cm, "LGA", pot)
    phi1 = float(score_peptides(cm, ["AGK"], pot)["score"][0])
    phi2 = float(score_peptides(cm, ["LGA"], pot)["score"][0])
    assert (d1 - d2) == pytest.approx(phi1 - phi2)


def test_reference_delta_zero_for_non_peptide_interface():
    cm, pot = _toy_contact_map(), _toy_potential()
    assert reference_delta(cm, "AGK", pot, interface="tcr_mhc") == 0.0


def test_ddg_sign_and_value():
    cm, pot = _toy_contact_map(), _toy_potential()
    # native "AGK": (A,A)=1.0 + (L,K)=3.0 = 4.0
    # mutant "AGA": (A,A)=1.0 + (L,A)=0.5 = 1.5  -> ddG = 4.0 - 1.5 = 2.5 (destabilising)
    assert ddg(cm, "AGK", "AGA", pot) == pytest.approx(2.5)


def test_alanine_scan_one_row_per_position():
    cm, pot = _toy_contact_map(), _toy_potential()
    native = "AGK"
    scan = alanine_scan(cm, native, pot)
    assert scan.columns == ["pos", "wt_aa", "ddG"]
    assert scan.height == len(native)
    assert scan["pos"].to_list() == [0, 1, 2]
    # wt_aa column reproduces the native peptide.
    assert scan["wt_aa"].to_list() == list(native)
    # pos0 is already 'A' -> mutating to Ala is a no-op -> ddG 0.
    assert scan.filter(pl.col("pos") == 0)["ddG"][0] == pytest.approx(0.0)
    # pos1 'G' has no TCR contact -> no-op -> ddG 0.
    assert scan.filter(pl.col("pos") == 1)["ddG"][0] == pytest.approx(0.0)
    # pos2 'K'->'A': contributes (L,K)=3.0 native vs (L,A)=0.5 -> ddG = 2.5.
    assert scan.filter(pl.col("pos") == 2)["ddG"][0] == pytest.approx(2.5)


def test_alanine_scan_position_matches_independent_ddg():
    cm, pot = _toy_contact_map(), _toy_potential()
    native = "AGK"
    scan = alanine_scan(cm, native, pot)
    for pos in range(len(native)):
        mutant = native[:pos] + "A" + native[pos + 1:]
        expected = ddg(cm, native, mutant, pot)
        assert scan.filter(pl.col("pos") == pos)["ddG"][0] == pytest.approx(expected)


def _toy_tcr_mhc_contact_map() -> ContactMap:
    # TCR 'A' contacts MHC within-region pos 0; TCR 'L' contacts MHC pos 2.
    # The peptide is NOT part of this interface, so a peptide mutation must not
    # change any score. residue.aa.to are MHC residues ('K', 'A'); were the
    # candidate peptide threaded onto pos.to (the bug), it would substitute these.
    contacts = pl.DataFrame(
        {
            "chain.type.from": ["TRA", "TRB"],
            "chain.type.to": ["MHCa", "MHCa"],
            "residue.aa.from": ["A", "L"],
            "residue.aa.to": ["K", "A"],
            "region.type.from": ["CDR3", "CDR3"],
            "residue.index.from": [10, 20],
            "residue.index.to": [0, 2],
            "region.start.from": [8, 18],
            "region.start.to": [0, 0],
            "pdb.id": ["toy", "toy"],
        }
    )
    return ContactMap(pdb_id="toy", contacts=contacts, peptide_length=3)


def test_ddg_tcr_mhc_is_zero():
    # A peptide mutation cannot affect the TCR-MHC interface (no peptide on it).
    cm, pot = _toy_tcr_mhc_contact_map(), _toy_potential()
    assert ddg(cm, "AGK", "AGA", pot, interface="tcr_mhc") == 0.0
    # Same when the peptide is unchanged.
    assert ddg(cm, "AGK", "AGK", pot, interface="tcr_mhc") == 0.0


def test_alanine_scan_tcr_mhc_all_zero():
    # Every per-position ΔΔG must be exactly 0 for an interface without the peptide.
    cm, pot = _toy_tcr_mhc_contact_map(), _toy_potential()
    native = "AGK"
    scan = alanine_scan(cm, native, pot, interface="tcr_mhc")
    assert scan.columns == ["pos", "wt_aa", "ddG"]
    assert scan.height == len(native)
    assert scan["pos"].to_list() == [0, 1, 2]
    assert scan["wt_aa"].to_list() == list(native)
    assert scan["ddG"].to_list() == [0.0, 0.0, 0.0]


def test_neoantigen_ddg_tcr_mhc_all_zero():
    cm, pot = _toy_tcr_mhc_contact_map(), _toy_potential()
    df = neoantigen_ddg(cm, "AGK", ["AGA", "KKK"], pot, interface="tcr_mhc")
    assert df["ddG"].to_list() == [0.0, 0.0]


def test_neoantigen_ddg():
    cm, pot = _toy_contact_map(), _toy_potential()
    native = "AGK"
    mutants = ["AGA", "AGK"]
    df = neoantigen_ddg(cm, native, mutants, pot)
    assert df.columns == ["native", "mutant", "ddG"]
    assert df["native"].to_list() == [native, native]
    assert df["mutant"].to_list() == mutants
    assert df.filter(pl.col("mutant") == "AGA")["ddG"][0] == pytest.approx(2.5)
    assert df.filter(pl.col("mutant") == "AGK")["ddG"][0] == pytest.approx(0.0)


def test_ddg_contact_weights_scale_the_virtual_path():
    """`weights` is what lets a contact PROBABILITY replace the map's 0/1 indicator, which is how
    a Potts `p_model` scan enters. Two contracts: all-ones is byte-identical to the default, and a
    uniform w scales ddG by exactly w, because the energy is linear in the per-contact weight."""
    import numpy as np

    cm, pot = _toy_contact_map(), _toy_potential()
    native, mutant = "AGK", "AGA"
    n_rows = cm.interface("tcr_peptide").height

    base = ddg(cm, native, mutant, pot)
    assert base != 0.0, "the toy mutation has to move the energy for this to test anything"
    assert ddg(cm, native, mutant, pot, weights=np.ones(n_rows)) == pytest.approx(base)
    for w in (0.5, 2.0):
        got = ddg(cm, native, mutant, pot, weights=np.full(n_rows, w))
        assert got == pytest.approx(w * base), (w, got, base)


def test_response_matrix_tcr_weights_scale_only_the_receptor_interface():
    """The same linearity contract as `ddg`, one level up, and confined to one interface.

    `response_matrix` is what the CPL reconstruction calls, so a Potts `p_model` scan enters the
    published panel through here. Three contracts: `None` and all-ones are byte-identical to the
    default, a uniform w scales the TCR:peptide term by exactly w, and the presentation term does
    not move at all -- the shipped Potts model is fitted on TCR:peptide and has no business
    reweighting the groove.
    """
    import numpy as np

    from tcren import response_matrix

    cm, pot = _toy_contact_map(), _toy_potential()
    n_rows = cm.interface("tcr_peptide").height

    base = response_matrix(cm, "AGK", tcr_potential=pot, mhc_potential=pot)
    assert np.any(base.phi_tcr != 0.0), "the toy map has to score for this to test anything"

    ones = response_matrix(cm, "AGK", tcr_potential=pot, mhc_potential=pot,
                           tcr_weights=np.ones(n_rows))
    assert np.allclose(ones.phi_tcr, base.phi_tcr, equal_nan=True)
    assert np.allclose(ones.phi_mhc, base.phi_mhc, equal_nan=True)

    for w in (0.5, 2.0):
        got = response_matrix(cm, "AGK", tcr_potential=pot, mhc_potential=pot,
                              tcr_weights=np.full(n_rows, w))
        assert np.allclose(got.phi_tcr, w * base.phi_tcr, equal_nan=True), w
        assert np.allclose(got.phi_mhc, base.phi_mhc, equal_nan=True), w


def test_complex_interface_sums_both_peptide_bearing_interfaces():
    """`interface="complex"` is TCR:peptide + peptide:MHC, each with its own potential.

    This is the contract that makes a LIBRARY ranking match a response-matrix cell. A cell has
    summed both interfaces since `response_matrix` existed, but a whole peptide could only be
    scored one interface at a time, so a combinatorial-library ROC silently saw the receptor term
    alone -- blind to a destroyed MHC anchor, which is the commonest reason a library peptide is
    inactive. Note the two effects are NOT separable in a library that varies every position;
    reporting them apart is the most the score can do.
    """
    cm, pot = _toy_complex_map(), _toy_complex_potential()
    tcr = ddg(cm, "AGK", "AAA", pot, interface="tcr_peptide")
    mhc = ddg(cm, "AGK", "AAA", pot, interface="peptide_mhc")
    both = ddg(cm, "AGK", "AAA", pot, interface="complex", mhc_potential=pot)
    assert mhc != pytest.approx(0.0), "the presentation term has to move for this to test anything"
    assert both == pytest.approx(tcr + mhc)
    assert both != pytest.approx(tcr), "scoring the receptor alone must not equal the complex"


def test_complex_interface_defaults_the_presentation_potential_to_mj():
    """Leaving `mhc_potential` unset uses Miyazawa-Jernigan, matching `cpl.response_matrix`."""
    from tcren.potential import mj

    cm, pot = _toy_complex_map(), _toy_complex_potential()
    assert ddg(cm, "AGK", "AAA", pot, interface="complex") == pytest.approx(
        ddg(cm, "AGK", "AAA", pot, interface="complex", mhc_potential=mj()))


def test_complex_weights_reach_only_the_receptor_channel():
    """`weights` reweights TCR:peptide alone, exactly as `response_matrix`'s `tcr_weights` does."""
    import numpy as np

    cm, pot = _toy_complex_map(), _toy_complex_potential()
    n = cm.interface("tcr_peptide").height
    tcr = ddg(cm, "AGK", "AAA", pot, interface="tcr_peptide")
    mhc = ddg(cm, "AGK", "AAA", pot, interface="peptide_mhc")
    got = ddg(cm, "AGK", "AAA", pot, interface="complex", mhc_potential=pot,
              weights=np.full(n, 2.0))
    assert got == pytest.approx(2.0 * tcr + mhc)


def test_reference_delta_complex_is_the_whole_complex_dphi():
    """The poly-alanine ΔΦ the CPL per-clone table reports is the sum of the two channels."""
    cm, pot = _toy_complex_map(), _toy_complex_potential()
    tcr = reference_delta(cm, "AGK", pot, interface="tcr_peptide")
    mhc = reference_delta(cm, "AGK", pot, interface="peptide_mhc")
    assert reference_delta(cm, "AGK", pot, interface="complex",
                           mhc_potential=pot) == pytest.approx(tcr + mhc)
