"""DFIRE reference states and the corrections they imply for a contact potential.

The interesting properties are the ones a sign error or a mis-binned reference would break:
that the radial energy is zero at the cutoff by construction, that the rotation term is a
divergence and so cannot be positive, that an isotropic pair earns nothing from it, and that
adding the corrections to a potential is a pure per-cell shift.
"""

from __future__ import annotations

import numpy as np
import polars as pl
import pytest

from tcren.potential import Potential, apply_corrections, corrections, radial_potential
from tcren.potential.dfire import ALPHA, RC, _rotation_term, select_scope


def _geom(rows: list[dict]) -> pl.DataFrame:
    """A pair-geometry frame from partial rows; unset fields take neutral defaults."""
    base = {"chain.type.from": "TRA", "chain.type.to": "PEPTIDE",
            "residue.aa.from": "A", "residue.aa.to": "L",
            "dist": 4.0, "cos.from": 0.0, "cos.to": 0.0, "contact": True}
    return pl.DataFrame([base | r for r in rows])


def _isotropic(aa_from: str, aa_to: str, per_cell: int) -> list[dict]:
    """Contacts spread evenly over the nine orientation cells — the reference distribution."""
    centres = (-2 / 3, 0.0, 2 / 3)
    return [{"residue.aa.from": aa_from, "residue.aa.to": aa_to, "cos.from": a, "cos.to": b}
            for a in centres for b in centres for _ in range(per_cell)]


# --- the radial reference -------------------------------------------------------------


def test_energy_is_zero_in_the_reference_bin():
    """``u(a, b, r_c) = 0`` by construction: the last bin *is* the reference."""
    rng = np.random.default_rng(0)
    g = _geom([{"dist": d} for d in rng.uniform(2.0, RC, 4000)])
    last = radial_potential(g).filter(
        (pl.col("residue.aa.from") == "A") & (pl.col("residue.aa.to") == "L")
    ).sort("r").tail(1)
    assert last["u"].item() == pytest.approx(0.0, abs=1e-12)


def test_a_population_following_the_reference_has_flat_zero_energy():
    """Sampling :math:`p(r) \\propto r^{\\alpha}` must return no structure at any distance."""
    rng = np.random.default_rng(1)
    r = RC * rng.uniform(0, 1, 400_000) ** (1 / (ALPHA + 1))
    u = (
        radial_potential(_geom([{"dist": float(d)} for d in r]))
        .filter((pl.col("residue.aa.from") == "A") & (pl.col("residue.aa.to") == "L"))
        .filter(pl.col("r") > 4.0)["u"].to_numpy()
    )
    assert np.abs(u).max() < 0.15


def test_a_pair_that_packs_close_is_favourable_at_short_range():
    """Excess occupancy below the reference curve must come out as negative energy."""
    rng = np.random.default_rng(2)
    far = rng.uniform(2.0, RC, 20_000)
    near = rng.uniform(3.0, 5.0, 20_000)
    rad = radial_potential(_geom([{"dist": float(d)} for d in np.concatenate([far, near])]))
    short = rad.filter(
        (pl.col("residue.aa.from") == "A") & (pl.col("residue.aa.to") == "L")
        & (pl.col("r") > 3.0) & (pl.col("r") < 5.0)
    )
    assert short["u"].max() < 0.0


# --- the rotation term ----------------------------------------------------------------


def test_rotation_term_is_never_positive():
    """It is minus a Kullback-Leibler divergence, so a positive value is a sign error."""
    rng = np.random.default_rng(3)
    rows = [{"residue.aa.from": a, "residue.aa.to": b,
             "cos.from": float(rng.uniform(-1, 1)), "cos.to": float(rng.uniform(-1, 1))}
            for a in "ALVK" for b in "ALVK" for _ in range(400)]
    assert _rotation_term(_geom(rows))["C_rot"].max() <= 0.0


def test_an_isotropic_pair_earns_nothing():
    got = _rotation_term(_geom(_isotropic("A", "L", 300)))["C_rot"].item()  # 2,700 contacts
    assert got == pytest.approx(0.0, abs=1e-9)


def test_a_pair_locked_in_one_orientation_earns_ln_nine():
    """Nine cells collapsed onto one is ``ln 9`` nats, less the Miller-Madow term."""
    n = 5000
    rows = [{"cos.from": 0.9, "cos.to": 0.9}] * n  # far above the count floor
    got = _rotation_term(_geom(rows))["C_rot"].item()
    assert got == pytest.approx(-(np.log(9) - 8 / (2 * n)), abs=1e-9)


def test_thin_cells_are_not_credited_at_all():
    """Nine samples over nine cells look structured, and Miller-Madow alone still believes it.

    At that count the null's own 99th percentile is 0.60 nats, so the count floor — not the
    bias term — is what keeps a thin cell out.
    """
    rng = np.random.default_rng(4)
    rows = [{"cos.from": float(rng.uniform(-1, 1)), "cos.to": float(rng.uniform(-1, 1))}
            for _ in range(9)]
    assert _rotation_term(_geom(rows))["C_rot"].item() == pytest.approx(0.0, abs=1e-9)


def test_the_count_floor_is_what_suppresses_a_thin_cell():
    """Same data, floor lifted: the raw estimator does credit it, which is the point."""
    rng = np.random.default_rng(4)
    rows = [{"cos.from": float(rng.uniform(-1, 1)), "cos.to": float(rng.uniform(-1, 1))}
            for _ in range(9)]
    assert _rotation_term(_geom(rows), min_oriented=0)["C_rot"].item() < -0.2


def test_glycine_has_no_direction_and_is_excluded():
    """Cα→Cβ is undefined for Gly, so its cosines are not finite and it contributes nothing."""
    rows = _isotropic("A", "L", 40) + [
        {"residue.aa.from": "G", "residue.aa.to": "L", "cos.from": np.nan, "cos.to": 0.5}
    ] * 50
    out = _rotation_term(_geom(rows))
    assert set(zip(out["residue.aa.from"], out["residue.aa.to"])) == {("A", "L")}


# --- the decomposition and its application ---------------------------------------------


def test_only_contacts_enter_the_corrections():
    """Pairs inside the radial cutoff but not in contact set the reference, not the energy."""
    rows = _isotropic("A", "L", 40) + [
        {"residue.aa.from": "V", "residue.aa.to": "K", "dist": 12.0, "contact": False}
    ] * 500
    t = corrections(_geom(rows)).table
    assert t.filter((pl.col("residue.aa.from") == "V") & (pl.col("residue.aa.to") == "K"))[
        "n_contacts"].item() == 0


def test_uncovered_cells_are_left_alone_by_the_correction():
    """A cell the decomposition never saw must come back unchanged, not shifted on no data."""
    pot = Potential(name="p", matrix=pl.DataFrame(
        {"residue.aa.from": ["A", "W"], "residue.aa.to": ["L", "W"], "value": [1.0, 2.0]}),
        alphabet=("A", "L", "W"))
    dec = corrections(_geom(_isotropic("A", "L", 60)))
    out = apply_corrections(pot, dec)
    assert out.matrix.filter(pl.col("residue.aa.from") == "W")["value"].item() == 2.0


def test_applying_no_terms_is_the_identity():
    pot = Potential(name="p", matrix=pl.DataFrame(
        {"residue.aa.from": ["A"], "residue.aa.to": ["L"], "value": [1.25]}),
        alphabet=("A", "L"))
    dec = corrections(_geom(_isotropic("A", "L", 60)))
    assert apply_corrections(pot, dec, terms=()).value("A", "L") == pytest.approx(1.25)


def test_the_two_terms_add_independently():
    pot = Potential(name="p", matrix=pl.DataFrame(
        {"residue.aa.from": ["A"], "residue.aa.to": ["L"], "value": [0.0]}),
        alphabet=("A", "L"))
    rng = np.random.default_rng(5)
    rows = [{"dist": float(rng.uniform(3, 5)), "cos.from": float(rng.uniform(-1, 1)),
             "cos.to": float(rng.uniform(0.4, 1))} for _ in range(3000)]
    dec = corrections(_geom(rows))
    d = apply_corrections(pot, dec, terms=("dist",)).value("A", "L")
    r = apply_corrections(pot, dec, terms=("rot",)).value("A", "L")
    both = apply_corrections(pot, dec, terms=("dist", "rot")).value("A", "L")
    assert both == pytest.approx(d + r)


def test_unknown_term_is_rejected():
    pot = Potential(name="p", matrix=pl.DataFrame(
        {"residue.aa.from": ["A"], "residue.aa.to": ["L"], "value": [0.0]}), alphabet=("A", "L"))
    with pytest.raises(ValueError, match="unknown correction"):
        apply_corrections(pot, corrections(_geom(_isotropic("A", "L", 20))), terms=("charge",))


def test_dfire2_is_the_orientation_free_term_plus_the_rotation_one():
    rng = np.random.default_rng(6)
    rows = [{"dist": float(rng.uniform(3, 6)), "cos.from": float(rng.uniform(-1, 1)),
             "cos.to": float(rng.uniform(-1, 1))} for _ in range(2000)]
    dec = corrections(_geom(rows))
    row = dec.table.filter((pl.col("residue.aa.from") == "A") & (pl.col("residue.aa.to") == "L"))
    assert dec.dfire2().value("A", "L") == pytest.approx(row["E0"].item() + row["C_rot"].item())


# --- scope selection --------------------------------------------------------------------


def test_scope_orients_a_pair_recorded_the_other_way_round():
    """A peptide→TCR row must be counted as TCR→peptide, not silently as its own cell."""
    g = _geom([
        {"chain.type.from": "TRA", "chain.type.to": "PEPTIDE",
         "residue.aa.from": "W", "residue.aa.to": "D"},
        {"chain.type.from": "PEPTIDE", "chain.type.to": "TRB",
         "residue.aa.from": "D", "residue.aa.to": "W"},
    ])
    out = select_scope(g, "tcr_peptide")
    assert out.height == 2
    assert set(out["residue.aa.from"]) == {"W"}


def test_scope_drops_the_other_interfaces():
    g = _geom([
        {"chain.type.from": "TRA", "chain.type.to": "PEPTIDE"},
        {"chain.type.from": "PEPTIDE", "chain.type.to": "MHCa"},
    ])
    assert select_scope(g, "tcr_peptide").height == 1
    assert select_scope(g, "peptide_mhc").height == 1
    assert select_scope(g, "all").height == 2


def test_unknown_scope_is_rejected():
    with pytest.raises(ValueError, match="scope must be"):
        select_scope(_geom([{}]), "tcr_tcr")
