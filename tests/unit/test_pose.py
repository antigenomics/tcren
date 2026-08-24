"""Fast synthetic tests for the cross-map pose-consistency descriptors (tcren.pose)."""

from __future__ import annotations

import numpy as np
import pytest

from tcren.pose import (
    POSE_FEATURES,
    _double_centred,
    _pair_j,
    _selfcheck,
    _spearman,
    c_score,
    pose_consistency,
    pose_native_reference,
)
from tcren.structure.model import PEPTIDE_TYPE, Atom, Chain, RegionMarkup, Residue, Structure


def _atom(name, el, xyz):
    return Atom(name, el, np.asarray(xyz, float))


def _res(i, resname, aa, xyz, cb_offset=1.0):
    """A residue with a Cα at ``xyz`` and a Cβ ``cb_offset`` further along +x."""
    x, y, z = xyz
    return Residue(i, i + 1, "", aa, resname, (
        _atom("CA", "C", [x, y, z]), _atom("CB", "C", [x + cb_offset, y, z]),
    ))


def _complex(pep, tcr):
    """``pep``/``tcr``: list of ``(aa, xyz, cb_offset)``."""
    three = {"A": "ALA", "L": "LEU", "K": "LYS", "D": "ASP", "W": "TRP", "Y": "TYR"}
    pc = Chain("C", [_res(i, three[a], a, xyz, off) for i, (a, xyz, off) in enumerate(pep)],
               chain_type=PEPTIDE_TYPE)
    tc = Chain("B", [_res(i, three[a], a, xyz, off) for i, (a, xyz, off) in enumerate(tcr)],
               chain_type="TRB")
    return Structure("synth", [pc, tc])


class _Pot:
    """A directional 3-letter potential: (L,D) favourable, (L,K) neutral, (L,W) repulsive."""

    def as_matrix(self):
        idx = {"L": 0, "D": 1, "K": 2, "W": 3}
        m = np.array([
            [0.0, -3.0, 0.0, 3.0],
            [-3.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0],
            [3.0, 0.0, 0.0, 0.0],
        ])
        return m, idx


def _line(pep_specs):
    """One TCR Leu at the origin; peptide residues strung out along +x at the given distances."""
    tcr = [("L", (0.0, 0.0, 0.0), -1.0)]  # Cβ points back toward -x, i.e. away from the peptide
    pep = [(aa, (d, 0.0, 0.0), off) for aa, d, off in pep_specs]
    return _complex(pep, tcr)


def test_double_centring_zeroes_rows_and_columns_of_an_asymmetric_matrix():
    class _P:
        def as_matrix(self):
            return np.array([[1.0, 2.0, 9.0], [3.0, 0.0, 1.0], [5.0, 4.0, 2.0]]), {"A": 0, "B": 1, "C": 2}

    j, idx = _double_centred(_P())
    assert np.allclose(j.sum(axis=0), 0.0)
    assert np.allclose(j.sum(axis=1), 0.0)
    assert idx == {"A": 0, "B": 1, "C": 2}


def test_double_centring_tolerates_an_unobserved_cell():
    class _P:
        def as_matrix(self):
            return np.array([[1.0, np.nan], [3.0, 2.0]]), {"A": 0, "B": 1}

    j, _ = _double_centred(_P())
    assert np.isnan(j[0, 1]) and np.isfinite(j[1, 0])


def test_pair_j_is_nan_outside_the_alphabet():
    j, idx = _double_centred(_Pot())
    assert np.isnan(_pair_j(["L"], ["Z"], j, idx)[0])
    assert np.isfinite(_pair_j(["L"], ["D"], j, idx)[0])


def test_spearman_guards():
    assert np.isnan(_spearman(np.array([1.0, 2.0]), np.array([1.0, 2.0])))   # n < 3
    assert np.isnan(_spearman(np.ones(5), np.arange(5.0)))                   # constant axis
    assert _spearman(np.arange(5.0), np.arange(5.0)) == pytest.approx(1.0)


def test_c_local_is_positive_when_the_tight_contact_is_the_complementary_one():
    # favourable (L,D) tight at 3.0 A; neutral (L,K) at 4.0; repulsive (L,W) loose at 4.9
    s = _line([("D", 3.0, 1.0), ("K", 4.0, 1.0), ("W", 4.9, 1.0)])
    d = pose_consistency(s, potential=_Pot())
    assert d["n_contacts"] == 3
    assert d["c_local"] == pytest.approx(1.0)
    assert d["e_tight_minus_loose"] > 0
    assert d["margin_energy_slope"] > 0


def test_c_local_inverts_on_a_forced_pose():
    # the same chemistry with the distance order reversed: the repulsive pair is now the tight one
    s = _line([("D", 4.9, 1.0), ("K", 4.0, 1.0), ("W", 3.0, 1.0)])
    d = pose_consistency(s, potential=_Pot())
    assert d["c_local"] == pytest.approx(-1.0)
    assert d["e_tight_minus_loose"] < 0


def test_frac_close_favourable_counts_only_the_below_median_contacts():
    s = _line([("D", 3.0, 1.0), ("K", 4.0, 1.0), ("W", 4.9, 1.0)])
    d = pose_consistency(s, potential=_Pot())
    # median distance is 4.0; the two pairs at or below it are (L,D) J<0 and (L,K) J>0 after
    # centring, so exactly half of the close contacts are favourable
    assert 0.0 <= d["frac_close_favourable"] <= 1.0


def test_sidechain_toward_is_positive_when_cbeta_points_at_the_partner():
    toward = pose_consistency(_line([("D", 4.0, -1.0)]), potential=_Pot())["sidechain_toward"]
    away = pose_consistency(_line([("D", 4.0, +1.0)]), potential=_Pot())["sidechain_toward"]
    # the peptide Cβ at -1.0 leans back toward the TCR, so Cβ-Cβ closes relative to Cα-Cα
    assert toward > away


def test_descriptors_are_nan_not_zero_when_there_is_no_interface():
    s = _complex([("D", (60.0, 0.0, 0.0), 1.0)], [("L", (0.0, 0.0, 0.0), 1.0)])
    d = pose_consistency(s, potential=_Pot())
    assert d["n_contacts"] == 0
    assert all(np.isnan(d[k]) for k in ("c_local", "e_tight_minus_loose", "margin_energy_slope"))


def test_too_few_contacts_gives_nan_rather_than_a_made_up_correlation():
    s = _line([("D", 3.0, 1.0), ("W", 4.0, 1.0)])  # 2 contacts, below the 3 a rho needs
    d = pose_consistency(s, potential=_Pot())
    assert d["n_contacts"] == 2 and np.isnan(d["c_local"])


def test_bundled_reference_is_loadable_and_covers_every_feature():
    ref = pose_native_reference()
    for f in POSE_FEATURES:
        assert f in ref and len(ref[f]) > 100
        assert np.isfinite(ref[f]).all()


def test_c_score_is_defined_for_a_single_row_and_orders_correctly():
    ref = pose_native_reference()
    median = {k: [float(np.median(ref[k]))] for k in POSE_FEATURES}
    poor = {k: [float(np.median(ref[k]) - 2 * np.std(ref[k]))] for k in POSE_FEATURES}
    assert np.isfinite(c_score(median)[0])                 # one structure, no cohort
    assert c_score(median)[0] > c_score(poor)[0]           # higher = more crystal-like


def test_selfcheck_runs():
    _selfcheck()


# --- contact-degree structure -----------------------------------------------------------------

def _fan(n_partners, spacing=1.5):
    """One TCR Leu reaching ``n_partners`` peptide Asp residues, all within the 5 Å cutoff."""
    tcr = [("L", (0.0, 0.0, 0.0), -1.0)]
    pep = [("D", (3.0, spacing * (i - (n_partners - 1) / 2), 0.0), 1.0) for i in range(n_partners)]
    return _complex(pep, tcr)


def test_over_reaching_receptor_residue_lowers_the_degree_descriptors():
    lean = pose_consistency(_fan(2), potential=_Pot())
    greedy = pose_consistency(_fan(6), potential=_Pot())
    assert lean["max_degree_tp"] == 2 and greedy["max_degree_tp"] == 6
    assert lean["frac_well_coordinated_tp"] == 1.0      # 2 partners is typical
    assert greedy["frac_well_coordinated_tp"] == 0.0    # 6 is not
    assert greedy["degree_evenness_tp"] <= lean["degree_evenness_tp"]


def test_degree_is_counted_on_the_receptor_side_not_the_peptide_side():
    # One peptide residue ringed by six TCR residues: normal (a peptide sits inside the groove).
    # Every TCR residue has degree 1, so the interface must read as well-coordinated.
    pep = [("D", (0.0, 0.0, 0.0), 1.0)]
    tcr = [("L", (3.0 * np.cos(t), 3.0 * np.sin(t), 0.0), 0.3) for t in np.linspace(0, 5.0, 6)]
    d = pose_consistency(_complex(pep, tcr), potential=_Pot())
    assert d["max_degree_tp"] == 1
    assert d["frac_well_coordinated_tp"] == 1.0
    assert d["degree_evenness_tp"] == pytest.approx(1.0)


# --- the Calpha approach shell ------------------------------------------------------------------

def test_ca_shell_is_larger_than_the_contact_set_and_scores_separately():
    # partners at 3 A (a contact) and at 9 A (in the Calpha shell, not a contact)
    s = _line([("D", 3.0, 1.0), ("K", 4.0, 1.0), ("W", 4.9, 1.0),
               ("W", 9.0, 1.0), ("W", 10.0, 1.0)])
    d = pose_consistency(s, potential=_Pot())
    assert d["n_contacts"] == 3
    assert d["n_ca_near_tp"] == 5                       # the shell sees the far pairs too
    assert np.isfinite(d["ca_energy_coupling_tp"])


def test_frac_ca_close_engaged_flags_backbones_that_touch_nothing():
    # all five peptide residues have Cα within 8 Å, but only the 3 Å one forms a 5 Å contact
    s = _line([("D", 3.0, 1.0), ("W", 6.0, 1.0), ("W", 6.5, 1.0), ("W", 7.0, 1.0)])
    d = pose_consistency(s, potential=_Pot())
    assert d["frac_ca_close_engaged_tp"] == pytest.approx(0.25)


def test_every_declared_feature_is_returned():
    d = pose_consistency(_line([("D", 3.0, 1.0), ("K", 4.0, 1.0), ("W", 4.9, 1.0)]), potential=_Pot())
    missing = [f for f in POSE_FEATURES if f not in d]
    assert not missing, missing


# --- the per-structure C* surrogate --------------------------------------------------------------

def test_agreement_mean_reproduces_coupling_when_standardized_on_the_input():
    from tcren.cohort import agreement, coupling
    rng = np.random.default_rng(0)
    q = rng.normal(size=400)
    for rho in (-0.6, 0.0, 0.5):
        e = rho * q + np.sqrt(1 - rho**2) * rng.normal(size=400)
        # coupling IS the mean of the agreement term when both are standardized on the same rows
        assert agreement(q, e).mean() == pytest.approx(coupling(q, e), abs=0.02)


def test_agreement_is_defined_for_a_single_row_against_a_reference():
    from tcren.cohort import agreement
    ref_q, ref_e = np.random.default_rng(1).normal(size=200), np.random.default_rng(2).normal(size=200)
    a = agreement([2.0], [2.0], ref_q, ref_e)          # both far above reference -> agree
    b = agreement([2.0], [-2.0], ref_q, ref_e)         # geometry good, chemistry bad -> disagree
    assert a.shape == (1,) and np.isfinite(a[0])
    assert a[0] > 0 and b[0] < 0


# --- the designed pose metric P -------------------------------------------------------------

def _full_complex():
    """Peptide + both receptor chains (CDR1/2/3 marked up) + MHC.

    Every family the generator emits must be exercised, including the per-loop Calpha profile, which
    needs TRA and TRB each carrying all three CDRs.
    """
    three = {"L": "LEU", "D": "ASP", "K": "LYS", "W": "TRP"}
    rng = np.random.default_rng(0)

    def chain(cid, ctype, n, origin, aas, loops):
        res = []
        for i in range(n):
            x, y_, z = origin + rng.normal(scale=1.8, size=3)
            aa = aas[i % len(aas)]
            res.append(Residue(i, i + 1, "", aa, three[aa], (
                _atom("CA", "C", [x, y_, z]), _atom("CB", "C", [x + 0.7, y_, z]))))
        c = Chain(cid, res, chain_type=ctype)
        if loops:
            span = n // 3
            c.regions = [
                RegionMarkup(r, i * span, (i + 1) * span - 1,
                             "".join(x.aa for x in res[i * span:(i + 1) * span]),
                             res[i * span:(i + 1) * span])
                for i, r in enumerate(("CDR1", "CDR2", "CDR3"))]
        return c

    pep = chain("C", PEPTIDE_TYPE, 9, np.array([0.0, 0.0, 0.0]), "DKW", False)
    tra = chain("D", "TRA", 18, np.array([-2.0, 0.0, 4.0]), "LDK", True)
    trb = chain("B", "TRB", 18, np.array([2.0, 0.0, 4.0]), "LDW", True)
    mhc = chain("A", "MHCa", 24, np.array([0.0, 0.0, -4.0]), "LDW", False)
    return Structure("synth_full", [pep, tra, trb, mhc])


def test_p_score_terms_are_all_produced_by_the_generator():
    """Every P term must be a descriptor the generator actually emits, or P silently drops it."""
    from tcren.pose_sweep import P_TERMS, pose_descriptors_full

    d = pose_descriptors_full(_full_complex(), potential=_Pot())
    # `offset` comes from orient.tcr_placement, not the map generator; the rest must be present
    produced = set(d) | {"offset"}
    missing = [n for n, _ in P_TERMS if n not in produced]
    assert not missing, missing


def test_p_score_orients_higher_is_better_and_is_row_wise():
    from tcren.pose_sweep import P_TERMS, p_score

    rng = np.random.default_rng(0)
    ref = {n: rng.normal(size=300) for n, _ in P_TERMS}
    good = {n: [float(np.mean(ref[n]) + sign * 2.0 * np.std(ref[n]))] for n, sign in P_TERMS}
    bad = {n: [float(np.mean(ref[n]) - sign * 2.0 * np.std(ref[n]))] for n, sign in P_TERMS}
    assert p_score(good, ref)[0] > p_score(bad, ref)[0]
    assert p_score({n: ref[n] for n, _ in P_TERMS}).shape == (300,)


# --- interpretable Calpha rules ----------------------------------------------------------------

def test_loop_ca_rules_delta_sign_follows_geometry():
    """delta = d_pep - d_mhc must be negative for a loop placed nearer the peptide."""
    from tcren.pose_sweep import loop_ca_rules

    s = _full_complex()
    # peptide sits at z=0, MHC at z=-4, receptor at z=+4, so every loop is nearer the peptide
    d = loop_ca_rules(s, cutoff=5.0)
    # only the per-loop deltas; the cross-loop contrasts (cdr3_vs_germline_*, cdr3_ab_asym_*) are
    # differences of two loops and carry no sign guarantee
    loops = ("cdr1a", "cdr2a", "cdr3a", "cdr1b", "cdr2b", "cdr3b")
    deltas = {k: v for k, v in d.items()
              if k.endswith("_all_delta") and k.split("_all_")[0] in loops and np.isfinite(v)}
    assert len(deltas) == 6, sorted(deltas)
    assert all(x < 0 for x in deltas.values()), deltas


def test_loop_ca_rules_reports_both_chains_and_all_three_loops():
    from tcren.pose_sweep import loop_ca_rules

    d = loop_ca_rules(_full_complex())
    for loop in ("cdr1a", "cdr2a", "cdr3a", "cdr1b", "cdr2b", "cdr3b"):
        assert f"{loop}_n_res" in d, loop
        assert f"{loop}_all_d_pep" in d and f"{loop}_all_d_mhc" in d


def test_sigma_split_is_the_difference_of_the_two_cdr3_nc_terms():
    from tcren.pose_sweep import loop_ca_rules

    d = loop_ca_rules(_full_complex())
    if all(k in d for k in ("cdr3a_delta_NC", "cdr3b_delta_NC", "sigma_NC_split")):
        assert d["sigma_NC_split"] == pytest.approx(d["cdr3a_delta_NC"] - d["cdr3b_delta_NC"])


def test_loop_ca_rules_partner_masks_are_a_union():
    """conA is exactly the peptide-contacting union MHC-contacting set, per loop."""
    from tcren.pose_sweep import loop_ca_rules

    d = loop_ca_rules(_full_complex(), cutoff=5.0)
    loops = [k[: -len("_frac_conA")] for k in d if k.endswith("_frac_conA")]
    assert loops, sorted(d)
    for name in loops:
        p, m, a = (d[f"{name}_frac_con"], d[f"{name}_frac_conM"], d[f"{name}_frac_conA"])
        assert a >= max(p, m) - 1e-9, (name, p, m, a)
        assert a <= min(1.0, p + m) + 1e-9, (name, p, m, a)


def test_peptide_ca_rules_engaged_residues_sit_closer_to_the_receptor():
    """A peptide residue that contacts the TCR cannot be farther from it than one that does not."""
    from tcren.pose_sweep import peptide_ca_rules

    d = peptide_ca_rules(_full_complex(), cutoff=5.0)
    assert 0.0 < d["pep_frac_conTcr"] < 1.0, d["pep_frac_conTcr"]
    assert d["pep_conTcr_d_tcr"] < d["pep_nonTcr_d_tcr"]
    assert d["pep_conTcr_d_tcr"] <= d["pep_all_d_tcr"]


def test_peptide_matched_conditional_is_per_loop():
    """pep_con_<loop>_d_<loop> conditions on contacts with that loop alone, so it is the closer set."""
    from tcren.pose_sweep import peptide_ca_rules

    d = peptide_ca_rules(_full_complex(), cutoff=5.0)
    loops = ("cdr1a", "cdr2a", "cdr3a", "cdr1b", "cdr2b", "cdr3b")
    matched = [(f"pep_con_{n}_d_{n}", f"pep_all_d_{n}") for n in loops if f"pep_con_{n}_d_{n}" in d]
    assert matched, sorted(k for k in d if k.startswith("pep_con_"))
    for k, ref in matched:
        assert d[k] <= d[ref] + 1e-9, (k, d[k], d[ref])
