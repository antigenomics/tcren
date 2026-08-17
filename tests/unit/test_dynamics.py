"""Unit tests for peptide conformational stability (:mod:`tcren.dynamics`)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from tcren.dynamics import backbone_torsions, peptide_stability, stability_table

PDB_DIR = Path(__file__).resolve().parents[1] / "assets" / "pdb"


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


# --- the torsion tree -----------------------------------------------------------------------------
def _chain(n_res):
    """A toy peptide's flat atom list: N, CA, C, O, CB per residue."""
    return [(r, n) for r in range(n_res) for n in ("N", "CA", "C", "O", "CB")]


def test_two_torsions_per_residue():
    tors = backbone_torsions(_chain(5))
    assert len(tors) == 10                                   # phi and psi for each of five residues


def test_phi_moves_the_side_chain_but_psi_does_not():
    """A torsional rotation splits the chain at the bond. Cβ hangs off Cα, so it is on the moving
    side of N–Cα and on the fixed side of Cα–C."""
    atoms = _chain(3)
    tors = backbone_torsions(atoms)
    phi, psi = tors[0], tors[1]                              # residue 0

    names = [f"{n}{r}" for r, n in atoms]
    assert (names[phi[0]], names[phi[1]]) == ("N0", "CA0")
    assert (names[psi[0]], names[psi[1]]) == ("CA0", "C0")
    assert "CB0" in [names[i] for i in phi[2]]
    assert "CB0" not in [names[i] for i in psi[2]]
    assert "O0" in [names[i] for i in psi[2]]


def test_every_torsion_moves_all_later_residues():
    atoms = _chain(4)
    res_of = np.array([r for r, _ in atoms])
    for start, _end, mov in backbone_torsions(atoms):
        this_res = res_of[start]
        later = set(np.flatnonzero(res_of > this_res).tolist())
        assert later <= set(mov.tolist())


def test_a_residue_missing_a_backbone_atom_defines_no_torsion():
    atoms = [(0, "N"), (0, "CA")] + [(1, n) for n in ("N", "CA", "C", "O")]
    tors = backbone_torsions(atoms)
    assert all(np.asarray([r for r, _ in atoms])[s] == 1 for s, _e, _m in tors)


def test_the_last_residue_still_has_torsions():
    """Its phi/psi move only its own atoms, but they are real degrees of freedom."""
    tors = backbone_torsions(_chain(2))
    assert len(tors) == 4


# --- the sampler ----------------------------------------------------------------------------------
@pytest.mark.slow
def test_stability_is_deterministic_given_a_seed(annotated):
    a = peptide_stability(annotated, n_steps=600, seed=7)
    b = peptide_stability(annotated, n_steps=600, seed=7)
    assert a.rmsf == b.rmsf and a.drift == b.drift and a.energy == b.energy


@pytest.mark.slow
def test_stability_fields_are_sane(annotated):
    st = peptide_stability(annotated, n_steps=800, seed=0)
    assert st.peptide == "LLFGYPVYV"
    assert st.rmsf > 0 and st.drift >= 0
    assert 0.0 < st.accept_rate <= 1.0
    assert st.n_samples == 800 - 800 // 8
    assert st.energy <= st.energy_start                      # MC keeps the best pose it ever saw
    assert st.energy_gap >= 0


@pytest.mark.slow
def test_zero_steps_leaves_the_pose_alone(annotated):
    st = peptide_stability(annotated, n_steps=0, seed=0)
    assert st.n_samples == 0
    assert st.rmsf == 0.0 and st.drift == 0.0
    assert st.energy == st.energy_start


@pytest.mark.slow
def test_a_colder_sampler_moves_less(annotated):
    """Temperature has to do what temperature does, or the readout means nothing."""
    cold = np.mean([peptide_stability(annotated, n_steps=1500, temperature=1.0, seed=s).rmsf
                    for s in range(3)])
    hot = np.mean([peptide_stability(annotated, n_steps=1500, temperature=40.0, seed=s).rmsf
                   for s in range(3)])
    assert cold < hot


@pytest.mark.slow
def test_anchor_restraint_holds_the_peptide_down(annotated):
    """Without it the peptide can leave the groove, which is not the question being asked."""
    pinned = np.mean([peptide_stability(annotated, n_steps=1500, temperature=20.0,
                                        anchor_weight=5.0, seed=s).drift for s in range(3)])
    loose = np.mean([peptide_stability(annotated, n_steps=1500, temperature=20.0,
                                       anchor_weight=0.0, seed=s).drift for s in range(3)])
    assert pinned < loose


@pytest.mark.slow
def test_dropping_the_intra_term_changes_the_energy_not_the_geometry_bookkeeping(annotated):
    with_intra = peptide_stability(annotated, n_steps=800, intra_weight=1.0, seed=0)
    without = peptide_stability(annotated, n_steps=800, intra_weight=0.0, seed=0)
    assert with_intra.energy_start != without.energy_start   # the term is actually in the sum
    assert with_intra.n_samples == without.n_samples
    assert with_intra.peptide == without.peptide


@pytest.mark.slow
def test_stability_table_reports_the_paired_difference(annotated):
    tab = stability_table([annotated], n_steps=800, seed=0)
    assert tab.height == 1
    assert {"rmsf_intra1", "rmsf_intra0", "delta_rmsf", "delta_drift"} <= set(tab.columns)
    row = tab.row(0, named=True)
    assert row["delta_rmsf"] == pytest.approx(row["rmsf_intra0"] - row["rmsf_intra1"])


def test_an_untyped_structure_is_rejected():
    from tcren.structure import parse_structure

    s = parse_structure(PDB_DIR / "1ao7.pdb")               # never chain-typed
    with pytest.raises(ValueError, match="no peptide chain"):
        peptide_stability(s)
