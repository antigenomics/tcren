"""Unit tests for rotamer-averaged contacts (:mod:`tcren.rotamers`)."""

from __future__ import annotations

import copy
from pathlib import Path

import numpy as np
import pytest

from tcren.rotamers import (
    N_CHI,
    _rotate,
    chi_axes,
    contact_probabilities,
    residue_rotamers,
    soft_energy,
)

PDB_DIR = Path(__file__).resolve().parents[1] / "assets" / "pdb"
KEYS = ["chain.id.from", "residue.index.from", "chain.id.to", "residue.index.to"]


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


# --- the torsion geometry -------------------------------------------------------------------------
def test_rotate_is_a_rigid_motion_of_the_moving_set():
    rng = np.random.default_rng(0)
    coords = rng.normal(size=(6, 3))
    moving = np.array([3, 4, 5])
    out = _rotate(coords, coords[0], coords[1], moving, 73.0)
    assert np.allclose(out[:3], coords[:3])                       # the fixed part does not move
    d_in = np.linalg.norm(coords[moving][:, None] - coords[moving][None], axis=-1)
    d_out = np.linalg.norm(out[moving][:, None] - out[moving][None], axis=-1)
    assert np.allclose(d_in, d_out)                               # internal geometry preserved
    # distance to the axis end is preserved too
    assert np.allclose(np.linalg.norm(coords[moving] - coords[1], axis=1),
                       np.linalg.norm(out[moving] - out[1], axis=1))


def test_a_full_turn_is_the_identity():
    rng = np.random.default_rng(1)
    coords = rng.normal(size=(5, 3))
    out = _rotate(coords, coords[0], coords[1], np.array([2, 3, 4]), 360.0)
    assert np.allclose(out, coords, atol=1e-9)


@pytest.mark.slow
def test_chi_axes_match_the_known_torsion_counts(annotated):
    """Rotatable-bond count per residue must equal the standard chi count."""
    seen = set()
    for chain in annotated.chains:
        for res in chain.residues:
            if res.aa in seen or res.aa not in N_CHI:
                continue
            n = len(chi_axes(res))
            if n < N_CHI[res.aa]:
                continue                                          # truncated side chain in the model
            seen.add(res.aa)
            assert n == N_CHI[res.aa], f"{res.aa}: {n} axes, expected {N_CHI[res.aa]}"
    assert len(seen) > 10, "too few residue types seen to be a real check"


@pytest.mark.slow
def test_chi1_axis_is_ca_cb_and_moves_everything_past_cb(annotated):
    for res in _peptide(annotated).residues:
        axes = chi_axes(res)
        if not axes:
            assert N_CHI.get(res.aa, 0) == 0 or "CB" not in [a.name for a in res.atoms]
            continue
        start, end, moving = axes[0]
        names = [a.name for a in res.atoms]
        assert names[start] == "CA" and names[end] == "CB"
        assert "CB" not in [names[i] for i in moving]
        assert all(names[i] not in ("N", "CA", "C", "O") for i in moving)


@pytest.mark.slow
def test_rotamer_count_is_three_to_the_number_of_sampled_chis(annotated):
    for res in _peptide(annotated).residues:
        n_chi = min(len(chi_axes(res)), 2)
        assert len(residue_rotamers(res, max_chi=2)) == 3 ** n_chi
        assert len(residue_rotamers(res, max_chi=1)) == 3 ** min(n_chi, 1)


@pytest.mark.slow
def test_the_input_conformation_is_always_the_first_rotamer(annotated):
    """A caller that keeps the best rotamer can never do worse than the pose it was given."""
    for res in _peptide(annotated).residues:
        rot = residue_rotamers(res, max_chi=2)
        assert np.allclose(rot[0], [a.coord for a in res.atoms])


@pytest.mark.slow
def test_rotamers_keep_the_backbone_fixed(annotated):
    for res in _peptide(annotated).residues:
        rot = residue_rotamers(res, max_chi=2)
        bb = [i for i, a in enumerate(res.atoms) if a.name in ("N", "CA", "C", "O", "CB")]
        assert np.allclose(rot[:, bb, :], rot[0][bb][None], atol=1e-9)


# --- contact probabilities ------------------------------------------------------------------------
@pytest.mark.slow
def test_contact_probabilities_are_probabilities(annotated):
    p = contact_probabilities(annotated, "tcr_peptide")
    assert p.height > 0
    v = p["p"].to_numpy()
    assert (v > 0).all() and (v <= 1.0 + 1e-12).all()
    assert set(KEYS) <= set(p.columns)


@pytest.mark.slow
def test_soft_map_covers_the_hard_one(annotated):
    """Every hard contact must keep some probability; the soft map may add pairs the pose missed."""
    from tcren.contactmap import ContactMap

    hard = {tuple(r) for r in
            ContactMap.from_structure(annotated).interface("tcr_peptide").select(KEYS).rows()}
    soft = {tuple(r) for r in contact_probabilities(annotated, "tcr_peptide").select(KEYS).rows()}
    assert hard <= soft
    assert len(soft) > len(hard)


@pytest.mark.slow
def test_zero_temperature_limit_keeps_the_native_pose_contacts(annotated):
    """As T -> 0 the weights collapse onto the best rotamer, so probabilities become 0/1."""
    p = contact_probabilities(annotated, "tcr_peptide", temperature=1e-3)
    v = p["p"].to_numpy()
    assert ((v < 1e-6) | (v > 1 - 1e-6)).mean() > 0.9


@pytest.mark.slow
def test_averaging_survives_a_wrong_rotamer_far_better_than_a_hard_map(annotated):
    """The measured claim: mean |dPhi| under a deliberately wrong chi1 falls ~10x."""
    from tcren.contactmap import ContactMap
    from tcren.pipeline import _interface_energy
    from tcren.potential import tcren
    from tcren.structure.model import Atom, Chain, Residue, Structure

    def perturb(structure):
        chains = []
        for chain in structure.chains:
            if chain.chain_type != "PEPTIDE":
                chains.append(chain)
                continue
            residues = []
            for res in chain.residues:
                axes = chi_axes(res)
                if not axes:
                    residues.append(res)
                    continue
                xyz = np.asarray([a.coord for a in res.atoms], float)
                a, b, moving = axes[0]
                xyz = _rotate(xyz, xyz[a], xyz[b], moving, 120.0)
                residues.append(Residue(
                    res.seq_index, res.pdb_index, res.insertion_code, res.aa, res.resname,
                    tuple(Atom(at.name, at.element, xyz[k]) for k, at in enumerate(res.atoms))))
            chains.append(Chain(chain.chain_id, residues, chain_type=chain.chain_type,
                                chain_supertype=chain.chain_supertype, regions=chain.regions))
        return Structure(structure.pdb_id, chains, complex_species=structure.complex_species,
                         cell_type=structure.cell_type)

    wrong = perturb(copy.deepcopy(annotated))
    pot = tcren()
    hard = abs(_interface_energy(ContactMap.from_structure(wrong).interface("tcr_peptide"), pot)
               - _interface_energy(ContactMap.from_structure(annotated).interface("tcr_peptide"),
                                   pot))
    soft = abs(soft_energy(wrong, pot) - soft_energy(annotated, pot))
    assert soft < hard, f"rotamer averaging did not help: soft {soft:.3f} vs hard {hard:.3f}"


@pytest.mark.slow
def test_unknown_interface_raises(annotated):
    with pytest.raises(ValueError, match="unknown interface"):
        contact_probabilities(annotated, "tcr_tcr")


# --- the native repack kernel ---------------------------------------------------------------------
@pytest.mark.slow
def test_repack_leaves_a_crystal_where_it_found_it(annotated):
    """The input conformer is index 0 of every residue's enumeration, so a crystal side chain the
    potential already prefers must come back unmoved."""
    from tcren.rotamers import repack

    out, report = repack(annotated)
    before = {(r.seq_index, a.name): a.coord
              for c in annotated.chains if c.chain_type == "PEPTIDE"
              for r in c.residues for a in r.atoms}
    after = {(r.seq_index, a.name): a.coord
             for c in out.chains if c.chain_type == "PEPTIDE" for r in c.residues for a in r.atoms}
    shift = np.array([np.linalg.norm(before[k] - after[k]) for k in before])
    assert shift.mean() < 0.5, f"crystal side chains moved by {shift.mean():.2f} A on average"
    assert (shift > 0.5).sum() <= 0.15 * len(shift)
    assert report.height == len(_peptide(annotated).residues)


@pytest.mark.slow
def test_repack_agrees_with_the_python_enumeration(annotated):
    """The C++ kernel must reproduce the energy the Python prototype computes, exactly."""
    from scipy.spatial import cKDTree

    from tcren import _relax
    from tcren.refine import _dope
    from tcren.rotamers import repack, residue_rotamers

    _out, report = repack(annotated)
    table, amap, x0, dx, nb = _dope()

    xyz, cls, own = [], [], []
    for chain in annotated.chains:
        for res in chain.residues:
            for a in res.atoms:
                if a.element == "H":
                    continue
                xyz.append(a.coord)
                cls.append(amap.get(f"{res.resname}:{a.name}", -1))
                own.append((chain.chain_id, res.seq_index))
    xyz = np.asarray(xyz)
    cls = np.asarray(cls, np.int32)
    tree = cKDTree(xyz)

    pep = next(c for c in annotated.chains if c.chain_type == "PEPTIDE")
    for res in pep.residues:
        key = (pep.chain_id, res.seq_index)
        mine = np.array([i for i, k in enumerate(own) if k == key])
        hits = tree.query_ball_point(xyz[mine], 12.0)
        near = np.unique(np.concatenate([np.asarray(h, int) for h in hits if h]))
        near = np.array([j for j in near if own[j] != key])
        heavy = [i for i, a in enumerate(res.atoms) if a.element != "H"]
        best = min(
            float(_relax.interface_energy(
                np.ascontiguousarray(c[heavy]), cls[mine], np.ascontiguousarray(xyz[near]),
                np.ascontiguousarray(cls[near]), table, table.shape[0], nb, x0, dx))
            for c in residue_rotamers(res, max_chi=2))
        got = float(report.filter(report["residue.index"] == res.seq_index)["energy"][0])
        assert got == pytest.approx(best, abs=1e-9)


@pytest.mark.slow
def test_repack_recovers_a_wrong_rotamer(annotated):
    """The claim the packer exists for: a discrete re-sample crosses the torsional barrier a local
    minimiser cannot. Measured on five crystals, side-chain RMSD 4.13 A -> 2.36 A."""
    from tcren.docking._transform import apply_rigid, kabsch
    from tcren.docking.align import _matched_anchors
    from tcren.rotamers import repack
    from tcren.structure.model import Atom, Chain, Residue, Structure

    def perturb(structure):
        chains = []
        for chain in structure.chains:
            if chain.chain_type != "PEPTIDE":
                chains.append(chain)
                continue
            residues = []
            for res in chain.residues:
                axes = chi_axes(res)
                if not axes:
                    residues.append(res)
                    continue
                xyz = np.asarray([a.coord for a in res.atoms], float)
                a, b, mov = axes[0]
                xyz = _rotate(xyz, xyz[a], xyz[b], mov, 120.0)
                residues.append(Residue(
                    res.seq_index, res.pdb_index, res.insertion_code, res.aa, res.resname,
                    tuple(Atom(at.name, at.element, xyz[k]) for k, at in enumerate(res.atoms))))
            chains.append(Chain(chain.chain_id, residues, chain_type=chain.chain_type,
                                chain_supertype=chain.chain_supertype, regions=chain.regions))
        return Structure(structure.pdb_id, chains, complex_species=structure.complex_species,
                         cell_type=structure.cell_type)

    def sc_rmsd(model, ref):
        mob, rp = _matched_anchors(model, ref)
        rot, tran, _ = kabsch(mob, rp)
        bb = {"N", "CA", "C", "O"}

        def atoms(s):
            c = next(c for c in s.chains if c.chain_type == "PEPTIDE")
            return {(r.seq_index, a.name): a.coord for r in c.residues for a in r.atoms
                    if a.element != "H" and a.name not in bb}

        m, r = atoms(model), atoms(ref)
        keys = sorted(m.keys() & r.keys())
        M = apply_rigid(np.array([m[k] for k in keys]), rot, tran)
        return float(np.sqrt(((M - np.array([r[k] for k in keys])) ** 2).sum(1).mean()))

    wrong = perturb(copy.deepcopy(annotated))
    fixed, _ = repack(wrong)
    before, after = sc_rmsd(wrong, annotated), sc_rmsd(fixed, annotated)
    assert before > 2.0, "the perturbation should actually move the side chains"
    assert after < 0.75 * before, f"repack did not recover: {before:.2f} -> {after:.2f} A"


@pytest.mark.slow
def test_repack_report_is_one_row_per_repacked_residue(annotated):
    from tcren.rotamers import repack

    _out, report = repack(annotated)
    pep = next(c for c in annotated.chains if c.chain_type == "PEPTIDE")
    assert report.height == len(pep.residues)
    assert set(report["chain.id"].unique()) == {pep.chain_id}
    assert (report["n_conformers"] >= 1).all()
    assert ((report["p_best"] > 0) & (report["p_best"] <= 1.0 + 1e-12)).all()
    # Gly and Pro have no rotatable chi, so exactly one conformer.
    for row in report.iter_rows(named=True):
        if row["residue.aa"] in ("G", "A", "P"):
            assert row["n_conformers"] == 1


@pytest.mark.slow
def test_repack_only_touches_the_selected_chains(annotated):
    """The default repacks the peptide; the TCR and MHC must come back untouched."""
    from tcren.rotamers import repack

    out, _ = repack(annotated)
    for before, after in zip(annotated.chains, out.chains):
        if before.chain_type == "PEPTIDE":
            continue
        moved = max(float(np.linalg.norm(a.coord - b.coord))
                    for r0, r1 in zip(before.residues, after.residues)
                    for a, b in zip(r0.atoms, r1.atoms))
        assert moved == 0.0, f"chain {before.chain_id} moved"


def test_repack_rejects_an_empty_selection(annotated):
    from tcren.rotamers import repack

    with pytest.raises(ValueError, match="no residues selected"):
        repack(annotated, chains=("NOT_A_CHAIN_TYPE",))
