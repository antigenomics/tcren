"""The coordinate-space geometry model: a draw is a chain, not a mean of distance matrices.

The algebra is pinned on synthetic data with a known answer, independently of any crystal, and the
contracts that keep the rest of `tcren.potts` working are pinned on the fixtures.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from tcren.geometry import BOND_CA, orientation_cosine, virtual_cb
from tcren.potts.chain import (BLOCKS, PEPTIDE_BLOCK, ChainPrior, fit_chain_prior,
                               predicted_sites_chain, sample_geometries)
from tcren.potts.sites import ORIENT_COLUMN, SITE_COLUMNS

PDB_DIR = Path(__file__).resolve().parents[1] / "assets" / "pdb"
pytestmark = pytest.mark.skipif(not PDB_DIR.exists(), reason="pdb assets missing")

LOOPS = {"TRA:CDR1": "DRGSQS", "TRA:CDR2": "IYSNGD", "TRA:CDR3": "AVTTDSWGKLQ",
         "TRB:CDR1": "MNHEY", "TRB:CDR2": "SVGAGI", "TRB:CDR3": "ASRPGLAGGRPEQY"}
PEPTIDE = "LLFGYPVYV"


@pytest.fixture(scope="module")
def structures():
    pytest.importorskip("arda")
    from tcren.annotation.batch import iter_annotated_set
    from tcren.mhc.regions import annotate_mhc_batch

    loaded = list(iter_annotated_set(PDB_DIR))
    annotate_mhc_batch(loaded)
    return loaded


@pytest.fixture(scope="module")
def prior(structures):
    return fit_chain_prior(structures, grid=12, rank=4)


def test_the_factor_model_recovers_a_known_subspace():
    """No crystals: the PPCA algebra is pinned against data generated from a known W."""
    rng = np.random.default_rng(0)
    dim, rank, n = 40, 3, 4000
    w_true = rng.normal(0.0, 1.0, (dim, rank))
    noise = 0.1
    x = rng.normal(size=(n, rank)) @ w_true.T + rng.normal(0.0, noise, (n, dim))
    cov = x.T @ x / (n - 1)
    values, vectors = np.linalg.eigh(cov)
    values, vectors = values[::-1], vectors[:, ::-1]
    sigma2 = values[rank:].mean()
    w_hat = vectors[:, :rank] * np.sqrt(np.maximum(values[:rank] - sigma2, 0.0))
    # Compare the subspaces, not the loadings: W is identified only up to rotation.
    true_basis = np.linalg.qr(w_true)[0]
    hat_basis = np.linalg.qr(w_hat)[0]
    principal = np.linalg.svd(true_basis.T @ hat_basis, compute_uv=False)
    assert np.allclose(principal, 1.0, atol=1e-2), principal
    assert sigma2 == pytest.approx(noise ** 2, rel=0.2)


def test_every_block_is_fitted(prior, structures):
    assert set(prior.classes()) <= {"MHCI", "MHCII"}
    for cls in prior.classes():
        assert np.asarray(prior.w[cls]).shape[0] == prior.dim
        assert len(prior.psi[cls]) == prior.dim
        for block in BLOCKS:
            key = f"{block}|{cls}"
            assert np.asarray(prior.mu[key]).shape == (len(prior.bins[key]), prior.grid * 3)
            assert (np.asarray(prior.n_eff[key]) > 0).all()


def test_a_draw_is_a_chain(prior):
    """The whole point: every drawn trace sits exactly on the virtual bond."""
    lengths = {b: len(LOOPS[b]) for b in BLOCKS if b != PEPTIDE_BLOCK} | {
        PEPTIDE_BLOCK: len(PEPTIDE)}
    traces = sample_geometries(prior, lengths, draws=16, seed=0)
    for block, trace in traces.items():
        assert trace.shape == (16, lengths[block], 3)
        gaps = np.linalg.norm(np.diff(trace, axis=-2), axis=-1)
        assert np.allclose(gaps, BOND_CA, atol=2e-7), f"{block} is not a chain"


def test_an_unprojected_draw_is_not_a_chain(prior):
    """If it were, the projection would be doing nothing and could be dropped."""
    lengths = {b: len(LOOPS[b]) for b in BLOCKS if b != PEPTIDE_BLOCK} | {
        PEPTIDE_BLOCK: len(PEPTIDE)}
    raw = sample_geometries(prior, lengths, draws=16, seed=0, project=False)
    gaps = np.linalg.norm(np.diff(raw["TRB:CDR3"], axis=-2), axis=-1)
    assert np.abs(gaps - BOND_CA).max() > 0.05


def test_a_length_never_seen_still_reads_back(prior):
    """Length conditions the loop block alone, so an unseen (L, M) pair needs no crystal."""
    lengths = {b: 6 for b in BLOCKS if b != PEPTIDE_BLOCK} | {PEPTIDE_BLOCK: 20}
    lengths["TRB:CDR3"] = 22
    traces = sample_geometries(prior, lengths, draws=4, seed=0)
    assert traces["TRB:CDR3"].shape == (4, 22, 3)
    assert np.isfinite(np.concatenate([t.ravel() for t in traces.values()])).all()


def test_support_falls_away_from_the_crystals(prior):
    """`support` is what tells a caller it is extrapolating."""
    near = prior.support("TRB:CDR3", 12, "MHCI")
    far = prior.support("TRB:CDR3", 40, "MHCI")
    assert near > far


def test_a_missing_length_is_refused(prior):
    with pytest.raises(KeyError, match="missing"):
        sample_geometries(prior, {"TRA:CDR1": 6}, draws=2)


def test_predicted_sites_satisfy_the_potts_schema(prior):
    """The contract step 2 already holds: `site_codes` and `eta` must accept this unchanged."""
    from tcren.potts import PottsModel
    from tcren.potts.sites import eta, site_codes

    sites = predicted_sites_chain(prior, LOOPS, PEPTIDE, draws=4, seed=0)
    assert set(SITE_COLUMNS).issubset(sites.columns)
    assert ORIENT_COLUMN in sites.columns
    model = PottsModel.bundled()
    codes, _, _ = site_codes(sites.select(SITE_COLUMNS), model)
    assert np.isfinite(eta(codes, model)).all()


def test_predicted_sites_never_claim_an_observation(prior):
    """`sigma` is a placeholder here — a fit pointed at this table would be fitting zeros."""
    sites = predicted_sites_chain(prior, LOOPS, PEPTIDE, draws=4, seed=0)
    assert (sites["sigma"] == 0.0).all()
    assert sites["d_heavy"].is_null().all()
    assert sites["d_ca"].max() <= 15.0
    assert sites["pdb.id"].n_unique() == 4          # one id per draw, the shape scoring threads on


def test_a_seed_reproduces_a_draw(prior):
    lengths = {b: len(LOOPS[b]) for b in BLOCKS if b != PEPTIDE_BLOCK} | {
        PEPTIDE_BLOCK: len(PEPTIDE)}
    a = sample_geometries(prior, lengths, draws=8, seed=7)
    b = sample_geometries(prior, lengths, draws=8, seed=7)
    assert all(np.array_equal(a[k], b[k]) for k in a)


def test_excluding_a_structure_changes_the_fit(structures):
    held = fit_chain_prior(structures, grid=12, rank=4, exclude={structures[0].pdb_id})
    full = fit_chain_prior(structures, grid=12, rank=4)
    assert held.n_structures != full.n_structures or held.mu != full.mu


def test_the_npz_round_trips(prior, tmp_path):
    prior.to_npz(tmp_path / "chain.npz")
    back = ChainPrior.from_npz(tmp_path / "chain.npz")
    assert back.grid == prior.grid and back.bins == prior.bins
    for cls in prior.classes():
        assert np.allclose(np.asarray(back.w[cls]), np.asarray(prior.w[cls]), atol=1e-5)
        assert np.allclose(np.asarray(back.psi[cls]), np.asarray(prior.psi[cls]), atol=1e-5)


def test_a_set_with_no_complete_complex_is_refused(structures):
    from tcren.structure.model import PEPTIDE_TYPE

    stripped = [structures[0].__class__(
        s.pdb_id, [c for c in s.chains if c.chain_type != PEPTIDE_TYPE],
        s.complex_species, s.cell_type) for s in structures]
    with pytest.raises(ValueError, match="all six CDR loops and a peptide"):
        fit_chain_prior(stripped, grid=12, rank=4)


def test_orientation_cosine_signs_are_the_way_round():
    """+1 is a side chain aimed at the partner. If this inverts, `g_orient` reads backwards."""
    ca = np.array([[0.0, 0.0, 0.0], [3.8, 0.0, 0.0], [6.0, 3.1, 0.0]])
    cb = virtual_cb(ca)
    direction = cb[1] - ca[1]
    toward = ca[1] + 10.0 * direction / np.linalg.norm(direction)
    away = ca[1] - 10.0 * direction / np.linalg.norm(direction)
    cos = orientation_cosine(ca, np.stack([toward, away]))
    assert cos[1, 0] == pytest.approx(1.0, abs=1e-6)
    assert cos[1, 1] == pytest.approx(-1.0, abs=1e-6)


def test_a_collinear_trace_does_not_produce_nan():
    straight = np.array([[0.0, 0, 0], [3.8, 0, 0], [7.6, 0, 0], [11.4, 0, 0]])
    cb = virtual_cb(straight)
    assert np.isfinite(cb).all()
    assert np.allclose(np.linalg.norm(cb - straight, axis=1), np.linalg.norm(cb[0] - straight[0]))


def test_the_orientation_block_is_opt_in_on_both_sides(prior):
    """A model without `g_orient`, or a frame without `cos_cb`, must score exactly as before."""
    from tcren.potts import PottsModel
    from tcren.potts.sites import eta, site_codes

    model = PottsModel.bundled()
    assert model.g_orient == []
    sites = predicted_sites_chain(prior, LOOPS, PEPTIDE, draws=2, seed=0)
    plain, _, _ = site_codes(sites.select(SITE_COLUMNS), model)
    withcos, _, _ = site_codes(sites, model)
    assert len(plain) == len(withcos) == 6            # the block is absent from both
    assert np.allclose(eta(plain, model), eta(withcos, model))


def test_a_model_with_g_orient_needs_the_column(prior):
    from tcren.potts import PottsModel
    from tcren.potts.model import N_ORIENT
    from tcren.potts.sites import eta, site_codes

    model = PottsModel.bundled()
    model.g_orient = [0.1] * N_ORIENT
    sites = predicted_sites_chain(prior, LOOPS, PEPTIDE, draws=2, seed=0)
    codes, sizes, _ = site_codes(sites, model)
    assert len(codes) == 7 and sizes[-1] == N_ORIENT
    assert np.isfinite(eta(codes, model)).all()
    # Without the column the frame yields six blocks and the model seven -- refused, not silent.
    plain, _, _ = site_codes(sites.select(SITE_COLUMNS), model)
    with pytest.raises(ValueError, match="blocks"):
        eta(plain, model)
