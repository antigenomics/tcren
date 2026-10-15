r"""A geometry drawn from a fitted joint, in coordinates rather than in distances.

:mod:`tcren.potts.geom` models the receptor:peptide Cα *distance surface* and averages it over
crystals. Two things follow from that choice, and both are why its placement fails --- Jensen--
Shannon 0.4200 against the crystal's own 0.7052, contacts put 2.3 Å (CDR3) to 5.8 Å (CDR2) too far
out over 7,732 held-out contacts.

**A mean of distance matrices is not a distance matrix any chain can realise.** Consecutive Cα atoms
are a fixed virtual bond apart --- measured 3.7985 Å, s.d. 0.0668, over 642 CDR gaps --- and nothing
in a per-cell average enforces that. Recovering coordinates from a distance matrix is a Euclidean
distance-matrix completion and is not convex; the tractable relaxation
(:math:`|d(i{+}1,j) - d(i,j)| \le b`) is necessary and not sufficient.

**And the cells are not independent.** A loop moves as a body, so the distances within one map
covary; a per-cell mean and s.d. cannot express that.

This module works where both are free. ``Canonical2026`` puts every complex in one frame, so loop Cα
coordinates are directly comparable across crystals, and the model is a joint Gaussian over them:

.. math::  x \sim \mathcal{N}\big(\mu(\text{lengths}),\ W W^{\mathsf T} + \Psi\big)

over the concatenated graphon-knot coordinates of all six CDR loops and the peptide. A draw is a set
of 3D traces, :func:`tcren.geometry.enforce_bond` puts each one exactly on the virtual bond, and the
distance matrix is then *derived* rather than modelled --- so every geometry this emits is one some
chain could adopt, with no relaxation anywhere.

**Length stops being a joint problem.** The distance surface had to be conditioned on the *pair*
:math:`(L, M)`, and a CDR3 of 16 against an 11-mer may be a combination no crystal solved --- 359
crystals cover 62 of the 131 ``(L, M)`` combinations the receptors and peptides between them imply.
Here :math:`L` conditions the loop block and :math:`M` the peptide block, separately, and their
distances follow from the geometry. The pair never has to have been observed.

What is fitted from what
------------------------

The receptor block comes from the 359 TCRs and the peptide block from the 359 peptides.
:mod:`tcren.potts.cross` crossing every TCR with every pMHC *is* the product of those two marginals,
so it adds no coordinate information --- what it adds is the empirical
:math:`P(d_{C\alpha} \mid \text{region}, i, L, j, M)` at length pairs no crystal has, which is what
:func:`distance_agreement` checks a draw against.
"""
from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from pathlib import Path

import numpy as np
import polars as pl

from ..geometry import (BOND_CA, _read_weights, enforce_bond, graphon_weights,
                        orientation_cosine)
from .cross import receptor_residues, partner_residues
from .geom import REGIONS_MODELLED, _CHAINS
from .sites import ORIENT_COLUMN, SITE_COLUMNS

__all__ = ["BLOCKS", "BUNDLED_CHAIN", "ChainPrior", "OUT_OF_REACH_ETA", "PEPTIDE_BLOCK",
           "annealed_contact_map", "distance_agreement", "fit_chain_prior",
           "predicted_sites_chain", "sample_geometries"]

#: One-body log-odds given to a site a particular draw put beyond the radius. ``exp(-60)`` is
#: numerically zero, so the site is present in the vector and cannot be occupied --- which is what
#: makes every draw a distribution over the **same** configuration space.
OUT_OF_REACH_ETA: float = -60.0

#: The peptide's block in the joint vector.
PEPTIDE_BLOCK: str = "PEPTIDE"

#: Blocks of the joint, in order. Six receptor loops and the peptide, concatenated into one vector
#: so the covariance carries what moves them **together** --- above all the docking angle, a
#: whole-receptor rigid-body degree of freedom that no per-loop model can represent. That is the one
#: CDR2 needs: CDR2 is germline and barely varies in length, so length conditioning buys it nothing.
BLOCKS: tuple[str, ...] = tuple(f"{c}:{r}" for c in _CHAINS for r in REGIONS_MODELLED) \
    + (PEPTIDE_BLOCK,)

#: The chain prior shipped in the wheel.
BUNDLED_CHAIN: str = "geom_tcr_peptide_chain"

#: Registration convention for the knots. ``"endpoint"`` puts residue 0 on knot 0 and residue
#: ``n-1`` on the last knot, so loop **termini align across structures of different length** --- the
#: operation an average over crystals actually performs --- and its round trip is *exact*, which a
#: coordinate model needs and the distance surface never had. See :func:`tcren.geometry.graphon_weights`.
_CONVENTION = "endpoint"


def _trace(structure, block: str) -> np.ndarray | None:
    """The ``(n, 3)`` Cα trace of one block, or ``None`` if it is absent or too short."""
    if block == PEPTIDE_BLOCK:
        par = partner_residues(structure)
        if len(par) < 2:
            return None
        return par.select("ca_x", "ca_y", "ca_z").to_numpy().astype(float)
    chain, region = block.split(":")
    rec = receptor_residues(structure, regions=(region,), chains=(chain,))
    if len(rec) < 2:
        return None
    return rec.select("ca_x", "ca_y", "ca_z").to_numpy().astype(float)


def _to_knots(trace: np.ndarray, grid: int) -> np.ndarray:
    """``(grid * 3,)`` knot coordinates of an ``(n, 3)`` trace."""
    return (graphon_weights(len(trace), grid, convention=_CONVENTION) @ trace).ravel()


def _from_knots(knots: np.ndarray, length: int, grid: int) -> np.ndarray:
    """``(..., length, 3)`` trace from ``(..., grid * 3)`` knots."""
    shaped = knots.reshape(*knots.shape[:-1], grid, 3)
    return np.einsum("lg,...gk->...lk", _read_weights(length, grid, _CONVENTION), shaped)


def _kernel(lengths: np.ndarray, bins: np.ndarray, bandwidth: float) -> np.ndarray:
    """``(len(bins), len(lengths))`` Gaussian weights over integer lengths."""
    w = np.exp(-0.5 * ((bins[:, None] - lengths[None, :]) / bandwidth) ** 2)
    return w


@dataclass
class ChainPrior:
    """A joint Gaussian over graphon-knot Cα coordinates, per MHC class.

    Attributes:
        grid: Knots per block, ``G``. The joint vector is ``len(BLOCKS) * G * 3`` long.
        rank: Factors in the low-rank covariance.
        mu: ``key -> (n_bins, G*3)`` kernel-smoothed conditional mean, keyed
            ``"{block}|{class}"``; ``bins`` gives the length each row is the mean at.
        bins: ``key -> (n_bins,)`` the integer lengths ``mu`` is tabulated at.
        n_eff: ``key -> (n_bins,)`` kernel-weighted structures behind each row of ``mu``. A caller
            reading a length far from any crystal can see that it is extrapolating.
        w: ``class -> (D, rank)`` factor loadings over the **whole** joint vector, so a draw moves
            all six loops and the peptide together.
        psi: ``class -> (D,)`` diagonal residual variance.
        bond: Virtual Cα–Cα bond the draws are projected onto.
        bandwidth: Gaussian kernel width over lengths, in residues.
        n_structures: Structures the fit saw, per class.
        notes: Free text.
    """

    grid: int
    rank: int
    mu: dict[str, list[list[float]]]
    bins: dict[str, list[int]]
    n_eff: dict[str, list[float]]
    w: dict[str, list[list[float]]]
    psi: dict[str, list[float]]
    bond: float = BOND_CA
    bandwidth: float = 1.5
    n_structures: dict[str, int] = field(default_factory=dict)
    notes: str = ""

    @property
    def dim(self) -> int:
        """Length of the joint vector, ``len(BLOCKS) * grid * 3``."""
        return len(BLOCKS) * self.grid * 3

    def classes(self) -> list[str]:
        return sorted(self.w)

    def _resolve(self, mhc_class: str) -> str:
        if mhc_class in self.w:
            return mhc_class
        if not self.w:
            raise ValueError("this prior carries no class")
        return self.classes()[0]

    def block_mean(self, block: str, length: int, mhc_class: str) -> np.ndarray:
        """``(G*3,)`` conditional mean of one block at one length, interpolated between bins."""
        key = f"{block}|{self._resolve(mhc_class)}"
        bins = np.asarray(self.bins[key], float)
        table = np.asarray(self.mu[key], float)
        return np.stack([np.interp(length, bins, table[:, d]) for d in range(table.shape[1])])

    def support(self, block: str, length: int, mhc_class: str) -> float:
        """Kernel-weighted structures behind :meth:`block_mean` --- how far this is a guess."""
        key = f"{block}|{self._resolve(mhc_class)}"
        return float(np.interp(length, np.asarray(self.bins[key], float),
                               np.asarray(self.n_eff[key], float)))

    def to_npz(self, path: str | Path) -> None:
        """Everything, arrays included."""
        meta = {k: v for k, v in asdict(self).items() if k not in ("mu", "w", "psi")}
        arrays = {f"mu::{k}": np.asarray(v, np.float32) for k, v in self.mu.items()}
        arrays |= {f"w::{k}": np.asarray(v, np.float32) for k, v in self.w.items()}
        arrays |= {f"psi::{k}": np.asarray(v, np.float32) for k, v in self.psi.items()}
        np.savez_compressed(path, meta=np.array(json.dumps(meta)), **arrays)

    @classmethod
    def from_npz(cls, path: str | Path) -> "ChainPrior":
        with np.load(path, allow_pickle=False) as z:
            meta = json.loads(str(z["meta"]))
            parts: dict[str, dict] = {"mu": {}, "w": {}, "psi": {}}
            for name in z.files:
                if "::" in name:
                    kind, key = name.split("::", 1)
                    parts[kind][key] = z[name].astype(float).tolist()
        return cls(**meta, **parts)

    @classmethod
    def bundled(cls, name: str = BUNDLED_CHAIN) -> "ChainPrior":
        """A chain prior shipped in the wheel, by name."""
        import io
        from importlib.resources import files

        res = files("tcren.data").joinpath(f"{name}.npz")
        if not res.is_file():
            raise FileNotFoundError(f"no bundled chain prior named {name!r}")
        return cls.from_npz(io.BytesIO(res.read_bytes()))


def fit_chain_prior(structures, *, grid: int = 24, rank: int = 12, bandwidth: float = 1.5,
                    exclude: set[str] = frozenset(), notes: str = "") -> ChainPrior:
    """Fit the joint over knot coordinates from a co-framed structure set.

    Args:
        structures: Chain-typed, MHC-annotated structures **in one canonical frame** --- the
            coordinates are compared across structures directly, so an unoriented set fits noise.
            ``data/Canonical2026`` is the set this is for.
        grid: Knots per block. 24 keeps the coordinate round trip under 0.67 Å at every loop length
            a crystal carries (0.12 Å at length 4); 12 costs up to 1.39 Å.
        rank: Factors. The joint is ``len(BLOCKS) * grid * 3`` long against a few hundred
            structures, so the covariance is estimated low-rank plus diagonal rather than in full.
        bandwidth: Gaussian kernel over integer lengths, residues. A length with no crystal is
            filled from its neighbours, and :meth:`ChainPrior.support` says how well.
        exclude: ``pdb.id`` values to leave out --- how a held-out fit is built.
        notes: Free text stored on the prior.

    Returns:
        A :class:`ChainPrior`.

    Raises:
        ValueError: if no structure contributed every block.
    """
    drop = {x.lower() for x in exclude}
    per_class: dict[str, list[tuple[np.ndarray, dict[str, int]]]] = {}
    for structure in structures:
        if str(structure.pdb_id).lower() in drop:
            continue
        traces = {b: _trace(structure, b) for b in BLOCKS}
        if any(t is None for t in traces.values()):
            continue                     # a partial vector cannot enter a joint covariance
        cls = "MHCII" if any(getattr(c, "mhc_class", None) == "MHCII"
                             for c in getattr(structure, "mhc_calls", []) or []) else "MHCI"
        vector = np.concatenate([_to_knots(traces[b], grid) for b in BLOCKS])
        per_class.setdefault(cls, []).append((vector, {b: len(traces[b]) for b in BLOCKS}))
    if not per_class:
        raise ValueError("no structure carried all six CDR loops and a peptide — were the chains "
                         "typed and the MHC annotated before this was called?")

    width = grid * 3
    mu: dict[str, list] = {}
    bins: dict[str, list] = {}
    n_eff: dict[str, list] = {}
    w_out: dict[str, list] = {}
    psi_out: dict[str, list] = {}
    n_struct: dict[str, int] = {}
    for cls, items in per_class.items():
        x = np.stack([v for v, _ in items])                       # (S, D)
        n_struct[cls] = len(x)
        centre = np.empty_like(x)
        for b, block in enumerate(BLOCKS):
            lengths = np.array([lens[block] for _, lens in items], float)
            grid_bins = np.arange(int(lengths.min()), int(lengths.max()) + 1)
            weights = _kernel(lengths, grid_bins.astype(float), bandwidth)   # (n_bins, S)
            block_x = x[:, b * width:(b + 1) * width]
            table = (weights @ block_x) / weights.sum(axis=1, keepdims=True)
            key = f"{block}|{cls}"
            mu[key] = table.tolist()
            bins[key] = grid_bins.tolist()
            n_eff[key] = weights.sum(axis=1).tolist()
            # Centre each structure on the mean at ITS OWN length, so the residual carries only the
            # variation the length does not explain -- which is what the covariance should describe.
            at = np.searchsorted(grid_bins, lengths.astype(int))
            centre[:, b * width:(b + 1) * width] = block_x - table[at]
        # Probabilistic PCA in closed form: the top `rank` eigenvectors of the residual covariance,
        # with the discarded tail becoming the isotropic noise floor.
        cov = centre.T @ centre / max(len(x) - 1, 1)
        values, vectors = np.linalg.eigh(cov)
        values, vectors = values[::-1], vectors[:, ::-1]
        k = min(rank, len(values) - 1, len(x) - 1)
        sigma2 = max(float(values[k:].mean()), 1e-9)
        loadings = vectors[:, :k] * np.sqrt(np.maximum(values[:k] - sigma2, 0.0))
        w_out[cls] = loadings.tolist()
        psi_out[cls] = np.maximum(np.diag(cov) - (loadings ** 2).sum(axis=1), 1e-6).tolist()
    return ChainPrior(grid=grid, rank=rank, mu=mu, bins=bins, n_eff=n_eff, w=w_out, psi=psi_out,
                      bandwidth=bandwidth, n_structures=n_struct, notes=notes)


def sample_geometries(prior: ChainPrior, lengths: dict[str, int], *, mhc_class: str = "MHCI",
                      draws: int = 256, seed: int = 0,
                      project: bool = True) -> dict[str, np.ndarray]:
    """Draw whole geometries from the joint and put each one on the virtual bond.

    Args:
        prior: A fitted :class:`ChainPrior`.
        lengths: ``block -> residues`` for every entry of :data:`BLOCKS`.
        mhc_class: ``"MHCI"`` or ``"MHCII"``; falls back across classes when one is unseen.
        draws: Geometries to draw.
        seed: RNG seed.
        project: Put every drawn trace exactly on ``prior.bond``. Leave it on --- an unprojected
            draw is a Gaussian sample, not a chain.

    Returns:
        ``block -> (draws, length, 3)`` Cα traces, all in the canonical frame and mutually
        consistent, because one draw of the joint produced all of them.

    Raises:
        KeyError: if ``lengths`` is missing a block.
    """
    missing = [b for b in BLOCKS if b not in lengths]
    if missing:
        raise KeyError(f"lengths is missing {missing}")
    cls = prior._resolve(mhc_class)
    width = prior.grid * 3
    mean = np.concatenate([prior.block_mean(b, lengths[b], cls) for b in BLOCKS])
    loadings = np.asarray(prior.w[cls], float)
    psi = np.asarray(prior.psi[cls], float)
    rng = np.random.default_rng(seed)
    factors = rng.standard_normal((draws, loadings.shape[1]))
    noise = rng.standard_normal((draws, len(psi))) * np.sqrt(psi)
    x = mean[None, :] + factors @ loadings.T + noise                 # (draws, D)
    out = {}
    for b, block in enumerate(BLOCKS):
        trace = _from_knots(x[:, b * width:(b + 1) * width], lengths[block], prior.grid)
        out[block] = enforce_bond(trace, prior.bond) if project else trace
    return out


def _role(peptide: str, mhc_class: str) -> list[str]:
    from ..refine.anchors import predict_anchors

    anchors = set(predict_anchors(peptide).anchors)
    return ["anchor" if i in anchors else "tcr_facing" for i in range(len(peptide))]


def predicted_sites_chain(prior: ChainPrior, loops: dict[str, str], peptide: str, *,
                          mhc_class: str = "MHCI", radius: float = 15.0, draws: int = 256,
                          seed: int = 0, project: bool = True,
                          pdb_id: str = "predicted") -> pl.DataFrame:
    """A site table from sequence alone, one row per ``(draw, receptor residue, peptide residue)``.

    Each draw is tagged as its own ``pdb.id`` --- the shape :func:`tcren.potts.score.contact_probabilities`
    already threads over, so a whole draw set is scored in one call rather than one call per draw.

    Args:
        prior: A fitted :class:`ChainPrior`.
        loops: ``"{chain}:{region}" -> sequence`` for the six CDR loops.
        peptide: Peptide sequence.
        mhc_class: ``"MHCI"`` or ``"MHCII"``.
        radius: Keep pairs at or under this Cα distance, Å.
        draws: Geometries to draw.
        seed: RNG seed.
        project: See :func:`sample_geometries`.
        pdb_id: Prefix for the per-draw ids, which are ``f"{pdb_id}#{k}"``.

    Returns:
        A frame carrying :data:`tcren.potts.sites.SITE_COLUMNS` plus ``cos_cb`` and ``draw``. ``d_heavy`` is null
        and ``sigma`` is 0 throughout: nothing here observed a contact, so this table places a pair
        and must never be fitted from.

    Raises:
        KeyError: if ``loops`` is missing one of the six.
    """
    missing = [b for b in BLOCKS if b != PEPTIDE_BLOCK and b not in loops]
    if missing:
        raise KeyError(f"loops is missing {missing}")
    lengths = {b: len(loops[b]) for b in BLOCKS if b != PEPTIDE_BLOCK} | {
        PEPTIDE_BLOCK: len(peptide)}
    traces = sample_geometries(prior, lengths, mhc_class=mhc_class, draws=draws, seed=seed,
                              project=project)
    pep = traces[PEPTIDE_BLOCK]                                   # (draws, M, 3)
    roles = _role(peptide, mhc_class)
    frames = []
    for block in BLOCKS:
        if block == PEPTIDE_BLOCK:
            continue
        chain, region = block.split(":")
        loop = traces[block]                                      # (draws, L, 3)
        d = np.linalg.norm(loop[:, :, None, :] - pep[:, None, :, :], axis=-1)   # (draws, L, M)
        cos = orientation_cosine(loop, pep)
        keep = d <= radius
        if not keep.any():
            continue
        k, i, j = np.nonzero(keep)
        frames.append(pl.DataFrame({
            "pdb.id": np.char.add(f"{pdb_id}#", k.astype(str)),
            "aa.rec": np.asarray(list(loops[block]))[i],
            "chain.rec": np.repeat(chain, len(k)),
            "region.rec": np.repeat(region, len(k)),
            "pos.rec": i.astype(np.int64),
            "aa.par": np.asarray(list(peptide))[j],
            "pos.par": j.astype(np.int64),
            "par.len": np.full(len(k), len(peptide), dtype=np.int64),
            "role.par": np.asarray(roles)[j],
            "cls": np.repeat("peptide", len(k)),
            "d_heavy": np.full(len(k), None, dtype=object),
            "d_ca": d[keep].astype(float),
            "sigma": np.zeros(len(k)),
            ORIENT_COLUMN: cos[keep].astype(float),
            "draw": k.astype(np.int64),
        }))
    if not frames:
        raise ValueError(f"no receptor:peptide pair within {radius} Å in any of {draws} draws")
    return (pl.concat(frames).select([*SITE_COLUMNS, ORIENT_COLUMN, "draw"])
            .sort("pdb.id", "chain.rec", "region.rec", "pos.rec", "pos.par"))


def distance_agreement(sites: pl.DataFrame, cross: pl.DataFrame) -> pl.DataFrame:
    """Compare a draw set's Cα distances against the empirical distribution the crystals imply.

    The draws come from a joint over *coordinates*; :mod:`tcren.potts.cross` holds the distances
    every real TCR makes against every real pMHC. If the two agree cell by cell, the coordinate
    model has reproduced the distance distribution without ever being shown a length pair.

    Args:
        sites: Output of :func:`predicted_sites_chain`.
        cross: A cross table, filtered to the ``(chain, region, len.rec, par.len)`` being scored.

    Returns:
        One row per ``(chain.rec, region.rec, pos.rec, pos.par)`` with ``drawn_mean``, ``drawn_sd``,
        ``crystal_mean``, ``crystal_sd``, ``delta`` and the two supporting counts.
    """
    keys = ["chain.rec", "region.rec", "pos.rec", "pos.par"]
    drawn = sites.group_by(keys).agg(pl.col("d_ca").mean().alias("drawn_mean"),
                                     pl.col("d_ca").std().alias("drawn_sd"),
                                     pl.len().alias("draws"))
    crystal = cross.group_by(keys).agg(pl.col("d_ca").mean().alias("crystal_mean"),
                                       pl.col("d_ca").std().alias("crystal_sd"),
                                       pl.len().alias("pairings"))
    return (drawn.join(crystal, on=keys, how="inner")
            .with_columns((pl.col("drawn_mean") - pl.col("crystal_mean")).alias("delta"))
            .sort(keys))


def _union_sites(prior: ChainPrior, loops: dict[str, str], peptide: str, mhc_class: str,
                 radius: float, draws: int, seed: int, project: bool):
    """``(frame, d, cos, reach)`` --- the union cell set and, per draw, its geometry.

    ``frame`` has one row per union cell; ``d``, ``cos`` and ``reach`` are ``(draws, n_cells)``.
    ``cos`` is :func:`tcren.geometry.orientation_cosine` off the drawn backbone, so a drawn geometry
    can be asked the same question a crystal is: does this residue turn toward the peptide?
    """
    lengths = {b: len(loops[b]) for b in BLOCKS if b != PEPTIDE_BLOCK} | {
        PEPTIDE_BLOCK: len(peptide)}
    traces = sample_geometries(prior, lengths, mhc_class=mhc_class, draws=draws, seed=seed,
                              project=project)
    pep = traces[PEPTIDE_BLOCK]
    roles = _role(peptide, mhc_class)
    rows, dists, coses = [], [], []
    for block in BLOCKS:
        if block == PEPTIDE_BLOCK:
            continue
        chain, region = block.split(":")
        loop = traces[block]
        d = np.linalg.norm(loop[:, :, None, :] - pep[:, None, :, :], axis=-1)   # (draws, L, M)
        cos = orientation_cosine(loop, pep)                                     # (draws, L, M)
        union = (d <= radius).any(axis=0)                                       # (L, M)
        if not union.any():
            continue
        i, j = np.nonzero(union)
        rows.append(pl.DataFrame({
            "aa.rec": np.asarray(list(loops[block]))[i],
            "chain.rec": np.repeat(chain, len(i)),
            "region.rec": np.repeat(region, len(i)),
            "pos.rec": i.astype(np.int64),
            "aa.par": np.asarray(list(peptide))[j],
            "pos.par": j.astype(np.int64),
            "par.len": np.full(len(i), len(peptide), dtype=np.int64),
            "role.par": np.asarray(roles)[j],
            "cls": np.repeat("peptide", len(i)),
        }))
        dists.append(d[:, i, j])                                                # (draws, n_ij)
        coses.append(cos[:, i, j])
    if not rows:
        raise ValueError(f"no receptor:peptide pair within {radius} Å in any of {draws} draws")
    frame = pl.concat(rows)
    d = np.concatenate(dists, axis=1)
    return frame, d, np.concatenate(coses, axis=1), d <= radius


def annealed_contact_map(prior: ChainPrior, loops: dict[str, str], peptide: str, model, *,
                         mhc_class: str = "MHCI", radius: float = 15.0, draws: int = 256,
                         beta: float = 0.0, seed: int = 0, project: bool = True,
                         particles: int = 32, steps: int = 128,
                         chains: int = 64, burn: int = 100, gibbs_draws: int = 100,
                         thin: int = 3, pdb_id: str = "predicted") -> pl.DataFrame:
    r"""Draw geometries, score each under the Hamiltonian, and weight them by their own ``log Z``.

    .. math::  w_g = \operatorname{softmax}_g(\beta \log Z_g), \qquad
               \bar p_a = \sum_g w_g\, p_a^{(g)}

    ``beta = 0`` is the flat average over draws and **is the default, because it measured best**.
    ``beta = 1`` is the marginal of the joint over (contact map, geometry) under a flat prior on
    geometry --- the draws already carry :math:`P(g)`, so reweighting by :math:`Z_g` is the
    posterior.

    **Every draw is scored on the same site set, and that is not a detail.** ``log Z`` is extensive
    in the number of sites in reach, so weighting draws that each carry their own site count would
    select for the largest envelope. The union of every draw's reachable cells is built once and a
    draw that put a cell out of reach gets :data:`OUT_OF_REACH_ETA` there, so every draw is a
    distribution over one :math:`\{0,1\}^n` and the weights compare like with like.

    **Even so, annealing costs rather than pays** (2026-09-07, 10 crystals held out by structure,
    96 draws, 22.6 observed contacts each):

    ======  ========  ========  ========  ========
    beta    E[n]      ratio     js        w_eff
    ======  ========  ========  ========  ========
    0       22.0      0.90      **0.4342**  96.0
    0.02    26.8      1.10      0.4278    69.6
    0.05    36.4      1.46      0.4042    18.8
    0.2     46.5      1.82      0.3601    1.9
    1       47.1      1.84      0.3538    1.1
    ======  ========  ========  ========  ========

    The common site set removes the site-count confound and the direction survives it, so the
    mechanism is the **contact** count: ``log Z`` is a free energy and rises with occupancy, so
    selecting on it inflates the map. This reproduces the recorded failure of Boltzmann-reweighting
    the fixed crystal members (MAE 6.01 → 8.75 at λ = 0.5) with the alternative explanation ruled
    out. ``beta`` is kept because it is the parameter the reading is defined in terms of, and
    ``w_eff`` is reported so a degenerate anneal is visible rather than silent.

    Args:
        prior: A fitted :class:`ChainPrior`.
        loops: ``"{chain}:{region}" -> sequence`` for the six CDR loops.
        peptide: Peptide sequence.
        model: A fitted :class:`tcren.potts.PottsModel`.
        mhc_class: ``"MHCI"`` or ``"MHCII"``.
        radius: Reach bound, Å.
        draws: Geometries drawn.
        beta: Annealing strength. 0 (the default, and the measured best) averages over draws; 1 is
            the joint marginal; large selects one draw.
        seed: RNG seed.
        project: See :func:`sample_geometries`.
        particles, steps: Annealed-importance-sampling settings for ``log Z``.
        chains, burn, gibbs_draws, thin: Block-Gibbs settings for the per-draw marginals.
        pdb_id: ``pdb.id`` written on every row.

    Returns:
        One row per union cell, carrying :data:`tcren.potts.sites.SITE_COLUMNS` plus ``d_ca_sd``
        (spread of the drawn distance), ``p_avail`` (the weighted fraction of draws that put the
        cell in reach), ``p_model`` (the weighted marginal), ``log_z_mean`` and ``w_eff`` --- the
        effective number of draws carrying the weight, so a degenerate anneal is visible rather
        than silent.
    """
    from .kernel import bucket_edges, colour, coupling_matrix, edges
    from .model import N_ORIENT, OBIN
    from .sample import ais_log_z, gibbs
    from .sites import eta as eta_of, site_codes

    frame, d, cos, reach = _union_sites(prior, loops, peptide, mhc_class, radius, draws, seed,
                                        project)
    n = len(frame)
    # One ordering for every draw: the codes and the coupling graph depend only on the cell set and
    # the sequences, so they are built once and only `eta`'s distance block moves per draw.
    ordered = frame.with_columns(pl.lit(pdb_id).alias("pdb.id"),
                                 pl.lit(0.0).alias("d_ca"), pl.lit(0.0).alias("sigma"),
                                 pl.lit(0.0).alias(ORIENT_COLUMN),
                                 pl.lit(None, dtype=pl.Float64).alias("d_heavy"))
    codes, _, q = site_codes(ordered.select([*SITE_COLUMNS, ORIENT_COLUMN]), model)
    edge_lists = edges(q, joint=model.joint, far=model.far)
    kv = np.asarray(model.kernel, float)
    if edge_lists and kv.size:
        E, _ = bucket_edges(edge_lists, np.zeros(1, np.int64), kv)
        a_matrix = coupling_matrix(n, E[:, 0], E[:, 1], E[:, 2], kv)
        colours = colour(n, E[:, 0], E[:, 1])
    else:
        a_matrix, colours = np.zeros((n, n)), [np.arange(n)]

    nb = int(np.floor(model.radius / model.dbin)) + 1
    g_dist = np.asarray(model.g_dist, float)
    g_orient = np.asarray(model.g_orient, float) if model.g_orient else None
    # Everything that does not move with the geometry, evaluated once. The distance block and --
    # when the model carries one -- the orientation block are the only per-draw terms.
    base = eta_of(codes, model) - g_dist[codes[2]]
    if g_orient is not None:
        base = base - g_orient[codes[6]]
    rng = np.random.default_rng(seed)
    p_draw = np.empty((draws, n))
    log_z = np.empty(draws)
    for g in range(draws):
        bins = np.clip(np.floor(d[g] / model.dbin).astype(np.int64), 0, nb - 1)
        field = base + g_dist[bins]
        if g_orient is not None:
            field = field + g_orient[np.clip(np.floor((cos[g] + 1.0) / OBIN).astype(np.int64),
                                             0, N_ORIENT - 1)]
        eta_g = np.where(reach[g], field, OUT_OF_REACH_ETA)
        p_draw[g], _ = gibbs(eta_g, a_matrix, colours, rng, chains=chains, burn=burn,
                             draws=gibbs_draws, thin=thin)
        log_z[g], _ = ais_log_z(eta_g, a_matrix, colours, rng, particles=particles, steps=steps)
    shifted = beta * (log_z - log_z.max())
    weights = np.exp(shifted) / np.exp(shifted).sum()
    w_eff = float(1.0 / (weights ** 2).sum())
    return q.drop("sid", "loop", "pchain").with_columns(
        pl.Series("d_ca", (weights[:, None] * d).sum(0)),
        pl.Series(ORIENT_COLUMN, (weights[:, None] * cos).sum(0)),
        pl.Series("d_ca_sd", d.std(axis=0)),
        pl.Series("p_avail", (weights[:, None] * reach).sum(0)),
        pl.Series("p_model", (weights[:, None] * p_draw).sum(0)),
        pl.lit(float(log_z.mean())).alias("log_z_mean"),
        pl.lit(w_eff).alias("w_eff"),
        pl.lit(float(beta)).alias("beta"))
