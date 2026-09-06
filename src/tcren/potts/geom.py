r"""A length-generalising prior over the receptor:peptide Cα distance matrix.

:func:`tcren.potts.available_pairs` reads :math:`d_{\mathrm{C}\alpha}` off a structure. That makes
every downstream read-out conditional on having a template of the right shape — and lengths are not
evenly covered, so a CDR3 of length 17 against an 11-mer may have no template at all. This module
models the distance matrix instead:

.. math::  P\big(d_{\mathrm{C}\alpha} \mid \text{loop}, i, L,\ j, M,\ \text{MHC class}\big)

and emits a site table with the same schema :func:`available_pairs` emits, so every function in
:mod:`tcren.potts` accepts it unchanged.

**Length generalisation is the normalised coordinate, not the parameterisation.** A loop of length
:math:`L` against a peptide of length :math:`M` is resampled onto :math:`(u, v) = (i/(L+1),
j/(M+1))` — the graphon coordinates :func:`tcren.contactmap.registered_map` already uses — so
:math:`L = 10` and :math:`L = 17` land on one surface that can be read back at any length. What is
averaged is a *shape*; what is predicted is that shape stretched to the lengths asked for.

**Feasibility is checked against the closed-form loop physics, not against the fit.** A predicted
geometry has to be one an Ω-loop of that length can actually adopt, and
:func:`tcren.geometry.reach_max` answers that in closed form from the discrete-Frenet virtual-bond
chain — the maximum apex reach of a length-``L`` loop pinned at neck separation ``R``. Predictions
that exceed it are flagged rather than scored.

**Three loops, not one.** On the 17 bundled crystals CDR3 carries 69.3% of TCR:peptide contacts and
CDR1 a further 26.5%; CDR1+CDR2+CDR3 is 98.4%. A CDR3-only prior caps recall below 0.70 by
construction, so all three are modelled. CDR1 and CDR2 are germline and barely vary in length —
they need the surface, not the length generalisation.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from pathlib import Path

import numpy as np
import polars as pl

from ..geometry import BOND_CA, reach_max
from .sites import SITE_COLUMNS

__all__ = ["GeomPrior", "REGIONS_MODELLED", "fit_geom_prior", "loops_from_structure",
           "predicted_distances", "predicted_sites", "region_peptide_dmap"]

#: Receptor regions the prior carries a surface for. Together 98.4% of TCR:peptide contacts on the
#: bundled crystals; CDR3 alone is 69.3%.
REGIONS_MODELLED: tuple[str, ...] = ("CDR1", "CDR2", "CDR3")

_CHAINS: tuple[str, ...] = ("TRA", "TRB")


def region_peptide_dmap(structure, chain_type: str, region: str) -> np.ndarray | None:
    """``(L, M)`` Cα distance matrix between one receptor region and the peptide.

    Generalises :func:`tcren.contactmap._cdr3_peptide_dmaps` to any annotated region, which is what
    lets the prior reach CDR1 and CDR2 rather than capping its own recall at CDR3's 69.3%.
    """
    from ..contacts.geometry import _chain_ca
    from ..structure.model import PEPTIDE_TYPE

    pep = _chain_ca(structure, (PEPTIDE_TYPE,))
    if pep is None or len(pep) < 2:
        return None
    for c in structure.chains:
        if getattr(c, "chain_type", None) != chain_type:
            continue
        for reg in getattr(c, "regions", []) or []:
            if reg.region_type != region:
                continue
            pts = [r.ca for r in reg.residues if r.ca is not None]
            if len(pts) < 2:
                return None
            cas = np.asarray(pts, float)
            return np.linalg.norm(cas[:, None, :] - pep[None, :, :], axis=2)
    return None


def _resample(dmap: np.ndarray, grid: int) -> np.ndarray:
    """Bilinear resample of an ``(L, M)`` map onto ``grid x grid`` graphon coordinates."""
    from scipy.ndimage import map_coordinates
    u, v = np.meshgrid(np.linspace(0, dmap.shape[0] - 1, grid),
                       np.linspace(0, dmap.shape[1] - 1, grid), indexing="ij")
    return map_coordinates(dmap, [u.ravel(), v.ravel()], order=1,
                           mode="nearest").reshape(grid, grid)


def _key(chain: str, region: str, mhc_class: str) -> str:
    return f"{chain}:{region}|{mhc_class}"


@dataclass
class GeomPrior:
    """Mean and spread of the receptor:peptide Cα distance surface, in graphon coordinates.

    Attributes:
        grid: Surface resolution ``G`` per axis.
        mean / sd / n: One ``G x G`` surface per ``"{chain}:{region}|{class}"`` key, and how many
            structures went into it.
        neck_mean / neck_sd: Loop neck separation ``R`` (Å) per ``"{chain}:{region}"``, which
            :func:`tcren.geometry.reach_max` needs for the feasibility check.
        bond: Virtual Cα–Cα bond length the reachability bound uses.
        n_structures: Structures the fit saw.
    """

    grid: int
    mean: dict[str, list[list[float]]]
    sd: dict[str, list[list[float]]]
    n: dict[str, int]
    neck_mean: dict[str, float] = field(default_factory=dict)
    neck_sd: dict[str, float] = field(default_factory=dict)
    bond: float = BOND_CA
    n_structures: int = 0
    notes: str = ""

    def keys(self) -> list[str]:
        return sorted(self.mean)

    def surface(self, chain: str, region: str, mhc_class: str) -> tuple[np.ndarray, np.ndarray]:
        """``(mean, sd)`` surfaces for one key, falling back across MHC class when unseen."""
        k = _key(chain, region, mhc_class)
        if k not in self.mean:
            alt = [c for c in self.mean if c.startswith(f"{chain}:{region}|")]
            if not alt:
                raise KeyError(f"no surface for {chain}:{region} at any MHC class; "
                               f"the prior carries {self.keys()}")
            k = alt[0]
        return np.asarray(self.mean[k]), np.asarray(self.sd[k])

    def to_json(self, path: str | Path) -> None:
        Path(path).write_text(json.dumps(asdict(self)))

    @classmethod
    def from_json(cls, path: str | Path) -> "GeomPrior":
        return cls(**json.loads(Path(path).read_text()))


def fit_geom_prior(structures, *, grid: int = 12, regions=REGIONS_MODELLED,
                   notes: str = "") -> GeomPrior:
    """Average the receptor:peptide Cα surfaces of a structure set in graphon coordinates.

    Args:
        structures: An iterable of chain-typed structures — pass
            :func:`tcren.annotation.batch.iter_annotated_set` rather than annotating in a loop.
        grid: Surface resolution per axis.
        regions: Receptor regions to model.
        notes: Free text stored on the prior.

    Returns:
        A :class:`GeomPrior`.
    """
    acc: dict[str, list[np.ndarray]] = {}
    necks: dict[str, list[float]] = {}
    seen = 0
    for s in structures:
        cls = "MHCII" if any(getattr(c, "mhc_class", None) == "MHCII"
                             for c in getattr(s, "mhc_calls", []) or []) else "MHCI"
        used = False
        for chain in _CHAINS:
            for region in regions:
                d = region_peptide_dmap(s, chain, region)
                if d is None or min(d.shape) < 2:
                    continue
                acc.setdefault(_key(chain, region, cls), []).append(_resample(d, grid))
                used = True
        for chain in _CHAINS:
            for region in regions:
                pts = _region_ca(s, chain, region)
                if pts is not None and len(pts) >= 2:
                    necks.setdefault(f"{chain}:{region}", []).append(
                        float(np.linalg.norm(pts[-1] - pts[0])))
        seen += int(used)
    if not acc:
        raise ValueError("no structure yielded a receptor:peptide Ca surface — were the chains "
                         "typed (classify_chains) before this was called?")
    return GeomPrior(
        grid=grid,
        mean={k: np.mean(v, axis=0).tolist() for k, v in acc.items()},
        sd={k: np.std(v, axis=0).tolist() for k, v in acc.items()},
        n={k: len(v) for k, v in acc.items()},
        neck_mean={k: float(np.mean(v)) for k, v in necks.items()},
        neck_sd={k: float(np.std(v)) for k, v in necks.items()},
        n_structures=seen, notes=notes)


def _region_ca(structure, chain_type: str, region: str) -> np.ndarray | None:
    for c in structure.chains:
        if getattr(c, "chain_type", None) != chain_type:
            continue
        for reg in getattr(c, "regions", []) or []:
            if reg.region_type == region:
                pts = [r.ca for r in reg.residues if r.ca is not None]
                return np.asarray(pts, float) if len(pts) >= 2 else None
    return None


def predicted_distances(prior: GeomPrior, chain: str, region: str, *, loop_len: int,
                        pep_len: int, mhc_class: str = "MHCI") -> tuple[np.ndarray, np.ndarray]:
    r"""``(mean, sd)`` Cα distance matrices of shape ``(loop_len, pep_len)``.

    The surface is read back at :math:`(u, v) = (i/(L+1), j/(M+1))`, which is what carries a shape
    averaged over one set of lengths to a pair of lengths it never saw.
    """
    from scipy.ndimage import map_coordinates

    m, sd = prior.surface(chain, region, mhc_class)
    g = prior.grid
    u = (np.arange(1, loop_len + 1) / (loop_len + 1)) * (g - 1)
    v = (np.arange(1, pep_len + 1) / (pep_len + 1)) * (g - 1)
    uu, vv = np.meshgrid(u, v, indexing="ij")
    out = [map_coordinates(x, [uu.ravel(), vv.ravel()], order=1,
                           mode="nearest").reshape(loop_len, pep_len) for x in (m, sd)]
    return out[0], out[1]


def feasible(prior: GeomPrior, chain: str, region: str, loop_len: int,
             distance: float) -> bool:
    """Can a loop of ``loop_len`` residues reach ``distance`` from its own neck?

    Closed form, from the virtual-bond chain — :func:`tcren.geometry.reach_max`. This is the
    physics the fitted surface is checked against; it reads no fitted parameter but the neck.
    """
    neck = prior.neck_mean.get(f"{chain}:{region}")
    if neck is None:
        return True
    return distance <= reach_max(loop_len, neck, prior.bond) + neck / 2.0


def loops_from_structure(structure, regions=REGIONS_MODELLED) -> dict[str, str]:
    """``{"TRA:CDR3": "CAV...", ...}`` — the receptor loop sequences a prediction needs.

    Only the *resolved* residues are taken, so the sequence lines up with the Cα positions the
    surfaces were averaged over.
    """
    out: dict[str, str] = {}
    for c in structure.chains:
        ct = getattr(c, "chain_type", None)
        if ct not in _CHAINS:
            continue
        for reg in getattr(c, "regions", []) or []:
            if reg.region_type in regions:
                seq = "".join(r.aa for r in reg.residues if r.ca is not None and r.aa)
                if seq:
                    out[f"{ct}:{reg.region_type}"] = seq
    return out


def predicted_sites(prior: GeomPrior, loops: dict[str, str], peptide: str, *,
                    mhc_class: str = "MHCI", radius: float = 15.0,
                    pdb_id: str = "predicted") -> pl.DataFrame:
    r"""A site table built from sequence and lengths alone — no structure, no template.

    Emits the schema :func:`tcren.potts.available_pairs` emits, so
    :func:`tcren.potts.site_codes`, :func:`tcren.potts.eta` and every scoring function accept it
    unchanged. ``d_heavy`` is null and ``sigma`` is 0 throughout: nothing here observed a contact,
    and a caller that treats ``sigma`` as data on this frame is reading a placeholder as a fact.

    Args:
        prior: A fitted :class:`GeomPrior`.
        loops: ``{"TRA:CDR3": sequence, ...}`` — see :func:`loops_from_structure`.
        peptide: The peptide sequence.
        mhc_class: ``"MHCI"`` or ``"MHCII"``.
        radius: Availability radius, Å. A predicted pair beyond it is not a site, exactly as in
            :func:`available_pairs`.
        pdb_id: Identifier written into the frame.

    Returns:
        One row per predicted site, with an extra ``d_ca_sd`` column carrying the surface's own
        spread at that cell and ``feasible`` flagging cells beyond the loop's closed-form reach.
    """
    from ..refine.anchors import predict_anchors

    try:
        anchors = set(predict_anchors(peptide).anchors)
    except Exception:
        anchors = set()
    M = len(peptide)
    rows = []
    for key, seq in sorted(loops.items()):
        chain, region = key.split(":")
        L = len(seq)
        if L < 1 or M < 1:
            continue
        mean, sd = predicted_distances(prior, chain, region, loop_len=L, pep_len=M,
                                       mhc_class=mhc_class)
        for i in range(L):
            for j in range(M):
                d = float(mean[i, j])
                if d > radius:
                    continue
                rows.append({"pdb.id": pdb_id, "aa.rec": seq[i], "chain.rec": chain,
                             "region.rec": region, "pos.rec": i, "aa.par": peptide[j],
                             "pos.par": j, "par.len": M,
                             "role.par": "anchor" if j in anchors else "tcr_facing",
                             "cls": "peptide", "d_heavy": None, "d_ca": d, "sigma": 0.0,
                             "d_ca_sd": float(sd[i, j]),
                             "feasible": feasible(prior, chain, region, L, d)})
    if not rows:
        return pl.DataFrame(schema={**{c: pl.Utf8 for c in SITE_COLUMNS},
                                    "d_ca_sd": pl.Float64, "feasible": pl.Boolean})
    return pl.DataFrame(rows).select(list(SITE_COLUMNS) + ["d_ca_sd", "feasible"])
