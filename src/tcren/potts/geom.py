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

from ..geometry import BOND_CA, read_back, reach_max, resample_graphon
from .sites import SITE_COLUMNS

__all__ = ["BUNDLED_ENSEMBLE", "BUNDLED_PRIOR", "GeomPrior", "REGIONS_MODELLED", "feasible",
           "fit_geom_prior", "loops_from_structure", "predict_contact_map",
           "predicted_distances", "predicted_ensemble", "predicted_sites",
           "region_peptide_dmap"]

#: The Ca-geometry prior shipped in the wheel. Averaged over the same 359 alpha-beta Native2026
#: crystals `tcren2` is derived from, so a prediction and a potential see one structure set::
#:
#:     tcren potts geom -s "$TCREN_DATA_DIR/Native2026" -o geom_tcr_peptide.json --grid 12
BUNDLED_PRIOR: str = "geom_tcr_peptide"

#: The ENSEMBLE prior shipped in the wheel: the same crystals, kept one surface each instead of
#: collapsed to their mean. Averaging the surfaces and averaging over the surfaces are different
#: estimators, and only the second is a geometry -- see :func:`predict_contact_map`::
#:
#:     tcren potts geom -s "$TCREN_DATA_DIR/Native2026" -o geom_tcr_peptide_ens.npz --grid 32
BUNDLED_ENSEMBLE: str = "geom_tcr_peptide_ens"

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


def _resample(dmap: np.ndarray, grid: int, convention: str = "endpoint") -> np.ndarray:
    """Bilinear resample of an ``(L, M)`` map onto ``grid x grid`` graphon coordinates.

    Delegates to :func:`tcren.geometry.resample_graphon`. The default is ``"endpoint"``, which is
    what this function has always done and what the shipped priors were built under --- see
    :func:`tcren.geometry.graphon_weights` for why that is *not* adjoint to the read-back, and
    :func:`fit_geom_prior`'s ``convention`` argument for the one that is.
    """
    return resample_graphon(dmap, grid, convention=convention)


def _key(chain: str, region: str, mhc_class: str) -> str:
    return f"{chain}:{region}|{mhc_class}"


@dataclass
class GeomPrior:
    """Mean and spread of the receptor:peptide Cα distance surface, in graphon coordinates.

    Attributes:
        grid: Surface resolution ``G`` per axis.
        mean / sd / n: One ``G x G`` surface per ``"{chain}:{region}|{class}"`` key, and how many
            structures went into it.
        surfaces: One ``(S, G, G)`` stack per key when the fit kept them --- ONE registered
            surface per contributing crystal, not their mean. This is what
            :func:`predict_contact_map` averages over, and it is the whole difference between a
            geometry and an average of geometries: a contacting pair sits at the near end of the
            distribution being averaged, so the mean places it several angstroms too far out, and
            the distances within one map are not independent, since consecutive Cα atoms are a
            fixed virtual bond apart. Every member here is a chain that was actually solved, so
            both are respected by construction. ``None`` on a prior fitted without them.
        neck_mean / neck_sd: Loop neck separation ``R`` (Å) per ``"{chain}:{region}"``, which
            :func:`tcren.geometry.reach_max` needs for the feasibility check.
        bond: Virtual Cα–Cα bond length the reachability bound uses.
        n_structures: Structures the fit saw.
        convention: Graphon convention the surfaces were *registered* under. Read-back is always
            :math:`(i+1)/(L+1)`, so only ``"interior"`` makes the two adjoint --- see
            :func:`tcren.geometry.graphon_weights`. Defaults to ``"endpoint"`` so a prior written
            before the argument existed loads with the behaviour it was built under.
    """

    grid: int
    mean: dict[str, list[list[float]]]
    sd: dict[str, list[list[float]]]
    n: dict[str, int]
    surfaces: dict[str, np.ndarray] | None = None
    neck_mean: dict[str, float] = field(default_factory=dict)
    neck_sd: dict[str, float] = field(default_factory=dict)
    bond: float = BOND_CA
    n_structures: int = 0
    convention: str = "endpoint"
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

    def ensemble(self, chain: str, region: str, mhc_class: str) -> np.ndarray:
        """The ``(S, G, G)`` surface stack for one key, falling back across MHC class."""
        if not self.surfaces:
            raise ValueError("this prior kept no per-structure surfaces; refit with "
                             "fit_geom_prior(..., keep_surfaces=True) or load the bundled "
                             f"{BUNDLED_ENSEMBLE!r}")
        k = _key(chain, region, mhc_class)
        if k not in self.surfaces:
            alt = [c for c in self.surfaces if c.startswith(f"{chain}:{region}|")]
            if not alt:
                raise KeyError(f"no surface for {chain}:{region} at any MHC class; "
                               f"the prior carries {self.keys()}")
            k = alt[0]
        return np.asarray(self.surfaces[k], float)

    def to_json(self, path: str | Path) -> None:
        """Mean and spread only. The surface stack does not fit in JSON --- use :meth:`to_npz`."""
        d = {k: v for k, v in asdict(self).items() if k != "surfaces"}
        Path(path).write_text(json.dumps(d))

    @classmethod
    def from_json(cls, path: str | Path) -> "GeomPrior":
        return cls(**json.loads(Path(path).read_text()))

    def to_npz(self, path: str | Path) -> None:
        """Everything, surfaces included. Distances are stored to 0.01 Å as uint16 hundredths.

        A stack of 2,154 surfaces at grid 32 is 2.2 M numbers; as JSON that is tens of megabytes
        and as float64 it is 17 MB. Hundredths of an angstrom are four decimal places finer than
        the 0.5 Å bin the model's distance term reads, so the quantisation is below the model's
        own resolution, and the clip at 655 Å is beyond any Cα pair in a ternary complex.
        """
        meta = {k: v for k, v in asdict(self).items() if k != "surfaces"}
        arrays = {f"surface::{k}": np.clip(np.asarray(v) * 100.0, 0, 65535).astype(np.uint16)
                  for k, v in (self.surfaces or {}).items()}
        np.savez_compressed(path, meta=np.array(json.dumps(meta)), **arrays)

    @classmethod
    def from_npz(cls, path: str | Path) -> "GeomPrior":
        with np.load(path, allow_pickle=False) as z:
            meta = json.loads(str(z["meta"]))
            surf = {k.split("::", 1)[1]: z[k].astype(float) / 100.0
                    for k in z.files if k.startswith("surface::")}
        return cls(surfaces=surf or None, **meta)

    @classmethod
    def bundled(cls, name: str = BUNDLED_PRIOR) -> "GeomPrior":
        """A prior shipped in the wheel, by name. See :data:`BUNDLED_PRIOR`.

        Resolves ``{name}.npz`` first and ``{name}.json`` second, so an ensemble prior and a
        mean-only one are asked for the same way.
        """
        from importlib.resources import files
        for suffix in (".npz", ".json"):
            res = files("tcren.data").joinpath(f"{name}{suffix}")
            if res.is_file():
                if suffix == ".npz":
                    import io
                    return cls.from_npz(io.BytesIO(res.read_bytes()))
                return cls(**json.loads(res.read_text(encoding="utf-8")))
        raise FileNotFoundError(f"no bundled geometry prior named {name!r}")


def fit_geom_prior(structures, *, grid: int = 12, regions=REGIONS_MODELLED,
                   keep_surfaces: bool = False, convention: str = "endpoint",
                   notes: str = "") -> GeomPrior:
    """Average the receptor:peptide Cα surfaces of a structure set in graphon coordinates.

    Args:
        structures: An iterable of chain-typed structures — pass
            :func:`tcren.annotation.batch.iter_annotated_set` rather than annotating in a loop.
        grid: Surface resolution per axis.
        regions: Receptor regions to model.
        keep_surfaces: Keep one registered surface per contributing crystal rather than only
            their mean, which is what :func:`predict_contact_map` needs. Costs about 2 MB on
            disk at ``grid=32`` and is the difference between reproducing the contact count and
            reading a fifth of it --- see :attr:`GeomPrior.surfaces`.
        convention: Graphon convention to *register* under. ``"endpoint"`` reproduces the shipped
            priors; ``"interior"`` is adjoint to the read-back and is the one to prefer for a new
            fit --- under ``"endpoint"`` a round trip contracts a loop's footprint by
            :math:`(L-1)/(L+1)`, a fixed 1.29 Å at the termini of a 6-mer that no grid refinement
            removes. See :func:`tcren.geometry.graphon_weights`.
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
                acc.setdefault(_key(chain, region, cls), []).append(_resample(d, grid, convention))
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
        surfaces={k: np.stack(v) for k, v in acc.items()} if keep_surfaces else None,
        neck_mean={k: float(np.mean(v)) for k, v in necks.items()},
        neck_sd={k: float(np.std(v)) for k, v in necks.items()},
        n_structures=seen, convention=convention, notes=notes)


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
    m, sd = prior.surface(chain, region, mhc_class)
    return read_back(m, loop_len, pep_len), read_back(sd, loop_len, pep_len)


def predicted_ensemble(prior: GeomPrior, chain: str, region: str, *, loop_len: int,
                       pep_len: int, mhc_class: str = "MHCI") -> np.ndarray:
    r"""``(S, loop_len, pep_len)`` --- every member surface read back at the target lengths.

    One :func:`tcren.geometry.read_back` call over the whole ``(S, G, G)`` stack, not one per
    member: the member axis is carried through untouched, so it interpolates only in :math:`(u, v)`.
    Each slice is one crystal's geometry carried to these lengths, joint intact.
    """
    return read_back(prior.ensemble(chain, region, mhc_class), loop_len, pep_len)


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
    from .model import AA

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
            if seq[i] not in AA:
                continue        # a non-standard residue has no code in the model's alphabet
            for j in range(M):
                if peptide[j] not in AA:
                    continue
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


def predict_contact_map(prior: GeomPrior, loops: dict[str, str], peptide: str, model, *,
                        mhc_class: str = "MHCI", radius: float = 15.0,
                        pdb_id: str = "predicted", members: int | None = None,
                        **gibbs) -> pl.DataFrame:
    r"""The contact map of a receptor:peptide pair, averaged over GEOMETRIES rather than over
    surfaces --- from sequence and lengths alone.

    :func:`predicted_sites` reads the mean surface and scores one map. That is the wrong order of
    operations: averaging first and scoring second asks the model about a geometry no chain can
    adopt, since a mean cannot reproduce a close approach and consecutive Cα atoms are a fixed
    virtual bond apart. Here every member of :attr:`GeomPrior.surfaces` is read back at the target
    lengths, the map is scored under *that* geometry, and the contact probabilities are averaged:

    .. math:: \bar p_a = \frac{1}{S}\sum_{s} p_a\big(d^{(s)}\big)\,
              \mathbb{1}\big[d^{(s)}_a \le r\big]

    with :math:`S` the members, :math:`d^{(s)}` member :math:`s`'s distance matrix and :math:`r`
    the availability radius. A pair beyond :math:`r` in a member contributes nothing from it, so
    ``p_avail`` --- the fraction of members placing the pair in reach --- is reported beside the
    probability.

    Held out by structure over 364 Native2026 crystals this recovers the contact count the mean
    surface loses: mean :math:`\sum_a \bar p_a` of 24.62 against 21.02 contacts observed
    (ratio 1.17), where the mean surface reads 9.80 (ratio 0.47); median contact ROC-AUC 0.9325
    against 0.9242. Both the level and the ordering improve.

    Args:
        prior: A :class:`GeomPrior` carrying surfaces --- :data:`BUNDLED_ENSEMBLE`, or a fit with
            ``keep_surfaces=True``.
        loops: ``{"TRA:CDR3": sequence, ...}`` --- see :func:`loops_from_structure`.
        peptide: The peptide sequence.
        model: A fitted :class:`tcren.potts.PottsModel`.
        mhc_class: ``"MHCI"`` or ``"MHCII"``.
        radius: Availability radius, Å, applied to each member's own distance.
        pdb_id: Identifier written into the frame.
        members: Score only the first ``members`` surfaces. ``None`` (default) uses all of them;
            the cost is linear and about 0.16 s per member on a map of this size.
        **gibbs: Passed to :func:`tcren.potts.contact_probabilities`.

    Returns:
        One row per site on the union of the members' site sets, carrying ``d_ca`` (the member
        mean), ``d_ca_sd``, ``p_avail``, and the three averaged probabilities.
    """
    from .score import contact_probabilities

    rows = _ensemble_rows(prior, loops, peptide, mhc_class=mhc_class, radius=radius,
                          pdb_id=pdb_id, members=members)
    if rows.is_empty():
        return rows
    # ONE call with every member as its own `pdb.id`, never a loop over members: the sampler
    # already splits by structure and threads across them, and a member is exactly a structure
    # to it. S separate calls is the same arithmetic at a fraction of the throughput.
    scored = contact_probabilities(rows.drop("member"), model, **gibbs)
    key = ["chain.rec", "region.rec", "pos.rec", "pos.par"]
    n = rows["member"].n_unique()
    return (scored.with_columns(rows["member"])
            .group_by(key)
            .agg(pl.col("aa.rec").first(), pl.col("aa.par").first(), pl.col("par.len").first(),
                 pl.col("role.par").first(), pl.col("cls").first(),
                 pl.col("d_ca").mean().alias("d_ca"), pl.col("d_ca").std().alias("d_ca_sd"),
                 (pl.len() / n).alias("p_avail"),
                 # summed and divided by every member, not averaged over the members that kept
                 # the pair: a member placing it out of reach votes zero, it does not abstain
                 (pl.col("p_independent").sum() / n).alias("p_independent"),
                 (pl.col("p_model").sum() / n).alias("p_model"),
                 (pl.col("p_conditional").sum() / n).alias("p_conditional"),
                 pl.col("feasible").any())
            .with_columns(pl.lit(pdb_id).alias("pdb.id"), pl.lit(None, pl.Float64).alias("d_heavy"),
                          pl.lit(0.0).alias("sigma"))
            .select(list(SITE_COLUMNS) + ["d_ca_sd", "feasible", "p_avail",
                                          "p_independent", "p_model", "p_conditional"])
            .sort(["chain.rec", "region.rec", "pos.rec", "pos.par"]))


def _ensemble_rows(prior: GeomPrior, loops: dict[str, str], peptide: str, *,
                   mhc_class: str, radius: float, pdb_id: str,
                   members: int | None) -> pl.DataFrame:
    """One row per (member, in-reach pair), each member tagged as its own ``pdb.id``."""
    from ..refine.anchors import predict_anchors
    from .model import AA

    try:
        anchors = set(predict_anchors(peptide).anchors)
    except Exception:
        anchors = set()
    M = len(peptide)
    out = []
    for key, seq in sorted(loops.items()):
        chain, region = key.split(":")
        L = len(seq)
        if L < 1 or M < 1:
            continue
        d = predicted_ensemble(prior, chain, region, loop_len=L, pep_len=M, mhc_class=mhc_class)
        if members is not None:
            d = d[:members]
        cells = [(i, j) for i in range(L) if seq[i] in AA for j in range(M) if peptide[j] in AA]
        for k in range(d.shape[0]):
            for i, j in cells:
                x = float(d[k, i, j])
                if x > radius:
                    continue
                out.append({"member": k, "pdb.id": f"{pdb_id}#{k:04d}", "aa.rec": seq[i],
                            "chain.rec": chain, "region.rec": region, "pos.rec": i,
                            "aa.par": peptide[j], "pos.par": j, "par.len": M,
                            "role.par": "anchor" if j in anchors else "tcr_facing",
                            "cls": "peptide", "d_heavy": None, "d_ca": x, "sigma": 0.0,
                            "feasible": feasible(prior, chain, region, L, x)})
    if not out:
        return pl.DataFrame()
    return pl.DataFrame(out).select(["member"] + list(SITE_COLUMNS) + ["feasible"])
