"""Wide single-structure descriptor family over the interface maps --- the search space.

:mod:`tcren.pose` ships a curated handful. This module generates the whole family they were drawn
from, so a metric can be *designed* against a benchmark rather than guessed: every pairing of

* a **map** --- ``d1`` closest heavy-atom, ``d2`` Cbeta, ``d3`` Calpha;
* a **valuation** on that map --- distance, contact indicator, the double-centred pair term ``J``,
  the raw potential entry, heavy-atom multiplicity, residue chemistry class;
* a **reduction** --- moments, quantiles, rank correlations between two valuations, OLS slopes,
  tercile contrasts, fractions, per-residue degree structure, concentration/participation ratios,
  per-peptide-position and per-CDR-loop shares;
* a **scope** --- TCR:peptide or TCR:MHC, the contact set or the Calpha approach shell, and the
  CDR1/CDR2/CDR3 partition of the receptor side.

Everything is computed from ONE ``multi_contacts`` call plus the potential, so the whole family
costs about what a single contact map costs. All are single-structure by construction: no cohort,
no labels, no generator-reported quantity.

Naming is ``<what>_<scope>``: scope ``tp`` = TCR:peptide, ``tm`` = TCR:MHC, and a loop suffix where
the receptor side is partitioned.
"""

from __future__ import annotations

import numpy as np
import polars as pl

from .contacts.definitions import ContactDefinition, multi_contacts
from .contacts.table import residue_annotation
from .pose import _CA_CLOSE, _CB_CLOSE, _double_centred, _pair_j, _spearman
from .structure.model import MHC_TYPES, PEPTIDE_TYPE, RECEPTOR_TYPES, Structure

__all__ = ["pose_descriptors_full", "P_TERMS", "p_score"]

#: The designed pose metric ``P``: ``(descriptor, sign)`` pairs, sign oriented so that **higher is
#: more binder-like**. Selected by greedy forward selection on macro PR over the **22 well-powered
#: cohorts** of the balanced VDJdb binder benchmark (1,089 AlphaFold models, 523 real / 566 mock,
#: 7 alleles) --- the widest receptor-ranking surface available, and the deposit behind the
#: template-split table. Candidates were first de-duplicated at |Spearman| <= 0.8 so no two terms
#: are restatements of each other.
#:
#: Every term is single-structure and coordinate-only: contact-distance spread and shape, the
#: Cbeta-vs-Calpha agreement of the two maps, the slope of contact favourability against Cbeta
#: distance, receptor- and partner-side partner counts, interface size, MHC-side composition, and
#: the CDR footprint offset. No cohort, no label, no generator-reported quantity enters it.
P_TERMS = (
    ("d2_fav_tp_slope", +1.0),     # favourability rises as Cbeta closes
    ("d1_tm_sd", -1.0),            # tight, uniform TCR:MHC contact shell
    ("degT_tp_mean", +1.0),        # receptor residues engaged, on average
    ("d1_tp_sd", -1.0),            # uniform TCR:peptide contact distances
    ("d3_d2_tp_rho", -1.0),        # backbone and side-chain maps not merely redundant
    ("n_contacts_tp", +1.0),       # interface size
    ("degT_tm_frac_le3", -1.0),
    ("d2_tp_skew", +1.0),
    ("degP_tm_sd", -1.0),
    ("d1_tp_q90", -1.0),           # no long tail of marginal contacts
    ("fracP_tiny_tm", +1.0),
    ("d3_d2_tp_slope", -1.0),
    ("offset", -1.0),              # CDR footprint centred on the peptide
)

_REP = 18.0
_KEY = ["chain.id.from", "residue.index.from", "chain.id.to", "residue.index.to"]

# Residue chemistry classes; a contact is described by the pair of classes it joins.
_CLASS = {
    "hyd": set("AVILMFWCP"), "aro": set("FWYH"), "pos": set("KRH"),
    "neg": set("DE"), "pol": set("STNQY"), "tiny": set("GAS"),
}


def _cls_mask(aa: list[str], name: str) -> np.ndarray:
    m = _CLASS[name]
    return np.array([a in m for a in aa], bool)


def _moments(x: np.ndarray, tag: str, out: dict) -> None:
    """mean/sd/skew/kurtosis and three quantiles of one valuation."""
    x = x[np.isfinite(x)]
    if len(x) < 3:
        return
    mu, sd = float(x.mean()), float(x.std())
    out[f"{tag}_mean"], out[f"{tag}_sd"] = mu, sd
    if sd > 1e-12:
        z = (x - mu) / sd
        out[f"{tag}_skew"] = float((z ** 3).mean())
        out[f"{tag}_kurt"] = float((z ** 4).mean())
    for q in (10, 50, 90):
        out[f"{tag}_q{q}"] = float(np.percentile(x, q))


def _assoc(x: np.ndarray, v: np.ndarray, tag: str, out: dict) -> None:
    """Rank correlation, OLS slope and tercile contrast of ``v`` against ``x``."""
    ok = np.isfinite(x) & np.isfinite(v)
    if ok.sum() < 4:
        return
    xs, vs = x[ok], v[ok]
    out[f"{tag}_rho"] = _spearman(xs, vs)
    if np.std(xs) > 1e-12:
        out[f"{tag}_slope"] = float(np.polyfit(xs, vs, 1)[0])
    lo, hi = np.percentile(xs, [33.3, 66.7])
    a, b = vs[xs <= lo], vs[xs >= hi]
    if len(a) and len(b):
        out[f"{tag}_terc"] = float(a.mean() - b.mean())


def _degree(keys: np.ndarray, tag: str, out: dict) -> None:
    """Partner-count structure on one side of the interface."""
    if not len(keys):
        return
    _, deg = np.unique(keys, return_counts=True)
    deg = deg.astype(float)
    out[f"{tag}_max"] = float(deg.max())
    out[f"{tag}_mean"] = float(deg.mean())
    out[f"{tag}_sd"] = float(deg.std())
    out[f"{tag}_even"] = float(deg.sum() ** 2 / (len(deg) * (deg ** 2).sum()))
    out[f"{tag}_frac_le3"] = float((deg <= 3).mean())
    p = deg / deg.sum()
    out[f"{tag}_entropy"] = float(-(p * np.log(p + 1e-12)).sum() / np.log(len(p) + 1e-12)) if len(p) > 1 else 0.0


def _concentration(w: np.ndarray, tag: str, out: dict) -> None:
    """How much of a weight (|J|, contact count) sits on the few largest entries."""
    w = np.abs(w[np.isfinite(w)])
    if len(w) < 2 or w.sum() <= 0:
        return
    s = np.sort(w)[::-1]
    out[f"{tag}_pr"] = float(w.sum() ** 2 / (len(w) * (w ** 2).sum()))
    out[f"{tag}_top1"] = float(s[0] / w.sum())
    out[f"{tag}_top3"] = float(s[:3].sum() / w.sum())


def _layers(structure: Structure, cutoff: float, partner) -> pl.DataFrame:
    """d1/d2/d3 on one row per pair, with receptor region and partner position attached."""
    stacked = multi_contacts(structure, ContactDefinition(d1=cutoff, d2=_REP, d3=_REP))
    ann = residue_annotation(structure)
    ctype = dict(zip(ann["chain.id"].to_list(), ann["chain.type"].to_list()))
    tcr, par = list(RECEPTOR_TYPES), list(partner)
    stacked = stacked.with_columns(
        pl.col("chain.id.from").replace_strict(ctype, default=None).alias("tf"),
        pl.col("chain.id.to").replace_strict(ctype, default=None).alias("tt"),
    )
    fwd = pl.col("tf").is_in(tcr) & pl.col("tt").is_in(par)
    rev = pl.col("tf").is_in(par) & pl.col("tt").is_in(tcr)
    stacked = stacked.filter(fwd | rev).with_columns(
        pl.when(fwd).then(pl.col("residue.aa.from")).otherwise(pl.col("residue.aa.to")).alias("aa_t"),
        pl.when(fwd).then(pl.col("residue.aa.to")).otherwise(pl.col("residue.aa.from")).alias("aa_p"),
        pl.when(fwd).then(pl.col("chain.id.from")).otherwise(pl.col("chain.id.to")).alias("ct"),
        pl.when(fwd).then(pl.col("residue.index.from")).otherwise(pl.col("residue.index.to")).alias("rt"),
        pl.when(fwd).then(pl.col("chain.id.to")).otherwise(pl.col("chain.id.from")).alias("cp"),
        pl.when(fwd).then(pl.col("residue.index.to")).otherwise(pl.col("residue.index.from")).alias("rp"),
    )
    if stacked.is_empty():
        return stacked
    wide = (stacked.filter(pl.col("layer") == "d1")
            .select(*_KEY, "aa_t", "aa_p", "ct", "rt", "cp", "rp",
                    pl.col("dist").alias("d1"),
                    pl.col("n_atom_contacts").alias("nat") if "n_atom_contacts" in stacked.columns
                    else pl.lit(1).alias("nat")))
    for layer in ("d2", "d3"):
        wide = wide.join(stacked.filter(pl.col("layer") == layer)
                         .select(*_KEY, pl.col("dist").alias(layer)), on=_KEY, how="full", coalesce=True)
    # receptor region label, for the CDR partition
    reg = ann.select(pl.col("chain.id").alias("ct"), pl.col("residue.index").alias("rt"),
                     pl.col("region.type").alias("region"))
    return wide.join(reg, on=["ct", "rt"], how="left")


def _block(w: pl.DataFrame, jmat, jindex, cutoff: float, sfx: str, out: dict) -> None:
    """Every reduction, for one interface's layer frame."""
    con = w.filter(pl.col("d1").is_not_null())
    out[f"n_contacts{sfx}"] = float(con.height)
    shell8 = w.filter(pl.col("d3").is_not_null() & (pl.col("d3") <= _CA_CLOSE))
    out[f"n_ca8{sfx}"] = float(shell8.height)
    if shell8.height:
        out[f"frac_ca8_engaged{sfx}"] = float(shell8["d1"].is_not_null().mean())
    cb8 = w.filter(pl.col("d2").is_not_null() & (pl.col("d2") <= _CB_CLOSE))
    out[f"n_cb8{sfx}"] = float(cb8.height)
    if cb8.height:
        out[f"frac_cb8_engaged{sfx}"] = float(cb8["d1"].is_not_null().mean())
    if not con.height:
        return

    aa_t, aa_p = con["aa_t"].to_list(), con["aa_p"].to_list()
    j = _pair_j(aa_t, aa_p, jmat, jindex)
    fav = -j
    d1 = con["d1"].to_numpy(); d2 = con["d2"].to_numpy(); d3 = con["d3"].to_numpy()
    nat = con["nat"].to_numpy().astype(float)
    marg = cutoff - d1

    for nm, v in (("d1", d1), ("d2", d2), ("d3", d3), ("J", j), ("nat", nat), ("marg", marg)):
        _moments(v, f"{nm}{sfx}", out)
    # every distance map against the chemistry, and the maps against each other
    for xn, x in (("d1", d1), ("d2", d2), ("d3", d3)):
        _assoc(x, fav, f"{xn}_fav{sfx}", out)
    _assoc(d3, d1, f"d3_d1{sfx}", out)
    _assoc(d3, d2, f"d3_d2{sfx}", out)
    _assoc(d2, d1, f"d2_d1{sfx}", out)
    _assoc(nat, fav, f"nat_fav{sfx}", out)
    _assoc(d1, nat, f"d1_nat{sfx}", out)

    ok = np.isfinite(j)
    if ok.sum() >= 3:
        out[f"frac_fav{sfx}"] = float((j[ok] < 0).mean())
        med = np.median(d1[ok])
        near, far = j[ok][d1[ok] <= med], j[ok][d1[ok] > med]
        if len(near): out[f"frac_fav_near{sfx}"] = float((near < 0).mean())
        if len(far):  out[f"frac_fav_far{sfx}"] = float((far < 0).mean())
        out[f"J_sum{sfx}"] = float(j[ok].sum())
        out[f"J_per_contact{sfx}"] = float(j[ok].mean())
        _concentration(j[ok], f"Jconc{sfx}", out)
    _concentration(nat, f"natconc{sfx}", out)

    # side-chain lean: Cbeta closing relative to Calpha
    both = np.isfinite(d2) & np.isfinite(d3)
    if both.sum():
        out[f"lean_mean{sfx}"] = float((d3[both] - d2[both]).mean())
        out[f"lean_frac{sfx}"] = float((d2[both] < d3[both]).mean())

    _degree(con["ct"].to_numpy().astype(str) + ":" + con["rt"].to_numpy().astype(str),
            f"degT{sfx}", out)
    _degree(con["cp"].to_numpy().astype(str) + ":" + con["rp"].to_numpy().astype(str),
            f"degP{sfx}", out)

    # residue chemistry: which classes the interface is built from, and complementarity
    for c in _CLASS:
        out[f"fracT_{c}{sfx}"] = float(_cls_mask(aa_t, c).mean())
        out[f"fracP_{c}{sfx}"] = float(_cls_mask(aa_p, c).mean())
    pos_t, neg_t = _cls_mask(aa_t, "pos"), _cls_mask(aa_t, "neg")
    pos_p, neg_p = _cls_mask(aa_p, "pos"), _cls_mask(aa_p, "neg")
    out[f"salt_ok{sfx}"] = float((pos_t & neg_p).mean() + (neg_t & pos_p).mean())
    out[f"salt_bad{sfx}"] = float((pos_t & pos_p).mean() + (neg_t & neg_p).mean())
    hyd_t, hyd_p = _cls_mask(aa_t, "hyd"), _cls_mask(aa_p, "hyd")
    out[f"hyd_pair{sfx}"] = float((hyd_t & hyd_p).mean())
    out[f"hyd_mismatch{sfx}"] = float((hyd_t & ~hyd_p).mean() + (~hyd_t & hyd_p).mean())

    # how contact mass spreads over partner positions, and how uneven the energy is across them
    pk = con["cp"].to_numpy().astype(str) + ":" + con["rp"].to_numpy().astype(str)
    uq, inv = np.unique(pk, return_inverse=True)
    if len(uq) > 1 and ok.sum():
        share = np.bincount(inv, minlength=len(uq)).astype(float); share /= share.sum()
        out[f"posmass_top1{sfx}"] = float(share.max())
        out[f"posmass_entropy{sfx}"] = float(-(share*np.log(share+1e-12)).sum()/np.log(len(uq)))
        jj = np.where(np.isfinite(j), j, 0.0)
        per = np.bincount(inv, weights=jj, minlength=len(uq))
        out[f"posJ_sd{sfx}"] = float(per.std())
        out[f"posJ_min{sfx}"] = float(per.min())
        out[f"posJ_range{sfx}"] = float(per.max() - per.min())

    # the CDR partition of the receptor side
    regs = con["region"].to_list()
    for loop in ("CDR1", "CDR2", "CDR3"):
        m = np.array([r == loop for r in regs], bool)
        out[f"share_{loop}{sfx}"] = float(m.mean())
        if m.sum() and ok.sum():
            out[f"J_{loop}{sfx}"] = float(np.nansum(j[m]))


def pose_descriptors_full(structure: Structure, potential=None, cutoff: float = 5.0) -> dict:
    """The full single-structure descriptor family (see the module docstring).

    Args:
        structure: chain-typed, region-annotated complex.
        potential: residue-pair potential supplying ``J``; defaults to bundled TCRen2.
        cutoff: heavy-atom contact cutoff (A).

    Returns:
        A flat ``{name: value}`` dict; absent reductions are simply missing keys, never zeros.
    """
    if potential is None:
        from .potential import tcren2
        potential = tcren2()
    jmat, jindex = _double_centred(potential)
    out: dict[str, float] = {}
    for sfx, partner in (("_tp", (PEPTIDE_TYPE,)), ("_tm", MHC_TYPES)):
        w = _layers(structure, cutoff, partner)
        if w.height:
            _block(w, jmat, jindex, cutoff, sfx, out)
    return out


def p_score(table, reference=None) -> np.ndarray:
    """The designed pose metric ``P`` --- equal-weight mean of ``z(sign * descriptor)`` over
    :data:`P_TERMS`. Higher is more binder-like.

    Fit-free in the same sense as :func:`tcren.cohort.q_score`: the signs are fixed by
    :data:`P_TERMS` and nothing is regressed on labels at score time. The **subset** was chosen on
    the balanced VDJdb benchmark, so ``P`` is a designed metric, not a derived one --- see
    :data:`P_TERMS` for the surface it was designed on.

    Measured, macro over the 22 cohorts of that benchmark: ROC 0.721 / PR 0.731 / P@10 0.904,
    against ipTM 0.592 / 0.606 / 0.770 and ``S`` 0.571 / 0.575 / 0.719. On the template-free
    cohorts, where every other score falls to chance, ``P`` reads ROC 0.691 / PR 0.712 against
    ``S`` 0.502 / 0.519. On TCRvdb it composes rather than replaces:
    ``z(P) + z(S) + z(ipTM)`` reaches 0.815 / 0.841 / 0.975 against ``S`` 0.799 / 0.817 / 0.923.

    Args:
        table: mapping / polars / pandas frame carrying the :data:`P_TERMS` descriptor columns,
            e.g. rows built from :func:`pose_descriptors_full`.
        reference: cohort defining each term's mean and sd; ``None`` standardizes against the input.

    Returns:
        One value per row.
    """
    from .cohort import zscore

    def col(t, name):
        if hasattr(t, "columns") and not isinstance(t, dict):
            return np.asarray(t[name].to_numpy() if hasattr(t[name], "to_numpy") else t[name], float)
        return np.asarray(t[name], float)

    parts = []
    for name, sign in P_TERMS:
        ref = None if reference is None else sign * col(reference, name)
        parts.append(zscore(sign * col(table, name), ref))
    return np.mean(parts, axis=0)
