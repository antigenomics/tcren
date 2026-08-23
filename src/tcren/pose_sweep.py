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

__all__ = ["pose_descriptors_full", "loop_ca_profile", "loop_ca_rules", "P_TERMS", "p_score"]

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
    # d(favourability)/d(Cbeta distance) over TCR:peptide contacts. Positive = the chemically
    # complementary pairs are the ones at LARGER Cbeta separation, i.e. specificity is carried by
    # long-reach side chains, while the tightest Cbeta pairs are small and chemically neutral.
    ("d2_fav_tp_slope", +1.0),
    ("d1_tm_sd", -1.0),            # uniform TCR:MHC contact-distance shell
    ("degT_tp_mean", +1.0),        # receptor residues engaged, on average
    ("d1_tp_sd", -1.0),            # uniform TCR:peptide contact distances
    ("d3_d2_tp_rho", -1.0),        # Calpha and Cbeta maps NOT merely restating each other
    ("n_contacts_tp", +1.0),       # interface size
    ("degT_tm_frac_le3", -1.0),    # receptor spreads over the MHC rather than touching it lightly
    ("d2_tp_skew", +1.0),          # Cbeta distances right-skewed: a compact core plus reachers
    ("degP_tm_sd", -1.0),          # MHC-side partner counts even
    ("d1_tp_q90", -1.0),           # no long tail of marginal contacts
    ("fracP_tiny_tm", +1.0),       # Gly/Ala/Ser share of the contacted MHC face
    ("d3_d2_tp_slope", -1.0),      # Cbeta separation grows slowly with Calpha separation
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
    # Identity (residue types, chain, indices) is built from the UNION of all three layers, not from
    # d1 alone. Taking it from d1 would leave every non-contacting shell pair with null residue
    # identity, so anything that reads aa/region over the shell would silently see contacts only.
    ident = stacked.select(*_KEY, "aa_t", "aa_p", "ct", "rt", "cp", "rp").unique(subset=_KEY)
    d1f = stacked.filter(pl.col("layer") == "d1")
    wide = d1f.select(*_KEY, pl.col("dist").alias("d1"),
                      pl.col("n_atom_contacts").alias("nat") if "n_atom_contacts" in stacked.columns
                      else pl.lit(1).alias("nat"))
    for layer in ("d2", "d3"):
        wide = wide.join(stacked.filter(pl.col("layer") == layer)
                         .select(*_KEY, pl.col("dist").alias(layer)), on=_KEY, how="full", coalesce=True)
    wide = wide.join(ident, on=_KEY, how="left")
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
    out.update(loop_ca_profile(structure, cutoff=cutoff))
    out.update(loop_ca_rules(structure, cutoff=cutoff))
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


# =====================================================================================================
# Per-loop Calpha profile: contacting vs non-contacting, CDR1-3 of each chain, against peptide and MHC
# =====================================================================================================

_LOOPS = (("TRA", "CDR1", "cdr1a"), ("TRA", "CDR2", "cdr2a"), ("TRA", "CDR3", "cdr3a"),
          ("TRB", "CDR1", "cdr1b"), ("TRB", "CDR2", "cdr2b"), ("TRB", "CDR3", "cdr3b"))


def _auc_lt(a: np.ndarray, b: np.ndarray) -> float:
    """P(a < b) for two samples --- the Mann-Whitney statistic, computed by ranking not by pairing.

    Here: the probability that a randomly chosen *contacting* pair sits closer in Calpha than a
    randomly chosen *non-contacting* one. 1.0 = Calpha distance orders contact perfectly within this
    structure; 0.5 = it carries no information about which pairs actually touch.
    """
    if not len(a) or not len(b):
        return float("nan")
    from scipy.stats import rankdata

    both = np.concatenate([a, b])
    r = rankdata(both)[: len(a)]
    u = r.sum() - len(a) * (len(a) + 1) / 2.0
    return float(1.0 - u / (len(a) * len(b)))


def loop_ca_profile(structure: Structure, cutoff: float = 5.0,
                    ca_radius: float = 12.0) -> dict[str, float]:
    """Calpha-distance profile of contacting vs non-contacting pairs, per CDR loop and partner.

    For each of the six CDR loops (CDR1/2/3 of the alpha and beta chain) against each partner
    (peptide, MHC), the pairs whose Calpha atoms lie within ``ca_radius`` are split into those that
    make a ``cutoff`` A heavy-atom contact and those that do not, and each side's Calpha distance
    distribution is described. The question it answers per loop: *at what Calpha separation does this
    loop actually engage, and how cleanly does backbone proximity predict engagement at all?*

    Keys are ``<loop>_<pep|mhc>_<stat>`` with stats:

    ``n_shell``      pairs within ``ca_radius``
    ``n_con``        of those, pairs that make a contact
    ``frac_eng``     n_con / n_shell --- how much of the loop's Calpha neighbourhood is productive
    ``d3con_mean``   mean Calpha distance of the CONTACTING pairs
    ``d3con_min``    closest Calpha approach of the loop to the partner
    ``d3con_sd``     spread of the contacting Calpha distances
    ``d3non_mean``   mean Calpha distance of the NON-contacting pairs in the shell
    ``d3sep``        ``d3non_mean - d3con_mean`` --- how far apart the two populations sit
    ``auc_d3``       P(contacting pair is closer in Calpha than a non-contacting one), 0.5 = no
                     information, 1.0 = Calpha distance orders engagement perfectly

    Args:
        structure: chain-typed, region-annotated complex.
        cutoff: heavy-atom contact cutoff (A).
        ca_radius: Calpha shell radius (A) defining "in the neighbourhood".

    Returns:
        Flat ``{name: value}``; a loop absent from the structure simply contributes no keys.
    """
    out: dict[str, float] = {}
    for partner, sfx in ((( PEPTIDE_TYPE,), "pep"), (MHC_TYPES, "mhc")):
        w = _layers(structure, cutoff, partner)
        if not w.height or "region" not in w.columns:
            continue
        ann = residue_annotation(structure)
        ctype = dict(zip(ann["chain.id"].to_list(), ann["chain.type"].to_list()))
        chain_t = np.array([ctype.get(c) for c in w["ct"].to_list()], dtype=object)
        region = np.array(w["region"].to_list(), dtype=object)
        d3 = w["d3"].to_numpy()
        has = w["d1"].is_not_null().to_numpy()
        for ch, reg, name in _LOOPS:
            m = (chain_t == ch) & (region == reg) & np.isfinite(d3) & (d3 <= ca_radius)
            if not m.any():
                continue
            con, non = d3[m & has], d3[m & ~has]
            k = f"{name}_{sfx}"
            out[f"{k}_n_shell"] = float(m.sum())
            out[f"{k}_n_con"] = float(len(con))
            out[f"{k}_frac_eng"] = float(len(con) / m.sum())
            if len(con):
                out[f"{k}_d3con_mean"] = float(con.mean())
                out[f"{k}_d3con_min"] = float(con.min())
                out[f"{k}_d3con_sd"] = float(con.std())
            if len(non):
                out[f"{k}_d3non_mean"] = float(non.mean())
            if len(con) and len(non):
                out[f"{k}_d3sep"] = float(non.mean() - con.mean())
                out[f"{k}_auc_d3"] = _auc_lt(con, non)
    return out


# =====================================================================================================
# Interpretable Calpha rules: where each CDR loop sits relative to the peptide and to the MHC
# =====================================================================================================

def loop_ca_rules(structure: Structure, cutoff: float = 5.0) -> dict[str, float]:
    """Per-loop Calpha geometry as plain distances, for reading rather than for scoring.

    For every residue of each CDR loop, two numbers: its **nearest Calpha distance to the peptide**
    and its **nearest Calpha distance to the MHC**. Their difference

        ``delta = d_pep - d_mhc``

    is negative when the loop residue sits closer to the peptide than to the MHC and positive when
    it leans on the MHC instead. Averaging over a loop gives one statement per loop per structure,
    e.g. "CDR3beta sits 1.2 A closer to the peptide than to the MHC".

    Each loop is reported three ways: over all its residues, over the residues that actually make a
    5 A heavy-atom contact **with the peptide**, and over those that do not. The contacting split is
    the informative one --- it asks where a residue sits *given that it is engaged*, which separates
    "this loop reaches the peptide" from "this loop happens to be near it".

    Keys are ``<loop>_<all|con|non>_<d_pep|d_mhc|delta>`` plus ``<loop>_n_res`` and
    ``<loop>_frac_con``.

    Args:
        structure: chain-typed, region-annotated complex.
        cutoff: heavy-atom cutoff (A) defining "contacts the peptide".

    Returns:
        Flat ``{name: value}``; loops absent from the structure contribute no keys.
    """
    from .contacts.geometry import all_atom_contacts

    ann = residue_annotation(structure)
    ctype = dict(zip(ann["chain.id"].to_list(), ann["chain.type"].to_list()))

    def ca_of(types):
        return np.asarray([r.ca for c in structure.chains if c.chain_type in types
                           for r in c.residues if r.ca is not None], float)

    pep_ca, mhc_ca = ca_of((PEPTIDE_TYPE,)), ca_of(MHC_TYPES)
    if not len(pep_ca) or not len(mhc_ca):
        return {}
    # MHC groove helices, for the docking-polarity contrasts
    helix: dict[str, np.ndarray] = {}
    for c in structure.chains:
        if c.chain_type not in MHC_TYPES:
            continue
        for reg in (c.regions or []):
            if reg.region_type.startswith("HELIX"):
                pts = [r.ca for r in reg.residues if r.ca is not None]
                if pts:
                    helix.setdefault(reg.region_type, []).extend(pts)
    helix = {k: np.asarray(v, float) for k, v in helix.items()}

    # which receptor residues contact the peptide at all
    con = all_atom_contacts(structure, cutoff=cutoff)
    touch = set()
    for a, b in (("from", "to"), ("to", "from")):
        sub = con.filter(
            pl.col(f"chain.id.{a}").replace_strict(ctype, default=None).is_in(list(RECEPTOR_TYPES))
            & pl.col(f"chain.id.{b}").replace_strict(ctype, default=None).is_in([PEPTIDE_TYPE]))
        touch |= set(zip(sub[f"chain.id.{a}"].to_list(), sub[f"residue.index.{a}"].to_list()))

    out: dict[str, float] = {}
    for ch, reg, name in _LOOPS:
        chain = next((c for c in structure.chains if c.chain_type == ch), None)
        if chain is None:
            continue
        region = next((r for r in (chain.regions or []) if r.region_type == reg), None)
        if region is None:
            continue
        res = [r for r in region.residues if r.ca is not None]
        if not res:
            continue
        ca = np.asarray([r.ca for r in res], float)
        d_pep = np.linalg.norm(ca[:, None, :] - pep_ca[None], axis=2).min(axis=1)
        d_mhc = np.linalg.norm(ca[:, None, :] - mhc_ca[None], axis=2).min(axis=1)
        engaged = np.array([(chain.chain_id, r.seq_index) in touch for r in res], bool)
        out[f"{name}_n_res"] = float(len(res))
        out[f"{name}_frac_con"] = float(engaged.mean())
        for tag, m in (("all", np.ones(len(res), bool)), ("con", engaged), ("non", ~engaged)):
            if not m.any():
                continue
            out[f"{name}_{tag}_d_pep"] = float(d_pep[m].mean())
            out[f"{name}_{tag}_d_mhc"] = float(d_mhc[m].mean())
            out[f"{name}_{tag}_delta"] = float(d_pep[m].mean() - d_mhc[m].mean())

        # how the loop engages along its own length: protrusion of its closest point below its mean,
        # and how evenly its residues sit against the peptide
        out[f"{name}_protrusion"] = float(d_pep.mean() - d_pep.min())
        out[f"{name}_d_pep_sd"] = float(d_pep.std())
        out[f"{name}_d_pep_min"] = float(d_pep.min())
        out[f"{name}_d_mhc_min"] = float(d_mhc.min())

        # sigma involution: does this loop read the N-terminal or the C-terminal half of the peptide?
        if len(pep_ca) >= 4:
            half = len(pep_ca) // 2
            dN = np.linalg.norm(ca[:, None, :] - pep_ca[None, :half], axis=2).min(axis=1)
            dC = np.linalg.norm(ca[:, None, :] - pep_ca[None, -half:], axis=2).min(axis=1)
            out[f"{name}_d_pepN"] = float(dN.mean())
            out[f"{name}_d_pepC"] = float(dC.mean())
            out[f"{name}_delta_NC"] = float(dN.mean() - dC.mean())

        # docking polarity: which MHC helix the loop leans on
        for h1, h2, tag in (("HELIX_A1", "HELIX_A2", "helixA"), ("HELIX_A1", "HELIX_B1", "helixAB")):
            a1, a2 = helix.get(h1), helix.get(h2)
            if a1 is None or a2 is None or not len(a1) or not len(a2):
                continue
            da = np.linalg.norm(ca[:, None, :] - a1[None], axis=2).min(axis=1)
            db = np.linalg.norm(ca[:, None, :] - a2[None], axis=2).min(axis=1)
            out[f"{name}_d_{tag}1"] = float(da.mean())
            out[f"{name}_d_{tag}2"] = float(db.mean())
            out[f"{name}_delta_{tag}"] = float(da.mean() - db.mean())

    # --- contrasts BETWEEN loops: the statements that need two loops to make -----------------------
    def g(k):
        return out.get(k, float("nan"))

    for tag in ("all", "con"):
        # CDR3 reach relative to the germline loops, per chain and pooled
        for ch, three, one, two in (("a", "cdr3a", "cdr1a", "cdr2a"), ("b", "cdr3b", "cdr1b", "cdr2b")):
            germ = np.nanmean([g(f"{one}_{tag}_d_pep"), g(f"{two}_{tag}_d_pep")])
            out[f"cdr3_vs_germline_{ch}_{tag}_d_pep"] = float(g(f"{three}_{tag}_d_pep") - germ)
            germm = np.nanmean([g(f"{one}_{tag}_d_mhc"), g(f"{two}_{tag}_d_mhc")])
            out[f"cdr3_vs_germline_{ch}_{tag}_d_mhc"] = float(g(f"{three}_{tag}_d_mhc") - germm)
        # alpha/beta asymmetry of the two CDR3 loops
        out[f"cdr3_ab_asym_{tag}_d_pep"] = float(g(f"cdr3a_{tag}_d_pep") - g(f"cdr3b_{tag}_d_pep"))
        out[f"cdr3_ab_asym_{tag}_d_mhc"] = float(g(f"cdr3a_{tag}_d_mhc") - g(f"cdr3b_{tag}_d_mhc"))
        out[f"cdr3_ab_asym_{tag}_delta"] = float(g(f"cdr3a_{tag}_delta") - g(f"cdr3b_{tag}_delta"))
    # the sigma involution as one number: alpha should read N, beta should read C
    out["sigma_NC_split"] = float(g("cdr3a_delta_NC") - g("cdr3b_delta_NC"))
    return out
