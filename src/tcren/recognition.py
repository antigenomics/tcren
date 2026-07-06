"""Gaussian Bayesian-network classifier: real vs shuffled TCR-pMHC complexes.

A conditional-linear-Gaussian Bayes net. A DAG is learned (BIC hill-climbing) over the standardized interface
features on the *within-class-centred* data, so the edges capture genuine feature-feature dependence rather
than the class shift. The binary class ``y`` (real = 1 / shuffled = 0) and the MHC class are then added as
discrete parents of **every** feature node, shifting its conditional mean. Classification is the Gaussian
log-likelihood ratio ``log P(x | y=1) - log P(x | y=0)`` (plus the class-prior log-odds if not balanced).

Pure numpy (dep-light). Trained parameters serialise to gzipped JSON (:meth:`GaussianBNClassifier.save` /
:meth:`load`); :meth:`to_dot` renders the network with graphviz. Trained on the Shuffled2026 decoys from
:mod:`tcren.shuffle`.

This module also provides :class:`BayesianLogisticRecognizer` — a frozen distribution-aware Bayesian logistic
regression (fit externally with PyMC): each feature enters via its family's canonical form
(:func:`encode_features` — circular angles as cos/sin, bounded ratios as logit, counts/continuous linearly),
so unlike the Gaussian BN it does not mis-specify the count and angle features.
"""
from __future__ import annotations

import gzip
import json
import math
from pathlib import Path

import numpy as np

_EPS = 1e-9


def _bic_local(Z: np.ndarray, j: int, parents: list[int]) -> float:
    n = Z.shape[0]
    X = np.column_stack([np.ones(n), Z[:, parents]] if parents else [np.ones(n)])
    beta, *_ = np.linalg.lstsq(X, Z[:, j], rcond=None)
    rss = float(np.sum((Z[:, j] - X @ beta) ** 2)) or _EPS
    k = len(parents) + 2
    return -0.5 * n * math.log(rss / n) - 0.5 * k * math.log(n)


def _acyclic(parents: dict[int, set[int]], p: int) -> bool:
    indeg = {j: len(parents[j]) for j in range(p)}
    children = {j: [c for c in range(p) if j in parents[c]] for j in range(p)}
    q = [j for j in range(p) if indeg[j] == 0]
    seen = 0
    while q:
        u = q.pop()
        seen += 1
        for c in children[u]:
            indeg[c] -= 1
            if indeg[c] == 0:
                q.append(c)
    return seen == p


def _hill_climb(Z: np.ndarray, max_parents: int = 3) -> dict[int, list[int]]:
    """BIC-scored greedy structure search over the columns of ``Z`` (add/remove edges)."""
    p = Z.shape[1]
    parents = {j: set() for j in range(p)}
    score = {j: _bic_local(Z, j, []) for j in range(p)}
    improved = True
    while improved:
        improved = False
        best = (1e-6, None)
        for a in range(p):
            for b in range(p):
                if a == b:
                    continue
                if a in parents[b]:
                    cand = parents[b] - {a}
                elif len(parents[b]) < max_parents:
                    cand = parents[b] | {a}
                    trial = {k: set(v) for k, v in parents.items()}
                    trial[b] = cand
                    if not _acyclic(trial, p):
                        continue
                else:
                    continue
                delta = _bic_local(Z, b, sorted(cand)) - score[b]
                if delta > best[0]:
                    best = (delta, (b, cand))
        if best[1]:
            b, cand = best[1]
            parents[b] = cand
            score[b] = _bic_local(Z, b, sorted(cand))
            improved = True
    return {j: sorted(parents[j]) for j in range(p)}


class GaussianBNClassifier:
    """Conditional-linear-Gaussian BN classifier (see the module docstring)."""

    def __init__(self, feature_names: list[str], max_parents: int = 3):
        self.feature_names = list(feature_names)
        self.max_parents = max_parents

    # -- fit ---------------------------------------------------------------------------------------------
    def _standardize(self, X: np.ndarray) -> np.ndarray:
        X = np.asarray(X, float)
        X = np.where(np.isfinite(X), X, np.take(self.mu_, np.arange(X.shape[1]))[None, :])
        return (X - self.mu_) / self.sd_

    def fit(self, X, y, mhc_class=None) -> "GaussianBNClassifier":
        X = np.asarray(X, float)
        y = np.asarray(y, int)
        m = np.zeros(len(y)) if mhc_class is None else np.asarray(mhc_class, float)
        self.mu_ = np.nanmean(np.where(np.isfinite(X), X, np.nan), axis=0)
        self.sd_ = np.nanstd(np.where(np.isfinite(X), X, np.nan), axis=0) + _EPS
        Z = self._standardize(X)
        # structure on within-(y,m)-class-centred data: remove the class/covariate shift first
        Zc = Z.copy()
        for yv in (0, 1):
            for mv in np.unique(m):
                mask = (y == yv) & (m == mv)
                if mask.sum() > 1:
                    Zc[mask] -= Zc[mask].mean(axis=0)
        self.structure_ = _hill_climb(Zc, self.max_parents)
        # per-node conditional: x_j ~ N(b0 + b.parents + g.y + d.m, sigma^2)
        self.nodes_ = {}
        n = len(y)
        for j in range(len(self.feature_names)):
            pa = self.structure_[j]
            cov = np.column_stack([np.ones(n)] + ([Z[:, pa]] if pa else []) + [y.astype(float), m])
            beta, *_ = np.linalg.lstsq(cov, Z[:, j], rcond=None)
            resid = Z[:, j] - cov @ beta
            self.nodes_[j] = {"parents": pa, "beta": beta.tolist(),
                              "sigma": float(np.sqrt(np.mean(resid ** 2)) + _EPS)}
        self.prior_ = float(np.mean(y))
        return self

    # -- predict -----------------------------------------------------------------------------------------
    def _loglik(self, Z: np.ndarray, m: np.ndarray, yval: int) -> np.ndarray:
        n = Z.shape[0]
        ll = np.zeros(n)
        for j, nd in self.nodes_.items():
            pa = nd["parents"]
            beta = np.asarray(nd["beta"])
            sig = nd["sigma"]
            cov = np.column_stack([np.ones(n)] + ([Z[:, pa]] if pa else [])
                                  + [np.full(n, yval, float), m])
            mean = cov @ beta
            ll += -0.5 * math.log(2 * math.pi * sig ** 2) - 0.5 * ((Z[:, j] - mean) / sig) ** 2
        return ll

    def decision_function(self, X, mhc_class=None) -> np.ndarray:
        """Log-likelihood ratio ``log P(x|y=1) - log P(x|y=0)`` (balanced; add prior log-odds separately)."""
        Z = self._standardize(np.asarray(X, float))
        m = np.zeros(len(Z)) if mhc_class is None else np.asarray(mhc_class, float)
        return self._loglik(Z, m, 1) - self._loglik(Z, m, 0)

    def predict_proba(self, X, mhc_class=None, balanced: bool = True) -> np.ndarray:
        s = self.decision_function(X, mhc_class)
        if not balanced:
            s = s + math.log(self.prior_ / (1 - self.prior_ + _EPS))
        p = 1.0 / (1.0 + np.exp(-np.clip(s, -700, 700)))
        return np.column_stack([1 - p, p])

    # -- marginalization ---------------------------------------------------------------------------------
    def _joint_gaussian(self):
        """Reconstruct the class-conditional joint Gaussian from the DAG: shared covariance + the y/m means.

        Each standardized node is ``z_j = b0_j + sum b.parents + g_j y + d_j m + eps_j``. In matrix form
        ``z = (I-B)^{-1}(c + eps)`` for fixed (y, m), so the (homoscedastic) covariance is
        ``Sigma = (I-B)^{-1} diag(sigma^2) (I-B)^{-T}`` and the class-mean shift is ``(I-B)^{-1} g``.
        """
        p = len(self.feature_names)
        B = np.zeros((p, p)); c0 = np.zeros(p); g = np.zeros(p); d = np.zeros(p); D = np.zeros(p)
        for j, nd in self.nodes_.items():
            beta = np.asarray(nd["beta"]); pa = nd["parents"]
            c0[j] = beta[0]
            for k, a in enumerate(pa):
                B[j, a] = beta[1 + k]
            g[j] = beta[-2]; d[j] = beta[-1]; D[j] = nd["sigma"] ** 2
        IB = np.linalg.inv(np.eye(p) - B)
        return IB, c0, g, d, IB @ np.diag(D) @ IB.T

    def marginal_decision(self, X, keep, mhc_class=None) -> np.ndarray:
        """LLR ``log P(x_G|y=1) - log P(x_G|y=0)`` after **marginalizing out** every feature not in ``keep``.

        ``keep`` is a list of feature names (e.g. the geometry features, energy marginalised out). Because the
        covariance is shared across classes the marginal LLR is linear in the kept features.
        """
        idx = [self.feature_names.index(n) for n in keep]
        Z = self._standardize(np.asarray(X, float))
        m = np.zeros(len(Z)) if mhc_class is None else np.asarray(mhc_class, float)
        IB, c0, g, d, Sigma = self._joint_gaussian()
        Sinv = np.linalg.inv(Sigma[np.ix_(idx, idx)])
        shift = (IB @ g)[idx]                              # mu_1 - mu_0 on the kept block (m-independent)
        base = (IB @ (c0 + 0.5 * g))                       # the m=0 midpoint; add d*m per sample below
        dm = (IB @ d)[idx]
        Zk = Z[:, idx]
        mid = base[idx][None, :] + np.outer(m, dm)         # per-sample class midpoint on kept block
        return ((Zk - mid) @ Sinv) @ shift

    def marginal_proba(self, X, keep, mhc_class=None) -> np.ndarray:
        s = self.marginal_decision(X, keep, mhc_class)
        p = 1.0 / (1.0 + np.exp(-np.clip(s, -700, 700)))
        return np.column_stack([1 - p, p])

    # -- persistence + rendering -------------------------------------------------------------------------
    def to_dict(self) -> dict:
        return {"feature_names": self.feature_names, "max_parents": self.max_parents,
                "mu": self.mu_.tolist(), "sd": self.sd_.tolist(), "prior": self.prior_,
                "structure": {str(k): v for k, v in self.structure_.items()},
                "nodes": {str(k): v for k, v in self.nodes_.items()}}

    @classmethod
    def from_dict(cls, d: dict) -> "GaussianBNClassifier":
        obj = cls(d["feature_names"], d["max_parents"])
        obj.mu_ = np.asarray(d["mu"]); obj.sd_ = np.asarray(d["sd"]); obj.prior_ = d["prior"]
        obj.structure_ = {int(k): v for k, v in d["structure"].items()}
        obj.nodes_ = {int(k): v for k, v in d["nodes"].items()}
        return obj

    def save(self, path: str | Path) -> Path:
        path = Path(path)
        data = json.dumps(self.to_dict()).encode()
        (gzip.open(path, "wb") if str(path).endswith(".gz") else open(path, "wb")).write(data)
        return path

    @classmethod
    def load(cls, path: str | Path) -> "GaussianBNClassifier":
        path = Path(path)
        raw = (gzip.open(path, "rb") if str(path).endswith(".gz") else open(path, "rb")).read()
        return cls.from_dict(json.loads(raw))

    def to_dot(self, coef_threshold: float = 0.15) -> str:
        """Graphviz DAG: feature-feature edges (partial slopes) + class/MHC covariate edges above threshold."""
        names = self.feature_names
        lines = ["digraph BN {", '  rankdir=LR; node [shape=box, style=rounded, fontsize=9];',
                 '  y [shape=ellipse, style=filled, fillcolor="#ffd9d9", label="class (real/shuffled)"];',
                 '  mhc [shape=ellipse, style=filled, fillcolor="#d9e6ff", label="MHC class"];']
        for j, nm in enumerate(names):
            lines.append(f'  f{j} [label="{nm}"];')
        for j, nd in self.nodes_.items():
            beta = nd["beta"]
            pa = nd["parents"]
            for k, a in enumerate(pa):
                lines.append(f'  f{a} -> f{j} [label="{beta[k+1]:+.2f}", fontsize=7];')
            g, d = beta[-2], beta[-1]                      # y and mhc covariate slopes
            if abs(g) >= coef_threshold:
                lines.append(f'  y -> f{j} [color="#cc3333", label="{g:+.2f}", fontsize=7];')
            if abs(d) >= coef_threshold:
                lines.append(f'  mhc -> f{j} [color="#3355cc", label="{d:+.2f}", fontsize=7];')
        lines.append("}")
        return "\n".join(lines)


# ======================================================================================================
# Distribution-aware Bayesian logistic recognizer
# ======================================================================================================
# Encodings that respect each feature's natural distribution before the (linear) logistic predictor:
#   circular angle (von Mises) -> cos/sin ; bounded ratio (Beta) -> logit ; exact duplicate -> dropped.
# Counts (Poisson canonical) and continuous / unit-vector features already enter linearly, so a logistic
# regression -- unlike the Gaussian BN above -- does not mis-specify them.
_ENCODE = {"dock_torsion": "cos_sin", "chain_balance": "logit_half", "n_hbond": "drop"}
_HALF = 0.5


def encode_features(X, feature_names) -> tuple[np.ndarray, list[str]]:
    """Distribution-aware design matrix (pre-standardization).

    ``dock_torsion`` (circular, wraps) -> its von Mises sufficient statistics ``(cos, sin)``; ``chain_balance``
    ([0, 0.5] Beta) -> ``logit(2x)``; ``n_hbond`` dropped (exact duplicate of ``ct_tp_hydrogen_bond``);
    everything else (counts + continuous + unit-vector cos/sin components) enters linearly.

    Args:
        X: ``(n, len(feature_names))`` raw feature array.
        feature_names: column names of ``X``.

    Returns:
        ``(Z, encoded_names)`` — the encoded matrix and its column names.
    """
    X = np.asarray(X, float)
    idx = {n: i for i, n in enumerate(feature_names)}
    cols, names = [], []
    for n in feature_names:
        enc = _ENCODE.get(n, "linear")
        if enc == "drop":
            continue
        x = X[:, idx[n]]
        if enc == "cos_sin":
            cols += [np.cos(x), np.sin(x)]; names += [f"{n}_cos", f"{n}_sin"]
        elif enc == "logit_half":
            u = np.clip(x / _HALF, 1e-4, 1 - 1e-4)
            cols.append(np.log(u / (1 - u))); names.append(f"{n}_logit")
        else:
            cols.append(x); names.append(n)
    return np.column_stack(cols), names


class BayesianLogisticRecognizer:
    """Frozen distribution-aware Bayesian logistic (posterior-mean coefficients) — dep-light numpy predictor.

    Applies :func:`encode_features`, standardizes with the stored training statistics (nan -> train mean), and
    returns ``sigmoid(alpha + Z @ beta)``. Fit externally by PyMC (``appendix/logistic_stan/build.py``) and
    frozen here; serialises to gzipped JSON.
    """

    def __init__(self, feature_names, encoded_names, mean, sd, alpha, beta, prior: str = "normal"):
        self.feature_names = list(feature_names)
        self.encoded_names = list(encoded_names)
        self.mean = np.asarray(mean, float)
        self.sd = np.asarray(sd, float)
        self.alpha = float(alpha)
        self.beta = np.asarray(beta, float)
        self.prior = prior

    def _design(self, X) -> np.ndarray:
        Z, names = encode_features(X, self.feature_names)
        if names != self.encoded_names:
            raise ValueError("encoded feature names do not match the fitted model")
        Z = np.where(np.isfinite(Z), Z, self.mean[None, :])       # nan -> train mean
        return (Z - self.mean) / self.sd

    def decision_function(self, X) -> np.ndarray:
        return self.alpha + self._design(X) @ self.beta

    def predict_proba(self, X) -> np.ndarray:
        p = 1.0 / (1.0 + np.exp(-np.clip(self.decision_function(X), -700, 700)))
        return np.column_stack([1 - p, p])

    def to_dict(self) -> dict:
        return {"feature_names": self.feature_names, "encoded_names": self.encoded_names,
                "mean": self.mean.tolist(), "sd": self.sd.tolist(),
                "alpha": self.alpha, "beta": self.beta.tolist(), "prior": self.prior}

    @classmethod
    def from_dict(cls, d: dict) -> "BayesianLogisticRecognizer":
        return cls(d["feature_names"], d["encoded_names"], d["mean"], d["sd"], d["alpha"], d["beta"],
                   d.get("prior", "normal"))

    def save(self, path: str | Path) -> Path:
        path = Path(path)
        data = json.dumps(self.to_dict()).encode()
        (gzip.open(path, "wb") if str(path).endswith(".gz") else open(path, "wb")).write(data)
        return path

    @classmethod
    def load(cls, path: str | Path) -> "BayesianLogisticRecognizer":
        path = Path(path)
        raw = (gzip.open(path, "rb") if str(path).endswith(".gz") else open(path, "rb")).read()
        return cls.from_dict(json.loads(raw))


# ===================================================================================================
# Structure -> the 35-descriptor recognition vector the frozen recognizers consume, and P(real).
#
# Reproduces the extractor the shipped models were trained on (the manuscript's compute_features.py):
# docking geometry + per-interface TCRen/MJ energetics (F, poly-Ala ΔF) + contact-type tallies +
# biopython ΔSASA burial + MHC-class indicator. Heavy imports are function-local so that a bare
# ``import tcren`` (and ``import tcren.recognition``) stays dependency-light.
# ===================================================================================================

#: Feature names, in the order the frozen logistic recognizer's design matrix expects. The Gaussian BN
#: uses the same list minus ``mhc_class_bin`` (which it carries as a discrete node instead).
RECOGNITION_FEATURES = (
    "extent", "e_tcr_mhc", "chain_balance", "pitch", "crossing", "dock_d", "dock_torsion",
    "dock_tcr_uy", "dock_tcr_uz", "dock_mhc_uy", "dock_mhc_uz", "e_cdr12", "e_cdr3a", "e_cdr3b",
    "F_tcr_pep", "F_tcr_mhc", "F_pep_mhc", "dF_tcr_pep", "dF_pep_mhc", "n_contacts_tp",
    "n_pep_contacted", "n_contacts_tm", "ct_tp_salt_bridge", "ct_tm_salt_bridge",
    "ct_tp_hydrogen_bond", "ct_tm_hydrogen_bond", "ct_tp_aromatic", "ct_tm_aromatic",
    "ct_tp_hydrophobic", "ct_tm_hydrophobic", "ct_tp_other", "ct_tm_other", "n_hbond",
    "burial", "mhc_class_bin",
)
_CT_TYPES = ("salt_bridge", "hydrogen_bond", "aromatic", "hydrophobic", "other")
_TCR_TYPES = ("TRA", "TRB", "TRG", "TRD")


def _extent(cm) -> float:
    """Distinct TCR residues contacting the pMHC (interface size); default TCR-region selection."""
    import polars as pl
    df = pl.concat([cm.interface("tcr_peptide"), cm.interface("tcr_mhc")])
    nodes = set()
    if df.height:
        for a, i in zip(df["chain.id.from"].to_list(), df["residue.index.from"].to_list()):
            nodes.add((a, i))
    return float(len(nodes))


def _chain_balance(cm) -> float:
    """min(a,b)/(a+b) over TCR:peptide contacts by TCR chain (0.5 = both chains equal, 0 = one only)."""
    tp = cm.interface("tcr_peptide", tcr_regions="all")
    if tp.height == 0:
        return math.nan
    a = b = 0
    for t in tp["chain.type.from"].to_list():
        a += t == "TRA"
        b += t == "TRB"
    return min(a, b) / (a + b) if (a + b) else math.nan


def _burial(structure, tcr_ids, pmhc_ids) -> float:
    """Interface ΔSASA = SASA(TCR alone) + SASA(pMHC alone) − SASA(complex) via biopython ShrakeRupley
    (``n_points=100``), reproducing the training-time ``burial``. ΔSASA is an interface quantity, so the
    distal TCR constant domain cancels; computed on a temp PDB of the typed chains."""
    if not tcr_ids or not pmhc_ids:
        return math.nan
    import os
    import tempfile
    from copy import deepcopy

    from Bio.PDB import PDBParser
    from Bio.PDB.Model import Model
    from Bio.PDB.SASA import ShrakeRupley
    from Bio.PDB.Structure import Structure as BioStructure

    from .structure.io import write_pdb

    with tempfile.TemporaryDirectory() as td:
        path = os.path.join(td, "complex.pdb")
        write_pdb(structure, path)
        model = PDBParser(QUIET=True).get_structure("x", path)[0]      # parsed fully into memory
    sr = ShrakeRupley(n_points=100)

    def sasa_of(ids):
        m2 = Model(0)
        for ch in model:
            if ch.id in ids:
                m2.add(deepcopy(ch))
        s2 = BioStructure("t")
        s2.add(m2)
        sr.compute(s2, level="A")
        return sum(a.sasa for ch in m2 for res in ch if res.id[0] == " " for a in res.get_atoms())

    both = set(tcr_ids) | set(pmhc_ids)
    return float((sasa_of(set(tcr_ids)) + sasa_of(set(pmhc_ids))) - sasa_of(both))


def recognition_features(source, *, organism: str = "human", potential=None) -> dict[str, float]:
    """Extract the 35-descriptor recognition vector from a TCR–pMHC structure (path or parsed).

    Reproduces the feature set the shipped real-vs-shuffled recognizers were trained on
    (:data:`RECOGNITION_FEATURES`): docking geometry, per-interface TCRen/MJ energies (raw ``F`` and
    poly-alanine ``ΔF``), contact-type tallies, interface ΔSASA ``burial``, and the ``mhc_class_bin``
    indicator. The structure is chain-typed and MHC-annotated in place. Returns a dict keyed by
    :data:`RECOGNITION_FEATURES` (degenerate/undefined terms are ``NaN``).

    Feed the result to :func:`frozen_recognizers` (or :class:`BayesianLogisticRecognizer`) for
    ``P(real)`` — the probability the complex looks like a genuine TCR–pMHC recognition interface.
    """
    import polars as pl

    from .annotation import classify_chains
    from .contact_types import contact_type_counts
    from .contactmap import ContactMap
    from .ddg import reference_delta
    from .mhc import annotate_mhc
    from .oracle import _native_peptide
    from .orient.docking import docking_angles
    from .orient.tcrdock_geometry import docking_geometry
    from .pipeline import _interface_energy
    from .potential import mj as _mj
    from .potential import tcren as _tcren
    from .structure import Structure, import_structure

    s = source if isinstance(source, Structure) else import_structure(source)
    if all(c.chain_type is None for c in s.chains):
        classify_chains(s, organism=organism, autodetect_species=True)
    calls = annotate_mhc(s)
    tcren_pot = potential or _tcren()
    mj_pot = _mj()

    cm = ContactMap.from_structure(s)
    native = _native_peptide(s)
    row = {k: math.nan for k in RECOGNITION_FEATURES}

    try:                                                              # geometry (docking)
        da = docking_angles(s)
        row["pitch"], row["crossing"] = float(da.incident_angle), float(da.crossing_angle)
    except Exception:
        pass
    try:
        dg = docking_geometry(s)                                     # native TCRdock rigid-body params
        row.update(dock_d=float(dg.d), dock_torsion=float(dg.torsion),
                   dock_tcr_uy=float(dg.tcr_unit_y), dock_tcr_uz=float(dg.tcr_unit_z),
                   dock_mhc_uy=float(dg.mhc_unit_y), dock_mhc_uz=float(dg.mhc_unit_z))
    except Exception:
        pass

    tm = cm.interface("tcr_mhc", tcr_regions="all")                  # interface energetics
    row["F_tcr_mhc"] = row["e_tcr_mhc"] = float(_interface_energy(tm, mj_pot))
    tp = cm.interface("tcr_peptide", tcr_regions="all")
    reg, ch = pl.col("region.type.from"), pl.col("chain.type.from")
    row["F_tcr_pep"] = float(_interface_energy(tp, tcren_pot))
    row["F_pep_mhc"] = float(_interface_energy(cm.interface("peptide_mhc"), mj_pot))
    row["e_cdr12"] = float(_interface_energy(tp.filter(reg.is_in(["CDR1", "CDR2"])), tcren_pot))
    row["e_cdr3a"] = float(_interface_energy(tp.filter((reg == "CDR3") & (ch == "TRA")), tcren_pot))
    row["e_cdr3b"] = float(_interface_energy(tp.filter((reg == "CDR3") & (ch == "TRB")), tcren_pot))
    if native:
        try:
            row["dF_tcr_pep"] = float(reference_delta(cm, native, tcren_pot, interface="tcr_peptide"))
        except Exception:
            pass
        try:
            row["dF_pep_mhc"] = float(reference_delta(cm, native, mj_pot, interface="peptide_mhc"))
        except Exception:
            pass

    row["extent"] = _extent(cm)                                      # coverage
    row["chain_balance"] = _chain_balance(cm)
    row["n_contacts_tp"] = float(tp.height)
    row["n_pep_contacted"] = float(tp.select("residue.index.to").unique().height if tp.height else 0)
    row["n_contacts_tm"] = float(tm.height)

    ctp = contact_type_counts(cm, "tcr_peptide")                     # contact types
    ctm = contact_type_counts(cm, "tcr_mhc")
    for t in _CT_TYPES:
        row[f"ct_tp_{t}"] = float(ctp[f"pairs_{t}"])
        row[f"ct_tm_{t}"] = float(ctm[f"pairs_{t}"])
    row["n_hbond"] = float(ctp["pairs_hydrogen_bond"])

    tcr_ids = [c.chain_id for c in s.chains if c.chain_type in _TCR_TYPES]
    pmhc_ids = [c.chain_id for c in s.chains if c.chain_type is not None and c.chain_type not in _TCR_TYPES]
    row["burial"] = _burial(s, tcr_ids, pmhc_ids)
    row["mhc_class_bin"] = 1.0 if any(getattr(cl, "mhc_class", None) == "MHCII" for cl in calls) else 0.0
    return row


def frozen_recognizers():
    """Load the shipped real-vs-shuffled recognizers ``(logistic, bn)`` from ``tcren.data``.

    ``logistic`` is the headline distribution-aware :class:`BayesianLogisticRecognizer`
    (``shuffle_logistic.json.gz``); ``bn`` is the :class:`GaussianBNClassifier`
    (``shuffle_bn.json.gz``). Feed rows from :func:`recognition_features` to :func:`real_probability`.
    """
    from importlib import resources
    d = resources.files("tcren.data")
    lr = BayesianLogisticRecognizer.load(str(d.joinpath("shuffle_logistic.json.gz")))
    bn = GaussianBNClassifier.load(str(d.joinpath("shuffle_bn.json.gz")))
    return lr, bn


def real_probability(rows, *, recognizers=None) -> dict[str, np.ndarray]:
    """``P(real)`` for feature rows from :func:`recognition_features`.

    ``rows`` is a dict or a list of dicts keyed by :data:`RECOGNITION_FEATURES`. Returns
    ``{"logistic": p, "bn": p}`` — the headline logistic recognizer and the Gaussian BN, each an array
    of P(genuine TCR–pMHC interface). NaN features are imputed to the training mean by each model.
    """
    if isinstance(rows, dict):
        rows = [rows]
    lr, bn = recognizers or frozen_recognizers()
    Xlr = np.array([[r.get(k, np.nan) for k in lr.feature_names] for r in rows], float)
    Xbn = np.array([[r.get(k, np.nan) for k in bn.feature_names] for r in rows], float)
    m = np.array([r.get("mhc_class_bin", 0.0) for r in rows], float)
    return {"logistic": lr.predict_proba(Xlr)[:, 1], "bn": bn.predict_proba(Xbn, m)[:, 1]}
