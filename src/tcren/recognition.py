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
from functools import lru_cache
from pathlib import Path

import numpy as np

_EPS = 1e-9


def _sigmoid(x: np.ndarray) -> np.ndarray:
    """Numerically-safe logistic (clips the exponent to avoid overflow warnings)."""
    return 1.0 / (1.0 + np.exp(-np.clip(x, -700, 700)))


def _dump_json_gz(d: dict, path: str | Path) -> Path:
    """Serialise ``d`` to JSON at ``path`` (gzip-compressed iff the name ends in ``.gz``)."""
    path = Path(path)
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "wb") as fh:
        fh.write(json.dumps(d).encode())
    return path


def _load_json_gz(path: str | Path) -> dict:
    """Inverse of :func:`_dump_json_gz`."""
    path = Path(path)
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rb") as fh:
        return json.loads(fh.read())


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
        p = _sigmoid(s)
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
        p = _sigmoid(s)
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
        return _dump_json_gz(self.to_dict(), path)

    @classmethod
    def load(cls, path: str | Path) -> "GaussianBNClassifier":
        return cls.from_dict(_load_json_gz(path))

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
        p = _sigmoid(self.decision_function(X))
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
        return _dump_json_gz(self.to_dict(), path)

    @classmethod
    def load(cls, path: str | Path) -> "BayesianLogisticRecognizer":
        return cls.from_dict(_load_json_gz(path))


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

#: Interface-symmetry descriptors from per-loop TCR:peptide contact **counts** (not energies), emitted as
#: extra ``recognize`` output columns — **not** part of :data:`RECOGNITION_FEATURES` (the frozen models'
#: 35-vector is fixed). ``cdr3_dominance`` = CDR3(α+β) share of CDR contacts (higher = CDR3-dominated,
#: oriented positive); ``cdr3_ab_imbalance`` = ``|CDR3α−CDR3β|`` normalised (absolute); ``chain_cdr_imbalance``
#: = ``|α−β|`` whole-CDR normalised (absolute). See :func:`_interface_symmetry`.
INTERFACE_SYMMETRY_FEATURES = ("cdr3_dominance", "cdr3_ab_imbalance", "chain_cdr_imbalance")

#: CDR3-local frame features (18), the FramePose layer the whole-TCR :data:`RECOGNITION_FEATURES` miss.
#: Per loop, relative to the pMHC groove frame (u, w, n; origin = peptide Cα centroid):
#: ``reach`` = |loop centroid − origin|; ``o{u,w,n}`` = unit(centroid−origin)·(u,w,n) (where over the
#: groove the loop sits); ``a{u,w,n}`` = unit(Cα_N→Cα_C)·(u,w,n) (loop orientation over the groove);
#: ``topep`` = min Cα-Cα distance loop→peptide (engagement depth); ``ext`` = |Cα_C − Cα_N| (extension).
_CDR3_FRAME_KEYS = ("reach", "ou", "ow", "on", "au", "aw", "an", "topep", "ext")
CDR3_FRAME_FEATURES = tuple(f"{loop}_{k}" for loop in ("cdr3a", "cdr3b") for k in _CDR3_FRAME_KEYS)

#: Matrix-swap features (12): the same TCR:peptide contacts scored under TCRen vs the generic MJ
#: potential, per interface group. ``tcren_{g}``/``mj_{g}`` are the two energies and ``d_{g}`` their
#: difference (the recognition-specific component; generic packing cancels). ``g`` ∈ {tp, cdr12, cdr3a,
#: cdr3b}. Note ``tcren_tp``/``tcren_cdr12``/``tcren_cdr3a``/``tcren_cdr3b`` duplicate the core
#: ``F_tcr_pep``/``e_cdr12``/``e_cdr3a``/``e_cdr3b`` by construction (kept for full parity).
_MATRIX_SWAP_GROUPS = ("tp", "cdr12", "cdr3a", "cdr3b")
MATRIX_SWAP_FEATURES = tuple(f"{pre}_{g}" for g in _MATRIX_SWAP_GROUPS for pre in ("tcren", "mj", "d"))

#: The full feature vector: the 35 core recognition descriptors + the 18 CDR3-frame + 12 matrix-swap.
FULL_FEATURES = RECOGNITION_FEATURES + CDR3_FRAME_FEATURES + MATRIX_SWAP_FEATURES

#: Frozen "forced-pose" classifier: P(this pose is an AF-forced interface rather than a crystal-natural
#: one). A raw-feature logistic (no standardization) over interface *strain* — stretched CDR3 loops and
#: thin contacts. Trained ONLY on provenance (Canonical2026 crystals = 0 vs AF/TCRmodel2 models = 1;
#: n=2681, 268 crystal / 2413 forced), so it is independent of any binder label; 5-fold CV AUC 0.762.
#: High ``p_forced`` marks a "too-good-to-be-true" pose; the score grades crystal < AF-real < AF-decoy.
#:
#: .. note::
#:    For new work prefer the fit-free :func:`tcren.cohort.strain_z` (``S_strain``). It grades the
#:    same crystal < AF-real < AF-decoy provenance gradient by signed standardization of the strain
#:    terms, with no training set — so it is fully reproducible, unlike the coefficients below.
#:
#: .. warning::
#:    These coefficients are **frozen and not re-derivable** -- the n=2681 training set no longer
#:    exists. ``models/fit_frozen.py::forced_pose`` in the benchmark repo recovers the *procedure*
#:    (unstandardized L2 logistic, C=0.1, which reproduces the 0.762 CV above to within 0.001) but
#:    not the coefficients. Refitting on the surviving 1168-row fixture gives a **better** in-sample
#:    ROC (0.769 vs 0.745), which is how we know these were fit on different rows rather than
#:    overfit to what survives. Do not replace them with a refit without re-basing the benchmarks.
FORCED_POSE_MODEL = {
    "features": ("dock_d", "cdr3b_reach", "cdr3b_topep", "cdr3a_ext", "extent_per_ct", "chain_balance"),
    "coef": (-0.46517433874162056, 0.14437146872011086, -0.31411562068257676,
             -2.114810136001524, 1.198769596894963, -0.6237422800760706),
    "intercept": 26.11747560652168,
    "cv_auc": 0.762,
}


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


def _interface_symmetry(tp) -> dict[str, float]:
    """CDR3-dominance and TCR chain/loop imbalance from per-loop TCR:peptide contact **counts**.

    ``tp`` is the ``tcr_peptide`` interface table (``tcr_regions="all"``). Unlike ``e_cdr*`` (which are
    interface *energies*), these are pure contact-topology descriptors. Emitted as extra output columns
    (:data:`INTERFACE_SYMMETRY_FEATURES`), not part of :data:`RECOGNITION_FEATURES`.
    """
    import polars as pl
    reg, ch = pl.col("region.type.from"), pl.col("chain.type.from")
    h = lambda f: float(tp.filter(f).height)  # noqa: E731
    n12 = h(reg.is_in(["CDR1", "CDR2"]))                                   # germline CDR1/2 (both chains)
    n3a, n3b = h((reg == "CDR3") & (ch == "TRA")), h((reg == "CDR3") & (ch == "TRB"))
    nA = h(reg.is_in(["CDR1", "CDR2", "CDR3"]) & (ch == "TRA"))            # whole alpha CDRs
    nB = h(reg.is_in(["CDR1", "CDR2", "CDR3"]) & (ch == "TRB"))            # whole beta CDRs
    tot = n12 + n3a + n3b
    return {
        # CDR3 (a+b) share of CDR TCR:peptide contacts -- higher = CDR3-dominated (binder-like; oriented +)
        "cdr3_dominance": (n3a + n3b) / tot if tot else math.nan,
        # |CDR3a - CDR3b| normalised imbalance -- absolute magnitude (direction is tested, not assumed)
        "cdr3_ab_imbalance": abs(n3a - n3b) / (n3a + n3b) if (n3a + n3b) else math.nan,
        # |alpha - beta| whole-CDR contact imbalance, normalised -- absolute magnitude
        "chain_cdr_imbalance": abs(nA - nB) / (nA + nB) if (nA + nB) else math.nan,
    }


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


def _cdr3_frame_features(structure) -> dict[str, float]:
    """The 18 CDR3-local frame descriptors (:data:`CDR3_FRAME_FEATURES`) for a chain-typed structure.

    Both CDR3 loops are projected onto the pMHC groove frame (see :data:`CDR3_FRAME_FEATURES`). The
    structure must already be chain-typed (``classify_chains``) so its CDR3 regions are populated.
    Undefined terms (no groove frame, missing peptide, or a loop with < 3 Cα) are ``NaN``.
    """
    from .orient.docking import _chain_ca, _groove_frame

    out = {k: math.nan for k in CDR3_FRAME_FEATURES}
    try:
        u, w, n = _groove_frame(structure)
    except Exception:
        return out
    pep = _chain_ca(structure, ("PEPTIDE",))
    if len(pep) < 2:
        return out
    origin = pep.mean(axis=0)
    basis = np.stack([u, w, n])                                        # rows = groove basis
    for loop, ctype in (("cdr3a", "TRA"), ("cdr3b", "TRB")):
        cas = None
        for c in structure.chains:
            if c.chain_type != ctype:
                continue
            for reg in getattr(c, "regions", []) or []:
                if reg.region_type == "CDR3":
                    pts = [r.ca for r in reg.residues if r.ca is not None]
                    if len(pts) >= 3:
                        cas = np.asarray(pts)
                    break
        if cas is None:
            continue
        d = cas.mean(axis=0) - origin
        reach = float(np.linalg.norm(d))
        off = basis @ (d / (reach + 1e-9))
        av = cas[-1] - cas[0]
        ax = basis @ (av / (np.linalg.norm(av) + 1e-9))
        topep = float(np.linalg.norm(cas[:, None, :] - pep[None, :, :], axis=2).min())
        ext = float(np.linalg.norm(cas[-1] - cas[0]))
        for k, v in zip(_CDR3_FRAME_KEYS, (reach, *off, *ax, topep, ext)):
            out[f"{loop}_{k}"] = float(v)
    return out


def _matrix_swap_features(cm, tcren_pot, mj_pot) -> dict[str, float]:
    """The 12 matrix-swap descriptors (:data:`MATRIX_SWAP_FEATURES`) from a contact map.

    Scores the TCR:peptide contacts (whole interface + the CDR1/2, CDR3α, CDR3β groups) under both the
    TCRen and the generic MJ potential; the per-group difference ``d`` isolates the recognition-specific
    component (generic packing cancels since both potentials read the identical contacts).
    """
    import polars as pl

    from .pipeline import _interface_energy

    tp = cm.interface("tcr_peptide", tcr_regions="all")
    reg, ch = pl.col("region.type.from"), pl.col("chain.type.from")
    groups = {
        "tp": tp,
        "cdr12": tp.filter(reg.is_in(["CDR1", "CDR2"])),
        "cdr3a": tp.filter((reg == "CDR3") & (ch == "TRA")),
        "cdr3b": tp.filter((reg == "CDR3") & (ch == "TRB")),
    }
    out: dict[str, float] = {}
    for name, df in groups.items():
        et = float(_interface_energy(df, tcren_pot))
        em = float(_interface_energy(df, mj_pot))
        out[f"tcren_{name}"], out[f"mj_{name}"], out[f"d_{name}"] = et, em, et - em
    return out


def recognition_features(source, *, organism: str = "human", potential=None,
                         full: bool = False, annotate: bool = True) -> dict[str, float]:
    """Extract the 35-descriptor recognition vector from a TCR–pMHC structure (path or parsed).

    Reproduces the feature set the shipped real-vs-shuffled recognizers were trained on
    (:data:`RECOGNITION_FEATURES`): docking geometry, per-interface TCRen/MJ energies (raw ``F`` and
    poly-alanine ``ΔF``), contact-type tallies, interface ΔSASA ``burial``, and the ``mhc_class_bin``
    indicator. The structure is chain-typed and MHC-annotated in place. Returns a dict keyed by
    :data:`RECOGNITION_FEATURES` (degenerate/undefined terms are ``NaN``).

    Feed the result to :func:`frozen_recognizers` (or :class:`BayesianLogisticRecognizer`) for
    ``P(real)`` — the probability the complex looks like a genuine TCR–pMHC recognition interface.

    With ``full=True`` the row is extended with the 18 CDR3-frame (:data:`CDR3_FRAME_FEATURES`) and 12
    matrix-swap (:data:`MATRIX_SWAP_FEATURES`) descriptors — the complete :data:`FULL_FEATURES` vector.
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
    if annotate:                                                      # skip if pre-annotated (batch path)
        if all(c.chain_type is None for c in s.chains):
            classify_chains(s, organism=organism, autodetect_species=True)
        annotate_mhc(s)
    tcren_pot = potential or _tcren()
    mj_pot = _mj()

    cm = ContactMap.from_structure(s)
    native = _native_peptide(s)
    row = {k: math.nan for k in (FULL_FEATURES if full else RECOGNITION_FEATURES)}

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
    row["mhc_class_bin"] = 1.0 if any(getattr(c, "chain_supertype", None) == "MHCII"
                                      for c in s.chains) else 0.0

    if full:                                                          # FramePose CDR3 layer + matrix-swap
        row.update(_cdr3_frame_features(s))
        row.update(_matrix_swap_features(cm, tcren_pot, mj_pot))
    return row


def _stability_clash_columns(s) -> dict[str, float]:
    """Interface steric-clash + TCR:peptide contact-stability descriptors for the recognize table.

    Extra *output* columns, **not** part of :data:`RECOGNITION_FEATURES` or any fitted model: a
    coordinate-only read of forced-pose quality --- steric-clash burden (:func:`tcren.interface_clashes`)
    and contact fragility (:func:`tcren.contact_stability`). NaN where the structure lacks a peptide or
    receptor chain.
    """
    from .clashes import interface_clashes
    from .stability import contact_stability

    out: dict[str, float] = {}
    try:
        cl = interface_clashes(s)
        out["n_clashes"], out["clash_score"] = float(cl.n_clashes), float(cl.clash_score)
    except Exception:  # noqa: BLE001 - no peptide chain etc.
        out["n_clashes"] = out["clash_score"] = math.nan
    try:
        st = contact_stability(s)
        out["exp_lost"] = float(st.exp_lost)
        out["mean_margin"] = float(st.mean_margin)
        out["frac_robust"] = float(st.frac_robust)
    except Exception:  # noqa: BLE001 - no peptide/receptor chain etc.
        out["exp_lost"] = out["mean_margin"] = out["frac_robust"] = math.nan
    return out


def _symmetry_columns(s) -> dict[str, float]:
    """Interface-symmetry extra output columns (:data:`INTERFACE_SYMMETRY_FEATURES`) for the recognize
    table --- CDR3-dominance and α/β contact imbalance from a fresh contact map. NaN on failure."""
    from .contactmap import ContactMap
    try:
        cm = ContactMap.from_structure(s)
        return _interface_symmetry(cm.interface("tcr_peptide", tcr_regions="all"))
    except Exception:  # noqa: BLE001 - no peptide/receptor chain etc.
        return {k: math.nan for k in INTERFACE_SYMMETRY_FEATURES}


def recognition_table(items, *, organism: str = "human", full: bool = False, scores: bool = False,
                      with_p_real: bool = True, threads: int = 1, chunk: int = 64,
                      _cohort_scores: bool = True) -> list[dict]:
    """Batched feature (+score) extraction for a whole set of TCR–pMHC structures.

    ``items`` is an iterable of ``(id, structure-or-path)``. The set is annotated with a **single**
    arda call per organism (:func:`tcren.paper.helpers._batch_annotate`) and a **single** mmseqs MHC
    search (:func:`tcren.mhc.annotate_mhc_batch`) — the dataset-scale path that avoids the per-structure
    annotation cost — then :func:`recognition_features` (``full=``) is extracted for each. With
    ``with_p_real`` the ``p_real`` / ``p_real_bn`` recognizer columns are added; with ``scores`` the
    fit-free cohort scores ``q_bind`` / ``s_strain`` (**recommended**, see :mod:`tcren.cohort`) plus
    the fitted ``p_forced`` / ``p_bind`` (retained for reproducibility). Returns one row dict per
    structure (``complex.id`` + features [+ scores]); a structure that fails yields
    ``{"complex.id": id, "error": ...}`` so the batch stays resilient.

    ``threads`` > 1 splits the set into ``chunk``-sized batches and annotates/featurises them
    concurrently. Batching alone is not enough at cohort scale: arda's own search is single-threaded
    and its cost grows with the batch, so one 600-structure call is slower than ten 60-structure
    calls in parallel. mmseqs runs in a subprocess and releases the GIL, so plain threads suffice.
    """
    items = list(items)
    if threads > 1 and len(items) > chunk:
        from concurrent.futures import ThreadPoolExecutor
        batches = [items[i:i + chunk] for i in range(0, len(items), chunk)]
        with ThreadPoolExecutor(max_workers=threads) as ex:
            parts = list(ex.map(
                lambda b: recognition_table(b, organism=organism, full=full, scores=scores,
                                            with_p_real=with_p_real, threads=1,
                                            _cohort_scores=False), batches))
        rows = [r for p in parts for r in p]
        if scores:      # q_bind / s_strain are cohort-relative: computed once over the WHOLE set
            _add_cohort_scores(rows)
        return rows
    from .annotation import classify_chains
    from .annotation.arda_adapter import _import_arda
    from .mhc import annotate_mhc_batch
    from .paper.helpers import _batch_annotate
    from .structure import Structure, import_structure

    ids, structs, rows = [], [], []
    for id_, src in items:
        try:
            structs.append(src if isinstance(src, Structure) else import_structure(src))
            ids.append(id_)
        except Exception as exc:  # noqa: BLE001
            rows.append({"complex.id": id_, "error": f"{type(exc).__name__}: {str(exc)[:80]}"})

    if structs:                                                       # one arda + one mmseqs for the set
        recs = _batch_annotate(structs, _import_arda(), organisms=(organism, "mouse"))
        for i, s in enumerate(structs):
            try:
                classify_chains(s, organism=organism, autodetect_species=True,
                                precomputed_records=recs[i])
            except Exception:  # noqa: BLE001 - MHC-only / unannotatable chains stay unset
                pass
        annotate_mhc_batch(structs)

    recognizers = frozen_recognizers() if with_p_real else None
    if scores:
        from .binder import binder_features, binder_score

    for id_, s in zip(ids, structs):
        try:
            feats = recognition_features(s, organism=organism, full=full, annotate=False)
            row = {"complex.id": id_, **feats, **_stability_clash_columns(s), **_symmetry_columns(s)}
            if with_p_real:
                p = real_probability(feats, recognizers=recognizers)
                row["p_real"], row["p_real_bn"] = float(p["logistic"][0]), float(p["bn"][0])
            if scores:
                row["p_forced"] = forced_pose_score(feats)
                try:
                    row["p_bind"] = float(binder_score(binder_features(s)))
                except Exception:  # noqa: BLE001 - binder ext optional
                    row["p_bind"] = math.nan
            rows.append(row)
        except Exception as exc:  # noqa: BLE001
            rows.append({"complex.id": id_, "error": f"{type(exc).__name__}: {str(exc)[:80]}"})

    if scores and _cohort_scores:                                     # fit-free cohort scores (recommended)
        _add_cohort_scores(rows)
    return rows


def _add_cohort_scores(rows: list[dict]) -> None:
    """Append the fit-free cohort scores ``q_bind`` (:func:`tcren.cohort.q_score`) and ``s_strain``
    (:func:`tcren.cohort.strain_z`) in place. Cohort-relative, so they are computed over the whole
    batch at once and are the **recommended** binder / forced-pose scores (see :mod:`tcren.cohort`).
    Needs the ``full`` CDR3-frame features; NaN where a structure lacks them.
    """
    from . import cohort
    ok = [r for r in rows if "error" not in r]
    if len(ok) < 2:                                                   # cohort scores are undefined for <2
        for r in ok:
            r["q_bind"] = r["s_strain"] = math.nan
        return
    table = {k: [r.get(k, math.nan) for r in ok]
             for k in set().union(*(r.keys() for r in ok)) if k != "complex.id"}
    try:
        q, s = cohort.q_score(table), cohort.strain_z(table)
    except KeyError:                                                  # missing full features -> skip cleanly
        return
    for i, r in enumerate(ok):
        r["q_bind"], r["s_strain"] = float(q[i]), float(s[i])


@lru_cache(maxsize=None)
def frozen_recognizers():
    """Load the shipped real-vs-shuffled recognizers ``(logistic, bn)`` from ``tcren.data`` (cached).

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


def forced_pose_score(feats: dict[str, float]) -> float:
    """``P(forced)`` — probability a pose is an AF-forced interface, from :data:`FORCED_POSE_MODEL`.

    ``feats`` is a row from :func:`recognition_features` with ``full=True`` (it needs the CDR3-frame
    ``cdr3b_reach``/``cdr3b_topep``/``cdr3a_ext`` plus core ``dock_d``/``extent``/``n_contacts_tp``/
    ``chain_balance``). ``extent_per_ct`` is derived as ``extent / n_contacts_tp``. Returns ``NaN`` if
    any required feature is missing/undefined. High = "too good to be true" (see :data:`FORCED_POSE_MODEL`).
    """
    m = FORCED_POSE_MODEL
    nc = feats.get("n_contacts_tp", math.nan)
    derived = {"extent_per_ct": feats.get("extent", math.nan) / nc if nc else math.nan}
    z = m["intercept"]
    for name, w in zip(m["features"], m["coef"]):
        v = derived.get(name, feats.get(name, math.nan))
        if not (isinstance(v, (int, float)) and math.isfinite(v)):
            return math.nan
        z += w * float(v)
    return 1.0 / (1.0 + math.exp(-max(-700.0, min(700.0, z))))


def kit_score(p_bind, iptm) -> np.ndarray:
    """Synergistic AF × tcren binder score: ``z(p_bind) + z(iptm)`` over the scored cohort.

    Combines the intrinsic tcren binder score (:func:`tcren.binder.binder_score`, from
    ``recognize --scores``) with the AlphaFold/TCRmodel2 **ipTM** that ships free with every model. On the
    TCRvdb raw-label benchmark this fixed no-fit combination beats **either alone** at precision
    (macro-PR 0.847 vs ipTM 0.782 / p_bind 0.804; precision 0.969 at 10% recall vs ipTM 0.861; Δ macro-PR
    vs ipTM +0.065, 95% CI [+0.022, +0.100], P(Δ>0)=1.00). A leave-epitope-out logistic on the same two
    inputs gives the more conservative +0.041 [+0.005, +0.076] — a different estimator, not this score.
    Higher = more binder-like.

    Cohort-relative: ``z`` standardizes over the input arrays, so pass the **whole set** of AF models you
    are ranking (not one structure). NaNs are ignored by the mean/sd and propagate to their own entries.

    Args:
        p_bind: tcren binder scores for the cohort (``recognize --scores`` ``p_bind`` column).
        iptm: the matching AlphaFold ipTM values.

    Returns:
        The combined ranking score, one per structure.
    """
    p_bind = np.asarray(p_bind, float)
    iptm = np.asarray(iptm, float)

    def _z(a):
        return (a - np.nanmean(a)) / (np.nanstd(a) + 1e-9)

    return _z(p_bind) + _z(iptm)
