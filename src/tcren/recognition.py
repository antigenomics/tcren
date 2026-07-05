"""Gaussian Bayesian-network classifier: real vs shuffled TCR-pMHC complexes.

A conditional-linear-Gaussian Bayes net. A DAG is learned (BIC hill-climbing) over the standardized interface
features on the *within-class-centred* data, so the edges capture genuine feature-feature dependence rather
than the class shift. The binary class ``y`` (real = 1 / shuffled = 0) and the MHC class are then added as
discrete parents of **every** feature node, shifting its conditional mean. Classification is the Gaussian
log-likelihood ratio ``log P(x | y=1) - log P(x | y=0)`` (plus the class-prior log-odds if not balanced).

Pure numpy (dep-light). Trained parameters serialise to gzipped JSON (:meth:`GaussianBNClassifier.save` /
:meth:`load`); :meth:`to_dot` renders the network with graphviz. Trained on the Shuffled2026 decoys from
:mod:`tcren.shuffle`.
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
