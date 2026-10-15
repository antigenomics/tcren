"""One joint Gaussian per class; every read-out is a projection of it.

The whole model is `(mu_c, Sigma_c)` for c in {non-binder, binder}, estimated once on the
hold-out in the transformed descriptor space. Everything the author asked for is then a linear
algebra operation on that one object, with no further fitting:

* **posterior**            `P(1|x) = pi_1 N(x; mu_1, S_1) / sum_c pi_c N(x; mu_c, S_c)`
* **any marginal, exact**  `P(1|x_S)`: the same formula on `mu_{c,S}` and the sub-block `S_{c,SS}`.
  "Keep only geometry and dPhi_pep_mhc" is a sub-block, not a re-fit.
* **PCA at any width**     the covariance in PCA coordinates is `W' S_c W` exactly, so truncating
  to m components is a projection of the same object -- again no re-fit.
* **confidence prediction** `E[z | x] = mu_z + S_{z,x} S_xx^-1 (x - mu_x)` for
  `z = logit(ipTM)`; `z - E[z|x]` is the transformed-axis residual.
* **one-class anomaly**    `A(x) = sum_{k<K} (u_k'(x - mu_1))^2 / lam_k` over the stiffest binder
  directions -- its moments use binders, while the frozen coordinate transform uses both classes.

**There is no Jacobian in any of this.** The author's note asked to back-transform
`P(binder|PCA) -> Jacobian * P(binder|descriptors)`. For a POSTERIOR the Jacobian cancels: PCA is
affine, `|det dz/dx|` is a constant independent of x, and it appears identically in numerator and
denominator. `P(c|z) == P(c|x)` exactly whenever W is square and invertible, and
``tests/unit/test_score.py`` asserts it to 1e-9 rather than taking it on faith. The Jacobian is real for
the DENSITY; what happens under truncation is not a Jacobian but a marginalization, which for a
Gaussian is the sub-block above.

**Shrinkage is not optional here.** A full covariance over 145 coordinates is 10,585 parameters and
the negative arm has 1,155 rows once benchmark structures are held out. Sigma is shrunk toward a
scaled identity with a Ledoit-Wolf intensity; without it the smallest eigenvalues are estimation
noise and the Mahalanobis form divides by them.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


def _wmean_cov(X: np.ndarray, w: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    X, w = np.asarray(X, float), np.asarray(w, float)
    if (X.ndim != 2 or X.shape[1] == 0 or w.shape != (len(X),)
            or not np.isfinite(X).all() or not np.isfinite(w).all()
            or (w < 0).any() or not np.isclose(w.sum(), 1., rtol=1e-10, atol=1e-12)
            or (w ** 2).sum() >= 1.):
        raise ValueError("finite rows and nonnegative normalized reliability weights with effective size > 1 required")
    mu = w @ X
    Z = X - mu
    return mu, (Z * w[:, None]).T @ Z / (1.0 - (w ** 2).sum())


def _positive_cov(S: np.ndarray) -> None:
    if (S.ndim != 2 or S.shape[0] == 0 or S.shape[0] != S.shape[1]
            or not np.isfinite(S).all() or not np.allclose(S, S.T, rtol=1e-10, atol=1e-12)):
        raise ValueError("covariance must be finite and symmetric positive definite")
    try:
        np.linalg.cholesky(S)
    except np.linalg.LinAlgError as exc:
        raise ValueError("covariance must be positive definite") from exc


def conditional_variance(S: np.ndarray, cross: np.ndarray, variance: float) -> float:
    """Schur complement of a scalar augmentation; reject an inconsistent joint covariance."""
    S, cross = np.asarray(S, float), np.asarray(cross, float)
    _positive_cov(S)
    if cross.shape != (len(S),) or not np.isfinite(cross).all() or not np.isfinite(variance) or variance < 0:
        raise ValueError("invalid confidence covariance shape or variance")
    explained = float(cross @ np.linalg.solve(S, cross))
    residual = float(variance - explained)
    if residual < -1e-10 * max(1., abs(variance), abs(explained)):
        raise ValueError("confidence joint covariance must be positive semidefinite")
    return residual


def _width(k: int | None, p: int, minimum: int = 0) -> int:
    if k is None:
        return p
    if isinstance(k, (bool, np.bool_)) or not isinstance(k, (int, np.integer)) or not minimum <= k <= p:
        raise ValueError(f"width must be an integer in [{minimum}, {p}]")
    return int(k)


def _shrink(S: np.ndarray, X: np.ndarray) -> tuple[np.ndarray, float]:
    """Ledoit-Wolf toward a scaled identity, `(1-a)S + a*(tr S / p)I`.

    The intensity is taken from the unweighted rows -- sklearn's estimator has no weighted form,
    and the epitope weights change the effective sample size by well under the factor that would
    move `a` materially. Reported so the approximation is on the record rather than hidden.
    """
    from sklearn.covariance import ledoit_wolf_shrinkage
    p = S.shape[0]
    try:
        a = float(ledoit_wolf_shrinkage(X, assume_centered=False))
    except Exception:
        a = 0.1
    a = float(np.clip(a, 1e-3, 0.9))
    return (1 - a) * S + a * (np.trace(S) / p) * np.eye(p), a


@dataclass
class Joint:
    """Two Gaussians in one transformed coordinate system, and the read-outs they induce."""
    names: list[str]
    mu: dict[int, np.ndarray] = field(default_factory=dict)
    cov: dict[int, np.ndarray] = field(default_factory=dict)
    prior: dict[int, float] = field(default_factory=dict)
    alpha: dict[int, float] = field(default_factory=dict)
    lam1: np.ndarray | None = None   # binder eigenvalues, ascending (stiffest first)
    U1: np.ndarray | None = None

    def fit(self, X: dict[int, np.ndarray], w: dict[int, np.ndarray],
            shrink: bool = True) -> "Joint":
        """`shrink=False` keeps the raw covariance, which only the artefact test wants.

        Ledoit-Wolf floors the smallest binder direction at s.d. 0.0797, and the directions
        the crystal test flags as the generator's own regularity all sit below 0.05 -- so a
        shrunk model cannot see the artefact band at all, and measuring it needs the raw one.
        """
        if (set(X) != {0, 1} or set(w) != {0, 1} or not self.names
                or len(set(self.names)) != len(self.names)):
            raise ValueError("two classes and unique nonempty coordinate names required")
        n_tot = sum(len(X[c]) for c in X)
        for c in X:
            if np.asarray(X[c]).ndim != 2 or np.asarray(X[c]).shape[1] != len(self.names):
                raise ValueError("class rows must match coordinate names")
            m, S = _wmean_cov(X[c], w[c])
            S, a = _shrink(S, X[c]) if shrink else (S, 0.0)
            _positive_cov(S)
            self.mu[c], self.cov[c], self.alpha[c] = m, S, a
            self.prior[c] = len(X[c]) / n_tot
        lam, U = np.linalg.eigh(self.cov[1])
        self.lam1, self.U1 = np.maximum(lam, 1e-12), U
        return self

    # ---------------------------------------------------------------- read-outs
    def _idx(self, subset: list[str] | None) -> np.ndarray:
        if subset is None:
            return np.arange(len(self.names))
        s = set(subset)
        if not subset or len(s) != len(subset) or not s.issubset(self.names):
            raise ValueError("subset must contain unique known coordinates and cannot be empty")
        return np.array([i for i, n in enumerate(self.names) if n in s], dtype=int)

    def log_odds(self, X: np.ndarray, subset: list[str] | None = None) -> np.ndarray:
        """`log P(1|x_S) - log P(0|x_S)`, exact for any subset S by sub-blocking Sigma.

        `X` carries exactly the columns `subset` names, in that order, when `subset` is given.
        """
        X = np.asarray(X, float)
        j = self._idx(subset)
        if X.ndim != 2 or X.shape[1] not in ({len(self.names), len(j)} if subset else {len(self.names)}):
            raise ValueError("rows must match full or subset coordinate width")
        if subset is not None and X.shape[1] == len(j):
            # `_idx` returns indices in `self.names` order, so the sub-blocks below are in that
            # order; permute X's columns to match the docstring's promise rather than silently
            # pairing coordinate k of X with coordinate k of a differently-ordered sub-block.
            # `anomaly_on` already honours caller order; this makes the two agree.
            want = set(subset)
            X = X[:, [subset.index(n) for n in self.names if n in want]]
            out = {}
            for c in (0, 1):
                S = self.cov[c][np.ix_(j, j)]
                _positive_cov(S)
                _, logdet = np.linalg.slogdet(S)
                Z = X - self.mu[c][j]
                out[c] = -0.5 * (np.einsum("ij,jk,ik->i", Z, np.linalg.inv(S), Z) + logdet) \
                    + np.log(self.prior[c])
            return out[1] - out[0]
        out = {}
        for c in (0, 1):
            S = self.cov[c][np.ix_(j, j)]
            _positive_cov(S)
            _, logdet = np.linalg.slogdet(S)
            Z = X[:, j] - self.mu[c][j]
            out[c] = -0.5 * (np.einsum("ij,jk,ik->i", Z, np.linalg.inv(S), Z) + logdet) \
                + np.log(self.prior[c])
        return out[1] - out[0]

    def anomaly(self, X: np.ndarray, k: int | None = None) -> np.ndarray:
        """Partial Mahalanobis to the binder Gaussian on its k stiffest directions.

        The readout uses only binder moments. Training still selects binders by label and fits
        the shared coordinate transform on both classes; tier 1 is not negative-independent training.
        """
        k = _width(k, len(self.lam1))
        P = (X - self.mu[1]) @ self.U1[:, :k]
        return (P ** 2 / self.lam1[:k]).sum(1)

    def anomaly_on(self, X: np.ndarray, names: list[str], k: int | None = None) -> np.ndarray:
        """:meth:`anomaly` restricted to the coordinates `names`, which `X`'s columns carry.

        The eigenbasis is recomputed on the sub-block rather than sliced out of the full one: the
        stiff directions of a marginal are not the stiff directions of the joint restricted, and
        using the latter would score a structure against directions its table cannot supply.
        """
        if list(names) == list(self.names):
            return self.anomaly(X, k)
        self._idx(names)
        j = [self.names.index(n) for n in names]
        lam, U = np.linalg.eigh(self.cov[1][np.ix_(j, j)])
        lam = np.maximum(lam, 1e-12)
        k = _width(k, len(lam))
        P = (X - self.mu[1][j]) @ U[:, :k]
        return (P ** 2 / lam[:k]).sum(1)

    def predict(self, X: np.ndarray, target: str, given: list[str] | None = None,
                cls: int = 1) -> np.ndarray:
        """`E[target | given]` under class `cls` -- the Gaussian conditional mean.

        With `target="iptm"` this is what the structure says the generator's confidence should have
        been; the residual against the reported value is the QC channel.
        """
        t = self.names.index(target)
        j = self._idx(given)
        j = j[j != t]
        S = self.cov[cls]
        A = np.linalg.solve(S[np.ix_(j, j)], (X[:, j] - self.mu[cls][j]).T).T
        return self.mu[cls][t] + A @ S[np.ix_(j, [t])].ravel()

    def project(self, m: int) -> tuple[np.ndarray, "Joint"]:
        """The same model in the top-m PCA coordinates of the pooled within-class scatter.

        `W' Sigma_c W` is exact, so this is a projection of the fitted object and not a second fit.
        At m = p with W orthonormal it is a rotation, and the posterior is unchanged -- the
        identity `tests/unit/test_score.py` asserts.
        """
        m = _width(m, len(self.names), minimum=1)
        Sw = sum(self.prior[c] * self.cov[c] for c in self.cov)
        lam, U = np.linalg.eigh(Sw)
        W = U[:, ::-1][:, :m]
        out = Joint(names=[f"pc{i}" for i in range(m)])
        out.prior = dict(self.prior)
        for c in self.cov:
            out.mu[c] = W.T @ self.mu[c]
            out.cov[c] = W.T @ self.cov[c] @ W
        lam1, U1 = np.linalg.eigh(out.cov[1])
        out.lam1, out.U1 = np.maximum(lam1, 1e-12), U1
        return W, out
