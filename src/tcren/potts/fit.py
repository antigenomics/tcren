"""Fitting the coupled contact-map model by penalised pseudolikelihood.

The conditional ``P(sigma_a = 1 | sigma_{-a})`` is logistic in
``eta_a + sum_k K_k n_k(a)`` with ``n_k(a)`` the count of contacting neighbours in coupling class
``k``, so the coupled fit is an ordinary weighted-binomial GLM with a handful of extra integer
covariates, and it stays concave. That is Besag's pseudolikelihood; consistency for this model
class is Ravikumar, Wainwright & Lafferty (arXiv:1010.0311), and plmDCA (arXiv:1211.1281) is the
same recipe on Potts sequence variables. **No partition function is needed to fit** — only to score
(:func:`tcren.potts.ais_log_z`).

The design is one-hot per block plus ``J`` plus the neighbour counts, over-parametrised and
identified by an ℓ2 ridge, then projected to the zero-sum gauge. Penalise **then** project: an ℓ2
penalty silently picks its own gauge, which plmDCA flags and which is why the projection is a
separate step rather than a constraint.
"""

from __future__ import annotations

import numpy as np
import polars as pl
import scipy.sparse as sp

from .kernel import edges, neighbour_counts
from .model import PottsModel, kernel_names
from .sites import site_codes

#: ℓ2 ridge on every coefficient but the intercept — about one pseudo-observation of information.
DEFAULT_RIDGE = 1.0


def design(codes, sizes, counts: np.ndarray, fixed_coupling: np.ndarray | None = None):
    """Sparse row-level design ``[intercept | one-hot blocks | J | neighbour counts]``.

    ``fixed_coupling`` replaces the free 400-cell ``J`` block with a **single** column carrying the
    named potential's value at that cell, so its coefficient is the scale ``beta_matrix``.

    Returns:
        ``(X, block_slices, k0)`` with ``k0`` the first coupling-coefficient column.
    """
    rows = len(codes[0])
    ar = np.arange(rows)
    ri, ci, vi = [ar], [np.zeros(rows, np.int64)], [np.ones(rows)]
    slices, off = [], 1
    for ck, nk in zip(codes, sizes):
        ri.append(ar); ci.append(off + ck); vi.append(np.ones(rows))
        slices.append(slice(off, off + nk)); off += nk
    if fixed_coupling is None:
        ri.append(ar); ci.append(off + codes[0] * sizes[1] + codes[1]); vi.append(np.ones(rows))
        n_j = sizes[0] * sizes[1]
        slices.append(slice(off, off + n_j)); off += n_j
    else:
        ri.append(ar); ci.append(np.full(rows, off))
        vi.append(fixed_coupling[codes[0], codes[1]])
        slices.append(slice(off, off + 1)); off += 1
    k0 = off
    for k in range(counts.shape[1]):
        nz = np.nonzero(counts[:, k])[0]
        ri.append(nz); ci.append(np.full(len(nz), off)); vi.append(counts[nz, k])
        off += 1
    X = sp.csr_matrix((np.concatenate(vi), (np.concatenate(ri), np.concatenate(ci))),
                      shape=(rows, off))
    return X, slices, k0


def irls(X, y: np.ndarray, w: np.ndarray, ridge: float, *, free=(0,),
         tol: float = 1e-10, maxit: int = 200):
    """Penalised IRLS for a weighted Bernoulli GLM on a sparse design. Convex; Newton steps.

    Returns ``(coefficients, penalised Hessian)``. ``free`` names columns exempt from the ridge.
    """
    p = X.shape[1]
    lam = np.full(p, float(ridge))
    lam[list(free)] = 0.0
    b = np.zeros(p)
    for _ in range(maxit):
        mu = 1.0 / (1.0 + np.exp(-(X @ b)))
        H = (X.T @ sp.diags(w * mu * (1 - mu)) @ X).toarray() + np.diag(lam)
        step = np.linalg.solve(H, X.T @ (w * (y - mu)) - lam * b)
        b = b + step
        if np.max(np.abs(step)) < tol:
            break
    return b, H


def gauge(b: np.ndarray, slices, free_coupling: bool):
    """Project to the zero-sum (Ising) gauge, in place on a copy.

    Every one-body block is centred and ``J`` double-centred, the displaced level absorbed by the
    intercept, so ``eta`` is unchanged while as much as possible sits in the fields and as little
    as necessary in the couplings (Cocco et al., arXiv:1703.01222). A fixed-matrix coupling is a
    scale on an already double-centred matrix and is never re-centred.
    """
    b = b.copy()
    J = None
    if free_coupling:
        n = int(round(len(b[slices[-1]]) ** 0.5))
        J = b[slices[-1]].reshape(n, n)
        u, v, c = J.mean(1), J.mean(0), J.mean()
        J = J - u[:, None] - v[None, :] + c
        b[0] -= c
        b[slices[0]] += u
        b[slices[1]] += v
        b[slices[-1]] = J.ravel()
    for s in slices[:-1]:
        m = b[s].mean()
        b[s] -= m
        b[0] += m
    return b, J


def cluster_se(X, resid: np.ndarray, gid: np.ndarray, n_groups: int, H) -> np.ndarray:
    """Sandwich standard errors clustered on the structure.

    A structure contributes hundreds of correlated sites, so the model-based ``(X'WX)^-1``
    understates every standard error.
    """
    G = sp.csr_matrix((np.ones(len(gid)), (gid, np.arange(len(gid)))), shape=(n_groups, len(gid)))
    S = (G @ X.multiply(resid[:, None])).toarray()
    Hi = np.linalg.inv(H)
    return np.sqrt(np.clip(np.diag(Hi @ (S.T @ S) @ Hi), 0, None))


def fit_potts(sites: pl.DataFrame, *, radius: float = 15.0, cutoff: float = 5.0,
              couplings: bool = True, coupling_matrix: str | None = None,
              weights: dict[str, float] | None = None, ridge: float = DEFAULT_RIDGE,
              joint: bool | None = None, notes: str = "") -> PottsModel:
    """Fit the model to a table of available pairs.

    Args:
        sites: Rows from :func:`tcren.potts.available_pairs`, one or both partners concatenated.
        radius: Availability radius, Å. Sites beyond it are dropped.
        cutoff: Contact definition, Å — recorded on the model; ``sigma`` is already computed.
        couplings: Fit the ``sigma``–``sigma`` kernel. ``False`` gives the factorised model, whose
            partition function is then exact.
        coupling_matrix: Fix ``J`` to one scale on this bundled potential (``tcren2``, ``mj``, …)
            instead of fitting 400 free cells. Competing matrices fitted this way carry identical
            parameter counts, so their pseudo-log-likelihoods compare directly.
        weights: Per-structure weights, e.g. from
            :func:`tcren.potential.balanced_weights`, to down-weight redundancy. Missing ids get 1.
        ridge: ℓ2 penalty on every coefficient but the intercept.
        joint: Include the cross-class coupling family. Defaults to whether both partner classes
            are present.
        notes: Free text stored on the model.

    Returns:
        A :class:`PottsModel` in the zero-sum gauge.

    Example:
        >>> from tcren.potts import available_pairs, fit_potts    # doctest: +SKIP
        >>> pairs = pl.concat([available_pairs(s) for s in structures])   # doctest: +SKIP
        >>> model = fit_potts(pairs, notes="Native2026")          # doctest: +SKIP
        >>> model.to_json("potts.json")                           # doctest: +SKIP
    """
    from .model import centred_potential

    sites = sites.filter(pl.col("d_ca") <= radius).sort("pdb.id", maintain_order=True)
    if sites.is_empty():
        raise ValueError("no available pairs inside the radius")
    codes, sizes, q = site_codes(sites, radius=radius)
    if joint is None:
        joint = q["cls"].n_unique() > 1
    sigma = q["sigma"].to_numpy()
    w = (np.ones(q.height) if weights is None
         else q["pdb.id"].replace_strict(weights, default=1.0).to_numpy().astype(float))

    edge_lists = edges(q, joint=joint) if couplings else []
    counts = (neighbour_counts(edge_lists, sigma, q.height) if edge_lists
              else np.zeros((q.height, 0)))
    fixed = centred_potential(coupling_matrix) if coupling_matrix else None
    X, slices, k0 = design(codes, sizes, counts, fixed)
    b, H = irls(X, sigma, w, ridge)
    bg, J = gauge(b, slices, fixed is None)

    mu = 1.0 / (1.0 + np.exp(-(X @ b)))
    upid, gid = np.unique(q["pdb.id"].to_numpy(), return_inverse=True)
    se = cluster_se(X, w * (sigma - mu), gid, len(upid), H)
    pll = float(np.sum(w * (sigma * np.log(np.clip(mu, 1e-300, 1.0))
                            + (1 - sigma) * np.log(np.clip(1 - mu, 1e-300, 1.0)))))

    return PottsModel(
        alpha=float(bg[0]),
        h_rec=bg[slices[0]].tolist(), h_par=bg[slices[1]].tolist(),
        g_dist=bg[slices[2]].tolist(), g_region=bg[slices[3]].tolist(),
        g_role=bg[slices[4]].tolist(), g_class=bg[slices[5]].tolist(),
        kernel=bg[k0:].tolist(), kernel_se=se[k0:].tolist(),
        coupling=(J.tolist() if J is not None else None),
        beta_matrix=(None if fixed is None else float(bg[slices[-1]][0])),
        coupling_matrix_name=coupling_matrix,
        radius=radius, cutoff=cutoff, joint=bool(joint),
        n_structures=int(len(upid)), n_sites=int(q.height), n_contacts=int(sigma.sum()),
        pseudo_loglik=pll, notes=notes)


def kernel_table(model: PottsModel) -> pl.DataFrame:
    """The coupling coefficients with their cluster-robust standard errors, strongest first.

    A coefficient is the **log-odds added to a site's field per contacting neighbour** of that
    class. On TCR:peptide crystals every axial class comes out positive and every off-axis class
    negative: a made contact recruits its own sequence neighbours onto the same partner residue and
    suppresses the diagonal one.
    """
    names = kernel_names(model.joint)
    se = list(model.kernel_se) or [float("nan")] * len(model.kernel)
    return pl.DataFrame({"class": names[:len(model.kernel)],
                         "K": model.kernel, "se": se[:len(model.kernel)]}).with_columns(
        (pl.col("K") / pl.col("se")).alias("z")
    ).sort(pl.col("z").abs().fill_nan(-1.0), descending=True)   # a class with no edges sorts last
