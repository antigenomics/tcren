"""Bounded evidence and vectorized average precision helpers."""

import numpy as np


def align30(position, length):
    """IMGT CDR3 (anchors excluded), padding after its third residue."""
    p, n = np.asarray(position), np.asarray(length)
    if np.any((n > 30) | (n < 1) | (p < 1) | (p > n)):
        raise ValueError("CDR3 cannot be represented in 30 positions")
    return p + np.where(p > 3, 30 - n, 0)



def bounded(x, scale, kind):
    """Nonnegative extensive evidence → bounded monotone dimensionless score."""
    x = np.asarray(x, dtype=float)
    if scale <= 0 or not np.isfinite(x).all() or np.any(x < 0):
        raise ValueError("Expected finite nonnegative evidence and positive scale")
    if kind == "rational":
        return x / (scale + x)
    if kind == "exp":
        return -np.expm1(-x / scale)
    if kind == "hill2":
        return (x / scale) ** 2 / (1 + (x / scale) ** 2)
    raise ValueError(kind)



def ap_columns(y, scores):
    """Exact noninterpolated AP for columns, including tied score blocks."""
    y = np.asarray(y)
    s = np.asarray(scores)
    if s.ndim == 1:
        s = s[:, None]
    order = np.argsort(-s, axis=0, kind="stable")
    ys, ss = y[order], np.take_along_axis(s, order, axis=0)
    tp = np.cumsum(ys, axis=0)
    ends = np.r_[ss[:-1] != ss[1:], np.ones((1, s.shape[1]), dtype=bool)]
    # Endpoint index for each tied row, obtained by reverse minimum accumulation.
    indices = np.arange(len(y))[:, None]
    endidx = np.minimum.accumulate(np.where(ends, indices, len(y))[::-1], axis=0)[::-1]
    precision = np.take_along_axis(tp, endidx, axis=0) / (endidx + 1)
    return (precision * ys).sum(0) / y.sum()

