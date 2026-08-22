"""Substitution-informed smoothing of a sparse pair-count matrix.

A 20x20 potential derived from a few hundred structures is thin, and unevenly so: on the 374
reference crystals the TCR:peptide contacts give a median of about 21 observations per cell, but
tryptophan, cysteine and methionine columns hold a handful and some cells hold none at all. Those
cells are then set by the pseudocount rather than by data, and split-half derivations of the whole
matrix agree at only Pearson *r* = 0.44.

The classical fix is not more pseudocount but a *better* one. Henikoff & Henikoff's
substitution-matrix pseudocounts, as used by PSI-BLAST, replace the flat prior with one that says
what an unobserved cell should look like given the cells that *were* observed and how
interchangeable the residues are. If Ile:Leu contacts are common, Val:Leu is not really unknown.

Two steps.

**A conditional substitution model from BLOSUM62.** The published matrix holds rounded log-odds
:math:`s_{ab} = 2\\log_2 (q_{ab} / p_a p_b)` against its own target frequencies :math:`q` and
background :math:`p`. Neither :math:`q` nor :math:`p` is published alongside the scores, but both
are recoverable: :math:`q_{ab} \\propto p_a p_b 2^{s_{ab}/2}` together with the marginal condition
:math:`\\sum_b q_{ab} = p_a` determines :math:`p` as a fixed point, which
:func:`blosum_background` iterates. No constant is taken on trust and nothing is transcribed. The
conditional :math:`P(a \\mid a') = q_{aa'} / \\sum_a q_{aa'}` is what spreads an observation onto
its neighbours.

The background must be the **matrix's**, not the data's. Anchoring it to the composition being
smoothed looks appealing and is self-defeating: a residue that is rare in the data then gets a
near-zero prior as well, so the cells with least evidence are also the ones the prior refuses to
fill. That is the opposite of the intent.

**Adaptive blending.** The prior for cell :math:`(a, b)` is the substitution-weighted average of
every observed cell,

.. math::
    g(a, b) = \\sum_{a', b'} P(a \\mid a')\\, P(b \\mid b')\\, f(a', b')

and it is mixed in with a weight set by how much that cell actually saw:

.. math::
    \\tilde f(a, b) \\propto \\frac{n(a, b)\\, f(a, b) + \\beta\\, g(a, b)}{n(a, b) + \\beta}

A cell with :math:`n \\gg \\beta` is left where it is; an empty one becomes its prior outright.
This is the Henikoff scheme with :math:`\\beta` as the single knob, and :math:`\\beta = 0` is the
identity.

The blend is per **cell**, not per row, and that is deliberate. A row-wise version preserves the
row marginals exactly, which is tidy, but it can only redistribute mass a row already has -- so a
residue with no observations at all stays empty, and those are precisely the cells the prior exists
to fill. The proportionality is not a defect either: the log-odds
:math:`-\\ln\\big(n_{ab} N / n_{a\\cdot} n_{\\cdot b}\\big)` is invariant to a global rescaling of
the counts, so the normalisation this module applies is for interpretability and changes no derived
energy.
"""

from __future__ import annotations

import numpy as np
import polars as pl

from .model import AA20
from .._provenance import not_in_tcren2

#: Default pseudocount weight. A cell is pulled halfway to its substitution prior at this many
#: observations, so it is the count below which a cell is treated as under-observed.
DEFAULT_BETA = 20.0


def _scores(name: str) -> np.ndarray:
    """The substitution matrix over :data:`~tcren.potential.model.AA20`, via Biopython."""
    from Bio.Align import substitution_matrices

    mat = substitution_matrices.load(name)
    return np.array([[float(mat[a, b]) for b in AA20] for a in AA20])


@not_in_tcren2('Substitution-matrix pseudocounts are under evaluation against TCRen2, not part of it. See docs/potentials.rst for the measurements.')
def blosum_background(name: str = "BLOSUM62") -> np.ndarray:
    """The background frequencies implied by a substitution matrix's own scores.

    With :math:`w_{ab} = 2^{s_{ab}/2}` and :math:`q_{ab} = p_a p_b w_{ab} / Z`, the marginal
    condition :math:`\\sum_b q_{ab} = p_a` reduces to :math:`\\sum_b w_{ab} p_b = Z` for every
    :math:`a` -- one linear system, :math:`W p \\propto \\mathbf{1}`. Recovering the background
    this way means the published scores are the only thing taken from outside.

    Checked against the BLOSUM62 background as usually quoted: this returns Ala 0.082 (0.074),
    Leu 0.090 (0.099), Trp 0.012 (0.013), Cys 0.022 (0.025). The residual is the rounding of the
    scores to integers, which is the only lossy step in the inversion.

    Returns:
        A length-20 array summing to 1, indexed like :data:`~tcren.potential.model.AA20`.

    Raises:
        ValueError: if the solution is not strictly positive, i.e. the matrix does not admit a
            valid background under this model.

    Example:
        >>> import numpy as np
        >>> p = blosum_background()
        >>> bool(np.isclose(p.sum(), 1.0)) and bool((p > 0).all())
        True
    """
    p = np.linalg.solve(np.exp2(_scores(name) / 2.0), np.ones(20))
    if (p <= 0).any():
        raise ValueError(f"{name}: implied background is not strictly positive")
    return p / p.sum()


@not_in_tcren2('As blosum_background.')
def blosum_conditional(name: str = "BLOSUM62") -> np.ndarray:
    """``P(a | a')`` over :data:`~tcren.potential.model.AA20`, from a substitution matrix.

    Args:
        name: Any matrix Biopython can load (``BLOSUM62``, ``BLOSUM45``, ``PAM250``, ...).

    Returns:
        A ``(20, 20)`` array indexed like :data:`~tcren.potential.model.AA20`, whose **columns**
        sum to 1: entry ``[i, j]`` is the probability of residue ``i`` given an observation of
        residue ``j``.

    Example:
        >>> import numpy as np
        >>> p = blosum_conditional()
        >>> bool(np.allclose(p.sum(axis=0), 1.0))
        True
    """
    bg = blosum_background(name)
    q = np.outer(bg, bg) * np.exp2(_scores(name) / 2.0)
    q = (q + q.T) / 2.0  # the published scores are rounded, so the inverse is only near-symmetric
    return q / q.sum(axis=0, keepdims=True)


@not_in_tcren2('As blosum_background.')
def smooth_counts(
    counts: pl.DataFrame,
    beta: float = DEFAULT_BETA,
    matrix: str = "BLOSUM62",
    count_col: str = "count",
) -> pl.DataFrame:
    """Blend a pair-count table with its substitution-informed prior.

    Args:
        counts: Long table with ``residue.aa.from``, ``residue.aa.to`` and ``count_col``. Cells
            absent from the table are treated as zero and are still smoothed, so the result always
            covers the full 20x20 grid.
        beta: Pseudocount weight; the observation count at which a cell is pulled halfway to its
            prior. ``0`` returns the input counts on the full grid, unchanged.
        matrix: Substitution matrix to build the prior from.
        count_col: Name of the count column.

    Returns:
        The same schema, with smoothed counts, rescaled to the input's grand total. That rescaling
        is cosmetic -- the log-odds downstream is invariant to it -- and is applied so the numbers
        remain readable as counts.

    Example:
        >>> import polars as pl
        >>> c = pl.DataFrame({"residue.aa.from": ["I"], "residue.aa.to": ["F"], "count": [9.0]})
        >>> out = smooth_counts(c, beta=20.0)
        >>> bool(abs(out["count"].sum() - 9.0) < 1e-9) and out.height == 400
        True
    """
    if beta < 0:
        raise ValueError(f"beta must be non-negative, got {beta}")
    idx = {a: i for i, a in enumerate(AA20)}
    n = np.zeros((20, 20))
    for a, b, c in counts.select("residue.aa.from", "residue.aa.to", count_col).iter_rows():
        if a in idx and b in idx:
            n[idx[a], idx[b]] += float(c)

    total = n.sum()
    grid = pl.DataFrame([(a, b) for a in AA20 for b in AA20],
                        schema=["residue.aa.from", "residue.aa.to"], orient="row")
    if total <= 0 or beta == 0:
        out = n
    else:
        p = blosum_conditional(matrix)
        f = n / total
        g = p @ f @ p.T          # P(a|a') f(a',b') P(b|b')^T -- spread onto substitutable pairs
        out = (n * f + beta * g) / (n + beta)
        out *= total / out.sum()  # cosmetic: the log-odds downstream is scale-invariant

    return grid.with_columns(
        pl.Series(count_col, [out[idx[a], idx[b]] for a, b in
                              zip(grid["residue.aa.from"], grid["residue.aa.to"])])
    )


@not_in_tcren2('As blosum_background.')
def impute_thin_cells(
    counts: pl.DataFrame,
    min_count: int = 10,
    donors: int = 1,
    matrix: str = "BLOSUM62",
    count_col: str = "count",
) -> pl.DataFrame:
    """Rebuild under-observed cells from their nearest substitutable neighbours.

    The alternative to :func:`smooth_counts`, and a sharper instrument. Smoothing moves *every*
    cell toward a prior averaged over all 400, weighted by substitutability; this leaves
    well-observed cells exactly where they are and rebuilds only the thin ones, each from the
    single closest cell that has enough data. Nothing is blended into a cell that did not need it.

    What transfers is the donor's **enrichment**, not its count. With row and column marginals
    :math:`N_{a\\cdot}, N_{\\cdot b}` and grand total :math:`N`, the independence expectation for a
    cell is :math:`E_{ab} = N_{a\\cdot} N_{\\cdot b} / N` and its enrichment is
    :math:`\\rho_{ab} = n_{ab} / E_{ab}`. A thin cell takes the donor's :math:`\\rho` and keeps its
    own :math:`E`, so the imputed count is :math:`\\rho_{a'b'} E_{ab}`. Enrichment is what the
    log-odds :math:`-\\ln \\rho_{ab}` reads, so this is a statement about the energy and not about
    how often the residue pair happens to occur.

    Donors are ranked by :math:`s(a, a') + s(b, b')` in the substitution matrix over the cells
    holding at least ``min_count`` observations, and a cell is never its own donor. On the
    reference crystals this **subsumes the thin-row case**: cysteine on the TCR side carries 4
    contacts in total, so all 20 of its cells are thin and each is rebuilt from its own best
    partner rather than the whole row being copied from one neighbouring residue.

    Args:
        counts: Long table with ``residue.aa.from``, ``residue.aa.to`` and ``count_col``. Cells
            absent are zero, and zero is thin.
        min_count: Observations a cell needs to be left alone, and to be eligible as a donor.
            ``0`` disables imputation.
        donors: How many nearest donors to average the enrichment over. ``1`` is the nearest
            neighbour outright.
        matrix: Substitution matrix behind the neighbour ranking.
        count_col: Name of the count column.

    Returns:
        The same schema over the full 20x20 grid, rescaled to the input's grand total. As in
        :func:`smooth_counts` the rescaling is cosmetic.

    Raises:
        ValueError: if ``donors`` is not positive, or no cell clears ``min_count``.

    Example:
        >>> import polars as pl
        >>> c = pl.DataFrame({"residue.aa.from": ["I", "L", "C"], "residue.aa.to": ["F", "F", "F"],
        ...                   "count": [40.0, 60.0, 1.0]})
        >>> out = impute_thin_cells(c, min_count=10)
        >>> bool(out["count"].sum() > 0) and out.height == 400
        True
    """
    if donors < 1:
        raise ValueError(f"donors must be positive, got {donors}")
    idx = {a: i for i, a in enumerate(AA20)}
    n = np.zeros((20, 20))
    for a, b, c in counts.select("residue.aa.from", "residue.aa.to", count_col).iter_rows():
        if a in idx and b in idx:
            n[idx[a], idx[b]] += float(c)

    total = n.sum()
    out = n
    if total > 0 and min_count > 0:
        rows, cols = n.sum(axis=1), n.sum(axis=0)
        expected = np.outer(rows, cols) / total
        with np.errstate(divide="ignore", invalid="ignore"):
            rho = np.where(expected > 0, n / expected, 0.0)
        fat = n >= min_count
        if not fat.any():
            raise ValueError(f"no cell holds {min_count} observations; nothing can donate")
        s = _scores(matrix)
        fi, fj = np.nonzero(fat)
        out = n.copy()
        for i, j in zip(*np.nonzero(~fat)):
            near = np.argsort(-(s[i, fi] + s[j, fj]))[:donors]
            out[i, j] = float(np.mean(rho[fi[near], fj[near]])) * expected[i, j]
        out *= total / out.sum()

    grid = pl.DataFrame([(a, b) for a in AA20 for b in AA20],
                        schema=["residue.aa.from", "residue.aa.to"], orient="row")
    return grid.with_columns(
        pl.Series(count_col, [out[idx[a], idx[b]] for a, b in
                              zip(grid["residue.aa.from"], grid["residue.aa.to"])])
    )
