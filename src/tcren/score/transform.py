"""Descriptor -> a coordinate a Gaussian can live on.

Every descriptor gets ONE transform, chosen from its unit and its operator, and the choice is
monotone and variance-stabilising throughout. This is the first stage of :mod:`tcren.score`: the
descriptors :func:`tcren.recognition.recognition_table` emits live on nine different kinds of
scale -- bounded fractions, Poisson tallies, angles, energies in kT -- and a Gaussian cannot be
fitted across them until each is mapped to something with a comparable spread.

**Not a rank or quantile transform, and that is measured rather than argued.** Mapping each
marginal onto a uniform through the hold-out binder CDF took the per-cohort median ROC-AUC from
0.630 to 0.543 on the six templated cohorts and from 0.613 to 0.507 on the sixteen non-templated
ones. The reason is structural: the signal this pass is built on is that binders occupy a NARROW
range on certain axes while non-binders scatter, and flattening a marginal to uniform deletes
exactly that. Anything that equalises spread is the wrong tool here.

The unit string (`DETAIL[name][0]`) is the primary selector; :data:`OPERATOR` overrides it where
the unit is too coarse. Two cases where it is:
`ratio` holds six Pearson correlations on [-1,1] beside four unbounded means of products, and
`count` holds true Poisson tallies beside Hill numbers, which are exp(entropy) and want a log.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..descriptors.catalogue import OPERATOR
from ..recognition import DESCRIPTORS, DETAIL

#: operator name -> tuple of descriptors, inverted to descriptor -> operator
_OP = {n: op for op, names in OPERATOR.items() for n in names}

#: `pitch` reads the generator's confidence rather than the interface and is banned as a feature
#: by the library itself. The five `involves_tcr = False` columns are constant within an
#: epitope x allele cohort, so a receptor-ranking model reading one reaches the cohort label
#: without reading an interface -- excluded for the receptor tasks, restored for CPL where the
#: peptide is what varies and `Phi_pep_mhc` is a legitimate signal.
BANNED = ("pitch",)
NO_RECEPTOR = tuple(k for k, v in DESCRIPTORS.items() if not v[1])

#: Exactly determined by other emitted columns, so each contributes a null direction to any
#: covariance built over them. Every entry is also flagged in :data:`tcren.recognition.STATUS`,
#: and all fifteen were verified over 21,939 corpus structures at max relative residual 1.7e-11.
#: With them removed the catalogue is FULL RANK -- the list is closed, not partial.
DETERMINED = ("fp_chi_r7", "fp_chi_r8", "D1_cell", "J_cell", "offset", "n_loop_contacts",
              "neg_energy", "S_tot", "aniso", "couple_total", "crossing",
              # Dropped from the MODELLING BASIS only. The catalogue keeps emitting all four --
              # they are the familiar names a reader looks for, and Appendix A quotes two of them.
              "sc_gap_index", "dPhi_tcr_soft", "n_contacts_tm", "n_contacts_tp")

_HALDANE = 0.5  # the same correction `L_canon` uses, so a bounded 0 or 1 does not go infinite


def kind(name: str) -> str:
    """The transform class for one descriptor."""
    unit, op = DETAIL[name][0], _OP.get(name, "")
    if name == "mhc_class_bin":
        return "categorical"
    if op == "hill" and unit == "count":
        return "log"           # D1, D2, S are exp(H) and 1/sum p^2 -- log recovers the entropy
    if op == "correlation":
        return "fisher"        # six Pearson r on [-1, 1]
    if unit in ("deg", "rad") and name in ("dock_torsion", "crossing_signed"):
        return "cossin"        # genuinely circular: wraps at +-pi / +-180
    if unit == "fraction":
        return "logit"
    if unit == "signed fraction":
        return "atanh"
    if unit == "count":
        return "linear" if name.startswith("fp_chi") else "anscombe"
    if unit == "cosine":
        return "linear"        # already a bounded, near-linear direction cosine
    return "yeo"               # A, A^2, A^3, kT, log-odds, ratio, N/m, N, J, deg


def _finite(x: np.ndarray) -> np.ndarray:
    return x[np.isfinite(x)]


@dataclass
class Transformer:
    """Fit the per-column parameters on a reference set, then apply anywhere.

    Only ``yeo`` carries a fitted parameter (the Yeo-Johnson lambda, and the mean/sd after it);
    every other class is a fixed function. Yeo-Johnson rather than Box-Cox because `sc_gap_mean`
    (median -1.7 A over 60 crystals), `shift_u`, `shift_w` and `m_face_*` are legitimately
    negative and Box-Cox is undefined there.
    """
    names: list[str]
    lam: dict[str, float] = field(default_factory=dict)
    loc: dict[str, float] = field(default_factory=dict)
    scale: dict[str, float] = field(default_factory=dict)

    def out_names(self, names: list[str] | None = None) -> list[str]:
        """Coordinate names, one per descriptor except a circular one, which yields cos and sin."""
        out = []
        for n in (self.names if names is None else names):
            out.extend([f"{n}_cos", f"{n}_sin"] if kind(n) == "cossin" else [n])
        return out

    #: alias, for callers that read it as "what will `transform(X, names=...)` produce"
    transform_names = out_names

    def _raw(self, col: np.ndarray, n: str) -> np.ndarray:
        k = kind(n)
        with np.errstate(divide="ignore", invalid="ignore"):
            if k == "logit":
                p = np.clip(col, 0.0, 1.0)
                # Haldane-Anscombe on the count scale the fraction came from is unavailable here,
                # so the offset is applied to the probability directly -- the standard empirical
                # logit, and the only choice that keeps a hard 0 or 1 finite.
                eps = _HALDANE / 1000.0
                return np.log((p + eps) / (1 - p + eps))
            if k == "atanh":
                return np.arctanh(np.clip(col, -1 + 1e-6, 1 - 1e-6))
            if k == "fisher":
                return np.arctanh(np.clip(col, -1 + 1e-6, 1 - 1e-6))
            if k == "log":
                return np.log(np.maximum(col, 1e-12))
            if k == "anscombe":
                return 2.0 * np.sqrt(np.maximum(col, 0.0) + 0.375)
            if k in ("linear", "categorical"):
                return col.astype(float)
            if k == "yeo":
                return _yeo(col, self.lam.get(n, 1.0))
        return col.astype(float)

    def fit(self, X: np.ndarray) -> "Transformer":
        from scipy.stats import yeojohnson_normmax
        for j, n in enumerate(self.names):
            if kind(n) == "yeo":
                v = _finite(X[:, j])
                # a degenerate or tiny column has no lambda to find; 1.0 is the identity
                try:
                    self.lam[n] = float(yeojohnson_normmax(v)) if v.size > 50 and v.std() > 0 else 1.0
                except Exception:
                    self.lam[n] = 1.0
                self.lam[n] = float(np.clip(self.lam[n], -3.0, 3.0))
        Z = self._apply_raw(X)
        for j, n in enumerate(self.out_names()):
            v = _finite(Z[:, j])
            self.loc[n] = float(v.mean()) if v.size else 0.0
            s = float(v.std(ddof=1)) if v.size > 1 else 1.0
            self.scale[n] = s if s > 1e-12 else 1.0
        return self

    def _apply_raw(self, X: np.ndarray, names: list[str] | None = None) -> np.ndarray:
        cols = []
        for j, n in enumerate(self.names if names is None else names):
            c = X[:, j].astype(float)
            if kind(n) == "cossin":
                a = np.deg2rad(c) if DETAIL[n][0] == "deg" else c
                cols.extend([np.cos(a), np.sin(a)])
            else:
                cols.append(self._raw(c, n))
        return np.column_stack(cols)

    #: transformed coordinates are clipped to +-CLIP reference standard deviations.
    #: A power transform fitted on one population extrapolates on another, and the extrapolation is
    #: a power: one CPL structure (KMFLYQEEVE, clone mel8) has a standoff height of -2.51 A where
    #: the hold-out minimum is +6.59 A, and Yeo-Johnson's negative branch turns that into ``|z|`` =
    #: 5.1e4 while the 99.9th percentile of the same cohort is 8.5. That structure IS maximally
    #: anomalous and should score as such -- the clip preserves that -- but at 5e4 it dominates
    #: every sum it enters and takes an out-of-fold R^2 to -1.1e4. Ten reference sd is far outside
    #: anything real; the difference between 10 and 5e4 is arithmetic, not biology.
    CLIP = 10.0

    def transform(self, X: np.ndarray, *, names: list[str] | None = None,
                  count_clipped: bool = False):
        """Transformed and standardized against the fitted reference. Non-finite stays non-finite.

        `names` scores a SUBSET of the fitted descriptors -- the columns of `X`, in order. The
        parameters are per descriptor, so a subset is exact rather than an approximation; it is how
        a feature table that omits a whole family is still scored, by marginalization downstream.
        """
        Z = self._apply_raw(X, names=names)
        out = self.out_names(names)
        mu = np.array([self.loc[n] for n in out])
        sd = np.array([self.scale[n] for n in out])
        Z = (Z - mu) / sd
        Z = np.where(np.isfinite(Z), Z, np.nan)
        n_clip = int(np.nansum(np.abs(Z) > self.CLIP))
        Z = np.clip(Z, -self.CLIP, self.CLIP)
        return (Z, n_clip) if count_clipped else Z


def _yeo(x: np.ndarray, lam: float) -> np.ndarray:
    """Yeo-Johnson, defined on all of R. scipy's `yeojohnson` refuses a NaN, so this is elementwise."""
    x = x.astype(float)
    out = np.full_like(x, np.nan)
    ok = np.isfinite(x)
    v = x[ok]
    pos = v >= 0
    r = np.empty_like(v)
    if abs(lam) < 1e-8:
        r[pos] = np.log1p(v[pos])
    else:
        r[pos] = ((v[pos] + 1.0) ** lam - 1.0) / lam
    if abs(lam - 2.0) < 1e-8:
        r[~pos] = -np.log1p(-v[~pos])
    else:
        r[~pos] = -(((-v[~pos] + 1.0) ** (2.0 - lam)) - 1.0) / (2.0 - lam)
    out[ok] = r
    return out


def working_set(*, receptor_task: bool = True) -> list[str]:
    """The descriptors a model may read: catalogue minus determined, banned, and (optionally)
    the five computed without the receptor."""
    drop = set(DETERMINED) | set(BANNED)
    if receptor_task:
        drop |= set(NO_RECEPTOR)
    else:
        drop |= {"mhc_class_bin"}       # a stratifier, never a coordinate
    return [n for n in DESCRIPTORS if n not in drop]
