"""The TCRen2 score set: one frozen object, five read-outs, defined for a single structure.

Everything here is a projection of **one** fitted object -- a transform frozen on a hold-out, and a
Gaussian per class over the transformed descriptor coordinates. Nothing is estimated from the rows
being scored, so every score below is defined for one structure with no cohort, and a user can run
it on a single AlphaFold model.

============================  ====  ==========================================================
score                         tier  what is estimated
============================  ====  ==========================================================
:func:`peptide_score`         0     nothing; the direction is fixed by the potential
:func:`pose_score`            1     a covariance over hold-out binders -- no negative, no label
:func:`confidence_residual`   1     the same covariance, read as a conditional mean
:func:`binder_score`          2     class means and covariances, from hold-out binder labels
:func:`channel_scores`        2     the same object, marginalized to one descriptor family
============================  ====  ==========================================================

**The channels are what makes a number explainable.** A marginal of a Gaussian is a sub-block of
its covariance -- exact, closed form, no re-fit -- so asking *which part of the structure says
this complex is real* costs nothing beyond an index:

``placement``
    where the receptor sits in the groove frame.
``interface``
    how much interface it makes, and of what chemistry.
``shape``
    the footprint's shape, free of its size.
``energetics``
    the contact chemistry, in kT.
``mechanics``
    the interface read as a network of breakable springs.

Those five are the same partition the fit-free :func:`tcren.cohort.q_score`,
:func:`tcren.reliability.t_score` and the energy block were built on by hand; the difference is
that the weights inside a channel, and between a channel and the rest, come from one covariance
rather than from a choice. The fit-free scores are still shipped and still reported -- ``S`` leads
the functionally validated receptor screen on its own, and it *composes* with
:func:`binder_score` rather than being replaced by it.

**Provenance, and why this is not the read-out that was withdrawn.** ``P_native`` was removed from
the project because its coefficients were frozen against a training set no reader could
reconstruct. These are frozen against a named one: :func:`holdout_manifest` returns the structure
ids and labels the fit used, they ship inside the wheel, and ``tcren fit-holdout`` regenerates
:data:`MODEL_FILE` from them. The model also records the SHA-256 digest of the descriptor
catalogue it was fitted under, and assessing a feature table written under a different catalogue
raises rather than silently mixing generations.

Example::

    tcren features -s models/ -o feats.tsv
    tcren assess --features feats.tsv -o scores.tsv

or in process::

    import polars as pl
    from tcren.score import score_table
    scores = score_table(pl.read_csv("feats.tsv", separator="\\t"))
"""
from __future__ import annotations

import json
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import numpy as np

from ..recognition import DESCRIPTORS
from .model import Joint
from .transform import Transformer, kind, working_set

__all__ = ["CHANNELS", "CHANNEL_OF", "MODEL_FILE", "ScoreModel", "Transformer", "binder_score",
           "channel_scores", "confidence_residual", "holdout_manifest", "holdout_model",
           "kind", "peptide_score", "pose_score", "residue_deltas", "score_table",
           "working_set"]

#: descriptor family -> the plain-language channel a reader is given. ``potts`` joins
#: ``energetics`` because both are contact energies in kT and splitting them helps nobody.
CHANNEL_OF = {"placement": "placement", "interface": "interface", "topology": "shape",
              "energetics": "energetics", "potts": "energetics", "kinetics": "mechanics"}
CHANNELS = ("placement", "interface", "shape", "energetics", "mechanics")

#: `tcren features` keys its rows on `complex.id`; several benchmark tables use `pdb.id`. Both are
#: carried through untouched, and neither is required.
ID_COLUMNS = ("complex.id", "pdb.id")

#: the frozen fit, and the structures it was fitted on. Both ship inside the wheel.
MODEL_FILE = Path(__file__).resolve().parents[1] / "data" / "holdout_model.npz"
MANIFEST_FILE = Path(__file__).resolve().parents[1] / "data" / "holdout_manifest.csv.gz"

#: ipTM enters as a coordinate of the BINDER Gaussian only, on the logit scale. Every hold-out
#: positive carries it and one whole negative arm does not, so a two-class joint over it would
#: learn "no ipTM implies non-binder" -- a property of the deposit, not of the interface.
CONF_COORD = "iptm_logit"


def _logit(p: np.ndarray, lo: float = 1e-3) -> np.ndarray:
    q = np.clip(np.asarray(p, float), lo, 1 - lo)
    return np.log(q / (1 - q))


@dataclass
class ScoreModel:
    """The frozen hold-out fit. Load it with :func:`holdout_model`, never construct it by hand."""

    transform: Transformer
    joint: Joint
    #: the coordinates a receptor-ranking read-out may use: the five descriptors computed without
    #: the receptor are constant across every structure of one epitope on one allele, so a model
    #: reading them reaches the cohort label without reading an interface.
    receptor_names: tuple[str, ...]
    #: binder-class augmentation: E[logit ipTM | x] comes from this row of the binder covariance.
    conf_mu: float
    conf_cov: np.ndarray
    conf_var: float
    catalogue_digest: str
    n_pos: int
    n_neg: int
    n_epitopes: int
    tcren_version: str

    # ------------------------------------------------------------------ coordinates
    #: below this many usable coordinates a covariance read-out is not worth reporting.
    MIN_COORDINATES = 20

    def coordinates(self, table) -> tuple[np.ndarray, np.ndarray, list[str]]:
        """``(complete-case mask, coordinates, the coordinate names they carry)``.

        Two kinds of absence, handled differently and both deliberately.

        A **column** the table does not have is *marginalized out*. That is exact: a marginal of a
        Gaussian is a sub-block of its covariance, so a table without the ``potts`` family is
        scored by the model restricted to what it does have, not by a model with a hole in it.
        This is what lets ``tcren features -i placement,interface,topology`` still produce a
        binder score.

        A **row** missing one of the columns that are present is dropped rather than imputed: the
        whole construction is a covariance, and a filled-in value is a fabricated correlation.
        """
        have = _columns(table)
        cols = [c for c in self.transform.names if c in have]
        if len(cols) < self.MIN_COORDINATES:
            raise KeyError(
                f"only {len(cols)} of the {len(self.transform.names)} descriptors this model "
                f"reads are in the table. Regenerate with `tcren features -s <pdbs> "
                f"-i placement,interface,topology,energetics,potts`.")
        X = np.column_stack([np.asarray(_column(table, c), float) for c in cols])
        ok = np.isfinite(X).all(1)
        names = self.transform.transform_names(cols)
        return ok, self.transform.transform(X[ok], names=cols), names

    def available(self, table) -> list[str]:
        """The model coordinates this table can supply, in model order."""
        return self.coordinates(table)[2]

    def _subset(self, receptor: bool, have: list[str]) -> list[str]:
        keep = set(have)
        if receptor:
            keep &= set(self.receptor_names)
        return [n for n in self.joint.names if n in keep]

    # ------------------------------------------------------------------ the read-outs
    def pose_score(self, table) -> np.ndarray:
        """Tier 1. How far this interface sits from the manifold real binders occupy.

        A partial Mahalanobis to the binder Gaussian. **No negative and no binder label enter it**,
        which is the same standing as the shipped ``q_score`` and ``t_score``. Higher is more
        plausible: the sign is flipped from the raw distance so that every score in this module
        reads "higher is better".
        """
        ok, Z, have = self.coordinates(table)
        return _place(ok, -self.joint.anomaly_on(Z, have))

    def binder_score(self, table, *, receptor: bool = True) -> np.ndarray:
        """Tier 2. Log-odds that this complex is a genuine recognition interface.

        ``receptor=True`` (the default) marginalizes out the five descriptors computed without the
        receptor. Set it False only when the *peptide* is what varies across the structures being
        compared, as in a combinatorial peptide library, where the presentation energy is signal
        rather than the cohort's name.
        """
        ok, Z, have = self.coordinates(table)
        sub = self._subset(receptor, have)
        return _place(ok, self.joint.log_odds(_align(Z, have, sub), subset=sub))

    def channel_scores(self, table, *, receptor: bool = True) -> dict[str, np.ndarray]:
        """Tier 2, one log-odds per channel: **which part of the structure says so.**

        Each is an exact marginal of the same covariance, so the five are on one scale and can be
        read against each other. They do not sum to :func:`binder_score`, and should not: the whole
        model also reads the correlations between channels, which is the part a per-channel view
        cannot show.
        """
        ok, Z, have = self.coordinates(table)
        keep = set(self._subset(receptor, have))
        out = {}
        for ch in CHANNELS:
            sub = [n for n in _channel_names(self.joint.names, ch) if n in keep]
            if len(sub) >= 3:
                out[ch] = _place(ok, self.joint.log_odds(_align(Z, have, sub), subset=sub))
        return out

    def confidence_residual(self, table, iptm) -> np.ndarray:
        """Tier 1. Reported confidence minus what the coordinates say it should have been.

        ``E[logit ipTM | x]`` under the binder Gaussian, subtracted from the reported value. A
        large positive residual is a model the generator is more certain of than its own geometry
        and chemistry warrant, which is the failure mode a confident non-binder presents as. No
        binder label enters it.
        """
        ok, Z, have = self.coordinates(table)
        j = [self.joint.names.index(n) for n in have]
        mu, S = self.joint.mu[1][j], self.joint.cov[1][np.ix_(j, j)]
        pred = self.conf_mu + np.linalg.solve(S, (Z - mu).T).T @ self.conf_cov[j]
        return _place(ok, _logit(np.asarray(iptm, float))[ok] - pred)


# ---------------------------------------------------------------------- helpers
def _columns(table) -> set[str]:
    return set(getattr(table, "columns", []) or [])


def _column(table, name):
    return table[name].to_numpy() if hasattr(table[name], "to_numpy") else table[name]


def _place(ok: np.ndarray, v: np.ndarray) -> np.ndarray:
    """Scatter a score computed on complete-case rows back to the table's full length."""
    out = np.full(len(ok), np.nan)
    out[ok] = v
    return out


def _align(Z: np.ndarray, have: list[str], want: list[str]) -> np.ndarray:
    """Re-index the columns of `Z` (labelled `have`) onto `want`, which `log_odds` expects."""
    pos = {n: i for i, n in enumerate(have)}
    return Z[:, [pos[n] for n in want]]


def _channel_names(names, channel: str) -> list[str]:
    """Transformed coordinate names belonging to one channel.

    Keyed on the transformed names, because a circular descriptor becomes a ``(cos, sin)`` pair and
    a channel that named the raw descriptor would silently drop both halves of it.
    """
    out = []
    for n in names:
        base = n[:-4] if n.endswith(("_cos", "_sin")) else n
        if CHANNEL_OF.get(DESCRIPTORS.get(base, (None,))[0]) == channel:
            out.append(n)
    return out


@lru_cache(maxsize=1)
def holdout_model(path: str | None = None) -> ScoreModel:
    """The frozen hold-out fit that ships with the package.

    Raises a message naming ``tcren fit-holdout`` rather than a bare ``FileNotFoundError``, because
    a source checkout that has never run the fitter is the common way to arrive here.
    """
    p = Path(path) if path else MODEL_FILE
    if not p.exists():
        raise FileNotFoundError(
            f"{p} is not present. It ships with the wheel; in a source checkout, regenerate it "
            f"with `tcren fit-holdout --features <hold-out feature table> "
            f"--manifest {MANIFEST_FILE}`.")
    z = np.load(p, allow_pickle=False)
    meta = json.loads(str(z["meta"]))
    tr = Transformer(names=list(meta["descriptors"]))
    tr.lam = {k: float(v) for k, v in meta["lam"].items()}
    tr.loc = {k: float(v) for k, v in meta["loc"].items()}
    tr.scale = {k: float(v) for k, v in meta["scale"].items()}
    j = Joint(names=list(meta["coordinates"]))
    j.mu = {0: z["mu0"], 1: z["mu1"]}
    j.cov = {0: z["cov0"], 1: z["cov1"]}
    j.prior = {0: float(meta["prior"][0]), 1: float(meta["prior"][1])}
    j.alpha = {0: float(meta["alpha"][0]), 1: float(meta["alpha"][1])}
    lam, U = np.linalg.eigh(j.cov[1])
    j.lam1, j.U1 = np.maximum(lam, 1e-12), U
    return ScoreModel(transform=tr, joint=j,
                      receptor_names=tuple(meta["receptor_coordinates"]),
                      conf_mu=float(meta["conf_mu"]), conf_cov=z["conf_cov"],
                      conf_var=float(meta["conf_var"]),
                      catalogue_digest=meta["catalogue_digest"],
                      n_pos=int(meta["n_pos"]), n_neg=int(meta["n_neg"]),
                      n_epitopes=int(meta["n_epitopes"]),
                      tcren_version=meta["tcren_version"])


def holdout_manifest():
    """The structures the shipped model was fitted on: id, dataset, epitope and binder label.

    This is what makes the frozen coefficients reproducible rather than merely stated. The
    descriptors themselves are not shipped -- 8,292 rows by 147 columns is 19 MB and a wheel is
    not the place for it -- but they are recomputable: ``tcren fetch-data`` brings down the
    structure sets these ids name, ``tcren features`` writes the table, and ``tcren fit-holdout``
    returns the shipped model from the two.
    """
    import gzip

    import polars as pl
    with gzip.open(MANIFEST_FILE, "rt") as fh:
        return pl.read_csv(fh.read().encode())


# ---------------------------------------------------------------------- tier 0
def peptide_score(table) -> np.ndarray:
    """Tier 0, and nothing is fitted anywhere in it.

    The poly-alanine-referenced recognition energy over the two interfaces the peptide is part of,
    sign-flipped so higher is better. This is the instrument for **peptide ranking against a fixed
    receptor**, and it is not the instrument for receptor ranking: on a receptor benchmark it reads
    below chance, which is a property of the reference frame rather than a defect.
    """
    cols = _columns(table)
    need = {"dPhi_tcr_pep", "dPhi_pep_mhc"}
    if not need <= cols:
        raise KeyError(f"peptide_score needs {sorted(need)}; missing {sorted(need - cols)}")
    return -(np.asarray(_column(table, "dPhi_tcr_pep"), float)
             + np.asarray(_column(table, "dPhi_pep_mhc"), float))


# ---------------------------------------------------------------------- module-level sugar
def pose_score(table, *, model: ScoreModel | None = None) -> np.ndarray:
    """See :meth:`ScoreModel.pose_score`."""
    return (model or holdout_model()).pose_score(table)


def binder_score(table, *, receptor: bool = True, model: ScoreModel | None = None) -> np.ndarray:
    """See :meth:`ScoreModel.binder_score`."""
    return (model or holdout_model()).binder_score(table, receptor=receptor)


def channel_scores(table, *, receptor: bool = True,
                   model: ScoreModel | None = None) -> dict[str, np.ndarray]:
    """See :meth:`ScoreModel.channel_scores`."""
    return (model or holdout_model()).channel_scores(table, receptor=receptor)


def confidence_residual(table, iptm, *, model: ScoreModel | None = None) -> np.ndarray:
    """See :meth:`ScoreModel.confidence_residual`."""
    return (model or holdout_model()).confidence_residual(table, iptm)


def score_table(table, *, receptor: bool = True, iptm=None, model: ScoreModel | None = None):
    """Every read-out for a ``tcren features`` table, as one polars frame.

    Adds ``binder_iptm`` when ``iptm`` is supplied: the naive-Bayes sum of two log-odds, which
    needs no coefficient because both terms are already on that scale and stays defined for a
    single structure. On the functionally validated receptor screen it reads 0.771 against the
    posterior's 0.768 and ipTM's 0.795, and it is the recommended read when a confidence is
    available.
    """
    import polars as pl

    m = model or holdout_model()
    out = {}
    for c in ID_COLUMNS:
        if c in _columns(table):
            out[c] = _column(table, c)
            break
    out["pose_score"] = m.pose_score(table)
    out["binder_score"] = m.binder_score(table, receptor=receptor)
    for ch, v in m.channel_scores(table, receptor=receptor).items():
        out[f"channel_{ch}"] = v
    try:
        out["peptide_score"] = peptide_score(table)
    except KeyError:
        pass
    if iptm is not None:
        out["confidence_residual"] = m.confidence_residual(table, iptm)
        out["binder_iptm"] = out["binder_score"] + _logit(np.asarray(iptm, float))
    return pl.DataFrame(out)


from .explain import residue_deltas  # noqa: E402  (see __all__)
