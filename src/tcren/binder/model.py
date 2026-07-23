"""Frozen binder/non-binder classifier over the native interface descriptors.

A 5-feature standardized logistic regression (StandardScaler -> LogisticRegression) frozen to fixed
coefficients — no sklearn at inference, just ``P = sigmoid(b + Σ wᵢ (fᵢ − μᵢ)/σᵢ)``. The features are
AF-orthogonal interface geometry + the CDR1/2-vs-CDR3α TCRen potential term; shape complementarity is
deliberately omitted (it adds only ~0.006 — not worth the molecular-surface kernel).

**Reported performance is on raw labels** (``padj < 1e-5``, no label cleaning): TCRvdb macro ROC-AUC
**0.796** / pooled **0.810**, against AlphaFold/TCRmodel2 ipTM 0.794 / 0.793. Those are the numbers to
quote.

Label denoising (TCRNET motif-cluster consistency) is a **separate algorithm** and is not part of this
package's evaluation. Numbers computed on denoised labels measure tcren *and* that algorithm together,
so they are not reported here.

Caveats, both real:
  * The coefficients were **fit on a denoised subset**, so the training labels went through that
    separate algorithm even though the reported evaluation does not. Re-fitting on raw labels is an
    open item — see ``scripts/binder_validate.py``.
  * They are frozen from a 2-epitope (HLA-A*02:01: GLCTLVAML, YLQPRTFLL) training set;
    cross-allele/epitope generalization is untested.
"""

from __future__ import annotations

import math

#: Feature order the coefficients below are aligned to.
FEATURES = ("pm_cov_ntcr", "chain_balance", "n_hbond", "dSASA", "pp_combo")

BINDER_MODEL = {
    "features": FEATURES,
    "mu": (26.6498, 0.3277, 7.3387, 1951.2735, 0.0574),
    "sigma": (4.1887, 0.1506, 4.6099, 287.7851, 1.3329),
    "w": (0.9686, 1.0221, 0.5189, 1.1133, 1.0624),
    "b": -0.8275,
    # Frozen per-dataset z-constants for the potential term (raw CDR-sum mean/sd), from the
    # training set (denoised; see docstring caveat): pp_combo = z(ΣJ_CDR12) − z(ΣJ_CDR3α).
    "pp_z": {"cdr12": (0.2856, 0.6969), "cdr3a": (0.0971, 0.8955)},
    # Raw-label performance -- the reportable numbers. `denoised_auc` retained only as a
    # provenance record of the fit; denoising is a separate algorithm (see module docstring).
    "macro_auc_raw": 0.796,
    "pooled_auc_raw": 0.810,
    "denoised_auc": 0.928,
}


def binder_score(feats: dict[str, float]) -> float:
    """P(binder) from the 5 native descriptors (keys = :data:`FEATURES`)."""
    m = BINDER_MODEL
    z = m["b"]
    for f, mu, sg, w in zip(m["features"], m["mu"], m["sigma"], m["w"]):
        z += w * (float(feats[f]) - mu) / sg
    return 1.0 / (1.0 + math.exp(-z))


def _demo() -> None:
    """Self-check: a strong, well-buried interface scores high; a weak one low."""
    strong = {"pm_cov_ntcr": 34, "chain_balance": 0.45, "n_hbond": 14,
              "dSASA": 2400, "pp_combo": 1.6}
    weak = {"pm_cov_ntcr": 18, "chain_balance": 0.10, "n_hbond": 2,
            "dSASA": 1500, "pp_combo": -1.2}
    ps, pw = binder_score(strong), binder_score(weak)
    assert 0.0 <= pw < 0.5 < ps <= 1.0, (pw, ps)
    print(f"binder_score demo: strong={ps:.3f}  weak={pw:.3f}  OK")


if __name__ == "__main__":
    _demo()
