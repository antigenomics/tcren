"""Frozen binder/non-binder classifier over the native interface descriptors.

A 5-feature standardized logistic regression (StandardScaler -> LogisticRegression) frozen to fixed
coefficients — no sklearn at inference, just ``P = sigmoid(b + Σ wᵢ (fᵢ − μᵢ)/σᵢ)``. Fit on the
denoised (tcrnet-consistent) TCRvdb set; denoised in-sample AUC = 0.928, beating AlphaFold/TCRmodel2
confidence (0.872 denoised). The features are AF-orthogonal interface geometry + the CDR1/2-vs-CDR3α
TCRen potential term; shape complementarity is deliberately omitted (it adds only ~0.006 — not worth
the molecular-surface kernel). See scripts/binder_validate.py for the marginal-over-AF validation.

Caveat: coefficients are frozen from a 2-epitope (HLA-A*02:01: GLCTLVAML, YLQPRTFLL) training set;
cross-allele/epitope generalization is untested. Re-fit via scripts/binder_validate.py for new data.
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
    # denoised training set: pp_combo = z(ΣJ_CDR12) − z(ΣJ_CDR3α).
    "pp_z": {"cdr12": (0.2856, 0.6969), "cdr3a": (0.0971, 0.8955)},
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
