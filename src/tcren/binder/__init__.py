"""TCR binder/non-binder classification from AF-orthogonal interface geometry.

The shipped model (:func:`binder_score`, :data:`BINDER_MODEL`) scores a TCR-pMHC complex from native
interface descriptors (`tcren._geom`) plus the CDR1/2-vs-CDR3α TCRen potential term — signal that
beats AlphaFold/TCRmodel2 confidence for ranking candidate TCRs against a fixed pMHC (denoised AUC
0.928 vs AF 0.872). Feature extraction (:func:`binder_features`) is added once its native potential
term is validated; the frozen classifier is available now.
"""

from __future__ import annotations

from .model import BINDER_MODEL, FEATURES, binder_score

__all__ = ["binder_score", "BINDER_MODEL", "FEATURES"]


def __getattr__(name):  # lazy: features pulls in _geom + scoring, keep import light
    if name == "binder_features":
        from .features import binder_features
        return binder_features
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
