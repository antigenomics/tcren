"""Interface-sanity flag for a modelled complex.

The fitted binder score this package once carried was removed in 2.26.0 -- its coefficients were
frozen against a training set that no longer exists. What remains is the pre-energy check that a
TCR:pMHC interface is a plausible dock at all, which is a rule over contact count and docking
geometry rather than a model.
"""

from .noise import is_real_interface

__all__ = ["is_real_interface"]
