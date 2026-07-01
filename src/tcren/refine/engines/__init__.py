"""Peptide-modelling engine registry.

Engines share one interface (:class:`~tcren.refine.engines.base.Engine`): given a chain-typed
structure whose peptide is already threaded to the candidate sequence plus the predicted anchors, they
return a :class:`~tcren.refine.engines.base.ModelResult`. Two engines run out of the box (``dope``,
``ccd``); two are optional, gated on heavy installs (``openmm``, ``promod3``).
"""

from __future__ import annotations

from .base import Engine, EngineUnavailable, ModelResult
from .ccd import CcdEngine
from .dope import DopeEngine
from .openmm_engine import OpenMMEngine
from .promod3_engine import ProMod3Engine

ENGINES: dict[str, Engine] = {
    e.name: e
    for e in (DopeEngine(), CcdEngine(), OpenMMEngine(), ProMod3Engine())
}

__all__ = ["ENGINES", "Engine", "EngineUnavailable", "ModelResult",
           "DopeEngine", "CcdEngine", "OpenMMEngine", "ProMod3Engine",
           "get_engine", "available_engines"]


def get_engine(name: str) -> Engine:
    """Return the engine registered under ``name`` (raises KeyError with the valid choices)."""
    try:
        return ENGINES[name]
    except KeyError:
        raise KeyError(f"unknown engine {name!r}; choose from {sorted(ENGINES)}") from None


def available_engines() -> list[str]:
    """Names of engines whose backend can actually run in this environment."""
    return [name for name, e in ENGINES.items() if e.available()]
