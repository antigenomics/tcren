"""Common types for peptide-modelling engines.

Every engine takes a chain-typed :class:`~tcren.structure.model.Structure` whose peptide has already
been threaded to the candidate sequence (via :func:`tcren.refine.substitute_peptide`) plus the
predicted anchors, and returns a :class:`ModelResult`: the refined structure, an engine-specific
energy/score, and bookkeeping. Engines that need an unavailable dependency raise
:class:`EngineUnavailable` at call time (never at import), so ``import tcren`` always succeeds.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Protocol

from ...structure.model import Structure
from ..anchors import Decomposition


class EngineUnavailable(RuntimeError):
    """Raised when an engine's backend (OpenMM, ProMod3, …) is not installed."""


@dataclass(slots=True)
class ModelResult:
    """Outcome of (re)modelling a peptide into the groove."""

    structure: Structure
    energy: float  # engine-specific (DOPE energy, CCD anchor-RMSD, OpenMM potential, …)
    engine: str
    anchors: tuple[int, ...]
    iterations: int = 0
    info: dict = field(default_factory=dict)


class Engine(Protocol):
    """A peptide-modelling backend."""

    name: str

    def available(self) -> bool:
        """True if the backend can run in this environment (cheap import probe, no side effects)."""
        ...

    def run(self, structure: Structure, decomp: Decomposition, *, seed: int = 0) -> ModelResult:
        """Model the (already substituted) peptide; return the refined structure + score."""
        ...
