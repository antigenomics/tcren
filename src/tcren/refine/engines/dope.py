"""DOPE rigid-body Monte-Carlo engine (wraps the existing ``tcren._refine`` kernel).

This is the engine that already ships: a knowledge-based rigid-body refine of the peptide pose against
its TCR+MHC partners under the DOPE statistical potential. It does not use the anchors (the whole
peptide moves as a rigid body), but it is the always-available baseline every other engine is compared
against. See :func:`tcren.refine.refine_peptide`.
"""

from __future__ import annotations

from ...structure.model import Structure
from ..anchors import Decomposition
from .base import EngineUnavailable, ModelResult


class DopeEngine:
    name = "dope"

    def available(self) -> bool:
        try:
            from ... import _refine  # noqa: F401
        except ImportError:
            return False
        return True

    def run(self, structure: Structure, decomp: Decomposition, *, seed: int = 0,
            n_steps: int = 2000) -> ModelResult:
        if not self.available():
            raise EngineUnavailable("tcren._refine extension not built (run pip install -e .)")
        from .. import refine_peptide

        refined, energy = refine_peptide(structure, n_steps=n_steps, seed=seed)
        return ModelResult(refined, energy, self.name, decomp.anchors,
                           info={"potential": "DOPE", "mode": "rigid-body MC"})
