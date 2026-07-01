"""Open-source peptide (re)modelling: substitute → predict anchors → model in the groove.

:func:`model_peptide` is the single entry point that replaces the licensed FlexPepDock/MODELLER path.
It threads a candidate peptide onto the groove backbone (:func:`tcren.refine.substitute_peptide`),
predicts which residues anchor into the MHC (:func:`tcren.refine.anchors.predict_anchors`), and hands
off to a modelling engine (``dope`` and ``ccd`` run out of the box; ``openmm`` and ``promod3`` are
optional). The returned :class:`~tcren.refine.engines.base.ModelResult` carries the refined structure,
an engine score, and the anchor set used.
"""

from __future__ import annotations

from ..structure.model import Structure
from .anchors import native_peptide, predict_anchors
from .engines import ENGINES, ModelResult, get_engine
from .substitute import substitute_peptide

__all__ = ["model_peptide", "ModelResult"]


def model_peptide(structure: Structure, new_peptide: str | None = None, *,
                  engine: str = "ccd", seed: int = 0, **engine_kwargs) -> ModelResult:
    """Model ``new_peptide`` into ``structure``'s groove with the chosen engine.

    Args:
        structure: A chain-typed, MHC-annotated complex (peptide chain ``chain_type == 'PEPTIDE'``).
        new_peptide: Candidate one-letter sequence (equal length to the native peptide). ``None``
            re-models the native sequence in place — the self-reconstruction case used by the
            benchmark.
        engine: One of :data:`~tcren.refine.engines.ENGINES` (``dope``, ``ccd``, ``openmm``,
            ``promod3``).
        seed: RNG seed (engines that sample).
        **engine_kwargs: Forwarded to the engine's ``run`` (e.g. ``perturb``/``anchor_targets`` for
            ``ccd``, ``n_steps`` for ``dope``).

    Returns:
        A :class:`ModelResult` with the refined structure, engine score, and anchors used.
    """
    eng = get_engine(engine)
    seq = new_peptide if new_peptide is not None else native_peptide(structure)
    threaded = substitute_peptide(structure, seq)
    decomp = predict_anchors(seq, structure)
    return eng.run(threaded, decomp, seed=seed, **engine_kwargs)
