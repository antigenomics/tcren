"""Peptide register / forced-pose detection and correction.

A generated (AlphaFold / TCRmodel) TCR:pMHC complex can seat the peptide in the **wrong register** ---
anchors not in their MHC pockets, the peptide shifted along the groove --- which corrupts the
TCR-facing contacts the score reads. On the CPL benchmark this is the ila1 anomaly: the forced AF pose
drives the raw TCRen ranking below chance (ROC 0.35), and re-seating the peptide on the correctly
registered crystal recovers it (≈0.77).

:func:`check_register` flags such poses. It always reports the interface clash burden (via
:func:`tcren.clashes.interface_clashes`); when a correctly-registered ``reference`` (a crystal, or any
trusted pose of the same complex) is supplied it also measures the **anchor-Cα RMSD** in the MHC-groove
frame (:func:`tcren.refine.peptide_rmsd`) --- the reliable register signal, since a mis-registered
peptide's anchors land far from where they should sit. Note a heavy clash burden alone does *not* prove
a register error: AlphaFold peptide-swap models are routinely clashy, so register needs the reference.

:func:`fix_register` corrects a wrong-register model by re-threading its peptide sequence onto the
correctly-registered ``template`` backbone and re-refining through the open-source modelling path
(:func:`tcren.refine.model_peptide`) --- the FlexPepDock-functional replacement, no PyRosetta.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from ..clashes import ClashReport, interface_clashes
from ..structure.model import Structure
from .anchors import native_peptide, predict_anchors
from .rmsd import peptide_rmsd


@dataclass(frozen=True, slots=True)
class RegisterReport:
    """Register / forced-pose diagnostic for a generated TCR:pMHC complex.

    Attributes:
        clashes: Interface steric-clash burden of the model (always computed).
        anchors: 0-based peptide anchor indices used for the anchor RMSD.
        backbone_rmsd: Peptide backbone RMSD vs the reference (Å); ``nan`` without a reference.
        anchor_rmsd: Anchor-Cα RMSD vs the reference (Å) --- the register signal; ``nan`` without one.
        groove_rmsd: MHC-groove superposition residual (Å); ``nan`` without a reference.
        wrong_register: ``True``/``False`` when a reference is given (anchor RMSD over the cut),
            ``None`` when it cannot be determined (no reference).
        reason: Human-readable explanation of the verdict.
    """

    clashes: ClashReport
    anchors: tuple[int, ...]
    backbone_rmsd: float
    anchor_rmsd: float
    groove_rmsd: float
    wrong_register: bool | None
    reason: str

    @property
    def suspect(self) -> bool:
        """True if the pose is a wrong register, or (absent a reference) carries a severe clash."""
        return bool(self.wrong_register) or self.clashes.n_severe > 0


def check_register(
    model: Structure,
    reference: Structure | None = None,
    *,
    anchor_rmsd_cut: float = 2.0,
    tolerance: float = 0.4,
    top: int = 8,
) -> RegisterReport:
    """Assess whether ``model``'s peptide is clashing and/or in the wrong register.

    Args:
        model: The generated complex to check (chain-typed; peptide ``chain_type == 'PEPTIDE'``).
        reference: A correctly-registered pose of the *same* complex (crystal / trusted model). When
            given, the anchor-Cα RMSD in the MHC-groove frame decides ``wrong_register``.
        anchor_rmsd_cut: Anchor-Cα RMSD (Å) above which the register is called wrong.
        tolerance: vdW-overlap tolerance passed to :func:`tcren.clashes.interface_clashes`.
        top: Worst clashing pairs to keep in the report.

    Returns:
        A :class:`RegisterReport`.
    """
    clashes = interface_clashes(model, tolerance=tolerance, top=top)
    if reference is None:
        reason = (
            f"no reference supplied — register undetermined; interface clash burden: "
            f"{clashes.n_clashes} clashes ({clashes.n_severe} severe, max {clashes.max_overlap:.2f} Å)"
        )
        return RegisterReport(clashes, (), math.nan, math.nan, math.nan, None, reason)

    seq = native_peptide(model)
    anchors = predict_anchors(seq, model).anchors
    prm = peptide_rmsd(model, reference, anchors)
    signal = prm.anchor_ca_rmsd if not math.isnan(prm.anchor_ca_rmsd) else prm.backbone_rmsd
    wrong = signal > anchor_rmsd_cut
    which = "anchor-Cα" if not math.isnan(prm.anchor_ca_rmsd) else "backbone"
    reason = (
        f"{which} RMSD {signal:.2f} Å vs reference "
        f"({'>' if wrong else '≤'} {anchor_rmsd_cut} Å cut); "
        f"{clashes.n_severe} severe clashes"
    )
    return RegisterReport(
        clashes, anchors, prm.backbone_rmsd, prm.anchor_ca_rmsd, prm.groove_rmsd, wrong, reason,
    )


def fix_register(
    model: Structure,
    template: Structure,
    *,
    engine: str = "ccd",
    seed: int = 0,
    **engine_kwargs,
):
    """Correct a wrong-register model by re-threading its peptide onto ``template``'s register.

    The model's peptide *sequence* is threaded onto the ``template``'s correctly-registered backbone
    (backbone + Cβ preserved) and re-refined by the chosen open-source engine, so the returned complex
    carries the model's peptide in the template's groove pose. This is the FlexPepDock-functional
    correction path (:func:`tcren.refine.model_peptide`); no PyRosetta.

    Args:
        model: The wrong-register complex (source of the peptide *sequence*).
        template: A correctly-registered complex of equal peptide length (source of the *pose*), e.g.
            the crystal structure of the same clone.
        engine: Modelling engine (``ccd``/``dope`` run out of the box; ``openmm``/``promod3`` optional).
        seed: RNG seed for engines that sample.
        **engine_kwargs: Forwarded to the engine (e.g. ``n_steps`` for ``dope``).

    Returns:
        A :class:`tcren.refine.engines.base.ModelResult` with the corrected structure, engine score,
        and anchors used.

    Raises:
        ValueError: If the model and template peptides differ in length (equal length is required for
            backbone-preserving substitution).
    """
    from .model import model_peptide  # lazy: pulls the engine registry only when a fix is requested

    seq = native_peptide(model)
    tmpl = native_peptide(template)
    if len(seq) != len(tmpl):
        raise ValueError(
            f"peptide length mismatch: model {len(seq)} ({seq}) vs template {len(tmpl)} ({tmpl}); "
            "register correction re-threads onto the template backbone and needs equal length"
        )
    return model_peptide(template, seq, engine=engine, seed=seed, **engine_kwargs)
