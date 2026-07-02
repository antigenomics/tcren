"""One-call facade over the tcren pipeline for the paper notebooks.

:func:`summarize_structure` turns a single TCR-pMHC structure into a bundle of
ready-to-tabulate polars frames by composing the existing milestones:

* **S1+S2** — :func:`tcren.pipeline.run` (annotate → orient → contacts → per-interface
  scores), giving the ``scores``, ``markup`` and ``contacts`` frames;
* **S3** — :func:`tcren.scoring_rank.percentile_rank` of the structure's native peptide
  against a random pMHC background, giving the ``rank`` frame;
* **S4** — :func:`tcren.ddg.alanine_scan` of the native peptide, giving the ``ddg`` frame.

Nothing is re-derived here: the facade only orchestrates the milestone functions and
collects their outputs. The ``scores`` frame is therefore byte-identical to what
``run`` produces for the same structure and arguments.
"""

from __future__ import annotations

from pathlib import Path

import polars as pl

from .contactmap import ContactMap
from .ddg import alanine_scan
from .pipeline import run
from .potential import Potential
from .scoring_rank import percentile_rank
from .structure.io import import_structure
from .structure.model import PEPTIDE_TYPE, Structure

#: Columns of the ``ddg`` frame (matches :func:`tcren.ddg.alanine_scan`).
_DDG_SCHEMA = {"pos": pl.Int64, "wt_aa": pl.Utf8, "ddG": pl.Float64}


def _native_peptide(structure: Structure) -> str:
    """Return the one-letter sequence of the structure's peptide chain."""
    for chain in structure.chains:
        if chain.chain_type == PEPTIDE_TYPE:
            return chain.sequence()
    raise ValueError(f"no peptide chain found in {structure.pdb_id!r}")


def summarize_structure(
    structure: str | Path | Structure,
    *,
    organism: str = "human",
    superimpose: bool = True,
    potentials: dict[str, str | Potential | None] | None = None,
    tcr_regions: str = "all",
    background: int = 1000,
    seed: int = 0,
    alanine: bool = False,
    contact_weight: str = "residue",
) -> dict[str, pl.DataFrame]:
    """Summarise one TCR-pMHC structure into a bundle of tables.

    Composes the tcren milestones on a single structure: the full pipeline (S1+S2)
    for per-interface energies, region markup and contacts; the percentile rank of the
    native peptide against a random background (S3); and, optionally, the per-position
    alanine scan (S4).

    Args:
        structure: A structure file (any tcren-readable format) or a parsed
            :class:`~tcren.structure.model.Structure`.
        organism: Organism for TCR annotation.
        superimpose: Also orient onto the canonical database (adds ``rmsd`` to ``scores``).
        potentials: Optional per-interface potential override, forwarded to
            :func:`tcren.pipeline.run`. The TCR↔peptide potential is also used for the
            ``rank`` and ``ddg`` frames.
        tcr_regions: Which TCR regions to keep on the TCR side (``"all"`` default,
            ``"cdr"`` or ``"cdr+fr"``); forwarded to every milestone.
        background: Number of random background peptides for the percentile rank (S3).
        seed: Random seed for the background generation.
        alanine: When ``True``, compute the per-position alanine scan (S4) for the
            ``ddg`` frame. When ``False`` (default) the scan is skipped and ``ddg`` is an
            empty frame with the same schema.
        contact_weight: ``"residue"`` (default, legacy: each contacting residue pair
            contributes ``potential[a, b] x 1``) or ``"atomic"`` (each pair contributes
            ``potential[a, b] x n_atom_contacts``). Applies to all three interfaces of the
            ``scores`` frame and to the ``rank``/``ddg`` frames. ``"residue"`` keeps every
            output byte-identical to the legacy facade.

    Returns:
        Mapping with five polars frames:

        * ``scores`` — one row of per-interface energies (and ``rmsd`` if superimposed),
          identical to :func:`tcren.pipeline.run`'s scores.
        * ``rank`` — one row: the native peptide's energy and its ``rank_pct`` against the
          background.
        * ``ddg`` — the alanine scan (columns ``pos``/``wt_aa``/``ddG``); empty unless
          ``alanine=True``.
        * ``markup`` — the per-residue region-markup table.
        * ``contacts`` — the annotated residue-contact table.
    """
    s = structure if isinstance(structure, Structure) else import_structure(structure)

    # S1 + S2: full pipeline (parses nothing extra — reuses the parsed structure and
    # annotates its chains in place, which is what lets us read the native peptide below).
    result = run(
        s,
        organism=organism,
        superimpose=superimpose,
        potentials=potentials,
        tcr_regions=tcr_regions,
        contact_weight=contact_weight,
    )
    native = _native_peptide(s)

    # The TCR↔peptide potential drives the rank/ddg frames; resolve it exactly as the
    # pipeline does so an override is honoured consistently across all three frames.
    tcr_peptide_pot = _resolve_tcr_peptide_potential(potentials)

    cm = ContactMap.from_structure(s, count_atoms=(contact_weight == "atomic"))

    scores = pl.DataFrame(
        {"pdb.id": result.pdb_id, "rmsd": result.rmsd, **result.scores}
    )

    rank_row = percentile_rank(
        cm,
        native,
        tcr_peptide_pot,
        n_background=background,
        seed=seed,
        tcr_regions=tcr_regions,
        contact_weight=contact_weight,
    )
    rank = pl.DataFrame({"pdb.id": result.pdb_id, **rank_row})

    if alanine:
        ddg = alanine_scan(
            cm, native, tcr_peptide_pot, tcr_regions=tcr_regions,
            contact_weight=contact_weight,
        )
    else:
        ddg = pl.DataFrame(schema=_DDG_SCHEMA)

    return {
        "scores": scores,
        "rank": rank,
        "ddg": ddg,
        "markup": result.markup,
        "contacts": result.contacts,
    }


def _resolve_tcr_peptide_potential(
    potentials: dict[str, str | Potential | None] | None,
) -> Potential:
    """Resolve the TCR↔peptide potential the same way :func:`tcren.pipeline.run` does."""
    from .pipeline import _resolve_potentials

    return _resolve_potentials(potentials)["tcr_peptide"]
