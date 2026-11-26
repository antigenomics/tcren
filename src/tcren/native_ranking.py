"""Explicit within-unit ranking with the frozen native reference."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
import polars as pl
from scipy.stats import spearmanr

from .native_reference import VERSION, load_native_reference

if TYPE_CHECKING:
    from .contacts.language import ContactLanguageReference

METHOD = VERSION
_ROW = "__native_screen_row"
_FEATURE_TYPES = {
    "n_atom_contacts_tp": pl.Int64, "dPhi_peptide_sum": pl.Float64,
    "d_charge": pl.Float64, "coverage": pl.Float64, "interlock_frac": pl.Float64,
}


def _numeric(group: pl.DataFrame, name: str) -> pl.Series:
    if name in group.columns:
        return group[name].cast(pl.Float64)
    return pl.Series(name, [None] * group.height, dtype=pl.Float64)


def _z(values: np.ndarray) -> np.ndarray:
    if not np.isfinite(values).all():
        return np.full_like(values, np.nan)
    if np.all(values == values[0]):
        return np.zeros_like(values)
    if values.size == 2:
        # Population-standardized pairs are exactly +/-1. Preserve cancellation
        # when confidence and energy order the same two candidates oppositely.
        return np.array([-1., 1.]) if values[0] < values[1] else np.array([1., -1.])
    sd = values.std()
    return (values - values.mean()) / sd if sd > 0 else np.zeros_like(values)


def _size_coupling(counts: np.ndarray, values: np.ndarray) -> float | None:
    valid = np.isfinite(counts) & np.isfinite(values)
    x, y = counts[valid], values[valid]
    if x.size < 2 or np.ptp(x) == 0 or np.ptp(y) == 0:
        return None
    return float(spearmanr(x, y).statistic)


def _terms(group: pl.DataFrame, task: str, unit: str, native: dict, counts: np.ndarray):
    iptm = _numeric(group, "iptm").to_numpy()
    if np.any(np.isfinite(iptm) & ((iptm < 0) | (iptm > 1))):
        raise ValueError("ipTM must be in [0, 1]")
    p = np.clip(iptm, 1e-6, 1 - 1e-6)
    logit = np.where(np.isfinite(iptm), np.log(p / (1 - p)), np.nan)
    energy = _numeric(group, "dPhi_peptide_sum").to_numpy()
    anchor = language_source = None
    if task == "epitope":
        # The producer re-standardizes the two-term base before weighting it.
        centre = (_z(logit) - _z(energy)) / 2
        raw = {
            "confidence_energy": centre, "charge": -_numeric(group, "d_charge").to_numpy(),
            "coverage": _numeric(group, "coverage").to_numpy(),
        }
        weights = (0.8, 0.1, 0.1)
        missing = (~np.isfinite(iptm)).astype(int) + (~np.isfinite(energy))
        missing += ~np.isfinite(raw["charge"])
        missing += ~np.isfinite(raw["coverage"])
    else:
        peptides = group["peptide"].cast(pl.String).str.strip_chars().str.to_uppercase()
        if peptides.null_count() or peptides.n_unique() != 1 or not peptides[0]:
            raise ValueError(f"receptor ranking unit {unit!r} must contain one peptide")
        if "mhc" in group.columns and group["mhc"].drop_nulls().n_unique() > 1:
            raise ValueError(f"receptor ranking unit {unit!r} contains multiple MHC values")
        anchor = native["energy_by_peptide"].get(peptides[0])
        language = _numeric(group, "contact_language")
        if language.null_count() == group.height:
            language_source = "pi_atom"
            language = _numeric(group, language_source)
        else:
            language_source = "contact_language"
        plddt = _numeric(group, "plddt")
        size = native["size"]
        raw = {"iptm": logit, "plddt": plddt.fill_null(plddt.median()).to_numpy(),
               "language": language.to_numpy()}
        if anchor is not None:
            raw.update(interlock=_numeric(group, "interlock_frac").to_numpy(),
                       energy_shell=-np.abs(energy - anchor["centre"]))
        else:
            raw["retained_edges"] = _numeric(group, "kept_edge_frac").to_numpy()
        raw["size_shell"] = -np.abs(counts - size["centre"]) / size["spread"]
        weights = (1 / len(raw),) * len(raw)
        missing = np.count_nonzero(~np.isfinite(list(raw.values())), axis=0)
    terms = {f"term_{name}": _z(value) * weight
             for (name, value), weight in zip(raw.items(), weights, strict=True)}
    return terms, missing, anchor, language_source


def rank_native_screen(
    table: pl.DataFrame, *, task: str, ranking_unit: str, min_cohort: int = 30,
    reference_path: str | Path | None = None,
) -> pl.DataFrame:
    """Score each explicitly declared ranking unit and preserve every input row.

    ``workflow_score`` is the recorded formula on the supplied unit. ``screen_score``
    applies the cohort floor and the zero-peptide-contact pose policy.
    For ordering, use ascending ``screen_rank``; unscored poses have no rank.
    """
    if task not in ("tcr", "epitope"):
        raise ValueError("native ranking requires explicit task tcr or epitope")
    if not isinstance(min_cohort, int) or isinstance(min_cohort, bool) or min_cohort < 1:
        raise ValueError("min_cohort must be a positive integer")
    required = {ranking_unit, "iptm", "n_atom_contacts_tp"}
    if task == "tcr":
        required.add("peptide")
    missing = required - set(table.columns)
    if missing:
        raise ValueError(f"missing native ranking columns: {', '.join(sorted(missing))}")
    unit_keys = table[ranking_unit]
    units = unit_keys.cast(pl.String)
    if not table.height or unit_keys.null_count() or (units.str.strip_chars() == "").any():
        raise ValueError("native ranking requires nonempty ranking units")
    native = load_native_reference(reference_path)
    pieces = []
    for group in table.with_row_index(_ROW).partition_by(ranking_unit, maintain_order=True):
        unit = group[ranking_unit][0]
        counts = _numeric(group, "n_atom_contacts_tp").to_numpy()
        finite = np.isfinite(counts)
        if np.any(counts[finite] < 0) or np.any(counts[finite] != np.floor(counts[finite])):
            raise ValueError("n_atom_contacts_tp must contain nonnegative integer counts")
        contacts = finite & (counts > 0)
        zero = counts == 0
        terms, missing_terms, anchor, language_source = _terms(group, task, unit, native, counts)
        workflow = np.sum(list(terms.values()), axis=0)
        eligible = int(contacts.sum())
        screen = np.where(contacts & (eligible >= min_cohort), workflow, np.nan)
        screen[zero] = 0
        diagnostics = {
            "score_method": METHOD, "task": task, "ranking_unit": unit,
            "native_reference_sha256": native["reference_sha256"],
            "native_reference_available": (anchor is not None) if task == "tcr" else None,
            "native_energy_centre": anchor["centre"] if anchor is not None else None,
            "language_source": language_source,
            "unit_n_models": group.height, "unit_n_contact_bearing": eligible,
            "unit_n_zero_contact": int(zero.sum()), "unit_n_unknown_contact": int((~finite).sum()),
            "workflow_size_rho": _size_coupling(counts, workflow),
            **{f"{name}_size_rho": _size_coupling(counts, value)
               for name, value in terms.items()},
        }
        pieces.append(group.select(_ROW).with_columns(
            *(pl.lit(value).alias(key) for key, value in diagnostics.items()),
            *(pl.Series(key, value, nan_to_null=True) for key, value in terms.items()),
            pl.Series("n_missing_terms", missing_terms),
            pl.Series("workflow_score", workflow, nan_to_null=True),
            pl.Series("screen_score", screen, nan_to_null=True),
            pl.Series("screen_rank", np.where(contacts, screen, np.nan), nan_to_null=True)
            .rank(method="dense", descending=True),
            pl.Series("no_interface", zero),
        ).with_columns(
            pl.when("no_interface").then(pl.lit("no_peptide_contact"))
            .when(pl.Series(~finite)).then(pl.lit("unknown_contact_count"))
            .when(pl.lit(eligible < min_cohort)).then(pl.lit("small_cohort"))
            .when(pl.col("workflow_score").is_null()).then(pl.lit("incomplete_ranking_unit"))
            .otherwise(None).alias("score_reason"),
        ))
    results = pl.concat(pieces, how="diagonal_relaxed").sort(_ROW).drop(_ROW)
    return table.with_columns(results.get_columns())


def _native_epitope_row(structure) -> dict:
    from .contactmap import ContactMap
    from .pipeline import run, score_row
    from .topology.surface import surface_complementarity, surface_map

    row = {"complex.id": structure.pdb_id, **dict.fromkeys(_FEATURE_TYPES)}
    try:
        result = run(structure, typed=True, superimpose=False, reference_aa="A",
                     reference_mode="structural", atom_contact_count=True)
        row["dPhi_peptide_sum"] = score_row(result)["dPhi_peptide_sum"]
        contacts = ContactMap(structure.pdb_id, result.contacts).tcr_peptide()
        row["n_atom_contacts_tp"] = contacts["n_atom_contacts"].sum()
    except ValueError as exc:
        row["native_energy_error"] = str(exc)
    try:
        complementarity = surface_complementarity(
            surface_map(structure, side="pmhc"), surface_map(structure, side="tcr"))
        row.update({name: complementarity.get(name) for name in ("d_charge", "coverage", "interlock_frac")})
    except ValueError as exc:
        row["native_surface_error"] = str(exc)
    return row


def native_epitope_features(structures, *, threads: int = 0) -> pl.DataFrame:
    """Extract the measured peptide ranking terms from a batch of typed complexes.

    Annotate TCR chain roles and MHC in a single batch before calling. Join
    generator confidence and the ranking-unit key before ranking.
    """
    from .runtime import parallel_map

    batch = _annotated_batch(structures)
    return pl.DataFrame(parallel_map(_native_epitope_row, batch, threads), schema_overrides=_FEATURE_TYPES)


def _annotated_batch(structures) -> list:
    batch = list(structures)
    if not batch or any(s.mhc_calls is None for s in batch):
        raise ValueError("provide a nonempty structure batch with MHC annotation")
    if len({s.pdb_id for s in batch}) != len(batch):
        raise ValueError("native feature extraction requires unique structure IDs")
    return batch


def _native_receptor_row(structure) -> tuple:
    from .contacts.opportunities import native_contact_table, native_opportunities

    atoms = native_contact_table(structure).select(
        "pdb.id", "interface", "chain.id.from", "chain.id.to", "residue.index.from",
        "residue.index.to", "atom.from", "atom.to",
    )
    return _native_epitope_row(structure), native_opportunities(structure), atoms


def native_receptor_features(
    structures, reference: ContactLanguageReference, *, threads: int = 0,
) -> pl.DataFrame:
    """Extract all receptor-ranking inputs from one typed application batch."""
    from .contacts.language import contact_language
    from .contacts.mixture import atom_contact_mixture
    from .contacts.veto import contact_veto_features
    from .runtime import parallel_map

    batch = _annotated_batch(structures)
    rows, opportunities, atoms = zip(*parallel_map(_native_receptor_row, batch, threads), strict=True)
    table = pl.DataFrame(rows, schema_overrides=_FEATURE_TYPES)
    pairs = pl.concat(opportunities)
    terms = (
        contact_language(pairs, reference.models, threads=threads),
        contact_veto_features(pairs, reference.residue_effects),
        atom_contact_mixture(pl.concat(atoms), reference.atom_effects, table["complex.id"].to_list()),
    )
    for part in terms:
        table = table.join(part, left_on="complex.id", right_on="id", how="left",
                           validate="1:1", maintain_order="left")
    return table.with_columns(
        pl.lit(reference.sha256).alias("contact_reference_sha256"),
        pl.lit(len(batch)).alias("contact_language_batch_n"),
    )
