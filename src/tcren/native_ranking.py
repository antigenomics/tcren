"""Explicit within-unit ranking with the frozen native reference."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import spearmanr

from .native_reference import VERSION, load_native_reference

METHOD = VERSION
_ROW = "__native_screen_row"


def _numeric(group: pl.DataFrame, name: str) -> pl.Series:
    if name in group.columns:
        return group[name].cast(pl.Float64)
    return pl.Series(name, [None] * group.height, dtype=pl.Float64)


def _z(values: np.ndarray) -> np.ndarray:
    if not np.isfinite(values).all():
        return np.full_like(values, np.nan)
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
