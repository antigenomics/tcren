"""Apply the recorded contact-language model to directed residue opportunities."""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
import polars as pl
from scipy.sparse import csr_matrix

from ..runtime import resolve_threads, single_threaded
from .mixture import _mixture_fraction

_SWEEPS = ("sweep_max", "sweep_mean", "sweep_n3", "sweep_compact", "pep_gather_max")
_STATISTICS = (
    "mean", "vote", "top20", "bottom20", "spread", "max", "pi",
    "degree", "pi_degree", "gwd", "config", "esp", "pi_config",
)
_RESIDUES = ("A_engagement", "C_multiplicity", "D_identity")
_AROMATIC_ATOMS = ("CG", "CD1", "CD2", "CE1", "CE2", "CE3", "CZ", "CZ2", "CZ3", "CH2")


FORMAT = "tcren-contact-language-1"
VERSION = "native387-contact-language-v1"
_EFFECT_SCHEMA = {
    "interface": pl.String, "block": pl.String, "cell": pl.String,
    "effect": pl.Float64, "complement": pl.Float64,
}


@dataclass(frozen=True)
class ContactLanguageReference:
    """Content-identified effects and coefficients for receptor contact features."""

    sha256: str
    models: list[tuple[pl.DataFrame, dict[str, float]]]
    atom_effects: dict[str, float]
    residue_effects: dict[str, float]


def _finite_mapping(values: dict) -> bool:
    return (isinstance(values, dict) and bool(values)
            and all(isinstance(v, (int, float)) and not isinstance(v, bool)
                    and np.isfinite(v) for v in values.values()))


def load_contact_language_reference(
    path: str | Path, *, expected_sha256: str | None = None,
) -> ContactLanguageReference:
    """Load an explicit receptor reference and record its byte-level identity."""
    raw = Path(path).read_bytes()
    digest = sha256(raw).hexdigest()
    if expected_sha256 is not None and digest != expected_sha256:
        raise ValueError("contact-language reference SHA-256 does not match")
    data = json.loads(raw)
    if not isinstance(data, dict) or data.get("format") != FORMAT or data.get("version") != VERSION:
        raise ValueError("unsupported contact-language reference format or version")
    atoms, residues = data.get("atom_effects"), data.get("residue_effects")
    if not _finite_mapping(atoms) or not _finite_mapping(residues):
        raise ValueError("contact-language reference requires finite atom and residue effects")
    groups = data.get("models")
    if not isinstance(groups, list) or not groups:
        raise ValueError("contact-language reference must contain models")
    names, models = set(), []
    for group in groups:
        if not isinstance(group, dict):
            raise ValueError("invalid contact-language model")
        name = group.get("name")
        coefficients = group.get("coefficients")
        if not isinstance(name, str) or not name or name in names:
            raise ValueError("contact-language models must have unique names")
        if not _finite_mapping(coefficients) or "intercept" not in coefficients:
            raise ValueError("contact-language models require finite coefficients and an intercept")
        names.add(name)
        if not isinstance(group.get("effects"), list) or not group["effects"]:
            raise ValueError("contact-language models require cell effects")
        effects = pl.DataFrame(group["effects"], schema=_EFFECT_SCHEMA, orient="row")
        if (any(effects.null_count().row(0))
                or effects.select("interface", "block", "cell").is_duplicated().any()
                or not effects.select(pl.all_horizontal(pl.col("effect", "complement").is_finite()).all()).item()):
            raise ValueError("contact-language cell effects must be finite and unique")
        models.append((effects, coefficients))
    return ContactLanguageReference(digest, models, atoms, residues)


def _chemical_class(side: str) -> pl.Expr:
    atom = pl.col(f"atom.{side}")
    return (
        pl.when(atom == "N").then(pl.lit("bbN"))
        .when(atom.is_in(["O", "OXT"])).then(pl.lit("bbO"))
        .when(atom.str.starts_with("N")).then(pl.lit("scN"))
        .when(atom.str.starts_with("O")).then(pl.lit("scO"))
        .when(atom.str.starts_with("S")).then(pl.lit("scS"))
        .when(atom.is_in(_AROMATIC_ATOMS) & pl.col(f"aa.{side}").is_in(["F", "Y", "W", "H"]))
        .then(pl.lit("arC")).otherwise(pl.lit("alC"))
    )


def _opportunity_graph(frame: pl.DataFrame) -> pl.DataFrame:
    """Attach full-graph degree and shared-partner counts to contacting edges."""
    nodes = pl.concat([
        frame.select("pdb.id", pl.col(f"chain.id.{side}").alias("chain"),
                     pl.col(f"residue.index.{side}").alias("residue"))
        for side in ("from", "to")
    ]).unique().with_row_index("node")
    for side in ("from", "to"):
        frame = frame.join(
            nodes.rename({"node": f"node.{side}"}),
            left_on=["pdb.id", f"chain.id.{side}", f"residue.index.{side}"],
            right_on=["pdb.id", "chain", "residue"],
            how="inner",
        )
    edges = frame.filter("contacted").select("node.from", "node.to").unique()
    if not edges.height:
        return frame.with_columns(pl.lit(1).alias("ku"), pl.lit(1).alias("kv"), pl.lit(0).alias("esp"))
    u, v = (edges[f"node.{side}"].to_numpy() for side in ("from", "to"))
    endpoints = np.concatenate([u, v])
    adjacency = csr_matrix((np.ones(endpoints.size, dtype=np.int32),
                            (endpoints, np.concatenate([v, u]))),
                           shape=(nodes.height, nodes.height))
    degree = np.diff(adjacency.indptr)
    edges = edges.with_columns(
        pl.Series("ku", degree[u]), pl.Series("kv", degree[v]),
        pl.Series("esp", np.asarray((adjacency @ adjacency)[u, v]).reshape(-1)),
    )
    return frame.join(edges, on=["node.from", "node.to"], how="left").with_columns(
        pl.col("ku", "kv").fill_null(1), pl.col("esp").fill_null(0),
    )


def _sweeps(frame: pl.DataFrame, models: pl.DataFrame) -> pl.DataFrame:
    """Summarize the consecutive peptide positions each contacting receptor residue spans."""
    hits = (frame.filter(pl.col("contacted") & (pl.col("interface") == "tcr_peptide"))
            .select("pdb.id", "node.from", pl.col("residue.index.to").alias("position"))
            .unique().sort("pdb.id", "node.from", "position"))
    receptor = ["pdb.id", "node.from"]
    runs = hits.with_columns(
        (pl.col("position") - pl.int_range(pl.len()).over(receptor)).alias("run"),
    ).group_by(*receptor, "run").len().group_by(receptor).agg(pl.col("len").max().alias("longest"))
    compact = hits.group_by(receptor).agg((
        pl.len() / (pl.col("position").max() - pl.col("position").min() + 1)
    ).alias("compact"))
    sweeps = runs.join(compact, on=receptor).group_by("pdb.id").agg(
        pl.col("longest").max().alias("sweep_max"),
        pl.col("longest").mean().alias("sweep_mean"),
        (pl.col("longest") >= 3).sum().alias("sweep_n3"),
        pl.col("compact").mean().alias("sweep_compact"),
    )
    gather = hits.group_by("pdb.id", "position").len().group_by("pdb.id").agg(
        pl.col("len").max().alias("pep_gather_max"))
    return models.join(sweeps, on="pdb.id", how="left", maintain_order="left").join(
        gather, on="pdb.id", how="left", maintain_order="left",
    ).with_columns(pl.col(*_SWEEPS).fill_null(0))


def _segments(part: pl.DataFrame, likelihood: pl.DataFrame, n: int, block: str) -> dict:
    """Compute each requested summary over a batch using segmented native reductions."""
    part = part.join(likelihood, on="cell", how="inner")
    if not part.height:
        return {stat: np.full(n, 0.5 if stat.startswith("pi") else 0.0) for stat in _STATISTICS}
    segment = part["segment"].to_numpy()
    value = part["effect"].to_numpy()
    if block == "A_engagement":
        value = np.where(part["contacted"].to_numpy(), value, part["complement"].to_numpy())
    elif block == "C_multiplicity":
        weight = part["multiplicity"].to_numpy()
        value = value * weight / max(weight.mean(), 1e-9)
    counts = np.bincount(segment, minlength=n)
    den = np.maximum(counts, 1)
    found = counts > 0

    def reduce(values: np.ndarray) -> np.ndarray:
        return np.bincount(segment, weights=values, minlength=n)

    def average(weight: np.ndarray) -> np.ndarray:
        return reduce(weight * value) / np.maximum(reduce(weight), 1e-12)

    def mixture(values: np.ndarray) -> np.ndarray:
        return np.nan_to_num(_mixture_fraction(segment, values, n), nan=0.5)

    sorted_values = value[np.lexsort((value, segment))]
    starts = np.concatenate([[0], counts.cumsum()])
    cumulative = np.concatenate([[0.0], sorted_values.cumsum()])
    k = np.maximum(1, np.rint(0.2 * counts)).astype(np.int64)
    end = np.minimum(starts[:-1] + k, starts[1:])
    start = np.maximum(starts[1:] - k, starts[:-1])
    bottom = (cumulative[end] - cumulative[starts[:-1]]) / k
    top = (cumulative[starts[1:]] - cumulative[start]) / k
    degree = part["degree"].to_numpy()
    ku, kv = (part[c].to_numpy() for c in ("ku", "kv"))
    config = 1 / np.sqrt(ku.astype(float) * kv)
    mean_degree = reduce(degree) / den
    mean_config = reduce(config) / den
    return {
        "mean": (cumulative[starts[1:]] - cumulative[starts[:-1]]) / den,
        "vote": reduce(value > 0) / den,
        "top20": top, "bottom20": bottom, "spread": top - bottom,
        "max": np.where(found, sorted_values[np.maximum(starts[1:] - 1, 0)], 0.0),
        "pi": mixture(value), "degree": average(degree),
        "pi_degree": mixture(value * degree / np.maximum(mean_degree[segment], 1e-9)),
        "gwd": average(2.0 - 0.5 ** ku - 0.5 ** kv),
        "config": average(config), "esp": average(1.0 + part["esp"].to_numpy()),
        "pi_config": mixture(value * config / np.maximum(mean_config[segment], 1e-9)),
    }


def _block(frame: pl.DataFrame, interface: str, block: str) -> pl.DataFrame:
    name, _, receptor = interface.partition("|")
    selected = frame.filter(pl.col("interface") == name)
    if receptor:
        selected = selected.filter(pl.col("chain.type.from") == receptor)
    selected = selected.filter("within_shell" if block == "A_engagement" else "contacted")
    if block in _RESIDUES:
        cell = pl.concat_str("aa.from", "aa.to", separator="|")
    elif block == "B_chemistry":
        cell = pl.concat_str("atom.from", "atom.to", separator="|")
    elif block == "E_chemclass":
        cell = pl.concat_str(_chemical_class("from"), _chemical_class("to"), separator="|")
    elif block == "F_position":
        cell = pl.concat_str("aa.from", "residue.index.to", separator="@")
    elif block == "G_region":
        cell = pl.concat_str(pl.col("region.type.from").fill_null("none"),
                             "residue.index.to", separator="@")
    else:
        raise ValueError(f"unknown contact-language block {block!r}")
    return selected.with_columns(
        cell.alias("cell"),
        (pl.len().over("node.from") + pl.len().over("node.to")).alias("degree"),
        (pl.col("n_atom_pairs") / (
            pl.col("n_heavy_from").cast(pl.Float64) * pl.col("n_heavy_to")
        ).clip(lower_bound=1)).alias("multiplicity"),
    ).select("segment", "cell", "contacted", "degree", "ku", "kv", "esp", "multiplicity")


@single_threaded
def contact_language(
    opportunities: pl.DataFrame,
    models: Iterable[tuple[pl.DataFrame, Mapping[str, float]]],
    *, threads: int = 0,
) -> pl.DataFrame:
    """Apply frozen effects and linear coefficients, and average the model logits.

    Each model takes a cell table with interface, block, cell, effect and
    complement columns, and a coefficient mapping containing an intercept.
    Input comprises the whole application batch because multiplicity is
    normalized by its mean within a block. Rows without opportunities are
    absent so a caller can retain their missingness when joining scores.
    """
    threads = resolve_threads(threads)
    models = list(models)
    if not models:
        raise ValueError("at least one contact-language model is required")
    ids = opportunities.select("pdb.id").unique().sort("pdb.id").with_row_index("segment")
    schema = {"id": pl.String, "contact_language": pl.Float64, "n_language_models": pl.Int64}
    if not ids.height:
        return pl.DataFrame(schema=schema)
    frame = _opportunity_graph(opportunities.join(ids, on="pdb.id"))
    sweeps = _sweeps(frame, ids)
    blocks = {
        key for effects, _ in models
        for key in effects.select("interface", "block").unique().iter_rows()
    }
    common = {key: _block(frame, *key) for key in sorted(blocks)}

    def predict(model) -> np.ndarray:
        effects, coefficients = model
        if "intercept" not in coefficients or not np.isfinite(list(coefficients.values())).all():
            raise ValueError("contact-language coefficients must be finite and include an intercept")
        parts = effects.partition_by(["interface", "block"], as_dict=True)
        terms: dict[tuple[str, str], dict] = {}
        prediction = np.full(ids.height, coefficients["intercept"], dtype=float)
        for feature, coefficient in coefficients.items():
            if feature == "intercept":
                continue
            if feature in _SWEEPS:
                prediction += coefficient * sweeps[feature].to_numpy()
                continue
            fields = feature.rsplit("|", 2)
            if len(fields) != 3 or fields[-1] not in _STATISTICS:
                raise ValueError(f"invalid contact-language feature {feature!r}")
            interface, block, statistic = fields
            key = (interface, block)
            if key not in parts:
                continue
            if key not in terms:
                terms[key] = _segments(
                    common[key], parts[key].select("cell", "effect", "complement"),
                    ids.height, block,
                )
            prediction += coefficient * terms[key][statistic]
        return prediction

    count = min(threads, len(models))
    if count > 1:
        # Polars and NumPy release the GIL. Share the common blocks across
        # model workers using the runtime's one-thread native default.
        with ThreadPoolExecutor(max_workers=count) as pool:
            predictions = list(pool.map(predict, models))
    else:
        predictions = [predict(model) for model in models]
    return pl.DataFrame({
        "id": ids["pdb.id"], "contact_language": np.mean(predictions, axis=0),
        "n_language_models": [len(models)] * ids.height,
    }, schema=schema)
