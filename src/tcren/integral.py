"""Experimental fixed-reference integration of confidence, chemistry, shape and graph cycles.

Equal bounded contributions are a declared scoring convention, not binding probabilities.
No query-cohort standardization, assay fitting, covariance weights or adaptive signs.
"""
from __future__ import annotations

import hashlib
from importlib.resources import files
import json
from pathlib import Path

import numpy as np
import polars as pl
import typer

FORMAT = "integral-reference-1"
COMPONENTS = ("energy", "geometry", "topology")


def _potential_hash():
    return hashlib.sha256(files("tcren.data").joinpath("TCRen2_potential.csv").read_bytes()).hexdigest()


def _validated(table):
    required = {"complex.id", "mhc_class_bin", "Phi_tcr_pep", "n_contacts_tp",
                "sc_shape", "g_cyclo_frac"}
    if required - set(table.columns):
        raise ValueError(f"missing columns: {sorted(required - set(table.columns))}")
    ids = table["complex.id"]
    if ids.dtype != pl.String or ids.null_count() or ids.n_unique() != len(ids) or (ids == "").any():
        raise ValueError("nonempty unique string complex.id values required")
    values = {}
    for column in required - {"complex.id"} | ({"iptm"} & set(table.columns)):
        if not table.schema[column].is_numeric() and table.schema[column] != pl.Null:
            raise ValueError(f"{column} must be numeric")
        values[column] = np.asarray(table[column].cast(pl.Float64), dtype=float)
        if np.isinf(values[column]).any():
            raise ValueError(f"infinite {column}")
    for column, lo, hi in (("iptm", 0, 1), ("sc_shape", -1, 1), ("g_cyclo_frac", 0, 1)):
        x = values.get(column, np.array([]))
        if ((x < lo) | (x > hi)).any():
            raise ValueError(f"{column} must be in [{lo}, {hi}] or missing")
    cls = values["mhc_class_bin"]
    if not np.isin(cls, [0, 1]).all():
        raise ValueError("explicit finite MHC class (0=I, 1=II) required")
    count, phi = values["n_contacts_tp"], values["Phi_tcr_pep"]
    finite = np.isfinite(count)
    if ((count[finite] < 0) | (count[finite] != np.floor(count[finite]))).any():
        raise ValueError("contact counts must be nonnegative integers or missing")
    if ((count == 0) & np.isfinite(phi) & (phi != 0)).any():
        raise ValueError("nonzero energy with zero contacts")
    mean = np.full(len(table), np.nan)
    np.divide(-phi, count, out=mean, where=np.isfinite(phi) & finite & (count > 0))
    return {"energy": mean, "geometry": values["sc_shape"],
            "topology": values["g_cyclo_frac"],
            "iptm": values.get("iptm", np.full(len(table), np.nan))}, cls


def build_reference(table: pl.DataFrame, *, provenance: dict) -> dict:
    """Freeze finite native class-I component distributions; provenance must name input hash.

    Caller supplies verified native structures. Class filtering cannot establish native origin.
    Reference estimation reads no assay outcomes or generator confidence.
    """
    digest = provenance.get("input_sha256", "")
    if not isinstance(digest, str) or len(digest) != 64 or any(c not in "0123456789abcdef" for c in digest):
        raise ValueError("provenance requires input_sha256")
    values, cls = _validated(table)
    samples = {key: np.sort(values[key][(cls == 0) & np.isfinite(values[key])]).tolist()
               for key in COMPONENTS}
    if any(not x for x in samples.values()):
        raise ValueError("each component needs finite native class-I reference support")
    return {"format": FORMAT, "mhc_class": "I", "potential_sha256": _potential_hash(),
            "n_reference": int((cls == 0).sum()), "n_excluded_class_ii": int((cls == 1).sum()),
            "samples": samples, "provenance": dict(provenance)}


def integral_scores(table: pl.DataFrame, reference: dict) -> pl.DataFrame:
    """Append four bounded contributions, missingness and an equal-weight integral score.

    Missing components receive 0.5; rows with no observed component have null scores.
    Unrelated query rows never change the reference or any existing row's score.
    Requires the bundled TCRen2, all-TCR 5-Angstrom Phi and matching contact counts.
    """
    if reference.get("format") != FORMAT or reference.get("mhc_class") != "I":
        raise ValueError("unsupported integral reference format or MHC class")
    if reference.get("potential_sha256") != _potential_hash():
        raise ValueError("reference potential differs from bundled TCRen2")
    values, cls = _validated(table)
    if (cls != 0).any():
        raise ValueError("integral-reference-1 supports only MHC class I")
    columns = {"complex.id": table["complex.id"]}
    contributions, available = [], np.zeros(len(table), dtype=np.int8)
    for key in ("iptm", *COMPONENTS):
        x = values[key]
        missing = ~np.isfinite(x)
        if key != "iptm":
            r = np.asarray(reference.get("samples", {}).get(key, []), dtype=float)
            if r.ndim != 1 or not len(r) or not np.isfinite(r).all() or (np.diff(r) < 0).any():
                raise ValueError(f"invalid sorted reference samples for {key}")
            if key in {"geometry", "topology"} and ((r < (-1 if key == "geometry" else 0)) | (r > 1)).any():
                raise ValueError(f"out-of-domain reference samples for {key}")
            x = (np.searchsorted(r, x, side="left") + np.searchsorted(r, x, side="right")) / (2 * len(r))
        u = np.where(missing, .5, x)
        columns["integral_" + key] = u
        columns["missing_" + key] = missing
        contributions.append(u)
        available += ~missing
    score = np.mean(contributions, axis=0)
    columns["n_components"] = available
    columns["integral_score"] = pl.Series(score).set(pl.Series(available == 0), None)
    return pl.DataFrame(columns)


def _read(path):
    return pl.read_parquet(path) if path.suffix == ".parquet" else pl.read_csv(
        path, separator="," if path.suffix == ".csv" else "\t", infer_schema_length=None,
        schema_overrides={"complex.id": pl.String, "receptor_id": pl.String})


def _sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def reference_command(
    features: Path = typer.Option(..., "--features", exists=True, dir_okay=False),
    output: Path = typer.Option(..., "-o", "--output"),
):
    """Freeze a class-I reference from verified native features (no assay fitting).

    Input: features -i interface,topology,energetics; caller verifies native origin.
    """
    if output.exists():
        raise FileExistsError(output)
    provenance = {"input_sha256": _sha(features)}
    sidecar = Path(str(features) + ".provenance.json")
    if sidecar.exists():
        provenance["cli_provenance_sha256"] = _sha(sidecar)
    resource = build_reference(_read(features), provenance=provenance)
    with output.open("x") as stream:
        stream.write(json.dumps(resource, indent=2, allow_nan=False) + "\n")
    typer.echo(f"Frozen {resource['n_reference']} native class-I rows; reference {_sha(output)}")


def integrate_command(
    features: Path = typer.Option(..., "--features", exists=True, dir_okay=False),
    reference: Path = typer.Option(..., "--reference", exists=True, dir_okay=False),
    confidence: Path = typer.Option(None, "--confidence", exists=True, dir_okay=False),
    output: Path = typer.Option(..., "-o", "--output"),
):
    """Score class-I descriptors using a fixed reference and optional keyed ipTM table.

    Confidence table: complex.id, iptm. Missing entries remain flagged and neutral.
    Uses bundled TCRen2 at 5 Angstrom; caller must retain input producer provenance.
    """
    from .provenance import stamp
    table = _read(features)
    lineage = {"input_sha256": _sha(features), "reference_sha256": _sha(reference),
               "format": FORMAT, "potential_sha256": _potential_hash()}
    if confidence is not None:
        if "iptm" in table.columns:
            raise ValueError("iptm supplied both in features and confidence table")
        conf = _read(confidence).select("complex.id", "iptm")
        if conf["complex.id"].null_count() or conf["complex.id"].is_duplicated().any():
            raise ValueError("confidence IDs must be unique and nonmissing")
        if conf.join(table.select("complex.id"), on="complex.id", how="anti").height:
            raise ValueError("confidence table contains unexpected IDs")
        table = table.join(conf, on="complex.id", how="left", validate="1:1", maintain_order="left")
        lineage["confidence_sha256"] = _sha(confidence)
    result = integral_scores(table, json.loads(reference.read_text()))
    result.write_csv(output, separator="\t")
    stamp(output, command="integrate", extra=lineage)
    typer.echo(f"Wrote {result.height} integral scores with component missingness")
