"""Task-specific ranking scores over frozen component references.

Epitope ranking and TCR ranking are different tasks and do not share a score. This module owns
the **epitope-ranking** score, where the peptide varies while the receptor and MHC are held
fixed. The separate cohort screen in :mod:`tcren.screen` supports both tasks;
for TCR ranking it keeps ipTM primary and exposes complementary structural readouts.

The score is a fixed arithmetic combination, not a fit. Each component is standardized by a mean
and scale frozen on a held-out training pool whose labels were never read, so the constants are a
change of units rather than learned weights, and the standardized values are averaged with equal
weights. Signs are fixed a priori by what each component means.

**Nothing here normalizes against the query cohort.** A model's score depends on its own
descriptors and the frozen reference, never on the other rows submitted beside it, so scoring one
model alone and scoring it inside a batch give the same number.

The readout is dimensionless and bounded by nothing. It is a ranking score, not a binding free
energy and not a calibrated probability.
"""

from __future__ import annotations

import hashlib
import json
from importlib.resources import files
from pathlib import Path

import numpy as np
import polars as pl
import typer

FORMAT = "tcren-epitope-rank-1"
RESOURCE = "epitope_rank_reference.json"


def bundled_reference() -> dict:
    """The shipped frozen reference, with its own hash recorded for provenance."""
    raw = files("tcren.data").joinpath(RESOURCE).read_bytes()
    reference = json.loads(raw)
    reference["reference_sha256"] = hashlib.sha256(raw).hexdigest()
    return reference


def _validate(reference: dict) -> list[dict]:
    if reference.get("format") != FORMAT:
        raise ValueError(f"unsupported reference format: {reference.get('format')!r}")
    components = reference.get("components")
    if not isinstance(components, list) or not components:
        raise ValueError("reference carries no components")
    for c in components:
        missing = {"name", "column", "direction", "mean", "sd"} - set(c)
        if missing:
            raise ValueError(f"component {c.get('name')!r} missing {sorted(missing)}")
        if c["direction"] not in (1, -1):
            raise ValueError(f"component {c['name']!r} direction must be +1 or -1")
        if not np.isfinite(c["mean"]):
            raise ValueError(f"component {c['name']!r} has a nonfinite mean")
        if not np.isfinite(c["sd"]) or c["sd"] <= 0:
            raise ValueError(f"component {c['name']!r} has a nonpositive scale")
    return components


def epitope_rank(table: pl.DataFrame, reference: dict | None = None) -> pl.DataFrame:
    """Score candidate peptides for one receptor, one row per model.

    Args:
        table: Descriptor table keyed by a unique nonempty ``complex.id``. It must carry every
            component column named in the reference: ``dPhi_tcr_pep`` from
            ``scoring --delta --reference-mode structural`` and ``iptm`` from the generator.
        reference: A frozen reference dict; the bundled one by default.

    Returns:
        One row per input row, in input order, with a standardized contribution and a
        missingness flag per component, ``n_components`` observed, and ``epitope_rank_score``.
        A row with no observed component scores null rather than defaulting to the middle.
    """
    components = _validate(reference if reference is not None else bundled_reference())
    if "complex.id" not in table.columns:
        raise ValueError("complex.id column required")
    ids = table["complex.id"]
    if ids.dtype != pl.String or ids.null_count() or ids.n_unique() != len(ids) or (ids == "").any():
        raise ValueError("complex.id must be unique, nonempty strings")
    absent = [c["column"] for c in components if c["column"] not in table.columns]
    if absent:
        raise ValueError(f"missing component columns: {absent}")

    out = {"complex.id": ids}
    contributions, observed = [], np.zeros(table.height, dtype=np.int16)
    for c in components:
        x = np.asarray(table[c["column"]].cast(pl.Float64), dtype=float)
        if np.isinf(x).any():
            raise ValueError(f"infinite values in {c['column']}")
        seen = np.isfinite(x)
        z = c["direction"] * (x - c["mean"]) / c["sd"]
        out[f"z_{c['name']}"] = np.where(seen, z, np.nan)
        out[f"missing_{c['name']}"] = ~seen
        contributions.append(np.where(seen, z, 0.0))
        observed += seen
    total = np.sum(contributions, axis=0)
    out["n_components"] = observed
    # A row with nothing observed is *unavailable*, carried as a null rather than a NaN or a
    # neutral midpoint: imputing a value there would put an unmeasured row on the same scale as
    # a measured one.
    score = pl.Series(total / np.maximum(observed, 1))
    out["epitope_rank_score"] = score.set(pl.Series(observed == 0), None)
    return pl.DataFrame(out)


def rank_command(
    features: Path = typer.Option(..., "--features", exists=True, dir_okay=False,
                                  help="descriptor table: complex.id, dPhi_tcr_pep, iptm"),
    output: Path = typer.Option(..., "-o", "--output", help="per-model score table"),
    reference: Path = typer.Option(None, "--reference", exists=True, dir_okay=False,
                                   help="frozen reference JSON (default: the bundled one)"),
) -> None:
    """Rank candidate peptides for one receptor against a frozen component reference.

    Inputs are ``scoring --delta --reference-mode structural`` for ``dPhi_tcr_pep`` and the
    generator's ``iptm``. Nothing is fitted at scoring time and no query-cohort normalization is
    applied, so a model's score does not depend on the other rows supplied with it.
    """
    from .provenance import stamp

    table = pl.read_csv(features, separator="\t", infer_schema_length=None,
                        schema_overrides={"complex.id": pl.String})
    if reference is None:
        ref = bundled_reference()
    else:
        raw = Path(reference).read_bytes()
        ref = json.loads(raw)
        ref["reference_sha256"] = hashlib.sha256(raw).hexdigest()
    result = epitope_rank(table, ref)
    result.write_csv(output, separator="\t")
    stamp(output, command="rank-epitopes",
          extra={"format": ref["format"], "reference_sha256": ref["reference_sha256"]})
    typer.echo(f"Scored {result.height} models against reference {ref['reference_sha256'][:16]}")
