"""Small confidence gates for the existing three-component recognition score."""

from pathlib import Path
import json
import time

import numpy as np
import polars as pl
import typer

from .motif_mixture import FEATURES, simplex_weights, linear_score, _metrics, _impute
from .atom_formula import ap_columns
from .profile import digest
from ..runtime import parallel_map, resolve_threads
from ..structure.io import mean_bfactor

COHORT = ["dataset", "stratum", "epitope"]
GATES = ["constant", "complex-iptm", "complex-plddt", "cohort-iptm", "cohort-plddt"]


def _read_confidence_chunk(records):
    result = []
    for record in records:
        path = Path(record["source"])
        if digest(path) != record["sha256"]:
            raise ValueError("Raw coordinate hash changed")
        result.append(dict(id=record["id"], coordinate_plddt=mean_bfactor(path, atom="CA")))
    return result


def confidence_command(
    sources: Path = typer.Option(..., exists=True),
    metadata: list[Path] = typer.Option(..., exists=True),
    output: Path = typer.Option(..., "--output", "-o"),
    threads: int = typer.Option(0, "--threads", "-t"),
):
    """Read generator mean pLDDT; recover missing values from hash-verified model Cα."""
    if output.exists():
        raise ValueError("Output exists")
    records = json.loads(sources.read_text())
    meta = pl.concat(
        [
            pl.read_csv(p, separator="\t", infer_schema_length=10000).select("id", "iptm", "plddt")
            for p in metadata
        ]
    )
    selected = pl.DataFrame({"id": [r["id"] for r in records]}).join(
        meta, on="id", how="left", validate="1:1"
    )
    missing = set(selected.filter(pl.col("plddt").is_null() | pl.col("plddt").is_nan())["id"])
    todo = [r for r in records if r["id"] in missing]
    if todo:
        budget = resolve_threads(threads)
        chunks = [
            list(part)
            for part in np.array_split(np.array(todo, dtype=object), min(budget, len(todo)))
        ]
        recovered = pl.DataFrame(
            [r for chunk in parallel_map(_read_confidence_chunk, chunks, budget) for r in chunk]
        )
        selected = (
            selected.join(recovered, on="id", how="left", validate="1:1")
            .with_columns(
                pl.when(pl.col("plddt").is_not_null() & pl.col("plddt").is_finite())
                .then(pl.col("plddt"))
                .otherwise(pl.col("coordinate_plddt"))
                .alias("plddt"),
                pl.when(pl.col("id").is_in(missing))
                .then(pl.lit("coordinate_CA"))
                .otherwise(pl.lit("generator_metadata"))
                .alias("plddt_source"),
            )
            .drop("coordinate_plddt")
        )
    else:
        selected = selected.with_columns(pl.lit("generator_metadata").alias("plddt_source"))
    if not selected["plddt"].is_finite().all() or not selected["plddt"].is_between(0, 100).all():
        raise ValueError("Incomplete or invalid generator pLDDT")
    selected.write_csv(output, separator="\t")
    output.with_suffix(".sources.json").write_text(
        json.dumps(
            dict(
                inputs={str(p): digest(p) for p in [sources, *metadata]},
                models=selected.height,
                recovered=len(todo),
                interpretation="Generator pLDDT only; Cα residue mean for missing metadata; not experimental B-factors",
            ),
            indent=2,
        )
        + "\n"
    )


def confidence_gate(rows, gate):
    """Observed confidence, optionally averaged over the declared unlabelled cohort."""
    if gate not in GATES:
        raise ValueError("Unknown confidence gate")
    if gate == "constant":
        return np.ones(rows.height)
    column = "iptm" if gate.endswith("iptm") else "plddt"
    value = pl.col(column) / (100 if column == "plddt" else 1)
    if gate.startswith("cohort-"):
        value = value.mean().over(COHORT)
    q = rows.select(value.alias("q"))["q"].to_numpy()
    if not np.isfinite(q).all() or np.any((q < 0) | (q > 1)):
        raise ValueError("Finite confidence in [0,1] required")
    return q


def gated_score(x, q, weights):
    """S=I+q*(b*(B-I)+m*(M-I)); b,m>=0 and b+m<=1."""
    x, q = np.asarray(x), np.asarray(q)
    if q.shape != (len(x),) or not np.isfinite(q).all() or np.any((q < 0) | (q > 1)):
        raise ValueError("One bounded confidence per model required")
    transformed = x[:, 0, None] + q[:, None] * (x - x[:, 0, None])
    return linear_score(transformed, weights)


def odds_gate(q, alpha):
    """One log-odds offset, expressed by bounded alpha; no fitted slope/power."""
    q, alpha = np.asarray(q), np.asarray(alpha)
    if not np.isfinite(alpha).all() or np.any((alpha < 0) | (alpha > 1)):
        raise ValueError("Bounded odds parameter required")
    denominator = alpha * q + (1 - alpha) * (1 - q)
    return np.divide(
        alpha * q,
        denominator,
        out=np.zeros(np.broadcast_shapes(q.shape, alpha.shape)),
        where=denominator > 0,
    )


def odds_score(x, q, parameters):
    """S=(1-g)*I+g*((1-r)*B+r*M), g=alpha*q/(alpha*q+(1-alpha)*(1-q))."""
    x, q, p = np.asarray(x), np.asarray(q), np.asarray(parameters)
    scalar = p.ndim == 1
    p = np.atleast_2d(p)
    if p.shape[1] != 2 or not np.isfinite(p).all() or np.any((p < 0) | (p > 1)):
        raise ValueError("Bounded alpha and motif split required")
    linear_score(x, [1, 0, 0])
    if q.shape != (len(x),) or not np.isfinite(q).all() or np.any((q < 0) | (q > 1)):
        raise ValueError("One bounded confidence per model required")
    g = odds_gate(q[:, None], p[:, 0][None, :])
    r = p[:, 1][None, :]
    s = (1 - g) * x[:, 0, None] + g * ((1 - r) * x[:, 1, None] + r * x[:, 2, None])
    return s[:, 0] if scalar else s


def complementary_score(x, q, parameters):
    """Opposing expert weights: wB=b*q, wM=m*(1-q), wI=1-wB-wM."""
    x, q, p = np.asarray(x), np.asarray(q), np.asarray(parameters)
    scalar = p.ndim == 1
    p = np.atleast_2d(p)
    linear_score(x, [1, 0, 0])
    if q.shape != (len(x),) or not np.isfinite(q).all() or np.any((q < 0) | (q > 1)):
        raise ValueError("One bounded confidence per model required")
    if p.shape[1] != 2 or not np.isfinite(p).all() or np.any((p < 0) | (p > 1)):
        raise ValueError("Bounded expert proportions required")
    b = q[:, None] * p[:, 0][None, :]
    m = (1 - q[:, None]) * p[:, 1][None, :]
    s = (1 - b - m) * x[:, 0, None] + b * x[:, 1, None] + m * x[:, 2, None]
    return s[:, 0] if scalar else s


def _gate_chunk(task):
    x, q, y, cohort, w, odds, complementary = task
    scores = (
        complementary_score(x, q, w)
        if complementary
        else odds_score(x, q, w)
        if odds
        else gated_score(x, q, w)
    )
    return np.array([ap_columns(y[cohort == c], scores[cohort == c]) for c in np.unique(cohort)]).T


def command(
    features: Path = typer.Option(..., exists=True),
    metadata: list[Path] = typer.Option(..., exists=True),
    gate: str = typer.Option(...),
    output: Path = typer.Option(..., "--output", "-o"),
    threads: int = typer.Option(0, "--threads", "-t"),
    motifs: Path | None = typer.Option(None, exists=True),
    odds: bool = False,
    complementary: bool = False,
    rank_confidence: bool = False,
):
    """Fit two global gate proportions; preserve NLV/TTD/GIL floors from supplied scores."""
    output.mkdir(exist_ok=False)
    start = time.monotonic()
    if odds and complementary:
        raise ValueError("Choose one gate family")
    budget = resolve_threads(threads)
    rows = pl.read_parquet(features).sort("id")
    if motifs is not None:
        replacement = pl.read_csv(motifs, separator="\t").select("id", "enrichment")
        if set(replacement["id"]) != set(rows["id"]):
            raise ValueError("Replacement motif membership differs")
        rows = rows.rename({"enrichment": "previous_enrichment"}).join(
            replacement, on="id", validate="1:1", maintain_order="left"
        )
    meta = pl.concat(
        [
            pl.read_csv(p, separator="\t", infer_schema_length=10000).select("id", "iptm", "plddt")
            for p in metadata
        ]
    )
    meta = meta.rename({"iptm": "metadata_iptm"})
    rows = rows.join(meta, on="id", how="left", validate="1:1", maintain_order="left")
    old = rows["iptm"].to_numpy()
    new = rows["metadata_iptm"].to_numpy()
    valid = np.isfinite(old) & np.isfinite(new)
    np.testing.assert_allclose(old[valid], new[valid], atol=1e-12, rtol=0)
    labelled = rows.filter((pl.col("dataset") == "tcrvdb") & pl.col("y").is_in([0, 1]))
    _, mean = _impute(labelled.select(FEATURES).to_numpy(), np.arange(labelled.height))
    x = rows.select(FEATURES).to_numpy()
    missing = ~np.isfinite(x)
    x = np.where(missing, mean, x)
    rows = rows.with_columns(pl.Series("iptm", x[:, 0]))
    # Missing pLDDT remains explicit; impute from the same training cohort only.
    plddt_mean = labelled["plddt"].mean()
    rows = rows.with_columns(
        pl.col("plddt").is_null().alias("plddt_missing"), pl.col("plddt").fill_null(plddt_mean)
    )
    q = confidence_gate(rows, gate)
    cdf_reference = None
    if rank_confidence:
        if gate == "constant":
            raise ValueError("Constant gate has no confidence distribution")
        unlabelled = rows.with_columns(pl.Series("_q", q)).filter(pl.col("dataset") == "tcrvdb")
        if gate.startswith("cohort-"):
            cdf_reference = unlabelled.group_by(COHORT).agg(pl.col("_q").first())["_q"].to_numpy()
        else:
            cdf_reference = unlabelled["_q"].to_numpy()
        cdf_reference = np.sort(cdf_reference)
        q = np.searchsorted(cdf_reference, q, side="right") / len(cdf_reference)
    idx = np.flatnonzero(
        ((rows["dataset"] == "tcrvdb") & rows["y"].is_in([0, 1]).fill_null(False)).to_numpy()
    )
    train = rows[idx]
    y = train["y"].to_numpy()
    cohort = train["epitope"].to_numpy()
    names = np.unique(cohort)
    w = (
        np.array([(a / 100, r / 100) for a in range(101) for r in range(101)])
        if odds or complementary
        else simplex_weights()
    )
    tasks = [
        (x[idx], q[idx], y, cohort, part, odds, complementary)
        for part in np.array_split(w, min(budget, len(w)))
        if len(part)
    ]
    ap = np.concatenate(parallel_map(_gate_chunk, tasks, budget))
    macro = ap.mean(1)
    from sklearn.metrics import average_precision_score

    floors = {
        c: float(
            average_precision_score(
                train.filter(pl.col("epitope") == c)["y"],
                train.filter(pl.col("epitope") == c)["mixture_fullfit"],
            )
        )
        for c in ["NLVPMVATV", "TTDPSFLGRY", "GILGFVFTL"]
    }
    eligible = np.ones(len(w), dtype=bool)
    for c, floor in floors.items():
        eligible &= ap[:, np.flatnonzero(names == c)[0]] >= floor - 1e-12
    gil = train.filter(pl.col("epitope") == "GILGFVFTL")
    gil_target = float(average_precision_score(gil["y"], gil["prepared_scored"]))
    recovered = eligible & (ap[:, np.flatnonzero(names == "GILGFVFTL")[0]] >= gil_target - 1e-12)
    feasible = bool(eligible.any())
    selection = recovered if recovered.any() else eligible
    best = (
        int(np.argmax(np.where(selection, macro, -np.inf))) if feasible else int(np.argmax(macro))
    )
    fitted = w[best]
    score = (
        complementary_score(x, q, fitted)
        if complementary
        else odds_score(x, q, fitted)
        if odds
        else gated_score(x, q, fitted)
    )
    if complementary:
        wb = q * fitted[0]
        wm = (1 - q) * fitted[1]
        g = wb + wm
    elif odds:
        g = odds_gate(q, fitted[0])
        wb = g * (1 - fitted[1])
        wm = g * fitted[1]
    else:
        wb = q * fitted[1]
        wm = q * fitted[2]
        g = wb + wm
    rows = rows.with_columns(
        pl.col("mixture_fullfit").alias("previous_score"),
        pl.Series("mixture_fullfit", score),
        pl.lit(None, dtype=pl.Float64).alias("mixture_oof"),
        pl.Series("gate_q", q),
        pl.Series("w_I", 1 - g),
        pl.Series("w_B", wb),
        pl.Series("w_M", wm),
        *[pl.Series(f + "_scored", x[:, j]) for j, f in enumerate(FEATURES)],
        *[pl.Series(f + "_missing", missing[:, j]) for j, f in enumerate(FEATURES)],
    )
    rows.write_parquet(output / "predictions.parquet")
    _metrics(rows).write_csv(output / "metrics.tsv", separator="\t")
    rows.group_by(COHORT).agg(
        pl.len().alias("n"), pl.col("iptm", "plddt", "gate_q", "w_I", "w_B", "w_M").mean()
    ).sort(COHORT).write_csv(output / "cohort_weights.tsv", separator="\t")
    pl.DataFrame(
        w, schema=["b", "m"] if complementary else ["alpha", "r"] if odds else ["unused", "b", "m"]
    ).with_columns(
        pl.Series("macro_AP", macro),
        pl.Series("preserves_priority", eligible),
        pl.Series("recovers_GIL", recovered),
        *[pl.Series(c, ap[:, j]) for j, c in enumerate(names)],
    ).write_csv(output / "objectives.tsv", separator="\t")
    resource = dict(
        schema="tcren-confidence-motif-1",
        gate=gate,
        parameters=dict(b=float(fitted[0]), m=float(fitted[1]))
        if complementary
        else dict(alpha=float(fitted[0]), r=float(fitted[1]))
        if odds
        else dict(b=float(fitted[1]), m=float(fitted[2])),
        imputation=mean.tolist(),
        plddt_imputation=plddt_mean,
        formula="S=(1-b*q-m*(1-q))*ipTM+b*q*Betti+m*(1-q)*enrichment"
        if complementary
        else "S=(1-g)*ipTM+g*((1-r)*Betti+r*enrichment); g=alpha*q/(alpha*q+(1-alpha)*(1-q))"
        if odds
        else "S=(1-q*(b+m))*ipTM+q*b*Betti+q*m*enrichment",
        training_n=len(idx),
        training_macro_AP=float(macro[best]),
        floors=floors,
        preserves_priority=feasible,
        feasible_candidates=int(eligible.sum()),
        GIL_target=gil_target,
        recovered_candidates=int(recovered.sum()),
        recovers_GIL=bool(recovered[best]),
        cohort_keys=COHORT,
        rank_confidence=rank_confidence,
        confidence_cdf=cdf_reference.tolist() if cdf_reference is not None else None,
        interpretation="Full-fit TCRvdb optimization; cohort means include all supplied unlabelled feature rows; VDJdb labels not used",
        inputs={str(p): digest(p) for p in [features, *metadata, motifs] if p is not None},
        seconds=time.monotonic() - start,
    )
    (output / "weights.json").write_text(json.dumps(resource, indent=2) + "\n")


def residue_audit_command(
    residues: Path = typer.Option(..., exists=True),
    labels: Path = typer.Option(..., exists=True),
    output: Path = typer.Option(..., "--output", "-o"),
    annotations: list[Path] = typer.Option([], exists=True),
):
    """Describe named motif residues without fitting peptide-specific predictors or tests."""
    output.mkdir(exist_ok=False)
    detail = pl.read_parquet(residues).sort("id", "chain", "region", "residue_index")
    if not (detail["aa"].str.len_chars() == 1).all():
        raise ValueError("Residue audit expects single-residue features")
    if isinstance(annotations, list) and annotations:
        import gzip

        records = []
        for path in annotations:
            opener = gzip.open if path.suffix == ".gz" else open
            with opener(path, "rt") as handle:
                payload = json.load(handle)
            for model, organisms in payload.items():
                for chain, record in organisms.get("human", {}).items():
                    if record.get("locus") not in ["TRA", "TRB"]:
                        continue
                    end = record.get("v_sequence_end")
                    start = record.get("v_sequence_start")
                    records.append(
                        dict(
                            id=model,
                            chain_id=chain,
                            v_gene=record.get("v_call"),
                            v_start=int(float(start)) if start not in [None, ""] else None,
                            v_end=int(float(end)) if end not in [None, ""] else None,
                        )
                    )
        mapping = pl.DataFrame(records, schema_overrides={"v_start": pl.Int64, "v_end": pl.Int64})
        detail = detail.join(
            mapping, on=["id", "chain_id"], how="left", validate="m:1", maintain_order="left"
        )
        detail = detail.with_columns(
            ((pl.col("residue_index") + 1).is_between(pl.col("v_start"), pl.col("v_end")))
            .fill_null(False)
            .alias("V_aligned")
        )
    else:
        detail = detail.with_columns(pl.lit(False).alias("V_aligned"))
    detail = detail.with_columns(
        (pl.col("V_aligned") | pl.col("region").is_in(["CDR1", "CDR2"])).alias("V_aligned")
    )
    if "evidence_bits" in detail.columns:
        detail = detail.rename({"bits": "negative_tail_bits"}).with_columns(
            pl.col("evidence_bits").alias("bits")
        )
    scopes = ["id", "chain", "region"]
    detail = detail.with_columns(
        pl.col("aa").cum_count().over(scopes).alias("loop_position"),
        pl.col("aa").str.join("").over(scopes).alias("loop_sequence"),
        (
            (pl.col("aa") == "R") & (pl.col("aa").shift(-1).over(scopes) == "S")
            | (pl.col("aa") == "S") & (pl.col("aa").shift(1).over(scopes) == "R")
        )
        .fill_null(False)
        .alias("_rs"),
    )
    detail = detail.with_columns(
        (pl.col("_rs") | pl.col("_rs").shift(1).over(scopes) | pl.col("_rs").shift(-1).over(scopes))
        .fill_null(False)
        .alias("rs_neighborhood")
    )
    meta = pl.read_parquet(labels).select("id", "dataset", "stratum", "epitope", "y")
    if set(meta["id"]) != set(detail["id"]):
        raise ValueError("Residue/label membership differs")
    detail = detail.join(meta, on="id", validate="m:1").with_columns(
        (2 * np.log(2) * pl.col("bits").fill_null(0)).alias("fisher_contribution"),
        (
            (pl.col("chain") == "PEPTIDE")
            & pl.col("loop_sequence").str.contains("GIGI")
            & pl.col("loop_position").is_between(
                pl.col("loop_sequence").str.find("GIGI") + 1,
                pl.col("loop_sequence").str.find("GIGI") + 4,
            )
        )
        .fill_null(False)
        .alias("peptide_GIGI"),
        (
            (pl.col("epitope") == "TTDPSFLGRY")
            & (pl.col("chain") != "PEPTIDE")
            & (pl.col("aa") == "R")
        ).alias("TTD_R"),
        (
            (pl.col("epitope") == "GILGFVFTL")
            & (pl.col("chain") == "TRB")
            & (pl.col("region") == "CDR3")
            & pl.col("rs_neighborhood")
        ).alias("GIL_RS_flanks"),
        ((pl.col("epitope") == "GLCTLVAML") & pl.col("region").is_in(["CDR1", "CDR2"])).alias(
            "GLC_V_CDR12"
        ),
        (
            (pl.col("epitope") == "GLCTLVAML") & (pl.col("region") == "CDR3") & pl.col("V_aligned")
        ).alias("GLC_V_aligned_CDR3"),
    )
    detail = detail.with_columns(
        pl.when(pl.col("fisher_contribution").sum().over("id") > 0)
        .then(pl.col("fisher_contribution") / pl.col("fisher_contribution").sum().over("id"))
        .otherwise(0)
        .alias("fisher_share")
    )
    detail.drop("_rs").write_parquet(output / "residues.parquet")
    summary = detail.group_by("dataset", "epitope", "y", "chain", "region", "aa").agg(
        pl.col("id").n_unique().alias("models"),
        pl.len().alias("residues"),
        (pl.col("pairs") > 0).sum().alias("contacting_residues"),
        pl.col("cycles", "bits", "fisher_share").mean(),
    )
    summary.sort("dataset", "epitope", "y", "chain", "region", "aa").write_csv(
        output / "residue_types.tsv", separator="\t"
    )
    named = []
    for motif in ["TTD_R", "GIL_RS_flanks", "GLC_V_CDR12", "GLC_V_aligned_CDR3", "peptide_GIGI"]:
        named.append(
            detail.filter(motif)
            .group_by("dataset", "epitope", "y")
            .agg(
                pl.col("id").n_unique().alias("models"),
                pl.len().alias("residues"),
                (pl.col("pairs") > 0).sum().alias("contacting_residues"),
                pl.col("cycles", "bits", "fisher_share").mean(),
            )
            .with_columns(pl.lit(motif).alias("motif"))
        )
    pl.concat(named).sort("motif", "dataset", "epitope", "y").write_csv(
        output / "named_motifs.tsv", separator="\t"
    )


def pair_counts_command(
    residues: Path = typer.Option(..., exists=True),
    contacts: list[Path] = typer.Option(..., exists=True),
    output: Path = typer.Option(..., "--output", "-o"),
):
    """Two adjacent covalently connected residues as one focal contact subgraph."""
    from .residue_surprise import KEY, validate_counts
    import shutil

    if output.exists():
        raise ValueError("Output exists")
    rows = pl.read_parquet(residues).sort("id", "chain_id", "region", "residue_index")
    validate_counts(rows)
    if not rows["aa"].str.len_chars().eq(1).all():
        raise ValueError("Input must contain single residues")
    scopes = ["id", "chain_id", "region"]
    patch = (
        rows.with_columns(
            *[
                pl.col(c).shift(-1).over(scopes).alias("_next_" + c)
                for c in ["residue_index", "aa", "available", "complete"]
            ]
        )
        .filter(pl.col("_next_residue_index") == pl.col("residue_index") + 1)
        .with_columns(
            (pl.col("aa") + pl.col("_next_aa")).alias("aa"),
            (pl.col("available") + pl.col("_next_available")).alias("available"),
            (pl.col("complete") & pl.col("_next_complete")).alias("complete"),
        )
    )
    patch = patch.rename({"residue_index": "patch_start"})
    members = pl.concat(
        [
            patch.select("id", "chain_id", "patch_start", pl.col(c).alias("residue_index"))
            for c in ["patch_start", "_next_residue_index"]
        ]
    )
    c = pl.concat([pl.read_parquet(p) for p in contacts]).rename({"pdb.id": "id"})
    if set(c["id"]) - set(rows["id"]):
        raise ValueError("Unexpected contact models")
    directed = pl.concat(
        [
            c.select(
                "id",
                pl.col("chain.id." + side).alias("chain_id"),
                pl.col("residue.index." + side).alias("residue_index"),
                pl.col("atom." + side).alias("atom"),
                pl.col("chain.id." + other).alias("partner_chain"),
                pl.col("residue.index." + other).alias("partner_residue"),
            )
            for side, other in [("from", "to"), ("to", "from")]
        ]
    )
    edges = members.join(directed, on=KEY, how="inner")
    counts = edges.group_by("id", "chain_id", "patch_start").agg(
        pl.len().alias("pairs"),
        pl.struct("residue_index", "atom").n_unique().alias("occupied"),
        pl.struct("partner_chain", "partner_residue").n_unique().alias("partners"),
    )
    result = (
        patch.select(
            "id",
            "chain_id",
            "patch_start",
            "chain",
            "region",
            "aa",
            "available",
            "complete",
            "atom_scope",
        )
        .join(counts, on=["id", "chain_id", "patch_start"], how="left", validate="1:1")
        .with_columns(pl.col("pairs", "occupied", "partners").fill_null(0))
        .with_columns(
            (pl.col("pairs").cast(pl.Int64) - pl.col("partners").cast(pl.Int64)).alias("cycles")
        )
        .rename({"patch_start": "residue_index"})
        .sort(KEY)
    )
    validate_counts(result)
    result.write_parquet(output)
    identities = residues.with_suffix(".identities.tsv")
    if identities.exists():
        shutil.copy2(identities, output.with_suffix(".identities.tsv"))
    output.with_suffix(".sources.json").write_text(
        json.dumps(
            dict(
                inputs={str(p): digest(p) for p in [residues, *contacts]},
                models=result["id"].n_unique(),
                patches=result.height,
                graph="Adjacent observed residues within one CDR/peptide form a covalent focal patch; added cycles=distinct atom contacts-distinct partner residues; no peptide-specific motif bonus",
            ),
            indent=2,
        )
        + "\n"
    )


def score_command(
    features: Path = typer.Option(..., exists=True),
    resource: Path = typer.Option(..., exists=True),
    output: Path = typer.Option(..., "--output", "-o"),
):
    """Apply a frozen confidence gate and confidence CDF, without labels or fitting."""
    if output.exists():
        raise ValueError("Output exists")
    model = json.loads(resource.read_text())
    if model.get("schema") != "tcren-confidence-motif-1":
        raise ValueError("Unknown confidence resource")
    rows = (
        pl.read_parquet(features)
        if features.suffix == ".parquet"
        else pl.read_csv(features, separator="\t", infer_schema_length=10000)
    )
    apply_confidence_score(rows, model).write_csv(output, separator="\t")


def apply_confidence_score(rows, model):
    """Apply a frozen model to declared cohorts, preserving row order and missing flags."""
    if model.get("schema") != "tcren-confidence-motif-1":
        raise ValueError("Unknown confidence resource")
    x = rows.select(pl.col(FEATURES).cast(pl.Float64)).to_numpy()
    if np.isinf(x).any():
        raise ValueError("Infinite score components are not allowed")
    if "plddt" in rows.columns:
        rows = rows.with_columns(pl.col("plddt").cast(pl.Float64))
        observed_plddt = rows["plddt"].drop_nulls().to_numpy()
        if np.isinf(observed_plddt).any() or np.any((observed_plddt < 0) | (observed_plddt > 100)):
            raise ValueError("Generator plddt must lie in [0,100]")
    missing = ~np.isfinite(x)
    x = np.where(missing, model["imputation"], x)
    rows = rows.with_columns(pl.Series("iptm", x[:, 0]))
    if "plddt" in rows.columns:
        rows = rows.with_columns(
            pl.col("plddt").fill_nan(None).fill_null(model["plddt_imputation"])
        )
    q = confidence_gate(rows, model["gate"])
    if model.get("rank_confidence"):
        # Numerical ties at cohort means must survive row order and parallel
        # summation differences. Twelve decimals is comparison precision only.
        ref = np.round(np.asarray(model["confidence_cdf"]), 12)
        if (
            ref.ndim != 1
            or not len(ref)
            or not np.isfinite(ref).all()
            or np.any(np.diff(ref) < 0)
            or np.any((ref < 0) | (ref > 1))
        ):
            raise ValueError("A sorted bounded confidence CDF reference is required")
        q = np.searchsorted(ref, np.round(q, 12), side="right") / len(ref)
    p = model.get("parameters", model)
    if model["formula"].startswith("S=(1-b*q-m*(1-q))"):
        score = complementary_score(x, q, [p["b"], p["m"]])
        wb = p["b"] * q
        wm = p["m"] * (1 - q)
    elif "alpha" in p:
        score = odds_score(x, q, [p["alpha"], p["r"]])
        g = odds_gate(q, p["alpha"])
        wb = g * (1 - p["r"])
        wm = g * p["r"]
    else:
        score = gated_score(x, q, [1 - p["b"] - p["m"], p["b"], p["m"]])
        wb = p["b"] * q
        wm = p["m"] * q
    return rows.with_columns(
        *[pl.Series("missing_" + name, missing[:, i]) for i, name in enumerate(FEATURES)],
        pl.Series("score", score).set(pl.Series(missing.all(axis=1)), None),
        pl.Series("gate_q", q),
        pl.Series("w_I", 1 - wb - wm),
        pl.Series("w_B", wb),
        pl.Series("w_M", wm),
    )


def report_command(
    experiments: list[Path] = typer.Option(..., exists=True),
    output: Path = typer.Option(..., "--output", "-o"),
):
    """Assemble flat, complete round comparisons from library-produced metric tables."""
    output.mkdir(exist_ok=False)
    tables = []
    for path in experiments:
        table = pl.read_csv(path / "metrics.tsv", separator="\t").with_columns(
            pl.lit(path.parent.name + "/" + path.name).alias("experiment")
        )
        tables.append(table)
    all_metrics = pl.concat(tables)
    all_metrics.write_csv(output / "cohort_metrics.tsv", separator="\t")
    all_metrics.group_by("experiment", "population", "dataset", "stratum", "arm").agg(
        pl.len().alias("cohorts"), pl.col("n").sum(), pl.col("AP", "AUROC").mean()
    ).sort("experiment", "population", "dataset", "stratum", "arm").write_csv(
        output / "summary.tsv", separator="\t"
    )
