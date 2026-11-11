"""All-loop reciprocal contact multiplicity and empirical mock-tail information.

For internally connected residues, adding m focal-to-neighbour atom contacts to
q distinct neighbours increases graph beta1 by m-q. Intrinsic residue cycles
cancel. No contact between two neighbours or inter-residue covalent bond enters
this local graph. Reference tails are descriptive, not binding probabilities.
"""

from pathlib import Path
import json

import numpy as np
import polars as pl
import typer

from .saturation import EXPECTED
from .profile import BACKBONE, digest

KEY = ["id", "chain_id", "residue_index"]
CELL = ["chain", "region", "aa", "available"]
STATISTICS = ("occupied", "pairs", "cycles")


def residue_counts(atoms, contacts, atom_scope="all"):
    """Zero-inclusive observed heavy-atom counts, all six CDRs and peptide.

    Complete means canonical heavy atoms are present; no coordinate inference.
    The graph contraction is used only when focal and partner residues are
    complete. Missing atoms remain explicit rather than imputing contact edges.
    """
    if atom_scope not in ["all", "sidechain"]:
        raise ValueError("Atom scope must be all or sidechain")
    names = {
        "pdb.id": "id",
        "chain.id": "chain_id",
        "residue.index": "residue_index",
        "chain.type": "chain",
        "region.type": "region",
        "residue.aa": "aa",
    }
    inv = (
        atoms.select(*names, "atom.name")
        .rename(names)
        .filter(
            (
                pl.col("chain").is_in(["TRA", "TRB"])
                & pl.col("region").is_in(["CDR1", "CDR2", "CDR3"])
            )
            | (pl.col("chain") == "PEPTIDE")
        )
        .with_columns(
            pl.when(pl.col("chain") == "PEPTIDE")
            .then(pl.lit("PEPTIDE"))
            .otherwise(pl.col("region"))
            .alias("region")
        )
    )
    expected = pl.DataFrame(
        {
            "aa": list(EXPECTED),
            "expected": [sorted(set(x) | {"N", "CA", "C", "O"}) for x in EXPECTED.values()],
        }
    )
    inv = (
        inv.group_by(*KEY, "chain", "region", "aa")
        .agg(pl.col("atom.name").unique().alias("observed"))
        .join(expected, on="aa", how="left", validate="m:1")
        .with_columns(
            pl.col("observed").list.len().alias("available"),
            (
                pl.col("expected").is_not_null()
                & (pl.col("expected").list.set_difference(pl.col("observed")).list.len() == 0)
            )
            .fill_null(False)
            .alias("complete"),
        )
    )
    c = contacts.rename({"pdb.id": "id"})
    physical = [
        "id",
        "chain.id.from",
        "residue.index.from",
        "atom.from",
        "chain.id.to",
        "residue.index.to",
        "atom.to",
    ]
    if c.select(physical).is_duplicated().any():
        raise ValueError("Duplicate physical atom contacts")
    if (
        c["dist"].null_count()
        or not c["dist"].is_finite().all()
        or not c["dist"].is_between(0, 5).all()
    ):
        raise ValueError("Expected finite contacts within 5 A")
    ends = []
    for side, other in [("from", "to"), ("to", "from")]:
        ends.append(
            c.select(
                "id",
                pl.col("chain.id." + side).alias("chain_id"),
                pl.col("residue.index." + side).alias("residue_index"),
                pl.col("atom." + side).alias("atom"),
                pl.col("chain.id." + other).alias("partner_chain"),
                pl.col("residue.index." + other).alias("partner_residue"),
            )
        )
    directed = pl.concat(ends).join(
        inv.select(*KEY, pl.col("complete").alias("partner_complete")),
        left_on=["id", "partner_chain", "partner_residue"],
        right_on=KEY,
        how="left",
        validate="m:1",
    )
    if directed["partner_complete"].null_count():
        raise ValueError("Contacts must connect inventoried CDR and peptide residues")
    if atom_scope == "sidechain":
        directed = directed.filter(~pl.col("atom").is_in(BACKBONE))
        inv = inv.with_columns(
            pl.col("observed").list.set_difference(BACKBONE).list.len().alias("available")
        )
    counts = directed.group_by(KEY).agg(
        pl.len().alias("pairs"),
        pl.col("atom").n_unique().alias("occupied"),
        pl.struct("partner_chain", "partner_residue").n_unique().alias("partners"),
        pl.col("partner_complete").all().alias("partners_complete"),
    )
    return (
        inv.join(counts, on=KEY, how="left", validate="1:1")
        .with_columns(
            pl.col("pairs", "occupied", "partners").fill_null(0),
            pl.col("partners_complete").fill_null(True),
        )
        .with_columns(
            (pl.col("pairs").cast(pl.Int64) - pl.col("partners").cast(pl.Int64)).alias("cycles"),
            (pl.col("complete") & pl.col("partners_complete")).alias("complete"),
        )
        .with_columns(pl.lit(atom_scope).alias("atom_scope"))
        .drop("observed", "expected", "partners_complete")
        .sort(KEY)
    )


def validate_counts(rows):
    for name in STATISTICS:
        x = rows[name].to_numpy()
        if not np.isfinite(x).all() or np.any((x < 0) | (x != np.floor(x))):
            raise ValueError("Contact statistics must be nonnegative integers")


def fit_reference(rows):
    """Freeze actual integer residue histograms; no query-cohort fitting."""
    validate_counts(rows)
    usable = rows.filter("complete")
    if usable.is_empty():
        raise ValueError("No complete mock residues")
    if "weight" not in usable.columns:
        usable = usable.with_columns(pl.lit(1.0).alias("weight"))
    if not usable["weight"].is_finite().all() or (usable["weight"] <= 0).any():
        raise ValueError("Reference weights must be finite and positive")
    tables = {}
    for statistic in STATISTICS:
        tables[statistic] = (
            usable.group_by(*CELL, statistic)
            .agg(
                pl.len().alias("count"),
                pl.col("id").n_unique().alias("models"),
                pl.col("weight").sum().alias("mass"),
            )
            .sort(*CELL, statistic)
            .to_dicts()
        )
    scopes = rows["atom_scope"].unique().to_list() if "atom_scope" in rows.columns else ["all"]
    if len(scopes) != 1:
        raise ValueError("One atom scope is required per reference")
    return dict(
        schema="tcren-residue-surprise-1", atom_scope=scopes[0], cutoff_A=5, statistics=tables
    )


def score_residues(rows, resource, statistic, positive_reference=None, joint_residue=False):
    """Inclusive tail p=(1+#mock>=x)/(N+1); unsupported cells remain missing."""
    if resource.get("schema") != "tcren-residue-surprise-1" or statistic not in STATISTICS:
        raise ValueError("Unknown residue reference/statistic")
    validate_counts(rows)
    scopes = rows["atom_scope"].unique().to_list() if "atom_scope" in rows.columns else ["all"]
    if scopes != [resource.get("atom_scope", "all")]:
        raise ValueError("Query/reference atom scopes differ")
    rows = rows.with_columns(pl.col("available").cast(pl.Int64))
    hist = pl.DataFrame(resource["statistics"][statistic]).sort(*CELL, statistic)
    if "mass" not in hist.columns:
        hist = hist.with_columns(pl.col("count").cast(pl.Float64).alias("mass"))
    normalization = ["chain", "region"] if joint_residue else CELL
    hist = hist.with_columns(
        pl.col("count").sum().over(normalization).alias("reference_residues"),
        pl.col("mass").sum().over(normalization).alias("reference_mass"),
        pl.col("mass").reverse().cum_sum().reverse().over(CELL).alias("tail_count"),
    )
    # Sparse integer tails: a forward as-of match gives count at the first
    # observed value >= query; queries above the reference maximum have tail 0.
    total = hist.select(*normalization, "reference_residues", "reference_mass").unique()
    q = rows.with_row_index("_order").join(total, on=normalization, how="left", validate="m:1")
    q = q.with_columns(pl.col(statistic).cast(pl.Int64))
    hist = hist.with_columns(pl.col(statistic).cast(pl.Int64))
    q = (
        q.sort(*CELL, statistic)
        .join_asof(
            hist.select(*CELL, statistic, "tail_count"),
            on=statistic,
            by=CELL,
            strategy="forward",
            check_sortedness=False,
        )
        .sort("_order")
        .drop("_order")
    )
    result = q.with_columns(
        pl.when(pl.col("complete") & pl.col("reference_residues").is_not_null())
        .then(
            pl.when(pl.col(statistic) == 0)
            .then(1.0)
            .otherwise(
                ((1 + pl.col("tail_count").fill_null(0)) / (1 + pl.col("reference_mass"))).clip(
                    0, 1
                )
            )
        )
        .otherwise(None)
        .alias("tail_probability")
    ).with_columns(
        (-pl.col("tail_probability").log(2)).alias("bits"),
        (1 - pl.col("tail_probability")).alias("percentile"),
    )
    if positive_reference is not None:
        positive = score_residues(rows, positive_reference, statistic, joint_residue=joint_residue)
        result = result.with_columns(
            positive["tail_probability"].alias("positive_tail_probability"),
            positive["reference_residues"].alias("positive_reference_residues"),
            positive["reference_mass"].alias("positive_reference_mass"),
        )
        result = result.with_columns(
            (pl.col("positive_tail_probability") / pl.col("tail_probability"))
            .log(2)
            .alias("log_tail_ratio")
        )
        result = result.with_columns(
            pl.col("log_tail_ratio").clip(lower_bound=0).alias("enrichment")
        )
    return result


def model_identities(atoms):
    """Full observed receptor/epitope sequences, without clonotype collapse."""
    return (
        atoms.select("pdb.id", "chain.type", "chain.id", "residue.index", "residue.aa")
        .filter(pl.col("chain.type").is_in(["TRA", "TRB", "PEPTIDE"]))
        .unique()
        .sort("pdb.id", "chain.id", "residue.index")
        .group_by("pdb.id", "chain.type")
        .agg(pl.col("residue.aa").str.join("").alias("sequence"))
        .pivot(on="chain.type", index="pdb.id", values="sequence")
        .rename({"pdb.id": "id"})
    )


def counts_command(
    atoms: Path = typer.Option(..., exists=True),
    contacts: Path = typer.Option(..., exists=True),
    output: Path = typer.Option(..., "--output", "-o"),
    atom_scope: str = typer.Option("all"),
):
    """Reduce CLI atom inventory and CDR:peptide contacts into reciprocal residue tables."""
    if output.exists():
        raise ValueError("Output already exists")
    inventory = pl.read_parquet(atoms)
    identity = model_identities(inventory)
    accounting = (
        pl.read_csv(Path(str(atoms) + ".accounting.tsv"), separator="\t")
        .filter(pl.col("status") == "ok")
        .select(pl.col("pdb.id").alias("id"), "sha256")
    )
    identity = identity.join(accounting, on="id", how="left", validate="1:1")
    if identity["sha256"].null_count():
        raise ValueError("CLI inventory accounting must supply every coordinate hash")
    residue_counts(inventory, pl.read_parquet(contacts), atom_scope).write_parquet(output)
    identity.sort("id").write_csv(output.with_suffix(".identities.tsv"), separator="\t")


def fit_command(
    features: Path = typer.Option(..., exists=True),
    query_features: Path = typer.Option(..., exists=True),
    metadata: Path = typer.Option(..., exists=True),
    output: Path = typer.Option(..., "--output", "-o"),
    reference_label: int = typer.Option(0),
    balance_epitopes: bool = typer.Option(False),
):
    """Freeze empirical tails from an explicitly disjoint mock-residue inventory."""
    if output.exists():
        raise ValueError("Output already exists")
    mock, query = pl.read_parquet(features), pl.read_parquet(query_features)
    meta = pl.read_csv(metadata, separator="\t")
    if (
        reference_label not in [0, 1]
        or meta["y"].null_count()
        or meta["pdb.id"].n_unique() != meta.height
        or not (meta["y"] == reference_label).all()
        or set(meta["pdb.id"]) != set(mock["id"])
    ):
        raise ValueError("Metadata must match every reference model and its declared label 0 or 1")

    def identities(path, ids):
        table = pl.read_csv(path.with_suffix(".identities.tsv"), separator="\t")
        if table["id"].n_unique() != table.height or set(table["id"]) != set(ids):
            raise ValueError("Full receptor identities must match all feature models")
        if table.select("TRA", "TRB", "PEPTIDE", "sha256").null_count().sum_horizontal().item():
            raise ValueError("Full receptor/peptide sequences and coordinate hashes are required")
        return table

    qs = identities(query_features, query["id"])
    ms = identities(features, mock["id"])
    pairing = ["TRA", "TRB", "PEPTIDE"]
    overlap = ms.join(qs.select(pairing).unique(), on=pairing, how="semi")["id"]
    coordinate_overlap = ms.filter(pl.col("sha256").is_in(qs["sha256"].implode()))["id"]
    excluded = set(overlap) | set(coordinate_overlap) | (set(mock["id"]) & set(query["id"]))
    selected = ms.filter(~pl.col("id").is_in(excluded))
    selected = selected.with_columns(
        (
            pl.lit(selected.height) / selected["PEPTIDE"].n_unique() / pl.len().over("PEPTIDE")
            if balance_epitopes
            else pl.lit(1.0)
        ).alias("weight")
    )
    selected.write_csv(output.with_suffix(".models.tsv"), separator="\t")
    mock = mock.join(selected.select("id", "weight"), on="id", how="inner", validate="m:1")
    resource = fit_reference(mock)
    resource.update(
        reference_label=reference_label,
        models=mock["id"].n_unique(),
        excluded_models=len(excluded),
        downsampled_models=0,
        epitopes=selected["PEPTIDE"].n_unique(),
        query_sha256=digest(query_features),
        metadata_sha256=digest(metadata),
        identities_sha256=digest(features.with_suffix(".identities.tsv")),
        query_identities_sha256=digest(query_features.with_suffix(".identities.tsv")),
        exclusion="Exact full TRA+TRB+peptide pairing, coordinate content or model ID in benchmark only; shared TCRs or epitopes allowed",
        downsampling="none; all eligible source-labelled controls retained",
        weighting="equal total model weight per epitope, mean model weight one"
        if balance_epitopes
        else "one per model; no downweighting",
    )
    resource["input_sha256"] = digest(features)
    output.write_text(json.dumps(resource, indent=2) + "\n")


def score_command(
    features: Path = typer.Option(..., exists=True),
    reference: Path = typer.Option(..., exists=True),
    output: Path = typer.Option(..., "--output", "-o"),
    statistic: str = typer.Option("cycles"),
    positive_reference: Path | None = typer.Option(None, exists=True),
):
    """Export residue probabilities/bits and loop-balanced model surprise tables."""
    output.mkdir(exist_ok=False)
    positive = json.loads(positive_reference.read_text()) if positive_reference else None
    if positive and positive.get("reference_label") != 1:
        raise ValueError("Positive reference must declare label 1")
    detail = score_residues(
        pl.read_parquet(features), json.loads(reference.read_text()), statistic, positive
    )
    detail.write_parquet(output / "residues.parquet")
    raw_measures = ["bits", "percentile"] + (["log_tail_ratio", "enrichment"] if positive else [])
    loops = detail.group_by("id", "chain", "region").agg(
        pl.col(*raw_measures).fill_null(0).mean(),
        *[pl.col(m).count().alias(m + "_supported") for m in raw_measures],
        pl.len().alias("residues"),
        pl.col("bits").count().alias("supported"),
    )
    loops.write_csv(output / "loops.tsv", separator="\t")
    # Unsupported residues abstain with zero evidence in the full denominator;
    # detail stays null and support counts distinguish abstention from p=1.
    loop_counts = loops.filter(pl.col("chain") != "PEPTIDE").group_by("id").len()
    if (loop_counts["len"] != 6).any():
        raise ValueError("Every model requires all six CDR loops")
    sides = loops.with_columns(
        pl.when(pl.col("chain") == "PEPTIDE")
        .then(pl.lit("peptide"))
        .otherwise(pl.lit("tcr"))
        .alias("side")
    )
    measures = raw_measures
    sides = sides.group_by("id", "side").agg(
        *[
            pl.when(pl.col(m).null_count() == 0).then(pl.col(m).mean()).otherwise(None).alias(m)
            for m in measures
        ]
    )
    wide = sides.pivot(on="side", index="id", values=measures)
    for measure in measures:
        a, b = pl.col(measure + "_tcr"), pl.col(measure + "_peptide")
        wide = wide.with_columns(((a + b) / 2).alias(measure + "_mean"))
        if not measure.startswith("log_tail_ratio"):
            wide = wide.with_columns((a * b).sqrt().alias(measure + "_geometric"))
    wide.sort("id").write_csv(output / "scores.tsv", separator="\t")
    (output / "PROVENANCE.json").write_text(
        json.dumps(
            dict(
                statistic=statistic,
                features_sha256=digest(features),
                reference_sha256=digest(reference),
                positive_reference_sha256=digest(positive_reference)
                if positive_reference
                else None,
                interpretation="Descriptive mock-tail evidence; no independent-residue or binding-probability claim",
            ),
            indent=2,
        )
        + "\n"
    )


def evaluate_command(
    baseline: Path = typer.Option(..., exists=True),
    scores: list[Path] = typer.Option(..., exists=True),
    output: Path = typer.Option(..., "--output", "-o"),
    comparators: Path | None = typer.Option(None, exists=True),
    weights: list[float] = typer.Option([1.0], "--weight"),
    enrichment: bool = typer.Option(False),
):
    """Small explicit development comparison; all cohort/stratum deltas retained.

    One formula S=(prepared+w*sqrt(T*P))/(1+w), default w=1. T and P are
    mean residue rarity, or positive cognate/mock log enrichment with
    --enrichment. No arbitrary powers, coefficient grid or default promotion.
    Inspected cohorts remain internal development, not a new held-out test.
    """
    from sklearn.metrics import average_precision_score, roc_auc_score

    if not weights or not np.isfinite(weights).all() or np.any(np.asarray(weights) < 0):
        raise ValueError("Mixture weights must be finite and nonnegative")
    output.mkdir(exist_ok=False)
    rows = pl.read_parquet(baseline).sort("id")
    if rows["id"].n_unique() != rows.height:
        raise ValueError("Duplicate baseline IDs")
    arms = ["iptm", "prepared"]
    if comparators is not None:
        other = pl.read_parquet(comparators)
        names = [c for c in ["R5", "R7_interaction_p4"] if c in other.columns]
        check = rows.join(other.select("id", "Accepted"), on="id", how="inner", validate="1:1")
        np.testing.assert_allclose(check["prepared"], check["Accepted"], atol=1e-12, rtol=0)
        rows = rows.join(
            other.select("id", *names), on="id", how="left", validate="1:1", maintain_order="left"
        )
        arms.extend(names)
    comparator_count = len(arms)
    for path in scores:
        values = pl.read_csv(path, separator="\t")
        prefix = path.parent.name
        column = "enrichment_geometric" if enrichment else "percentile_geometric"
        if column not in values.columns:
            raise ValueError("Requested one-line score component is missing: " + column)
        columns = [column]
        if set(values["id"]) != set(rows["id"]):
            raise ValueError("Score/baseline membership differs")
        values = values.select("id", *columns).rename({c: prefix + "_" + c for c in columns})
        rows = rows.join(values, on="id", validate="1:1", maintain_order="left")
        for col in columns:
            feature = prefix + "_" + col
            for weight in weights:
                arm = feature + "_w" + str(weight)
                rows = rows.with_columns(
                    ((pl.col("prepared") + weight * pl.col(feature)) / (1 + weight)).alias(arm)
                )
                arms.append(arm)
    eligible = rows.filter(pl.col("y").is_in([0, 1]) & pl.all_horizontal(pl.col(*arms).is_finite()))
    records = []
    for population, subset in [
        ("all_labelled", eligible),
        ("observed_baselines", eligible.filter("matched")),
    ]:
        for key, part in subset.group_by("dataset", "stratum", "epitope"):
            y = part["y"].to_numpy()
            if len(set(y)) != 2:
                continue
            ref = average_precision_score(y, part["prepared"])
            for arm in arms:
                records.append(
                    dict(
                        zip(["dataset", "stratum", "epitope"], key),
                        arm=arm,
                        population=population,
                        n=part.height,
                        binders=int(y.sum()),
                        nonbinders=int(len(y) - y.sum()),
                        AP=float(average_precision_score(y, part[arm])),
                        AUROC=float(roc_auc_score(y, part[arm])),
                        delta_AP=float(average_precision_score(y, part[arm]) - ref),
                    )
                )
    metrics = pl.DataFrame(records).sort("dataset", "stratum", "epitope", "arm")
    metrics.write_csv(output / "metrics.tsv", separator="\t")
    rows.write_parquet(output / "predictions.parquet")
    metrics.group_by("population", "dataset", "stratum", "arm").agg(
        pl.len().alias("cohorts"),
        pl.col("n").sum(),
        pl.col("AP", "AUROC", "delta_AP").mean(),
        pl.col("delta_AP").min().alias("worst_delta_AP"),
    ).sort("dataset", "stratum", "arm").write_csv(output / "summary.tsv", separator="\t")
    (output / "PROVENANCE.json").write_text(
        json.dumps(
            dict(
                inputs={
                    str(p): digest(p)
                    for p in [baseline, *scores, *([comparators] if comparators else [])]
                },
                total=rows.height,
                retained=eligible.height,
                excluded=rows.height - eligible.height,
                candidates=len(arms) - comparator_count,
                weights=weights,
                interpretation="Label-informed internal development; no untouched test or promotion",
            ),
            indent=2,
        )
        + "\n"
    )
