"""Frozen low-order geometry compatibility and the declared contact penalties.

Reference observations supply the prior's actual sample count. A normal model
with reference prior p(mu,sigma^2) proportional to 1/sigma^2 gives a Student-t
predictive after the external observations, then the native update. Correlated
contacts are likelihood observations, not independent complexes. No labels enter.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import t
import typer

from .geometry import BACKBONE_ATOMS

VERSION = "contact-marginals-2"
ARMS = ("fixed", "matrix", "flip", "flip_favourable")
EMPTY = {"n": 0, "mean": 0., "m2": 0.}


def combine_stats(a, b):
    """Merge n, mean and sum of squared deviations without unstable raw squares."""
    n = a["n"] + b["n"]
    if not n:
        return dict(EMPTY)
    delta = b["mean"] - a["mean"]
    return {"n": n, "mean": a["mean"] + delta*b["n"]/n,
            "m2": a["m2"] + b["m2"] + delta*delta*a["n"]*b["n"]/n}


def predictive_compatibility(x, stats):
    """Two-sided transformed-coordinate predictive tail, in [0,1], no floor."""
    n = np.asarray(stats["n"], dtype=float)
    m2 = np.asarray(stats["m2"], dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        scale = np.sqrt(m2/(n-1)*(1+1/n))
        s = 2*t.sf(np.abs(np.asarray(x) - stats["mean"])/scale, n-1)
    return np.where((n > 1) & (m2 > 0), s, np.nan)


def penalty_contributions(energy, weight, matrix_max):
    """Author's four resulting contributions, preserving missing table entries."""
    e, u = np.asarray(energy, float), np.asarray(weight, float)
    if e.shape != u.shape or not np.isfinite(matrix_max) or matrix_max < 0:
        raise ValueError("equal shapes and finite nonnegative full-matrix maximum required")
    if np.isinf(e).any() or np.isinf(u).any() or ((u < 0) | (u > 1)).any():
        raise ValueError("finite-or-missing energies and weights in [0,1] required")
    return dict(fixed=e+u, matrix=e+u*matrix_max, flip=(1-2*u)*e,
                flip_favourable=e+2*u*np.maximum(-e, 0))


def marginal_samples(table: pl.DataFrame, *, fitting=False) -> pl.DataFrame:
    """Atom-symmetric distance and binder-to-target orientation marginals.

    Distance uses exact unordered amino-acid:atom identities; fallback chemistry
    pools by element and backbone/sidechain. Angles use residue-pair identities;
    omega retains binder,target residue order; theta/phi additionally identify
    the binder or target endpoint. Angular fallback retains the endpoint role.
    Native roles orient TCR→peptide/MHC and peptide→MHC. External tables with
    source_chain/source_resid identities and no roles declare from=binder.
    atanh(cos(angle)) discards torsion sign;
    it is an approximate marginal normal model, never an independent joint one.
    Exact transform boundaries and missing geometry stay unavailable, not clipped.
    """
    table = table.with_row_index("row")
    for side in ("from", "to"):
        table = table.with_columns(pl.col("chem."+side).str.split_exact(":", 1).struct.field("field_0").alias("aa."+side),
            pl.col("chem."+side).str.split_exact(":", 1).struct.field("field_1").alias("name."+side))
        table = table.with_columns(pl.concat_str(pl.col("name."+side).str.slice(0, 1),
            pl.when(pl.col("name."+side).is_in(list(BACKBONE_ATOMS))).then(pl.lit("bb")).otherwise(pl.lit("sc")), separator=":").alias("pool."+side))
    roles = {"role.from", "role.to"}
    external_ids = {"source_chain.from", "source_resid.from", "source_chain.to", "source_resid.to"}
    if roles <= set(table.columns):
        a, b = pl.col("role.from"), pl.col("role.to")
        forward = ((a == "tcr") & b.is_in(["peptide", "mhc"])) | ((a == "peptide") & (b == "mhc"))
        backward = ((b == "tcr") & a.is_in(["peptide", "mhc"])) | ((b == "peptide") & (a == "mhc"))
        if table.filter(~(forward | backward).fill_null(False)).height:
            raise ValueError("geometry requires biological binder/target roles: tcr→peptide/mhc or peptide→mhc")
        if "interface" in table.columns:
            expected = (pl.when(((a == "tcr") & (b == "peptide")) | ((b == "tcr") & (a == "peptide")))
                .then(pl.lit("tcr_peptide"))
                .when((a == "tcr") | (b == "tcr")).then(pl.lit("tcr_mhc"))
                .otherwise(pl.lit("peptide_mhc")))
            if table.filter((pl.col("interface") != expected).fill_null(True)).height:
                raise ValueError("geometry interface disagrees with biological binder/target roles")
    elif not (roles & set(table.columns)) and external_ids <= set(table.columns):
        forward = pl.lit(True)
    else:
        raise ValueError("geometry requires role.from/role.to or explicit external source_chain/source_resid identities")
    table = table.with_columns(
        pl.when(forward).then(pl.col("aa.from")).otherwise(pl.col("aa.to")).alias("aa.binder"),
        pl.when(forward).then(pl.col("aa.to")).otherwise(pl.col("aa.from")).alias("aa.target"),
        *[pl.when(forward).then(pl.col(f"{family}.{side}_rad"))
          .otherwise(pl.col(f"{family}.{other}_rad")).alias(f"{family}.{role}_rad")
          for family in ("theta", "phi")
          for role, side, other in (("binder", "from", "to"), ("target", "to", "from"))])
    def unordered(a, b):
        return pl.concat_str(pl.min_horizontal(a, b), pl.max_horizontal(a, b), separator="|")
    d = np.asarray(table["dist"], float)
    if np.isfinite(d).any() and ((d[np.isfinite(d)] <= 0) | (d[np.isfinite(d)] > 5)).any():
        raise ValueError("geometry observations require 0 < distance <= 5 Angstrom")
    with np.errstate(divide="ignore", invalid="ignore"):
        distance = np.log(d/(5-d))
    parts = [table.select("row",
        (pl.lit("distance|")+unordered("chem.from", "chem.to")).alias("key"),
        (pl.lit("distance|")+unordered("pool.from", "pool.to")).alias("pool"))
        .with_columns(pl.Series("x", distance))]
    if fitting:
        keys = (["pdb.id", "chain.id.from", "residue.index.from", "chain.id.to", "residue.index.to"]
                if "residue.index.from" in table.columns else
                ["pdb.id", "source_chain.from", "source_resid.from", "source_chain.to", "source_resid.to"])
        angles = ["omega_rad", "theta.binder_rad", "theta.target_rad", "phi.binder_rad", "phi.target_rad"]
        conflicts = table.group_by(keys).agg([pl.col(c).drop_nulls().n_unique().alias(c) for c in angles])
        if conflicts.select(pl.any_horizontal(pl.col(angles) > 1).any()).item():
            raise ValueError("conflicting recovered orientations for one physical residue edge")
        # Missing atom correspondence must not make reference fitting row-order dependent.
        table = table.group_by(keys).agg(pl.col("row").min(), pl.col("aa.binder").first(), pl.col("aa.target").first(),
            *[pl.col(c).drop_nulls().first() for c in angles])
    for family, side in (("omega", None), ("theta", "binder"), ("theta", "target"), ("phi", "binder"), ("phi", "target")):
        column = family + ("."+side if side else "") + "_rad"
        angle = np.asarray(table[column], float)
        with np.errstate(divide="ignore", invalid="ignore"):
            x = np.arctanh(np.cos(angle))
        key = pl.concat_str("aa.binder", "aa.target", separator="|")
        if side:
            key = pl.concat_str(key, pl.lit(side), separator="|")
        parts.append(table.select("row", (pl.lit(family+"|")+key).alias("key"),
            pl.lit(family+"|"+(side or "all")).alias("pool")).with_columns(pl.Series("x", x)))
    return pl.concat(parts).filter(pl.col("x").is_finite() & pl.col("key").is_not_null())


def _stats(samples, group):
    return samples.group_by(group).agg(pl.len().alias("n"), pl.col("x").mean().alias("mean"),
        (pl.col("x").var(ddof=0)*pl.len()).alias("m2"))


def fit_resource(external: pl.DataFrame, native: pl.DataFrame):
    """External-data posterior becomes the prior; native sufficient stats update it."""
    es, ns = marginal_samples(external, fitting=True), marginal_samples(native, fitting=True)
    ext = {r["key"]: r for r in _stats(es, ["key", "pool"]).to_dicts()}
    pools = {r["pool"]: r for r in _stats(es, "pool").to_dicts()}
    nat = {r["key"]: r for r in _stats(ns, ["key", "pool"]).to_dicts()}
    records = []
    for key in sorted(ext.keys() | nat.keys()):
        pool = (ext.get(key) or nat[key])["pool"]
        prior = ext.get(key, EMPTY)
        direct = prior["n"] > 1 and prior["m2"] > 0
        if not direct:
            prior = pools.get(pool, EMPTY)
        update = nat.get(key, EMPTY)
        post = combine_stats(prior, update)
        records.append(dict(key=key, pool=pool, prior_source="exact" if direct else "pooled_chemistry",
            prior_n=prior["n"], prior_mean=prior["mean"], prior_m2=prior["m2"],
            native_n=update["n"], native_mean=update["mean"], native_m2=update["m2"], **post))
    def frequencies(table):
        return table.select(pl.min_horizontal("chem.from", "chem.to").alias("a"),
                            pl.max_horizontal("chem.from", "chem.to").alias("b")).drop_nulls().group_by("a", "b").len()
    categories = frequencies(external).rename({"len": "external_n"}).join(
        frequencies(native).rename({"len": "native_n"}), on=["a", "b"], how="full", coalesce=True).fill_null(0)
    supported = categories.filter(pl.col("external_n") > 0)
    normalizer = supported.select((pl.col("external_n") + pl.col("native_n")).sum()).item()
    categories = categories.with_columns(pl.when(pl.col("external_n") > 0)
        .then((pl.col("external_n")+pl.col("native_n"))/normalizer).otherwise(None).alias("posterior_mass_conditional_external_support"))
    return dict(version=VERSION, cutoff_A=5., arms=list(ARMS),
        angular_schema="omega|binderAA|targetAA; theta/phi|binderAA|targetAA|binder-or-target; endpoint-specific angular pools",
        biological_direction="TCR→peptide/MHC; peptide→MHC; external source_chain/source_resid tables without roles assume from=binder, to=target",
        reduction="worst_supported_atom_and_marginal_per_physical_residue_pair",
        prior_strength="actual external observation count; reference prior p(mu,sigma^2) proportional to 1/sigma^2",
        transforms="log(d/(5-d)); atanh(cos(angle)); boundaries missing; signed torsion information discarded",
        compatibility="2*StudentT.sf(abs(x-mean)/sqrt(m2/(n-1)*(1+1/n)),df=n-1)",
        weight="1-min(available predictive compatibility); unavailable marginals reported; all missing yields missing",
        interpretation="Approximate marginal predictive compatibility, not P(natural), binding probability or free energy; correlated observations; no marginal independence assertion",
        records=records, pools=list(pools.values()), categories=categories.sort("a", "b").to_dicts(),
        categorical_model="Dirichlet alpha=external integer counts on observed external support; no invented positive mass for unseen categories; native outside-support mass unavailable; frequencies do not determine unnaturalness",
        reference_counts=dict(external_atom_pairs=external.height, native_atom_pairs=native.height,
            external_marginal_observations=es.height, native_marginal_observations=ns.height))


def marginal_compatibilities(table: pl.DataFrame, resource: dict) -> pl.DataFrame:
    """Expose the exact frozen predictive tails and parameter lookup for auditing."""
    if resource.get("version") != VERSION or resource.get("arms") != list(ARMS):
        raise ValueError("unsupported frozen contact resource: directed contact-marginals-2 required; refit historical v1 resources")
    samples = marginal_samples(table)
    params = pl.DataFrame(resource["records"]).select("key", "n", "mean", "m2", "prior_source")
    pools = pl.DataFrame(resource["pools"]).select("pool", "n", "mean", "m2").rename(
        {c: c+"_pool" for c in ("n", "mean", "m2")})
    joined = samples.join(params, on="key", how="left", validate="m:1").join(pools, on="pool", how="left", validate="m:1")
    joined = joined.with_columns(pl.when(pl.col("n").is_null()).then(pl.lit("scoring_pool")).otherwise(pl.col("prior_source")).alias("parameter_source"),
        *[pl.coalesce(c, c+"_pool").alias(c) for c in ("n", "mean", "m2")])
    s = predictive_compatibility(joined["x"].to_numpy(), {c: joined[c].to_numpy() for c in ("n", "mean", "m2")})
    return joined.with_columns(pl.Series("compatibility", s).fill_nan(None))


def atom_weights(table: pl.DataFrame, resource: dict) -> pl.DataFrame:
    """Apply immutable predictive parameters; no update or label-dependent choice."""
    joined = marginal_compatibilities(table, resource)
    s = joined["compatibility"].to_numpy()
    ok, rows = np.isfinite(s), joined["row"].to_numpy()
    compatibility, supported = np.full(table.height, np.inf), np.zeros(table.height, dtype=np.int64)
    np.minimum.at(compatibility, rows[ok], s[ok])
    np.add.at(supported, rows[ok], 1)
    weight = np.where(supported > 0, 1-compatibility, np.nan)
    return table.with_columns(pl.Series("unnaturalness", weight).fill_nan(None),
                              pl.Series("supported_marginals", supported))


def fit_command(
    external: Path = typer.Option(..., "--external", exists=True, dir_okay=False),
    native: Path = typer.Option(..., "--native", exists=True, dir_okay=False),
    out: Path = typer.Option(..., "-o", "--out"),
):
    """Freeze external/native geometry sufficient statistics, without benchmark data."""
    from ..exports import _sha256
    if out.exists():
        raise FileExistsError(out)
    columns = ["pdb.id", "chem.from", "chem.to", "dist", "omega_rad", "theta.from_rad", "theta.to_rad", "phi.from_rad", "phi.to_rad"]
    external_keys = ["source_chain.from", "source_resid.from", "source_chain.to", "source_resid.to"]
    ex = pl.read_csv(external, separator="\t", columns=columns+external_keys,
        schema_overrides={c: pl.String for c in external_keys}, infer_schema_length=10000)
    na = pl.read_parquet(native)
    resource = fit_resource(ex, na)
    resource["sources"] = {"external_sha256": _sha256(external), "native_sha256": _sha256(native)}
    resource["producer_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    out.write_text(json.dumps(resource, indent=2, allow_nan=False)+"\n")
    # Descriptive reference checks only: no label access, model selection or tests.
    qa = []
    for name, table in (("external", ex), ("native", na)):
        sample = marginal_samples(table, fitting=True)
        stats = pl.DataFrame(resource["records"]).select("key", "n", "mean", "m2")
        joined = sample.join(stats, on="key", how="left", validate="m:1")
        compatibility = predictive_compatibility(joined["x"].to_numpy(), {c: joined[c].to_numpy() for c in ("n", "mean", "m2")})
        checked = joined.with_columns(pl.Series("compatibility", compatibility).fill_nan(None),
            pl.col("key").str.split("|").list.first().alias("family"))
        for row in checked.group_by("family").agg(pl.len().alias("n_marginal_observations"),
            pl.col("compatibility").count().alias("n_supported"), pl.col("compatibility").mean().alias("mean_compatibility"),
            pl.col("compatibility").median().alias("median_compatibility"), pl.col("compatibility").min().alias("minimum_compatibility")).to_dicts():
            qa.append(dict(reference=name, **row))
    pl.DataFrame(qa).write_csv(str(out)+".reference_checks.tsv", separator="\t")
    print(json.dumps({"resource": str(out), "sha256": _sha256(out), "reference_counts": resource["reference_counts"]}))


def score_contact_penalties(contacts, observations, potential, resource):
    """Preserve legacy all-TCR energy; apply CDR geometry increments once per edge.

    Unknown potential entries retain the legacy skipped-entry total and explicit
    counts/null contact contributions. Framework contacts remain unchanged and
    explicitly unassessed. A scored CDR edge without geometry makes new totals
    unavailable. Zero-contact models retain zero totals.
    """
    keys = ["chain.id.from", "residue.index.from", "chain.id.to", "residue.index.to"]
    if contacts.select(keys).is_duplicated().any():
        raise ValueError("penalties require unique physical residue-pair contacts")
    atoms = atom_weights(observations, resource)
    swap = pl.col("role.from") != "tcr"
    atoms = atoms.with_columns([
        pl.when(swap).then(pl.col(field+".to")).otherwise(pl.col(field+".from")).alias("rec_"+field)
        for field in ("chain.id", "residue.index")] + [
        pl.when(swap).then(pl.col(field+".from")).otherwise(pl.col(field+".to")).alias("pep_"+field)
        for field in ("chain.id", "residue.index")])
    groups = ["rec_chain.id", "rec_residue.index", "pep_chain.id", "pep_residue.index"]
    weights = atoms.group_by(groups).agg(pl.col("unnaturalness").max(),
        pl.len().alias("geometry_atom_pairs"), (pl.col("supported_marginals") == 0).sum().alias("unsupported_atom_pairs"),
        pl.col("supported_marginals").max().alias("max_supported_marginals"))
    weights = weights.rename(dict(zip(groups, keys)))
    result = contacts.join(weights, on=keys, how="left", validate="1:1")
    cdr = pl.col("region.type.from").is_in(["CDR1", "CDR2", "CDR3"]).fill_null(False)
    result = result.with_columns(pl.when(cdr).then(pl.col("unnaturalness")).otherwise(0.).alias("unnaturalness"),
        pl.when(~cdr).then(pl.lit("unassessed_framework"))
        .when(pl.col("unnaturalness").is_null()).then(pl.lit("unavailable_cdr_geometry"))
        .otherwise(pl.lit("supported_cdr_geometry")).alias("geometry_status"))
    matrix, index = potential.as_matrix()
    finite = matrix[np.isfinite(matrix)]
    if not finite.size:
        raise ValueError("potential has no finite entries")
    maximum = float(np.abs(finite).max())
    ai = np.array([index.get(a, -1) for a in result["residue.aa.from"]], dtype=np.int64)
    bi = np.array([index.get(a, -1) for a in result["residue.aa.to"]], dtype=np.int64)
    e = np.full(result.height, np.nan)
    valid = (ai >= 0) & (bi >= 0)
    e[valid] = matrix[ai[valid], bi[valid]]
    contributions = penalty_contributions(e, result["unnaturalness"].to_numpy(), maximum)
    result = result.with_columns(pl.Series("original_contribution", e).fill_nan(None),
                                pl.lit(potential.name).alias("potential"))
    summary = dict(potential=potential.name, original=float(np.nansum(e)),
        n_contacts=result.height, missing_potential_entries=int((~np.isfinite(e)).sum()),
        missing_scored_geometry=int((np.isfinite(e) & ~np.isfinite(result["unnaturalness"].to_numpy())).sum()),
        matrix_max_abs=maximum, matrix_scale_degenerate=maximum == 0,
        unassessed_framework_contacts=result.filter(pl.col("geometry_status")=="unassessed_framework").height)
    for arm, value in contributions.items():
        result = result.with_columns(pl.Series(arm+"_increment", value-e).fill_nan(None),
                                    pl.Series(arm+"_contribution", value).fill_nan(None))
        summary[arm] = None if summary["missing_scored_geometry"] else float(np.nansum(value))
    from ..energetics.scoring import _interface_energy
    if not np.isclose(summary["original"], _interface_energy(contacts, potential), rtol=0, atol=1e-12):
        raise ValueError("original energy differs from the existing scorer")
    return result, summary


def comparison_command(
    labels: Path = typer.Option(..., "--labels", exists=True, dir_okay=False),
    scores: Path = typer.Option(..., "--scores", exists=True, dir_okay=False),
    out: Path = typer.Option(..., "-o", "--out"),
):
    """Compare all arms and native ipTM on one common finite population per table."""
    from ..associations import GROUPS, _unique, _metrics
    from ..exports import _sha256
    label = pl.read_csv(labels, separator="\t", infer_schema_length=None)
    scored = pl.read_csv(scores, separator="\t", infer_schema_length=None)
    _unique(label, ["id"], "labels")
    _unique(scored, ["id", "potential"], "scores")
    if not set(GROUPS+["id", "y", "iptm"]) <= set(label.columns):
        raise ValueError("labels require explicit task/class/template strata, binary y and native iptm")
    y = label["y"].to_numpy()
    if not np.isin(y, [0, 1]).all() or label.select(pl.any_horizontal(pl.col(GROUPS).is_null()).any()).item():
        raise ValueError("finite binary labels and nonnull cohort keys required")
    if set(label["mhc_class"]) - {"I", "II"}:
        raise ValueError("MHC class must be I or II")
    classes = scored.select("id", "mhc_class").unique().join(label.select("id", "mhc_class"), on="id", suffix="_label")
    if classes.filter(pl.col("mhc_class").is_null() | (pl.col("mhc_class") != pl.col("mhc_class_label"))).height:
        raise ValueError("MHC class differs between coordinates and benchmark labels")
    finite_iptm = label.filter(pl.col("iptm").is_finite())["iptm"].to_numpy()
    if np.any((finite_iptm < 0) | (finite_iptm > 1)):
        raise ValueError("finite native iptm must be in [0,1]")
    if scored.join(label.select("id"), on="id", how="anti").height:
        raise ValueError("unexpected scored rows")
    out.mkdir(parents=True, exist_ok=False)
    rows, members = [], []
    arms = ["original", "weighted"] if "weighted" in scored.columns else ["original", *ARMS]
    if "normalized" in scored.columns:
        arms.append("normalized")
    if "atom_count" in scored.columns:
        arms.append("atom_count")
    if "surface_context_sum" in scored.columns:
        arms.extend(["surface_raw_sum", "surface_context_sum"])
    for (potential,), current in scored.group_by("potential"):
        table = label.join(current.drop("mhc_class"), on="id", how="left", validate="1:1")
        table = table.with_columns(pl.all_horizontal(pl.col([*arms, "iptm"]).is_finite()).fill_null(False).alias("used"))
        members.append(table.select("id", *GROUPS, "used").with_columns(pl.lit(potential).alias("potential")))
        for key, part in table.group_by(GROUPS, maintain_order=True):
            context = dict(zip(GROUPS, key))
            chosen = part.filter(pl.col("used"))
            yy = chosen["y"].to_numpy()
            base = _metrics(yy, -chosen["original"].to_numpy())
            iptm = _metrics(yy, chosen["iptm"].to_numpy())
            for arm in [*arms, "iptm"]:
                metrics = _metrics(yy, chosen["iptm"].to_numpy() if arm == "iptm" else -chosen[arm].to_numpy())
                for metric, value in metrics.items():
                    original, confidence = base[metric], iptm[metric]
                    rows.append(dict(**context, potential=potential, arm=arm, metric=metric,
                        n_rows=part.height, n_used=chosen.height, n_excluded=part.height-chosen.height,
                        n_positive=int(yy.sum()), n_negative=int(chosen.height-yy.sum()),
                        value=value, original_value=original, iptm_value=confidence,
                        delta_original=None if value is None else value-original,
                        delta_iptm=None if value is None else value-confidence,
                        loss_original=None if value is None else value < original,
                        models_with_missing_potential_entries=chosen.filter(pl.col("missing_potential_entries")>0).height))
    pl.DataFrame(rows, infer_schema_length=None).sort(GROUPS+["potential", "metric", "arm"]).write_csv(out/"metrics.tsv", separator="\t")
    pl.concat(members).write_csv(out/"membership.tsv", separator="\t")
    (out/"provenance.json").write_text(json.dumps(dict(labels_sha256=_sha256(labels), scores_sha256=_sha256(scores),
        directions={**dict.fromkeys(arms, -1), "iptm": 1}, no_fitting=True,
        pairing="declared arms, original energy and ipTM finite on exactly the same IDs within each potential/cohort/class/template stratum"), indent=2)+"\n")
