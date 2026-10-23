"""Task-matched cohort screening with frozen references and explicit availability.

Peptide ranking fixes a receptor and varies epitopes; receptor ranking fixes a
pMHC and varies receptors. The structural alanine reference follows the varying
side. Topology is included only for the receptor-ranking task.

The bundled TCR readout is (6 * ipTM + cycle fraction) / 7. Epitope and legacy
references combine standardized components with weighted Stouffer. ipTM is caller-supplied and remains
the primary receptor-ranking comparator. The combined score is dimensionless,
not a calibrated probability or a physical binding free energy.

Models with no TCR:peptide contacts are retained with discard reasons and do not
count toward the cohort floor. Below the default 30-model floor, components are
reported without a combined score. Cohort flags report by default; explicit
responses preserve discarded and withheld rows. No query-cohort normalization
or fitting is performed. Equations and output semantics: docs/screening.rst.
"""

from __future__ import annotations

import hashlib
import json
from enum import Enum
from importlib.resources import files
from pathlib import Path

import numpy as np
import polars as pl

FORMAT = "tcren-screen-1"
#: Below this many models the combined score is withheld; see the module docstring.
MIN_COHORT = 30
#: Fallback confidence weight when a reference does not declare one. Each shipped reference
#: carries its own ``default_weight``, because the right weight is task-specific: ipTM is much
#: the strongest input for TCR ranking and the energy is stronger for epitope ranking.
DEFAULT_WEIGHT = 0.8


class CohortKind(str, Enum):
    """Which side varies across the cohort, and therefore which score applies."""

    EPITOPE_RANKING = "epitope"
    TCR_RANKING = "tcr"


RESOURCE = {CohortKind.TCR_RANKING: "screen_reference_tcr.json",
            CohortKind.EPITOPE_RANKING: "screen_reference_epitope.json"}


def bundled_reference(task: CohortKind) -> dict:
    """The shipped frozen reference for one task, with its own hash for provenance."""
    raw = files("tcren.data").joinpath(RESOURCE[CohortKind(task)]).read_bytes()
    reference = json.loads(raw)
    reference["reference_sha256"] = hashlib.sha256(raw).hexdigest()
    return reference


def template_coverage(peptide: str, metadata: Path) -> dict:
    """Epitope-level crystal availability from an explicit, hash-recorded metadata TSV.

    The caller supplies a native-crystal inventory (for example Native2026), not
    model metadata. Membership is case-insensitive; it does not establish matching
    MHC, construct correspondence or pose accuracy. No network lookup is performed.
    """
    raw = Path(metadata).read_bytes()
    table = pl.read_csv(raw, separator="\t", infer_schema_length=None)
    if "epitope" not in table.columns:
        raise ValueError("template metadata needs an epitope column")
    epitopes = set(table["epitope"].drop_nulls().cast(pl.String)
                   .str.strip_chars().str.to_uppercase().to_list()) - {""}
    peptide = peptide.strip().upper()
    if not peptide or not epitopes:
        raise ValueError("template metadata and peptide must be nonempty")
    return {"native_template": peptide in epitopes, "epitope": peptide,
            "n_epitopes": len(epitopes),
            "metadata_sha256": hashlib.sha256(raw).hexdigest()}


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
        fraction = c.get("structural_fraction")
        if fraction is not None and (
                c["name"] == "confidence" or not np.isfinite(fraction) or not 0 <= fraction <= 1):
            raise ValueError("structural_fraction must lie in [0, 1] on structural components")
        if "confidence_rho" in c:
            raise ValueError("confidence_rho is retired; remove it from the reference")
    method = reference.get("score_method", "stouffer")
    if method not in {"stouffer", "raw_weighted_mean"}:
        raise ValueError("unknown score_method")
    if method == "raw_weighted_mean":
        active = {c["name"]: c for c in components
                  if c["name"] == "confidence" or c.get("structural_fraction", 1) != 0}
        if (reference.get("task") != "tcr" or set(active) != {"confidence", "topology"}
                or active["confidence"]["column"] != "iptm"
                or active["topology"]["column"] != "g_cyclo_frac"
                or any(c["direction"] != 1 for c in active.values())
                or active["topology"].get("structural_fraction") != 1):
            raise ValueError("raw_weighted_mean requires raw ipTM and cycle fraction only")
    return components


def stouffer(z: np.ndarray, weights: np.ndarray) -> np.ndarray:
    """Weighted Stouffer combination, ``sum(w_i z_i) / sqrt(sum w_i^2)``, row-wise.

    Reduces to plain Stouffer when the weights are equal. Rows are combined over the columns that
    are finite, so a missing component is dropped from both sums rather than imputed; a row with
    nothing finite returns NaN.

    Args:
        z: ``(n_rows, n_components)`` standardized, direction-corrected values.
        weights: one non-negative weight per component.
    """
    z = np.atleast_2d(np.asarray(z, dtype=float))
    w = np.asarray(weights, dtype=float)
    if w.ndim != 1 or len(w) != z.shape[1]:
        raise ValueError("one weight per component is required")
    if (w < 0).any() or not np.isfinite(w).all() or w.sum() <= 0:
        raise ValueError("weights must be finite, non-negative and not all zero")
    seen = np.isfinite(z)
    num = np.nansum(np.where(seen, z, 0.0) * w, axis=1)
    den = np.sqrt(((w ** 2) * seen).sum(axis=1))
    return np.where(den > 0, num / np.where(den > 0, den, 1.0), np.nan)


def cohort_signature(peptides, receptors) -> tuple[CohortKind, dict]:
    """Decide which side is held fixed, or refuse and say what differs.

    Args:
        peptides: one peptide identity per structure, in input order.
        receptors: one receptor identity per structure, typically the ``(CDR3a, CDR3b)`` pair
            read off the annotated FR3/FR4 boundaries. Never reconstruct a junction as
            ``"C" + cdr3 + "F"``; the anchors are a lookup.

    Returns:
        The cohort kind and a summary dict with the distinct counts actually observed.
    """
    peptides, receptors = list(peptides), list(receptors)
    if len(peptides) != len(receptors):
        raise ValueError("one peptide and one receptor identity per structure is required")
    if not peptides:
        raise ValueError("no structures to screen")
    n_pep, n_tcr = len({*peptides}), len({*receptors})
    summary = {"n_structures": len(peptides), "n_distinct_peptides": n_pep,
               "n_distinct_receptors": n_tcr}
    if n_tcr == 1 and n_pep > 1:
        return CohortKind.EPITOPE_RANKING, summary
    if n_pep == 1 and n_tcr > 1:
        return CohortKind.TCR_RANKING, summary
    if n_pep == 1 and n_tcr == 1:
        raise ValueError(
            f"all {len(peptides)} structures share one peptide and one receptor, so nothing "
            "varies and there is nothing to rank")
    raise ValueError(
        f"not a cohort: {n_pep} distinct peptides and {n_tcr} distinct receptors across "
        f"{len(peptides)} structures. A cohort holds one side fixed — split the input by peptide "
        "or by receptor, or use `tcren features` for a descriptor table")


def identities(structures) -> tuple[list[str], list[tuple[str, str]]]:
    """Per structure, the peptide sequence and the ``(CDR3a, CDR3b)`` pair.

    The CDR3s are read off the **annotated** CDR3 regions, whose boundaries come from the
    FR3/FR4 markup. They are never reconstructed as ``"C" + cdr3 + "F"`` -- the junction anchors
    are a lookup, and that shortcut breaks on tryptophan anchors and on CDR3s ending in F.

    Raises:
        ValueError: if a structure carries no peptide chain or no CDR3 on either TCR chain,
            which means it cannot be placed in a cohort.
    """
    peptides, receptors = [], []
    for s in structures:
        peptide = next((c.sequence() for c in s.chains if c.chain_type == "PEPTIDE"), None)
        if not peptide:
            raise ValueError(f"{s.pdb_id}: no annotated peptide chain, so it has no cohort identity")
        loops = {}
        for chain in s.chains:
            if chain.chain_type not in ("TRA", "TRB"):
                continue
            cdr3 = next((r.sequence for r in chain.regions if r.region_type == "CDR3"), None)
            if cdr3:
                loops[chain.chain_type] = cdr3
        if not loops:
            raise ValueError(f"{s.pdb_id}: no annotated CDR3 on either TCR chain, so it has no "
                             "cohort identity")
        peptides.append(peptide)
        receptors.append((loops.get("TRA", ""), loops.get("TRB", "")))
    return peptides, receptors


#: What each flag means, in the order a reader should read them. Kept here so the CLI, the docs
#: and the tests cannot drift apart.
FLAGS = {
    "small_cohort":
        "fewer models than --min-cohort, so the combined score is withheld and only the "
        "components are reported",
    "forced_pose":
        "the cohort's mean interface strain is above the reference cut: the generator had to "
        "force these poses, and the contact energy is known to invert on forced poses",
    "low_confidence":
        "the cohort's mean ipTM sits well below the reference distribution of per-epitope means, "
        "the regime where the generator's own ranking is least reliable; read "
        "cohort_confidence_z for how far below",
    "no_interface":
        "models in this cohort have no TCR:peptide contact at all -- the receptor was not placed "
        "on the groove, so every structural component for them is meaningless",
    "confidence_disagrees":
        "ipTM and the structural block rank this cohort differently, so at most one of them can "
        "be right about it",
}


def cohort_flags(
    table: pl.DataFrame,
    reference: dict,
    *,
    min_cohort: int = MIN_COHORT,
    strain: np.ndarray | None = None,
    calibrated: bool = False,
) -> dict:
    """Label-free warnings about the cohort as a whole. **Advisory, never applied to the score.**

    Every flag reads the cohort's own descriptors and no outcome, so all of them are available on
    an unseen epitope. They mark a *regime*, not a verdict: a flagged cohort is one where the
    evidence says a ranking may be untrustworthy, and the right response is to look rather than to
    accept a number.

    They do not act on the score unless a caller asks. Applicability observations
    do not validate an automatic corrective response: see :data:`WITHHOLDING_FLAGS`.

    Args:
        table: the same descriptor table passed to :func:`screen_scores`.
        reference: the frozen reference, whose ``cohort_confidence`` moments and ``forced_pose``
            cut supply the thresholds.
        min_cohort: the floor below which a score is withheld.
        strain: optional per-model interface strain, from :func:`tcren.cohort.strain_z`. Without
            it the forced-pose flag is reported as unavailable rather than guessed.
        calibrated: whether ``strain`` was standardized against the same population the reference
            cut came from. It is not comparable otherwise, so the flag stays unavailable.

    Returns:
        A dict of flag name to bool (or ``None`` where the inputs for that flag were absent),
        alongside the statistics the flags were derived from.
    """
    conf_ref = reference.get("cohort_confidence", {})
    forced_ref = reference.get("forced_pose", {})
    out: dict[str, object] = {}
    out["n_models"] = int(table.height)
    out["small_cohort"] = table.height < min_cohort

    # A model with no TCR:peptide contact is not a TCR:pMHC complex: the receptor was docked
    # elsewhere. Its energy, shape and topology are all computed over an empty interface, so they
    # carry no information, and a cohort where many models are like this cannot be ranked.
    contacts = (np.asarray(table["n_contacts_tp"].cast(pl.Float64), float)
                if "n_contacts_tp" in table.columns else None)
    if contacts is not None and np.isfinite(contacts).any():
        n_zero = int((contacts == 0).sum())
        out["small_cohort"] = table.height - n_zero < min_cohort
        out["n_no_interface"] = n_zero
        out["frac_no_interface"] = float(n_zero / max(len(contacts), 1))
        out["no_interface"] = bool(n_zero > 0)
    else:
        out["n_no_interface"] = None
        out["frac_no_interface"] = float("nan")
        out["no_interface"] = None

    iptm = (np.asarray(table["iptm"].cast(pl.Float64), float)
            if "iptm" in table.columns else np.array([np.nan]))
    mean_iptm = float(np.nanmean(iptm)) if np.isfinite(iptm).any() else float("nan")
    out["mean_iptm"] = mean_iptm
    # A z against the reference distribution, not a side of a line: the distance is the point.
    # A cohort at -2.7 and one at -0.7 are both "below the cut" and are not the same situation.
    mu, sd = conf_ref.get("mean"), conf_ref.get("sd")
    z = ((mean_iptm - mu) / sd if np.isfinite(mean_iptm) and mu is not None and sd
         else float("nan"))
    out["cohort_confidence_z"] = z
    flag_z = conf_ref.get("flag_z")
    out["low_confidence"] = (None if not np.isfinite(z) or flag_z is None else bool(z < flag_z))

    # strain_z standardizes within whatever set it is given, so a cohort's mean is only
    # comparable to the reference cut when both were standardized against the same population.
    # Without that the flag is UNAVAILABLE rather than False: silently reading False would claim
    # a cohort is not forced on the strength of an incomparable number.
    mean_strain = float(np.nanmean(strain)) if strain is not None and np.isfinite(strain).any() \
        else float("nan")
    out["mean_strain"] = mean_strain
    cut = forced_ref.get("mean_strain_cut")
    out["forced_pose"] = (None if not np.isfinite(mean_strain) or cut is None or not calibrated
                          else bool(mean_strain > cut))

    # Do the two blocks of evidence even agree about this cohort? If they rank it oppositely, at
    # most one of them is right, and the reader should know before trusting either.
    scored = screen_scores(table, reference, min_cohort=1)
    conf = np.asarray(scored["z_confidence"].cast(pl.Float64), float)
    structural = [c["name"] for c in _validate(reference) if c["name"] != "confidence"]
    if structural and np.isfinite(conf).sum() > 2:
        block = np.nanmean(np.column_stack(
            [np.asarray(scored[f"z_{n}"].cast(pl.Float64), float) for n in structural]), axis=1)
        m = np.isfinite(conf) & np.isfinite(block)
        if m.sum() > 2 and np.ptp(conf[m]) > 0 and np.ptp(block[m]) > 0:
            from scipy.stats import spearmanr
            rho = float(spearmanr(conf[m], block[m]).statistic)
            out["confidence_structure_agreement"] = rho
            out["confidence_disagrees"] = bool(rho < 0.0)
        else:
            out["confidence_structure_agreement"] = float("nan")
            out["confidence_disagrees"] = None
    else:
        out["confidence_structure_agreement"] = float("nan")
        out["confidence_disagrees"] = None
    return out


#: No cohort flag automatically withholds a score. Responses remain explicit;
#: the no-interface discard and cohort floor apply independently on both tasks.
WITHHOLDING_FLAGS: dict[CohortKind, tuple[str, ...]] = {
    CohortKind.TCR_RANKING: (),
    CohortKind.EPITOPE_RANKING: (),
}


def withhold_reason(flags: dict, task: CohortKind) -> str | None:
    """The first withholding flag that fired for this task, as the reason text, or ``None``."""
    for name in WITHHOLDING_FLAGS[CohortKind(task)]:
        if flags.get(name) is True:
            return f"{name}: {FLAGS[name]}"
    return None


class FlagResponse(str, Enum):
    """Explicit responses to a low-confidence cohort; REPORT leaves scores unchanged.

    Alternatives change only the combined score. Components, missingness and
    original reasons for unavailable scores remain visible. Reversal is not a
    validated correction for modelling failure.
    """

    REPORT = "report"
    STRUCTURE = "structure"
    GEOMETRY = "geometry"
    WITHHOLD = "withhold"
    REVERSE = "reverse"


def apply_flag_response(
    scored: pl.DataFrame,
    table: pl.DataFrame,
    reference: dict,
    *,
    response: FlagResponse,
    reason: str,
) -> pl.DataFrame:
    """Rewrite ``screen_score`` according to ``response``, recording what was done and why.

    The components are never rewritten -- only the combined score -- so a reader can always see
    the evidence the response was applied to.
    """
    response = FlagResponse(response)
    if response is FlagResponse.REPORT:
        return scored
    note = f"{response.value} applied — {reason}"
    if response is FlagResponse.WITHHOLD:
        out = scored.with_columns(
            screen_score=pl.lit(None, dtype=pl.Float64),
            score_unavailable=pl.lit(note))
    elif response is FlagResponse.REVERSE:
        out = scored.with_columns(screen_score=-pl.col("screen_score"))
    elif response is FlagResponse.STRUCTURE:
        out = scored.with_columns(
            screen_score=screen_scores(table, reference, weight=0.0,
                                       min_cohort=1)["screen_score"])
    else:
        geometry = next((c["name"] for c in _validate(reference) if c["name"] == "geometry"), None)
        if geometry is None:
            raise ValueError("this reference carries no geometry component")
        out = scored.with_columns(screen_score=pl.col(f"z_{geometry}"))
    # A response cannot restore a discarded/withheld score or erase its original reason.
    unavailable = scored["screen_score"].is_null()
    return out.with_columns(
        pl.when(unavailable).then(None).otherwise(pl.col("screen_score")).alias("screen_score"),
        pl.when(unavailable).then(scored["score_unavailable"])
        .otherwise(pl.col("score_unavailable")).alias("score_unavailable"),
        flag_response=pl.lit(response.value),
                            flag_response_reason=pl.lit(reason))


def screen_scores(
    table: pl.DataFrame,
    reference: dict,
    *,
    weight: float | None = None,
    native_template: bool | None = None,
    min_cohort: int = MIN_COHORT,
    suppress: str | None = None,
) -> pl.DataFrame:
    """Score one cohort, one row per model, in input order.

    Args:
        table: descriptor table keyed by a unique nonempty ``complex.id``, carrying every
            component column the reference names.
        reference: a frozen reference dict, normally from :func:`bundled_reference`.
        weight: confidence weight in the reference combination, in ``[0, 1]``. ``None`` takes the
            reference's own ``default_weight``, which is task-specific. Structural weights
            are ``(1 - weight) * structural_fraction``; omitted fractions share equally.
        native_template: optional TCR-ranking applicability gate. False uses confidence
            alone; True uses the requested weight. None keeps historical behavior.
            Components remain visible regardless of the gate. This is crystal
            availability, not a pose-quality judgement.
        min_cohort: below this many **scorable** rows the combined score is withheld. Models
            with no TCR:peptide interface are discarded and do not count toward it.
        suppress: when set, withhold the combined score for every row and record this as the
            reason. Used to act on a cohort flag -- see :func:`withhold_reason`.

    Returns:
        Per model: each component's raw value, its standardized contribution and a missingness
        flag, ``n_components``, ``no_interface``, ``screen_score`` and ``score_unavailable``. A
        row scores null rather than taking a neutral default when it has no observed component,
        when its receptor was never modelled on the groove, or when the cohort has fewer than
        ``min_cohort`` scorable models. Discarded rows are kept in the output so that nothing
        disappears silently.
    """
    components = _validate(reference)
    if weight is None:
        weight = float(reference.get("default_weight", DEFAULT_WEIGHT))
    if not np.isfinite(weight) or not 0.0 <= weight <= 1.0:
        raise ValueError("weight must lie in [0, 1]")
    if native_template is not None:
        if not isinstance(native_template, bool) or reference.get("task") != "tcr":
            raise ValueError("native_template must be a boolean for a TCR-ranking reference")
        if not native_template:
            weight = 1.0
    if "complex.id" not in table.columns:
        raise ValueError("complex.id column required")
    ids = table["complex.id"]
    if ids.dtype != pl.String or ids.null_count() or ids.n_unique() != len(ids) or (ids == "").any():
        raise ValueError("complex.id must be unique, nonempty strings")
    absent = [c["column"] for c in components if c["column"] not in table.columns]
    if absent:
        raise ValueError(f"missing component columns: {absent}")

    raw, z_cols, weights = [], [], []
    for c in components:
        x = np.asarray(table[c["column"]].cast(pl.Float64), dtype=float)
        if np.isinf(x).any():
            raise ValueError(f"infinite values in {c['column']}")
        seen = np.isfinite(x)
        z = c["direction"] * (x - c["mean"]) / c["sd"]
        raw.append(x)
        z_cols.append(np.where(seen, z, np.nan))
        # Confidence carries `weight`; the reference declares each structural fraction.
        weights.append(weight if c["name"] == "confidence" else 0.0)
    structural = [i for i, c in enumerate(components) if c["name"] != "confidence"]
    for i in structural:
        weights[i] = (1.0 - weight) * components[i].get("structural_fraction", 1.0 / len(structural))

    out = {"complex.id": ids}
    for i, c in enumerate(components):
        out[c["column"]] = raw[i]
        out[f"z_{c['name']}"] = z_cols[i]
        out[f"missing_{c['name']}"] = ~np.isfinite(raw[i])

    Z = np.column_stack(z_cols)
    observed = np.isfinite(Z).sum(axis=1).astype(np.int16)
    if reference.get("score_method") == "raw_weighted_mean":
        values = np.column_stack(raw)
        active = np.asarray(weights) > 0
        if np.any(np.isfinite(values[:, active]) &
                  ((values[:, active] < 0) | (values[:, active] > 1))):
            raise ValueError("raw weighted components must lie in [0, 1]")
        seen = np.isfinite(values)
        denominator = (seen * weights).sum(axis=1)
        numerator = (np.where(seen, values, 0) * weights).sum(axis=1)
        score = np.divide(numerator, denominator, out=np.full(len(values), np.nan),
                          where=denominator > 0)
    else:
        score = stouffer(Z, np.asarray(weights))
    no_weighted = ~np.isfinite(score)
    # A model with no TCR:peptide contact is not a TCR:pMHC complex -- the receptor was not
    # placed on the groove. Every structural component for it is computed over an empty
    # interface, so it is DISCARDED rather than scored, and it does not count toward the cohort
    # size either: a cohort that only reaches the floor by including failed models has not
    # really reached it.
    if "n_contacts_tp" in table.columns:
        contacts = np.asarray(table["n_contacts_tp"].cast(pl.Float64), float)
        no_interface = np.isfinite(contacts) & (contacts == 0)
    else:
        no_interface = np.zeros(table.height, dtype=bool)
    n_scorable = int(table.height - no_interface.sum())
    small = n_scorable < min_cohort
    reason = np.where(
        no_interface, "no TCR:peptide interface — the receptor was not modelled on the groove",
        np.where(small,
                 f"cohort of {n_scorable} scorable models is smaller than the minimum "
                 f"of {min_cohort}",
                 np.where(observed == 0, "no component observed",
                          np.where(no_weighted, "no weighted component observed",
                                   suppress if suppress else ""))))
    out["n_components"] = observed
    out["no_interface"] = no_interface
    series = pl.Series(score)
    out["screen_score"] = series.set(
        pl.Series(no_weighted | small | no_interface | bool(suppress)), None)
    out["score_unavailable"] = pl.Series(reason).set(pl.Series(reason == ""), None)
    if native_template is not None:
        out["native_template"] = [native_template] * table.height
        out["confidence_weight"] = [weight] * table.height
    return pl.DataFrame(out)


def describe_reference(reference: dict) -> str:
    """One line per component: its name, column, direction and frozen scale."""
    return "\n".join(
        f"  {c['name']:11s} {c['column']:18s} direction {c['direction']:+d}  "
        f"mean {c['mean']:+.4f}  sd {c['sd']:.4f}  n {c.get('n', '?')}"
        for c in _validate(reference))


def read_reference(path: Path) -> dict:
    """Load a reference from disk, recording its hash the same way the bundled one does."""
    raw = Path(path).read_bytes()
    reference = json.loads(raw)
    reference["reference_sha256"] = hashlib.sha256(raw).hexdigest()
    return reference
