"""Size-conditioned contact action under the frozen native multinomial."""

from __future__ import annotations

import numpy as np
import polars as pl
from scipy.special import gammaln

from .contacts.opportunities import native_contact_table
from .native_reference import load_native_reference
from .structure.model import Structure

DRAWS = 2000
SEED = 20261123


def native_contact_observations(
    structure: Structure, *, alphabet: str, scope: str,
) -> pl.DataFrame:
    """Extract directed 5 Å observations with the Native387 anchor requirement."""
    if alphabet not in {"atom", "residue"} or scope not in {"complex", "tcr_peptide"}:
        raise ValueError(f"unsupported native contact profile {alphabet}/{scope}")
    contacts = native_contact_table(
        structure, atom_pairs=alphabet == "atom",
        interface="all" if scope == "complex" else "tcr_peptide",
    )
    prefix = "atom" if alphabet == "atom" else "residue.aa"
    cells = pl.concat_str(f"{prefix}.from", f"{prefix}.to", separator="|")
    if scope == "complex":
        cells = pl.concat_str("interface", cells, separator="/")
    return contacts.select(pl.col("pdb.id").alias("id"), cells.alias("cell"))


def _null_moments(q: np.ndarray, sizes: set[int]) -> dict[int, tuple[float, float, float, float]]:
    """Extend the same simulated contact maps across sizes, recording requested moments."""
    if not sizes:
        return {}
    rng = np.random.default_rng(SEED)
    logq = np.log(q)
    normalized_q = q / q.sum()
    counts = np.zeros((DRAWS, len(q)), dtype=np.int32)
    rows = np.arange(DRAWS)
    surprise = np.zeros(DRAWS)
    log_factorials = np.zeros(DRAWS)
    moments = {}
    for size in range(1, max(sizes) + 1):
        cell = rng.choice(len(q), size=DRAWS, p=normalized_q)
        counts[rows, cell] += 1
        log_factorials += np.log(counts[rows, cell])
        surprise -= logq[cell]
        if size in sizes:
            action = -gammaln(size + 1.0) + log_factorials + surprise
            mean, sd = float(action.mean()), float(action.std(ddof=1))
            if sd > 0:
                centered = (action - mean) / sd
                second = centered**2
                moments[size] = (mean, sd, np.mean(centered * second), np.var(second))
    return moments


def contact_typicality(
    observations: pl.DataFrame, ids: list[str], *, alphabet: str = "atom",
    scope: str = "tcr_peptide",
) -> pl.DataFrame:
    """Report action, its conditional z-score, and simulation MCSE for each ID.

    Supply one row per observation with string columns ``id`` and ``cell``. Cells
    use directional reference keys, prefixed with the interface in complex scope.
    Include IDs without observations to retain their rows. Unknown cells, empty
    contact sets and degenerate nulls carry a reason and a null z-score.
    """
    reference = load_native_reference()
    profile = f"{alphabet}/{scope}"
    if profile not in reference["contact_profiles"]:
        raise ValueError(f"unsupported native contact profile {profile!r}")
    probability = reference["contact_profiles"][profile]["probability"]
    cells = sorted(probability)
    q = np.asarray([probability[cell] for cell in cells])
    # The stored probabilities keep their full precision. Normalize only the
    # simulation input for NumPy's categorical sampler.
    ids_table = pl.DataFrame({"id": ids}, schema={"id": pl.String})
    if len(set(ids)) != len(ids) or observations.join(ids_table, on="id", how="anti").height:
        raise ValueError("declare unique structure IDs including every observation")
    vocabulary = pl.DataFrame({"cell": cells, "logq": np.log(q)})
    counts = (observations.group_by("id", "cell").len().sort("id", "cell")
              .join(vocabulary, on="cell", how="left", maintain_order="left"))
    counts = counts.with_columns(
        pl.Series("log_factorial", gammaln(counts["len"].to_numpy() + 1.0)))
    grouped = counts.group_by("id").agg(
        pl.col("len").cast(pl.Int64).sum().alias("contact_count"),
        pl.col("logq").is_null().any().alias("unknown"),
        (pl.col("log_factorial") - pl.col("len") * pl.col("logq")).sum().alias("partial"),
    )
    table = ids_table.join(grouped, on="id", how="left", maintain_order="left")
    table = table.with_columns(pl.col("contact_count").fill_null(0))
    sizes = set(table.filter(~pl.col("unknown"))["contact_count"].to_list())
    moments = _null_moments(q, sizes)
    actions, scores, errors, reasons = [], [], [], []
    for size, unknown, partial in table.select("contact_count", "unknown", "partial").iter_rows():
        action = score = error = reason = None
        if size == 0:
            reason = "no_contacts"
        elif unknown:
            reason = "unknown_cell"
        else:
            action = float(partial - gammaln(size + 1.0))
            if size in moments:
                mean, sd, third, variance_second = moments[size]
                score = (action - mean) / sd
                # Delta method for (A - mean) / SD, including covariance of
                # the estimated mean and variance from the same simulations.
                error = np.sqrt(max(0.0, 1 + score * third + score**2 * variance_second / 4)
                                / DRAWS)
            else:
                reason = "degenerate_null"
        actions.append(action)
        scores.append(score)
        errors.append(error)
        reasons.append(reason)
    return table.select("id", "contact_count").with_columns(
        pl.Series("contact_action", actions, dtype=pl.Float64),
        pl.Series("contact_typicality", scores, dtype=pl.Float64),
        pl.Series("contact_typicality_mcse", errors, dtype=pl.Float64),
        pl.Series("typicality_reason", reasons, dtype=pl.String),
        pl.lit(profile).alias("contact_profile"),
        pl.lit(reference["version"]).alias("reference_version"),
        pl.lit(reference["reference_sha256"]).alias("reference_sha256"),
    )
