"""Benchmark harness for candidate (re-derived) TCRen potentials.

A driver (see ``scripts/rederive_sweep.py``) derives one candidate TCRen CSV per
configuration and calls the functions here to score it on the manuscript's oracle
benchmarks. Every function takes a *candidate* potential — either an in-memory
contacts/markup pair (cognate-rank, fully self-contained) or a candidate CSV path that
is threaded into ``summarize_structure``/``alanine_scan`` as the TCR↔peptide potential —
and returns the headline metric(s) for that candidate.

Benchmarks
----------
* :func:`cognate_rank_auc` — leave-one-out cognate-epitope rank against ~1000
  anchor-preserving random decoys (notebook 02). Self-contained: needs only the cached
  ``contacts``/``markup`` and the non-redundant ``ids``.
* :func:`n199_r2` — refit ``sum_lj_coul ~ tcren + mj_hla_peptide + mj_cdr_hla`` over the
  ~187 regenerable structures of the manuscript ``lj_coul_tcren_mj.csv`` oracle, after
  regenerating ONLY the ``tcren`` column with the candidate (HC3 OLS).
* :func:`as_r2` — analogous refit of the ankylosing-spondylitis Fig-6 model
  ``Total_score_with_const ~ TCRen + MJ_HLA_peptide + MJ_TCR_HLA`` over the AS oracle.
* :func:`ddg_direct_r2` — ``ΔΔG ~ ΔTCRen(tcr_peptide)`` over the ATLAS structures, with
  ΔTCRen the summed per-position alanine scan on the **peptide-bearing** ``tcr_peptide``
  interface only (avoids the known ``tcr_mhc`` alanine-scan artifact).

Oracle CSVs live in the manuscript repo; their absolute paths are the module defaults.
``summarize_structure`` results are cached per candidate so the same candidate is never
re-scored twice within a sweep.
"""

from __future__ import annotations

import tempfile
from functools import lru_cache
from pathlib import Path

import numpy as np
import polars as pl

from tcren import summarize_structure
from tcren.ddg import alanine_scan
from tcren.paths import reference_structure_path
from tcren.potential import Potential, derive_tcren, derive_tcren_loo

# --- manuscript oracle CSVs (absolute, as supplied by the task) -----------------------
_ORACLE = Path("/Users/mikesh/vcs/manuscripts/2026-tcren2/data/oracle")
N199_CSV = _ORACLE / "plots" / "lj_coul_tcren_mj.csv"
AS_CSV = _ORACLE / "as" / "total_energy_summary.csv"
DDG_CSV = _ORACLE / "ddg_direct" / "ddg_direct_deltas.csv"

_AA = list("LFIMVWYCHAGPTSQNDERK")
_AIDX = {a: i for i, a in enumerate(_AA)}


# =====================================================================================
# (a) cognate-epitope rank — leave-one-out, anchor-preserving random decoys
# =====================================================================================
def _mat_from(df: pl.DataFrame, vcol: str) -> np.ndarray:
    """20x20 potential matrix (rows = ``residue.aa.from``, cols = ``residue.aa.to``)."""
    m = np.zeros((20, 20))
    for r in df.iter_rows(named=True):
        i, j = _AIDX.get(r["residue.aa.from"]), _AIDX.get(r["residue.aa.to"])
        if i is not None and j is not None:
            m[i, j] = r[vcol]
    return m


def cognate_rank_auc(
    contacts: pl.DataFrame,
    markup: pl.DataFrame,
    ids: list[str],
    n_decoy: int = 1000,
    seed: int = 0,
    contact_weight: str = "residue",
) -> dict:
    """Leave-one-out cognate-epitope rank vs anchor-preserving random decoys.

    For each non-redundant structure, derive the TCRen potential leaving that structure
    out (:func:`derive_tcren_loo`), thread the cognate peptide plus ``n_decoy`` random
    decoys (positions 2 and last pinned to the cognate anchors) onto the structure's
    TCR↔peptide contacts, and record the cognate's percentile rank (fraction of decoys
    scoring strictly lower). Lift of notebook ``02_benchmark_cognate_unrelated.ipynb``.

    Args:
        contacts: Cached TCR↔peptide contact table (from ``annotate_structure_set``);
            must carry ``pos.to``/``chain.type.from``/``residue.aa.from`` (and
            ``n_atom_contacts`` when ``contact_weight="atomic"``).
        markup: Cached per-structure markup with ``pdb.id``/``peptide``.
        ids: Non-redundant αβ ``pdb.id`` set (the LOO inclusion set).
        n_decoy: Number of random decoys per structure (default 1000).
        seed: RNG seed.
        contact_weight: ``"residue"`` (default, each contact counts once) or ``"atomic"``
            (each contact weighted by its ``n_atom_contacts`` heavy-atom-pair count).

    Returns:
        ``{"median_rank_pct", "rank_auc", "n"}`` — median cognate rank (%), rank-based
        AUC (``mean(1 - rank/100)``), and the number of scored structures.
    """
    if contact_weight not in ("residue", "atomic"):
        raise ValueError(f"contact_weight must be 'residue' or 'atomic', got {contact_weight!r}")
    if contact_weight == "atomic" and "n_atom_contacts" not in contacts.columns:
        raise ValueError(
            "contact_weight='atomic' needs an n_atom_contacts column on the contacts table"
        )
    rng = np.random.default_rng(seed)
    ids = sorted(ids)
    loo = derive_tcren_loo(contacts, ids)
    loo_mat = {
        pid: _mat_from(loo.filter(pl.col("pdb.id") == pid), "TCRen.LOO") for pid in ids
    }
    peptide_of = {
        r["pdb.id"]: r["peptide"]
        for r in markup.iter_rows(named=True)
        if r.get("peptide")
    }
    ab = contacts.filter(pl.col("chain.type.from").is_in(["TRA", "TRB"]))

    ranks: list[float] = []
    for pid in ids:
        cog = peptide_of.get(pid)
        if not cog or len(cog) < 4 or any(a not in _AIDX for a in cog):
            continue
        m = loo_mat.get(pid)
        if m is None:
            continue
        sub = ab.filter(pl.col("pdb.id") == pid)
        pos = np.array(sub["pos.to"].to_list())
        tcr = np.array([_AIDX.get(a, -1) for a in sub["residue.aa.from"].to_list()])
        if contact_weight == "atomic":
            wgt = np.array(sub["n_atom_contacts"].to_list(), dtype=np.float64)
        else:
            wgt = np.ones(len(pos), dtype=np.float64)
        keep = (pos >= 0) & (pos < len(cog)) & (tcr >= 0)
        pos, tcr, wgt = pos[keep], tcr[keep], wgt[keep]
        if len(pos) == 0:
            continue
        L = len(cog)
        cogv = np.array([_AIDX[a] for a in cog])
        dec = rng.integers(0, 20, size=(n_decoy, L))
        dec[:, 1] = cogv[1]
        dec[:, L - 1] = cogv[L - 1]
        allp = np.vstack([cogv[None, :], dec])  # (N+1, L), row 0 = cognate
        contact_aa = allp[:, pos]  # (N+1, K)
        sc = (m[tcr[None, :], contact_aa] * wgt[None, :]).sum(axis=1)  # (N+1,)
        ranks.append(float((sc[1:] < sc[0]).mean() * 100))

    if not ranks:
        return {"median_rank_pct": float("nan"), "rank_auc": float("nan"), "n": 0}
    return {
        "median_rank_pct": float(np.median(ranks)),
        "rank_auc": float(1 - np.mean(ranks) / 100),
        "n": len(ranks),
    }


# =====================================================================================
# regenerate the candidate's per-structure tcr_peptide energy (cached per candidate)
# =====================================================================================
@lru_cache(maxsize=None)
def _tcr_peptide_energy(
    structure_name: str, candidate_csv: str, contact_weight: str = "residue"
) -> float | None:
    """tcr_peptide energy of a structure under a candidate TCRen potential, or ``None``.

    ``None`` when the structure is not a deposited Native2026 PDB or annotation fails.
    Cached per ``(structure_name, candidate_csv, contact_weight)`` so a candidate is
    scored once. ``contact_weight="atomic"`` weights each contact by its heavy-atom-pair
    count (the substitution is byte-equal to the legacy energy under ``"residue"``).
    """
    try:
        pdb = reference_structure_path(structure_name)
    except FileNotFoundError:
        return None
    try:
        out = summarize_structure(
            pdb,
            superimpose=False,
            potentials={"tcr_peptide": candidate_csv},
            background=10,
            contact_weight=contact_weight,
        )
    except Exception:
        return None
    return float(out["scores"]["tcr_peptide"][0])


def _ols_r2(df: pl.DataFrame, response: str, predictors: list[str]) -> dict:
    """HC3-robust OLS; return R², n, and per-predictor sign + significance (p<0.05)."""
    import statsmodels.api as sm

    sub = df.select([response, *predictors]).drop_nulls()
    y = sub[response].to_numpy()
    X = sm.add_constant(sub.select(predictors).to_numpy())
    fit = sm.OLS(y, X).fit(cov_type="HC3")
    # params/pvalues are [const, *predictors]
    coefs = {
        name: {
            "coef": float(fit.params[i + 1]),
            "sign": "+" if fit.params[i + 1] >= 0 else "-",
            "p": float(fit.pvalues[i + 1]),
            "significant": bool(fit.pvalues[i + 1] < 0.05),
        }
        for i, name in enumerate(predictors)
    }
    return {"r2": float(fit.rsquared), "n": sub.height, "coefficients": coefs}


# =====================================================================================
# (b) n=199 ergodicity refit — regenerate ONLY the tcren column per candidate
# =====================================================================================
def _loo_potential_csv(
    contacts: pl.DataFrame,
    derivation_ids: list[str],
    leave_out: str,
    workdir: Path,
    **derive_kwargs,
) -> str:
    """Derive a leave-one-out TCRen potential (excluding ``leave_out``) and write it.

    Derives from ``derivation_ids`` minus ``leave_out`` via :func:`derive_tcren`, writes a
    wide ``residue.aa.from, residue.aa.to, TCRen`` CSV (the layout ``Potential.from_csv``
    reads) under ``workdir``, and returns its path.
    """
    include = [i for i in derivation_ids if i != leave_out]
    pot = derive_tcren(contacts, include=include, **derive_kwargs)
    out = pot.matrix.rename({"value": "TCRen"})
    path = workdir / f"tcren_loo_{leave_out}.csv"
    out.write_csv(str(path))
    return str(path)


def n199_r2(
    candidate_csv: str | Path | None = None,
    oracle_csv: str | Path = N199_CSV,
    *,
    contact_weight: str = "residue",
    loo: bool = False,
    contacts: pl.DataFrame | None = None,
    derivation_ids: list[str] | None = None,
    **derive_kwargs,
) -> dict:
    """Refit ``sum_lj_coul ~ tcren + mj_hla_peptide + mj_cdr_hla`` for a candidate.

    Two derivation modes share the same refit:

    * **non-LOO** (default) — regenerate the ``tcren`` column with a single fixed
      ``candidate_csv`` (the candidate TCRen as the ``tcr_peptide`` potential) over the
      regenerable structures of the ``lj_coul_tcren_mj.csv`` oracle.
    * **LOO** (``loo=True``) — for each regenerable structure, derive a TCRen potential
      from ``derivation_ids`` *excluding that structure* (``derive_tcren(contacts,
      include=derivation_ids \\ {id})``), score the structure with that held-out
      potential, then refit. Requires ``contacts`` and ``derivation_ids``.

    Either mode keeps the oracle's physical response and MJ predictors and refits HC3 OLS.

    Args:
        candidate_csv: Fixed candidate TCRen CSV (non-LOO mode). Ignored when ``loo``.
        oracle_csv: The ``lj_coul_tcren_mj.csv`` oracle (default :data:`N199_CSV`).
        contact_weight: ``"residue"`` (default) or ``"atomic"``; threaded into scoring so
            the regenerated ``tcren`` column can use atomic weighting.
        loo: When ``True`` derive a held-out potential per structure (leave-one-out).
        contacts: Annotated contact table fed to :func:`derive_tcren` (LOO only).
        derivation_ids: Derivation inclusion set; each structure is left out in turn
            (LOO only).
        **derive_kwargs: Extra keyword args forwarded to :func:`derive_tcren` (LOO only),
            e.g. ``variant``, ``weights``.

    Returns:
        ``{"r2", "n", "coefficients"}`` where ``coefficients`` maps each predictor to its
        ``coef``/``sign``/``p``/``significant``.
    """
    # Validate arguments before any filesystem access so misuse raises deterministically
    # and the unit tests do not depend on the manuscript oracle CSV (N199_CSV) being present.
    if loo:
        if contacts is None or derivation_ids is None:
            raise ValueError("loo=True requires both contacts and derivation_ids")
    elif candidate_csv is None:
        raise ValueError("non-LOO mode requires candidate_csv")

    oracle = pl.read_csv(oracle_csv)
    rows = []

    if loo:
        with tempfile.TemporaryDirectory() as tmp:
            workdir = Path(tmp)
            for r in oracle.iter_rows(named=True):
                name = r["structure_name"]
                csv = _loo_potential_csv(
                    contacts, derivation_ids, name, workdir, **derive_kwargs
                )
                e = _tcr_peptide_energy(name, csv, contact_weight)
                if e is None:
                    continue
                rows.append({
                    "sum_lj_coul": r["sum_lj_coul"],
                    "tcren": e,
                    "mj_hla_peptide": r["mj_hla_peptide"],
                    "mj_cdr_hla": r["mj_cdr_hla"],
                })
    else:
        candidate_csv = str(candidate_csv)
        for r in oracle.iter_rows(named=True):
            e = _tcr_peptide_energy(r["structure_name"], candidate_csv, contact_weight)
            if e is None:
                continue
            rows.append({
                "sum_lj_coul": r["sum_lj_coul"],
                "tcren": e,
                "mj_hla_peptide": r["mj_hla_peptide"],
                "mj_cdr_hla": r["mj_cdr_hla"],
            })

    df = pl.DataFrame(rows)
    return _ols_r2(
        df, "sum_lj_coul", ["tcren", "mj_hla_peptide", "mj_cdr_hla"]
    )


# =====================================================================================
# (c) AS Fig-6 refit — regenerate the TCRen column per candidate (n up to 14)
# =====================================================================================
def as_r2(candidate_csv: str | Path, oracle_csv: str | Path = AS_CSV) -> dict:
    """Refit ``Total_score_with_const ~ TCRen + MJ_HLA_peptide + MJ_TCR_HLA``.

    Regenerates the ``TCRen`` predictor with the candidate over the regenerable AS
    structures (only the deposited Native2026 PDBs of the 14-row oracle resolve; the
    FoldX ``_2`` variants and ``non*`` MD decoys are not in Native2026), keeps the
    oracle's response and MJ predictors, and refits HC3 OLS.

    Returns:
        ``{"r2", "n", "coefficients"}``. ``n`` reflects how many AS rows were
        regenerable.
    """
    candidate_csv = str(candidate_csv)
    oracle = pl.read_csv(oracle_csv)
    rows = []
    for r in oracle.iter_rows(named=True):
        e = _tcr_peptide_energy(r["Structure"], candidate_csv)
        if e is None:
            continue
        rows.append({
            "Total_score_with_const": r["Total_score_with_const"],
            "TCRen": e,
            "MJ_HLA_peptide": r["MJ_HLA_peptide"],
            "MJ_TCR_HLA": r["MJ_TCR_HLA"],
        })
    df = pl.DataFrame(rows)
    return _ols_r2(
        df, "Total_score_with_const", ["TCRen", "MJ_HLA_peptide", "MJ_TCR_HLA"]
    )


# =====================================================================================
# (d) direct ΔΔG — ΔΔG ~ ΔTCRen(tcr_peptide) via alanine scan over ATLAS structures
# =====================================================================================
@lru_cache(maxsize=None)
def _delta_tcren_tcr_peptide(structure_name: str, candidate_csv: str) -> float | None:
    """Summed per-position ΔTCRen of the alanine scan on the tcr_peptide interface.

    Mutates each peptide position to Ala and sums ``E(native) - E(Ala@pos)`` on the
    **peptide-bearing** ``tcr_peptide`` interface only (the ``tcr_mhc`` interface holds
    no peptide residues, so scanning it is the known artifact; we never touch it).
    Uses the structure's own native peptide. ``None`` if unresolved/annotation fails.
    """
    try:
        pdb = reference_structure_path(structure_name)
    except FileNotFoundError:
        return None
    try:
        from tcren.contactmap import ContactMap
        from tcren.annotation import classify_chains
        from tcren.structure import parse_structure
        from tcren.structure.model import PEPTIDE_TYPE

        s = parse_structure(pdb, pdb_id=structure_name)
        classify_chains(s, organism="human", autodetect_species=True)
        native = next((c.sequence() for c in s.chains if c.chain_type == PEPTIDE_TYPE), None)
        if not native:
            return None
        cm = ContactMap.from_structure(s)
        pot = Potential.from_csv(candidate_csv)
        scan = alanine_scan(cm, native, pot, interface="tcr_peptide")
    except Exception:
        return None
    return float(scan["ddG"].sum())


def ddg_direct_r2(candidate_csv: str | Path, oracle_csv: str | Path = DDG_CSV) -> dict:
    """Refit ``ΔΔG ~ ΔTCRen(tcr_peptide)`` for a candidate over the ATLAS structures.

    The predictor is the summed alanine-scan ΔTCRen on the ``tcr_peptide`` interface
    (peptide is the substituted side; the ``tcr_mhc`` artifact is avoided). The response
    is the oracle ``ddG`` (experimental ATLAS ΔΔG). HC3 OLS.

    Returns:
        ``{"r2", "n", "coefficients"}`` with the single ``d_tcr_pep`` predictor.
    """
    candidate_csv = str(candidate_csv)
    oracle = pl.read_csv(oracle_csv)
    rows = []
    for r in oracle.iter_rows(named=True):
        d = _delta_tcren_tcr_peptide(r["structure_name"], candidate_csv)
        if d is None:
            continue
        rows.append({"ddG": r["ddG"], "d_tcr_pep": d})
    df = pl.DataFrame(rows)
    return _ols_r2(df, "ddG", ["d_tcr_pep"])
