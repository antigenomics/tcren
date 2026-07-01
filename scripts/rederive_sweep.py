"""Re-derivation sweep driver.

For each ``(set, redundancy_t, variant, pseudocount)`` configuration: cluster the cached
contacts/markup of the set, derive a candidate TCRen potential, run the oracle
benchmarks in :mod:`bench_harness`, and append one metrics row to
``scratch/rederive_metrics.csv``.

Contacts are *t-independent*: ``annotate_structure_set`` is run once per set and cached
to disk (``scratch/contacts_<set>.csv`` / ``scratch/markup_<set>.csv``), so each config
only re-clusters + re-derives (sub-second). The structure-level benchmarks (n199/AS/ΔΔG)
regenerate the candidate's ``tcren`` column via ``summarize_structure`` and cache per
candidate CSV inside ``bench_harness``.

This driver is intentionally a thin orchestration layer; it does NOT touch the bundled
default potential. Usage::

    conda run -n tcren-nb python scripts/rederive_sweep.py \
        --set 2026 --t 6.0 --t none --variant classic --pseudocount 1

Defaults sweep the 2026 set at ``t in {None, 6.0}``, ``variant=classic``,
``pseudocount=1`` and run all four benchmarks.
"""

from __future__ import annotations

import argparse
import itertools
from pathlib import Path

import polars as pl

import bench_harness as bh
from tcren.paper import annotate_structure_set
from tcren.potential import alphabeta_ids, derive_tcren, nonredundant_ids

_REPO = Path(__file__).resolve().parents[1]
_SCRATCH = _REPO / "scratch"
_SETS = {"2022": _REPO / "data" / "Native2022", "2026": _REPO / "data" / "Native2026"}
_METRICS = _SCRATCH / "rederive_metrics.csv"


def _cached_contacts(set_name: str) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Contacts + markup for a set, assembled once and cached to ``scratch/``."""
    _SCRATCH.mkdir(exist_ok=True)
    cp, mp = _SCRATCH / f"contacts_{set_name}.csv", _SCRATCH / f"markup_{set_name}.csv"
    if cp.exists() and mp.exists():
        return pl.read_csv(cp), pl.read_csv(mp)
    contacts, markup = annotate_structure_set(_SETS[set_name])
    contacts.write_csv(cp)
    markup.write_csv(mp)
    return contacts, markup


def _parse_t(value: str) -> float | None:
    """Parse a ``--t`` value: ``none``/``off`` -> ``None`` (redundancy off), else float."""
    return None if value.lower() in ("none", "off") else float(value)


def run_config(
    set_name: str,
    contacts: pl.DataFrame,
    markup: pl.DataFrame,
    t: float | None,
    variant: str,
    pseudocount: int,
    benchmarks: set[str],
) -> dict:
    """Derive one candidate and run the requested benchmarks; return a metrics row."""
    ab = alphabeta_ids(contacts)
    include = nonredundant_ids(markup.filter(pl.col("pdb.id").is_in(ab)), t=t)

    pot = derive_tcren(contacts, include=include, variant=variant, pseudocount=pseudocount)
    tag = f"{set_name}_t{'off' if t is None else t}_{variant}_pc{pseudocount}"
    cand_csv = _SCRATCH / f"candidate_{tag}.csv"
    pot.to_csv(cand_csv)

    row: dict = {
        "set": set_name,
        "redundancy_t": "off" if t is None else t,
        "variant": variant,
        "pseudocount": pseudocount,
        "n_alphabeta": len(ab),
        "n_nonredundant": len(include),
        "candidate_csv": str(cand_csv),
    }

    if "cognate" in benchmarks:
        cr = bh.cognate_rank_auc(contacts, markup, include)
        row["cognate_median_rank_pct"] = cr["median_rank_pct"]
        row["cognate_rank_auc"] = cr["rank_auc"]
        row["cognate_n"] = cr["n"]

    if "n199" in benchmarks:
        r = bh.n199_r2(cand_csv)
        row["n199_r2"] = r["r2"]
        row["n199_n"] = r["n"]
        row["n199_tcren_sign"] = r["coefficients"]["tcren"]["sign"]
        row["n199_tcren_sig"] = r["coefficients"]["tcren"]["significant"]

    if "as" in benchmarks:
        r = bh.as_r2(cand_csv)
        row["as_r2"] = r["r2"]
        row["as_n"] = r["n"]
        row["as_tcren_sign"] = r["coefficients"]["TCRen"]["sign"]
        row["as_tcren_sig"] = r["coefficients"]["TCRen"]["significant"]

    if "ddg" in benchmarks:
        r = bh.ddg_direct_r2(cand_csv)
        row["ddg_r2"] = r["r2"]
        row["ddg_n"] = r["n"]
        row["ddg_tcren_sign"] = r["coefficients"]["d_tcr_pep"]["sign"]

    return row


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--set", dest="sets", action="append", choices=sorted(_SETS),
                    help="structure set(s); repeat (default: 2026)")
    ap.add_argument("--t", dest="ts", action="append",
                    help="redundancy cutoff(s); 'none' = off; repeat (default: none, 6.0)")
    ap.add_argument("--variant", dest="variants", action="append", choices=["classic", "am"],
                    help="derivation variant(s); repeat (default: classic)")
    ap.add_argument("--pseudocount", dest="pseudocounts", action="append", type=int,
                    help="pseudocount(s); repeat (default: 1)")
    ap.add_argument("--benchmarks", default="cognate,n199,as,ddg",
                    help="comma-separated subset of cognate,n199,as,ddg")
    args = ap.parse_args()

    sets = args.sets or ["2026"]
    ts = [_parse_t(t) for t in (args.ts or ["none", "6.0"])]
    variants = args.variants or ["classic"]
    pseudocounts = args.pseudocounts or [1]
    benchmarks = {b.strip() for b in args.benchmarks.split(",") if b.strip()}

    _SCRATCH.mkdir(exist_ok=True)
    metrics_rows: list[dict] = []
    for set_name in sets:
        contacts, markup = _cached_contacts(set_name)
        print(f"[{set_name}] {contacts.height} contacts over "
              f"{contacts['pdb.id'].n_unique()} structures")
        for t, variant, pc in itertools.product(ts, variants, pseudocounts):
            print(f"  config t={t} variant={variant} pc={pc} ...", flush=True)
            row = run_config(set_name, contacts, markup, t, variant, pc, benchmarks)
            metrics_rows.append(row)
            headline = {k: v for k, v in row.items() if k.endswith(("_r2", "_auc"))}
            print(f"    -> {headline}")

    out = pl.DataFrame(metrics_rows)
    # append to any existing metrics, aligning columns
    if _METRICS.exists():
        prev = pl.read_csv(_METRICS)
        out = pl.concat([prev, out], how="diagonal_relaxed")
    out.write_csv(_METRICS)
    print(f"wrote {len(metrics_rows)} new row(s) -> {_METRICS}")


if __name__ == "__main__":
    main()
