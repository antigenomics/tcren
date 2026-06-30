"""Full re-derivation sweep over the exact configs supplied by the orchestration task.

Loads the *pre-cached* contacts/markup at ``scratch/cache/contacts_<set>.csv`` and
``scratch/cache/markup_<set>.csv`` (NOT re-annotated). The ``union`` set is the concat of
Native2022 + Native2026 deduped by ``pdb.id`` (each id's rows taken from the first source
that has it). For each config:

    ids = nonredundant_ids(markup, t) ∩ alphabeta_ids(contacts)
    pot = derive_tcren(contacts, include=ids, variant=..., pseudocount=...)
    write scratch/cache/TCRen_<label>.csv
    bench: cognate_rank_auc (LOO), n199_r2 (+sign), as_r2, ddg_direct_r2

Also benchmarks the paper's legacy potential and the current v2 default through the same
harness for the baseline rows. Writes a JSON metrics array to stdout and to
``scratch/cache/rederive_full_metrics.json``.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parent))
import bench_harness as bh  # noqa: E402

from tcren.potential import (  # noqa: E402
    alphabeta_ids,
    derive_tcren,
    nonredundant_ids,
)

_REPO = Path(__file__).resolve().parents[1]
_CACHE = _REPO / "scratch" / "cache"

LEGACY_CSV = _REPO / "tests" / "assets" / "oracle" / "data" / "TCRen_potential.csv"
V2_DEFAULT_CSV = _REPO / "src" / "tcren" / "data" / "TCRen_potential.csv"

CONFIGS = [
    {"set": "Native2022", "t": 6, "variant": "classic", "pseudocount": 1, "label": "legacy-control"},
    {"set": "Native2026", "t": None, "variant": "classic", "pseudocount": 1, "label": "2026-off"},
    {"set": "Native2026", "t": 3, "variant": "classic", "pseudocount": 1, "label": "2026-t3"},
    {"set": "Native2026", "t": 6, "variant": "classic", "pseudocount": 1, "label": "2026-t6-current"},
    {"set": "Native2026", "t": 9, "variant": "classic", "pseudocount": 1, "label": "2026-t9"},
    {"set": "Native2026", "t": 12, "variant": "classic", "pseudocount": 1, "label": "2026-t12"},
    {"set": "union", "t": 6, "variant": "classic", "pseudocount": 1, "label": "union-t6"},
    {"set": "union", "t": None, "variant": "classic", "pseudocount": 1, "label": "union-off"},
    {"set": "Native2026", "t": 6, "variant": "classic", "pseudocount": 0.5, "label": "2026-t6-pc0.5"},
    {"set": "Native2026", "t": 6, "variant": "am", "pseudocount": 1, "label": "2026-t6-am"},
]


def load_set(name: str) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Return (contacts, markup) for a set name. ``union`` = concat dedup by pdb.id."""
    if name in ("Native2022", "Native2026"):
        c = pl.read_csv(_CACHE / f"contacts_{name}.csv")
        m = pl.read_csv(_CACHE / f"markup_{name}.csv")
        return c, m
    if name == "union":
        c22 = pl.read_csv(_CACHE / "contacts_Native2022.csv")
        c26 = pl.read_csv(_CACHE / "contacts_Native2026.csv")
        m22 = pl.read_csv(_CACHE / "markup_Native2022.csv")
        m26 = pl.read_csv(_CACHE / "markup_Native2026.csv")
        ids22 = set(c22["pdb.id"].unique().to_list())
        c = pl.concat([c22, c26.filter(~pl.col("pdb.id").is_in(list(ids22)))])
        mids22 = set(m22["pdb.id"].unique().to_list())
        m = pl.concat([m22, m26.filter(~pl.col("pdb.id").is_in(list(mids22)))])
        return c, m
    raise ValueError(name)


def bench_candidate(label: str, cand_csv: Path, contacts: pl.DataFrame,
                    markup: pl.DataFrame, include: list[str]) -> dict:
    """Run all four benchmarks on a candidate CSV; return a metrics row."""
    cr = bh.cognate_rank_auc(contacts, markup, include)
    n199 = bh.n199_r2(cand_csv)
    asr = bh.as_r2(cand_csv)
    ddg = bh.ddg_direct_r2(cand_csv)
    tc = n199["coefficients"]["tcren"]
    # n=199 sign check: physical sum_lj_coul is negative-good; expect tcren coef +
    # (tcren correlates with the LJ+Coul energy). Record the actual sign + significance.
    return {
        "label": label,
        "cognate_median_rank_pct": cr["median_rank_pct"],
        "cognate_rank_auc": cr["rank_auc"],
        "cognate_n": cr["n"],
        "n199_r2": n199["r2"],
        "n199_n": n199["n"],
        "n199_tcren_sign": tc["sign"],
        "n199_tcren_sig": tc["significant"],
        "as_r2": asr["r2"],
        "as_n": asr["n"],
        "as_tcren_sign": asr["coefficients"]["TCRen"]["sign"],
        "ddg_direct_r2": ddg["r2"],
        "ddg_n": ddg["n"],
        "ddg_sign": ddg["coefficients"]["d_tcr_pep"]["sign"],
    }


def main() -> None:
    results: list[dict] = []

    # ---- baselines: legacy table + current v2 default through the SAME harness ----
    for label, csv in [("baseline-legacy-table", LEGACY_CSV),
                       ("baseline-v2-default", V2_DEFAULT_CSV)]:
        print(f"[baseline] {label} <- {csv}", flush=True)
        n199 = bh.n199_r2(csv)
        asr = bh.as_r2(csv)
        ddg = bh.ddg_direct_r2(csv)
        results.append({
            "label": label, "set": "(table)", "t": None, "variant": "(table)",
            "pseudocount": None, "n_structures": None, "n_clusters": None,
            "cognate_median_rank_pct": None, "cognate_rank_auc": None, "cognate_n": None,
            "n199_r2": n199["r2"], "n199_n": n199["n"],
            "n199_tcren_sign": n199["coefficients"]["tcren"]["sign"],
            "n199_tcren_sig": n199["coefficients"]["tcren"]["significant"],
            "as_r2": asr["r2"], "as_n": asr["n"],
            "as_tcren_sign": asr["coefficients"]["TCRen"]["sign"],
            "ddg_direct_r2": ddg["r2"], "ddg_n": ddg["n"],
            "ddg_sign": ddg["coefficients"]["d_tcr_pep"]["sign"],
        })
        print(f"    n199_r2={n199['r2']:.4f} as_r2={asr['r2']:.4f} ddg_r2={ddg['r2']:.4f}",
              flush=True)

    # ---- candidate configs ----
    set_cache: dict[str, tuple[pl.DataFrame, pl.DataFrame]] = {}
    for cfg in CONFIGS:
        sname = cfg["set"]
        if sname not in set_cache:
            set_cache[sname] = load_set(sname)
        contacts, markup = set_cache[sname]

        ab = alphabeta_ids(contacts)
        nr = nonredundant_ids(markup.filter(pl.col("pdb.id").is_in(ab)), t=cfg["t"])
        include = sorted(set(nr) & set(ab))

        pot = derive_tcren(contacts, include=include,
                           variant=cfg["variant"], pseudocount=cfg["pseudocount"])
        cand_csv = _CACHE / f"TCRen_{cfg['label']}.csv"
        pot.to_csv(cand_csv)

        print(f"[{cfg['label']}] set={sname} t={cfg['t']} variant={cfg['variant']} "
              f"pc={cfg['pseudocount']} | n_ab={len(ab)} n_clusters={len(include)}",
              flush=True)

        row = bench_candidate(cfg["label"], cand_csv, contacts, markup, include)
        row.update({
            "set": sname, "t": cfg["t"], "variant": cfg["variant"],
            "pseudocount": cfg["pseudocount"],
            "n_structures": len(ab), "n_clusters": len(include),
            "candidate_csv": str(cand_csv),
        })
        results.append(row)
        print(f"    cog_auc={row['cognate_rank_auc']:.4f} med_rank={row['cognate_median_rank_pct']:.2f} "
              f"n199_r2={row['n199_r2']:.4f}({row['n199_tcren_sign']},sig={row['n199_tcren_sig']}) "
              f"as_r2={row['as_r2']:.4f}(n={row['as_n']}) ddg_r2={row['ddg_direct_r2']:.4f}(n={row['ddg_n']})",
              flush=True)

    out = _CACHE / "rederive_full_metrics.json"
    out.write_text(json.dumps(results, indent=2))
    print(f"\nwrote {len(results)} rows -> {out}", flush=True)
    print("===JSON_BEGIN===")
    print(json.dumps(results))
    print("===JSON_END===")


if __name__ == "__main__":
    main()
