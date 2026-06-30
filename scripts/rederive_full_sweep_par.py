"""Parallel full re-derivation sweep — same configs/benchmarks as rederive_full_sweep,
but the per-structure summarize_structure / alanine-scan passes (n199, ddg) are mapped
across a ProcessPoolExecutor. Cheap steps (cognate LOO, OLS) run in the main process.

The two heavy per-structure workers live in bench_harness at module level and are
picklable. We bypass bench_harness.n199_r2 / ddg_direct_r2's internal serial loops and
instead fan their worker calls out to the pool, then reuse bench_harness._ols_r2 for the
identical regression.
"""

from __future__ import annotations

import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parent))
import bench_harness as bh  # noqa: E402

from tcren.potential import alphabeta_ids, derive_tcren, nonredundant_ids  # noqa: E402

_REPO = Path(__file__).resolve().parents[1]
_CACHE = _REPO / "scratch" / "cache"
LEGACY_CSV = _REPO / "tests" / "assets" / "oracle" / "data" / "TCRen_potential.csv"
V2_DEFAULT_CSV = _REPO / "src" / "tcren" / "data" / "TCRen_potential.csv"
N_WORKERS = min(12, (os.cpu_count() or 4))

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


# --- module-level pool workers (picklable) ------------------------------------------
def _n199_energy(args):
    name, cand = args
    return name, bh._tcr_peptide_energy(name, cand)


def _ddg_delta(args):
    name, cand = args
    return name, bh._delta_tcren_tcr_peptide(name, cand)


def load_set(name: str):
    if name in ("Native2022", "Native2026"):
        return (pl.read_csv(_CACHE / f"contacts_{name}.csv"),
                pl.read_csv(_CACHE / f"markup_{name}.csv"))
    if name == "union":
        c22 = pl.read_csv(_CACHE / "contacts_Native2022.csv")
        c26 = pl.read_csv(_CACHE / "contacts_Native2026.csv")
        m22 = pl.read_csv(_CACHE / "markup_Native2022.csv")
        m26 = pl.read_csv(_CACHE / "markup_Native2026.csv")
        ids22 = set(c22["pdb.id"].unique().to_list())
        c = pl.concat([c22, c26.filter(~pl.col("pdb.id").is_in(list(ids22)))])
        mids = set(m22["pdb.id"].unique().to_list())
        m = pl.concat([m22, m26.filter(~pl.col("pdb.id").is_in(list(mids)))])
        return c, m
    raise ValueError(name)


def _n199_par(ex, cand_csv: str) -> dict:
    oracle = pl.read_csv(bh.N199_CSV)
    names = oracle["structure_name"].to_list()
    energy = dict(ex.map(_n199_energy, [(n, cand_csv) for n in names]))
    rows = []
    for r in oracle.iter_rows(named=True):
        e = energy.get(r["structure_name"])
        if e is None:
            continue
        rows.append({"sum_lj_coul": r["sum_lj_coul"], "tcren": e,
                     "mj_hla_peptide": r["mj_hla_peptide"], "mj_cdr_hla": r["mj_cdr_hla"]})
    return bh._ols_r2(pl.DataFrame(rows), "sum_lj_coul",
                      ["tcren", "mj_hla_peptide", "mj_cdr_hla"])


def _ddg_par(ex, cand_csv: str) -> dict:
    oracle = pl.read_csv(bh.DDG_CSV)
    names = oracle["structure_name"].to_list()
    delta = dict(ex.map(_ddg_delta, [(n, cand_csv) for n in names]))
    rows = []
    for r in oracle.iter_rows(named=True):
        d = delta.get(r["structure_name"])
        if d is None:
            continue
        rows.append({"ddG": r["ddG"], "d_tcr_pep": d})
    return bh._ols_r2(pl.DataFrame(rows), "ddG", ["d_tcr_pep"])


def main() -> None:
    results: list[dict] = []
    with ProcessPoolExecutor(max_workers=N_WORKERS) as ex:
        # ---- baselines ----
        for label, csv in [("baseline-legacy-table", str(LEGACY_CSV)),
                           ("baseline-v2-default", str(V2_DEFAULT_CSV))]:
            print(f"[baseline] {label}", flush=True)
            n199 = _n199_par(ex, csv)
            asr = bh.as_r2(csv)
            ddg = _ddg_par(ex, csv)
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
            print(f"    n199_r2={n199['r2']:.4f}(n={n199['n']}) as_r2={asr['r2']:.4f}(n={asr['n']}) "
                  f"ddg_r2={ddg['r2']:.4f}(n={ddg['n']})", flush=True)

        # ---- candidate configs ----
        set_cache: dict = {}
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
            cc = str(cand_csv)

            print(f"[{cfg['label']}] set={sname} t={cfg['t']} variant={cfg['variant']} "
                  f"pc={cfg['pseudocount']} | n_struct={len(ab)} n_clusters={len(include)}",
                  flush=True)

            cr = bh.cognate_rank_auc(contacts, markup, include)  # cheap, main proc
            n199 = _n199_par(ex, cc)
            asr = bh.as_r2(cc)
            ddg = _ddg_par(ex, cc)

            results.append({
                "label": cfg["label"], "set": sname, "t": cfg["t"], "variant": cfg["variant"],
                "pseudocount": cfg["pseudocount"],
                "n_structures": len(ab), "n_clusters": len(include),
                "cognate_median_rank_pct": cr["median_rank_pct"],
                "cognate_rank_auc": cr["rank_auc"], "cognate_n": cr["n"],
                "n199_r2": n199["r2"], "n199_n": n199["n"],
                "n199_tcren_sign": n199["coefficients"]["tcren"]["sign"],
                "n199_tcren_sig": n199["coefficients"]["tcren"]["significant"],
                "as_r2": asr["r2"], "as_n": asr["n"],
                "as_tcren_sign": asr["coefficients"]["TCRen"]["sign"],
                "ddg_direct_r2": ddg["r2"], "ddg_n": ddg["n"],
                "ddg_sign": ddg["coefficients"]["d_tcr_pep"]["sign"],
                "candidate_csv": cc,
            })
            print(f"    cog_auc={cr['rank_auc']:.4f} med_rank={cr['median_rank_pct']:.2f} "
                  f"n199_r2={n199['r2']:.4f}({n199['coefficients']['tcren']['sign']},"
                  f"sig={n199['coefficients']['tcren']['significant']},n={n199['n']}) "
                  f"as_r2={asr['r2']:.4f}(n={asr['n']}) ddg_r2={ddg['r2']:.4f}(n={ddg['n']})",
                  flush=True)

    out = _CACHE / "rederive_full_metrics.json"
    out.write_text(json.dumps(results, indent=2))
    print(f"\nwrote {len(results)} rows -> {out}", flush=True)
    print("===JSON_BEGIN===", flush=True)
    print(json.dumps(results), flush=True)
    print("===JSON_END===", flush=True)


if __name__ == "__main__":
    main()
