#!/usr/bin/env python3
"""Cross-peptide docking accuracy benchmark for the open-source fold engines.

The honest accuracy question (not self-reconstruction): take pMHC structure **A**, replace its whole
peptide with a *different* peptide **P_B** that binds the **same MHC allele** and whose native complex
**B** we also know, model P_B into A's groove, and measure RMSD of the modelled P_B pose to its
**true** pose in B (superposing A's MHC groove onto B's).

  A (MHC_A + P_A native) --substitute P_B--> thread P_B on A's backbone --engine--> modelled P_B
  compare to  B (MHC_B + P_B native)   [MHC-groove superposition A->B]   => cross-dock RMSD

Pairs are same-allele, same-length peptides (``substitute_peptide`` needs equal length; class-I A*02
9-mers dominate). A ``baseline`` row (threaded, no engine) shows how much error is just the
A-vs-B backbone difference before any refinement.

    python scripts/fold_crossdock_benchmark.py --index-limit 150 --max-pairs 20
    python scripts/fold_crossdock_benchmark.py --engines dope,openmm --max-pairs 10 --oracle
"""

from __future__ import annotations

import argparse
import re
import time
from itertools import permutations
from pathlib import Path

import polars as pl

from tcren.annotation import classify_chains
from tcren.mhc import annotate_mhc
from tcren.refine import (model_peptide, native_peptide, peptide_rmsd, predict_anchors,
                          substitute_peptide)
from tcren.refine.engines import available_engines
from tcren.refine.oracle_flexpep import flexpep_available, flexpep_refine
from tcren.structure import parse_structure

_ALLELE = re.compile(r"([A-Z]+[0-9]?\*[0-9]{2})")


def _annotated(path: Path):
    s = parse_structure(path, pdb_id=path.name.split(".")[0])
    classify_chains(s, organism="human")
    annotate_mhc(s)
    return s


def _allele_key(structure) -> str | None:
    for c in structure.chains:
        if c.chain_type == "MHCa" and c.allele_info:
            m = _ALLELE.search(str(c.allele_info))
            if m:
                return m.group(1)  # gene + 2-digit field, e.g. A*02 (groups A*02:01 variants)
    return None


def build_index(dataset: Path, limit: int | None, cache: Path) -> pl.DataFrame:
    if cache.exists():
        df = pl.read_csv(cache)
        print(f"loaded index cache: {cache} ({df.height} structures)")
        return df
    paths = sorted(dataset.glob("*.pdb.gz")) + sorted(dataset.glob("*.pdb"))
    if limit:
        paths = paths[:limit]
    rows = []
    for i, p in enumerate(paths):
        try:
            s = _annotated(p)
            pep = native_peptide(s)
            allele = _allele_key(s)
            if allele and pep:
                rows.append({"pdb": p.name.split(".")[0], "path": str(p), "allele": allele,
                             "pep": pep, "length": len(pep)})
        except Exception:
            pass
        if (i + 1) % 25 == 0:
            print(f"  annotated {i + 1}/{len(paths)} ...")
    df = pl.DataFrame(rows)
    cache.parent.mkdir(parents=True, exist_ok=True)
    df.write_csv(cache)
    print(f"built index: {df.height} annotatable structures -> {cache}")
    return df


def pairs_from_index(df: pl.DataFrame) -> list[tuple[dict, dict]]:
    """Ordered (A, B) pairs sharing allele + length but with DIFFERENT peptide sequences."""
    out = []
    for (allele, length), g in df.group_by(["allele", "length"]):
        recs = g.to_dicts()
        if len(recs) < 2:
            continue
        for a, b in permutations(recs, 2):
            if a["pep"] != b["pep"]:
                out.append((a, b))
    return out


def run(df: pl.DataFrame, engines: list[str], max_pairs: int, seed: int,
        oracle: bool, rosetta_bin: str | None) -> pl.DataFrame:
    all_pairs = pairs_from_index(df)
    all_pairs.sort(key=lambda ab: (ab[0]["allele"], ab[0]["length"], ab[0]["pdb"], ab[1]["pdb"]))
    pairs = all_pairs[:max_pairs]
    print(f"{len(all_pairs)} same-allele/same-length ordered pairs; running {len(pairs)}")
    use_oracle = oracle and flexpep_available(rosetta_bin)

    struct_cache: dict[str, object] = {}

    def get(rec):
        if rec["pdb"] not in struct_cache:
            struct_cache[rec["pdb"]] = _annotated(Path(rec["path"]))
        return struct_cache[rec["pdb"]]

    rows = []
    for a, b in pairs:
        try:
            sa, sb = get(a), get(b)
            pep_b = b["pep"]
            decomp_b = predict_anchors(pep_b, sb)
            # baseline: thread P_B onto A, no refinement, vs B native.
            threaded = substitute_peptide(sa, pep_b)
            base_rm = peptide_rmsd(threaded, sb, anchors=decomp_b.anchors)
        except Exception as exc:
            print(f"  skip {a['pdb']}->{b['pdb']}: {exc}")
            continue

        tag = {"A": a["pdb"], "B": b["pdb"], "allele": a["allele"], "length": a["length"],
               "pepA": a["pep"], "pepB": pep_b}
        rows.append({**tag, "method": "baseline(threaded)", "backbone_rmsd": base_rm.backbone_rmsd,
                     "ca_rmsd": base_rm.ca_rmsd, "anchor_ca_rmsd": base_rm.anchor_ca_rmsd,
                     "groove_rmsd": base_rm.groove_rmsd, "ms": 0.0, "status": "ok"})

        for engine in engines:
            t0 = time.perf_counter()
            try:
                res = model_peptide(sa, pep_b, engine=engine, seed=seed)
                rm = peptide_rmsd(res.structure, sb, anchors=decomp_b.anchors)
                rows.append({**tag, "method": engine, "backbone_rmsd": rm.backbone_rmsd,
                             "ca_rmsd": rm.ca_rmsd, "anchor_ca_rmsd": rm.anchor_ca_rmsd,
                             "groove_rmsd": rm.groove_rmsd, "ms": (time.perf_counter() - t0) * 1e3,
                             "status": "ok"})
            except Exception as exc:
                rows.append({**tag, "method": engine, "backbone_rmsd": None, "ca_rmsd": None,
                             "anchor_ca_rmsd": None, "groove_rmsd": None,
                             "ms": (time.perf_counter() - t0) * 1e3, "status": f"fail: {exc}"})

        if use_oracle:
            t0 = time.perf_counter()
            try:
                refined = flexpep_refine(threaded, seed=seed)
                rm = peptide_rmsd(refined, sb, anchors=decomp_b.anchors)
                rows.append({**tag, "method": "flexpep(oracle)", "backbone_rmsd": rm.backbone_rmsd,
                             "ca_rmsd": rm.ca_rmsd, "anchor_ca_rmsd": rm.anchor_ca_rmsd,
                             "groove_rmsd": rm.groove_rmsd, "ms": (time.perf_counter() - t0) * 1e3,
                             "status": "ok"})
            except Exception as exc:
                print(f"  oracle skip {a['pdb']}->{b['pdb']}: {exc}")

    return pl.DataFrame(rows)


def summarize(df: pl.DataFrame) -> None:
    if df.is_empty():
        print("no rows")
        return
    n_pairs = df.select(pl.struct("A", "B").n_unique()).item()
    print(f"\n=== cross-dock accuracy: peptide RMSD to native P_B (Å), over {n_pairs} pairs ===")
    ok = df.filter(pl.col("status") == "ok")
    agg = (ok.group_by("method")
             .agg(pl.len().alias("n"),
                  pl.col("backbone_rmsd").median().round(2).alias("bb_med"),
                  pl.col("backbone_rmsd").mean().round(2).alias("bb_mean"),
                  pl.col("ca_rmsd").median().round(2).alias("ca_med"),
                  pl.col("backbone_rmsd").quantile(0.75).round(2).alias("bb_p75"),
                  pl.col("ms").median().round(0).alias("ms_med"))
             .sort("bb_med"))
    print(agg)
    print("baseline(threaded) = P_B on A's backbone, no refinement; engines should beat it to add value.")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dataset", type=Path, default=Path("data/Native2026"))
    ap.add_argument("--index-limit", type=int, default=150, help="how many structures to annotate for pairing")
    ap.add_argument("--max-pairs", type=int, default=20)
    ap.add_argument("--engines", default=None, help="comma list; default = all runnable")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--oracle", action="store_true", help="add FlexPepDock oracle column (slow)")
    ap.add_argument("--rosetta-bin", default=None)
    ap.add_argument("--cache", type=Path, default=Path("scratch/native2026_index.csv"))
    ap.add_argument("--out", type=Path, default=Path("scratch/fold_crossdock.csv"))
    args = ap.parse_args()

    engines = args.engines.split(",") if args.engines else available_engines()
    t0 = time.perf_counter()
    df_index = build_index(args.dataset, args.index_limit, args.cache)
    df = run(df_index, engines, args.max_pairs, args.seed, args.oracle, args.rosetta_bin)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    df.write_csv(args.out)
    summarize(df)
    print(f"\nwrote {args.out}  ({len(df)} rows, {time.perf_counter() - t0:.1f}s)")


if __name__ == "__main__":
    main()
