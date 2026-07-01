#!/usr/bin/env python3
"""Peptide-modelling recovery benchmark for the open-source fold engines.

The honest question for a refiner is: *starting from a displaced pose, how close to the native
crystal peptide can it get?* So for each TCR-pMHC complex we

  1. thread the native peptide back onto the groove backbone (the native reference),
  2. apply a **rigid displacement** to the peptide (small random rotation about its centroid +
     translation) — the same displaced start handed to every engine,
  3. re-model from that displaced start with each engine, and
  4. measure peptide RMSD to the native crystal pose (superposed on the MHC groove).

This is *not* native-in/native-out (which would be circular): the displacement makes recovery a real
test. Two caveats are reported, not hidden:
  * ``dope`` is a rigid-body refiner restrained to its input, so it can recover a rigid displacement
    but its RMSD floor is set by the restraint — it is a *local* refiner.
  * ``ccd`` is driven to the **native anchor Cα** (the only targets available without de-novo pocket
    prediction), so its *anchor* RMSD is closure residual, not an accuracy claim; its honest metric is
    the **backbone** RMSD — can it reconstruct the loop given correct anchors + a displaced start?

The native pose is ground truth. A FlexPepDock binary (``--rosetta-bin`` / ``$ROSETTA_BIN``) adds the
*oracle* column = the RMSD FlexPepDock itself achieves, i.e. the accuracy ceiling the license-free
engines should approach. The full per-(structure,engine) table and the attempted/ok/failed/skipped
denominators are written and printed (no silent caps). This is the QC half of milestone S6.

    python scripts/fold_benchmark.py --limit 8                       # quick smoke
    RUN_BENCHMARK=1 python scripts/fold_benchmark.py                 # full sweep
    python scripts/fold_benchmark.py --limit 20 --rosetta-bin /path/to/FlexPepDocking
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np
import polars as pl

from tcren.annotation import classify_chains
from tcren.mhc import annotate_mhc
from tcren.refine import model_peptide, native_peptide, peptide_rmsd, predict_anchors, substitute_peptide
from tcren.refine.engines import available_engines
from tcren.refine.oracle_flexpep import flexpep_available, flexpep_refine
from tcren.structure import parse_structure
from tcren.structure.model import PEPTIDE_TYPE, Atom, Chain, Residue, Structure


def _annotated(path: Path) -> Structure:
    s = parse_structure(path, pdb_id=path.name.split(".")[0])
    classify_chains(s, organism="human")
    annotate_mhc(s)
    return s


def _rodrigues(axis: np.ndarray, angle: float) -> np.ndarray:
    axis = axis / (np.linalg.norm(axis) + 1e-12)
    x, y, z = axis
    c, s, C = np.cos(angle), np.sin(angle), 1.0 - np.cos(angle)
    return np.array([
        [c + x * x * C, x * y * C - z * s, x * z * C + y * s],
        [y * x * C + z * s, c + y * y * C, y * z * C - x * s],
        [z * x * C - y * s, z * y * C + x * s, c + z * z * C],
    ])


def perturb_peptide(structure: Structure, trans: float, rot_deg: float, seed: int) -> Structure:
    """Return a copy with the peptide rigidly displaced (rotation about centroid + translation).

    A rigid move is the fair shared start: a rigid-body refiner (dope) can in principle undo it, and a
    dihedral closer (ccd) sees a genuinely displaced trace. Deterministic in ``seed``.
    """
    rng = np.random.default_rng(seed)
    pep = next(c for c in structure.chains if c.chain_type == PEPTIDE_TYPE)
    coords = np.array([a.coord for r in pep.residues for a in r.atoms])
    centroid = coords.mean(0)
    R = _rodrigues(rng.standard_normal(3), rng.normal(0.0, np.radians(rot_deg)))
    t = rng.normal(0.0, trans, 3)

    new_res = []
    for res in pep.residues:
        atoms = tuple(Atom(a.name, a.element, (a.coord - centroid) @ R.T + centroid + t)
                      for a in res.atoms)
        new_res.append(Residue(res.seq_index, res.pdb_index, res.insertion_code,
                               res.aa, res.resname, atoms))
    new_pep = Chain(pep.chain_id, new_res, chain_type=pep.chain_type,
                    chain_supertype=pep.chain_supertype)
    chains = [new_pep if c is pep else c for c in structure.chains]
    return Structure(structure.pdb_id, chains, complex_species=structure.complex_species,
                     cell_type=structure.cell_type)


def _native_anchor_targets(native: Structure, seq: str, anchors: tuple[int, ...]) -> np.ndarray:
    """Native Cα at the anchor positions — the targets ccd is driven to."""
    threaded = substitute_peptide(native, seq)
    pep = next(c for c in threaded.chains if c.chain_type == PEPTIDE_TYPE)
    ca = [r.ca for r in pep.residues]
    return np.array([ca[i] for i in anchors if 0 <= i < len(ca) and ca[i] is not None])


def run(dataset: Path, limit: int | None, engines: list[str], trans: float, rot_deg: float,
        seed: int, rosetta_bin: str | None, oracle: bool) -> tuple[pl.DataFrame, dict]:
    paths = sorted(dataset.glob("*.pdb.gz")) + sorted(dataset.glob("*.pdb"))
    if limit:
        paths = paths[:limit]
    # The FlexPepDock oracle is expensive (minutes/structure), so it is opt-in (--oracle), not
    # auto-enabled just because PyRosetta happens to be importable.
    use_oracle = oracle and flexpep_available(rosetta_bin)
    print(f"{len(paths)} structures | engines={engines} | displace=({trans} Å, {rot_deg}°) | "
          f"oracle(FlexPepDock)={'on' if use_oracle else 'OFF (not installed)'}")

    counts = {"structures": len(paths), "skipped": 0}
    rows: list[dict] = []
    for path in paths:
        pdb_id = path.name.split(".")[0]
        try:
            s = _annotated(path)
            seq = native_peptide(s)
            decomp = predict_anchors(seq, s)
            targets = _native_anchor_targets(s, seq, decomp.anchors)
            displaced = perturb_peptide(s, trans, rot_deg, seed)
        except Exception as exc:
            print(f"  skip {pdb_id}: {exc}")
            counts["skipped"] += 1
            continue

        for engine in engines:
            kw = {"anchor_targets": targets, "perturb": 0.0} if engine == "ccd" else {}
            t0 = time.perf_counter()
            base = {"pdb": pdb_id, "engine": engine, "mhc_class": decomp.mhc_class,
                    "pep_len": len(seq), "n_anchors": len(decomp.anchors)}
            try:
                # Re-model from the SAME displaced start; native seq threaded onto it inside model_peptide.
                res = model_peptide(displaced, engine=engine, seed=seed, **kw)
                rm = peptide_rmsd(res.structure, s, anchors=res.anchors)
                rows.append({**base, "backbone_rmsd": rm.backbone_rmsd, "ca_rmsd": rm.ca_rmsd,
                             "anchor_ca_rmsd": rm.anchor_ca_rmsd, "groove_rmsd": rm.groove_rmsd,
                             "ms": (time.perf_counter() - t0) * 1e3, "status": "ok"})
            except Exception as exc:
                rows.append({**base, "backbone_rmsd": None, "ca_rmsd": None, "anchor_ca_rmsd": None,
                             "groove_rmsd": None, "ms": (time.perf_counter() - t0) * 1e3,
                             "status": f"fail: {exc}"})

        if use_oracle:
            t0 = time.perf_counter()
            base = {"pdb": pdb_id, "engine": "flexpep(oracle)", "mhc_class": decomp.mhc_class,
                    "pep_len": len(seq), "n_anchors": len(decomp.anchors)}
            try:
                refined = flexpep_refine(displaced, rosetta_bin=rosetta_bin)
                rm = peptide_rmsd(refined, s, anchors=decomp.anchors)
                rows.append({**base, "backbone_rmsd": rm.backbone_rmsd, "ca_rmsd": rm.ca_rmsd,
                             "anchor_ca_rmsd": rm.anchor_ca_rmsd, "groove_rmsd": rm.groove_rmsd,
                             "ms": (time.perf_counter() - t0) * 1e3, "status": "ok"})
            except Exception as exc:
                rows.append({**base, "backbone_rmsd": None, "ca_rmsd": None, "anchor_ca_rmsd": None,
                             "groove_rmsd": None, "ms": (time.perf_counter() - t0) * 1e3,
                             "status": f"fail: {exc}"})

    return pl.DataFrame(rows), counts


def summarize(df: pl.DataFrame, counts: dict) -> None:
    print(f"\nattempted: {counts['structures']} structures "
          f"({counts['skipped']} skipped: no peptide / untypable)")
    if df.is_empty():
        print("no rows")
        return
    print("\n=== per-engine recovery to native (Å); n_ok / n_fail surfaced ===")
    rep = []
    for engine in df["engine"].unique().sort().to_list():
        e = df.filter(pl.col("engine") == engine)
        ok = e.filter(pl.col("status") == "ok")
        rep.append({"engine": engine, "n_ok": ok.height, "n_fail": e.height - ok.height,
                    "bb_med": round(ok["backbone_rmsd"].median(), 3) if ok.height else None,
                    "ca_med": round(ok["ca_rmsd"].median(), 3) if ok.height else None,
                    "anchor_med": round(ok["anchor_ca_rmsd"].median(), 3) if ok.height else None,
                    "ms_med": round(ok["ms"].median(), 1) if ok.height else None})
    print(pl.DataFrame(rep))
    print("note: ccd anchor_med is closure-to-native-target residual (input-driven), not an accuracy "
          "claim; its accuracy metric is bb_med (loop reconstruction from anchors).")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dataset", type=Path, default=Path("data/Native2026"))
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--engines", default=None, help="comma list; default = all runnable")
    ap.add_argument("--trans", type=float, default=1.0, help="rigid displacement translation σ (Å)")
    ap.add_argument("--rot", type=float, default=15.0, help="rigid displacement rotation σ (degrees)")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--rosetta-bin", default=None, help="FlexPepDock binary for the oracle column")
    ap.add_argument("--oracle", action="store_true",
                    help="run the FlexPepDock oracle column (slow: minutes/structure)")
    ap.add_argument("--out", type=Path, default=Path("scratch/fold_benchmark.csv"))
    args = ap.parse_args()

    engines = args.engines.split(",") if args.engines else available_engines()
    t0 = time.perf_counter()
    df, counts = run(args.dataset, args.limit, engines, args.trans, args.rot, args.seed,
                     args.rosetta_bin, args.oracle)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    df.write_csv(args.out)
    summarize(df, counts)
    print(f"\nwrote {args.out}  ({len(df)} rows, {time.perf_counter() - t0:.1f}s)")


if __name__ == "__main__":
    main()
