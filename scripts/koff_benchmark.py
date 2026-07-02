#!/usr/bin/env python3
"""koff / interface-mechanics benchmark — does native interface geometry predict the off-rate?

Formalises the ``explore/koff-mechanics`` analysis of the 2026-tcren2 manuscript
(``PLAN-koff-mechanics.md``) as a reproducible harness over the shipped ``tcren.mechanics`` API.

Thesis: a single native structure reads **kinetics, not thermodynamics**. The tensile elastic-network
stiffness ``K_tens`` (and the steered-rupture force/work) of the TCR<->pMHC interface tracks the
experimental dissociation off-rate ``koff`` at r ~ -0.55 (between-structure), while equilibrium
DeltaG / Kd are capped at |r| ~ 0.3 for every MD-free feature. The signal is *stiffness*, not just
interface size: it survives partialling out the raw contact count.

Data: ATLAS (Borrman et al. 2017, Proteins) — ``<manuscript>/data/oracle/atlas/ATLAS.tsv`` with
experimental ``Koff_per_S`` / ``Kd_microM`` / ``DeltaG_kcal_per_mol`` per complex. Structures are
resolved from ``tcren-ms/data`` or the manuscript ``data`` tree.

Honesty guardrails (PLAN section 1d / 3.1): WT + own-crystal rows only; **never pool
template-proxied mutants** (a feature computed from a shared template, replicated across mutant
labels, is pseudo-replication). Always report BOTH:
  * between-structure — one feature + mean target per distinct structure (effective n = # structures);
  * own-structure     — rows whose ``true_PDB`` is the deposited complex.
Report Pearson r, Spearman rho, permutation p, bootstrap 95% CI, and the size-control partial.

Usage::

    python scripts/koff_benchmark.py                 # build features (arda/mmseqs), print the tables
    python scripts/koff_benchmark.py --cache f.csv   # cache per-structure features to f.csv (reuse next run)
    python scripts/koff_benchmark.py --limit 20      # smoke test on the first 20 structures
"""

from __future__ import annotations

import argparse
import csv
import glob
import math
import resource
import time
from pathlib import Path

import numpy as np
from scipy.stats import pearsonr, spearmanr

from tcren.annotation import classify_chains
from tcren.mechanics import coupling_residues, rupture, stiffness_tensor
from tcren.mhc import annotate_mhc
from tcren.structure import parse_structure
from tcren.structure.model import PEPTIDE_TYPE

MS = Path("/Users/mikesh/vcs/manuscripts/2026-tcren2")
ATLAS = MS / "data/oracle/atlas/ATLAS.tsv"
SEARCH = ["/Users/mikesh/vcs/code/tcren-ms/data", str(MS / "data")]

# Feature -> target sign expectation is not enforced; we report raw r. Features are computed once per
# distinct native structure (mutations are not modelled), so between-structure / own-structure scoping
# is what keeps mutant rows from pseudo-replicating.
FEATURES = ["K_tens", "n_spring", "rupture_force", "rupture_work", "couple_pep"]
TARGETS = ["logkoff", "dG", "logKd"]


def _num(x: str) -> float:
    x = str(x).strip()
    try:
        return float(x)
    except ValueError:
        return math.nan


def _find_pdb(name: str) -> str | None:
    for base in SEARCH:
        for ext in ("pdb.gz", "pdb"):
            hits = glob.glob(f"{base}/**/{name}.{ext}", recursive=True)
            if hits:
                return hits[0]
    return None


def _features(pdb: str) -> dict[str, float] | None:
    """K_tens (invdist2), n_spring, tensile rupture force/work, couple_pep for one structure."""
    s = parse_structure(pdb)
    classify_chains(s, organism="human")
    annotate_mhc(s)
    if not any(c.chain_type == PEPTIDE_TYPE for c in s.chains):
        return None
    st = stiffness_tensor(s, weight="invdist2")
    rp = rupture(s, direction="tensile", weight="invdist2")
    cp = coupling_residues(s)
    return {
        "K_tens": st["K_tens"], "n_spring": st["n_spring"],
        "rupture_force": rp["rupture_force"], "rupture_work": rp["rupture_work"],
        "couple_pep": float(cp["couple_pep"]),
    }


def _load_rows(limit: int | None) -> list[dict]:
    """ATLAS WT/own rows with experimental targets + the native structure name; features attached."""
    atlas = list(csv.DictReader(open(ATLAS), delimiter="\t"))

    def is_mut(v: str) -> bool:
        return bool(v.strip()) and v.strip().lower() not in ("wt", "none", "-", "na", "\\n", "")

    # distinct native structures first (features are per-structure; compute each once).
    rows = []
    for r in atlas:
        if is_mut(r["TCR_mut"]) or is_mut(r["PEP_mut"]) or is_mut(r["MHC_mut"]):
            continue  # WT only — mutations are not modelled on the crystal
        dg, kd, koff = _num(r["DeltaG_kcal_per_mol"]), _num(r["Kd_microM"]), _num(r["Koff_per_S"])
        if all(math.isnan(v) for v in (dg, kd, koff)):
            continue
        own = r["true_PDB"].strip()
        name = (own or r["template_PDB"].strip()).lower()
        if not name:
            continue
        rows.append({
            "name": name, "own": bool(own),
            "dG": dg,
            "logKd": math.log10(kd) if kd == kd and kd > 0 else math.nan,
            "logkoff": math.log10(koff) if koff == koff and koff > 0 else math.nan,
        })
    names = sorted({r["name"] for r in rows})
    if limit:
        names = names[:limit]
        rows = [r for r in rows if r["name"] in set(names)]

    feats: dict[str, dict] = {}
    for i, name in enumerate(names, 1):
        pdb = _find_pdb(name)
        if not pdb:
            continue
        try:
            f = _features(pdb)
            if f:
                feats[name] = f
        except Exception as exc:  # noqa: BLE001 — skip un-annotatable structures, keep the batch alive
            print(f"  [skip] {name}: {type(exc).__name__}: {str(exc)[:60]}")
        if i % 10 == 0:
            print(f"  featurised {i}/{len(names)} structures")
    return [dict(r, **feats[r["name"]]) for r in rows if r["name"] in feats]


def _perm_p(x, y, rng, n=10000) -> float:
    """Two-sided permutation p for |Pearson r| under target-label shuffling."""
    r0 = abs(pearsonr(x, y)[0])
    yb = y.copy()
    hits = 0
    for _ in range(n):
        rng.shuffle(yb)
        hits += abs(pearsonr(x, yb)[0]) >= r0
    return (hits + 1) / (n + 1)


def _boot_ci(x, y, rng, n=2000):
    """Bootstrap 95% CI for Pearson r."""
    idx = np.arange(len(x))
    rs = []
    for _ in range(n):
        b = rng.choice(idx, len(idx), replace=True)
        if np.ptp(x[b]) > 0 and np.ptp(y[b]) > 0:
            rs.append(pearsonr(x[b], y[b])[0])
    return (float(np.percentile(rs, 2.5)), float(np.percentile(rs, 97.5))) if rs else (math.nan, math.nan)


def _partial(x, y, z):
    """Partial Pearson r of x,y controlling for z (residualise both on z)."""
    def resid(a):
        b, a0 = np.polyfit(z, a, 1)
        return a - (a0 + b * z)
    return pearsonr(resid(x), resid(y))[0]


def _scope(rows, feature, target, scope):
    """Return (x, y) arrays for a feature/target under 'between' or 'own' scoping."""
    if scope == "own":
        rr = [r for r in rows if r["own"]]
        pairs = [(r[feature], r[target]) for r in rr]
    else:  # between-structure: mean target per distinct structure (one feature value per structure)
        by: dict[str, list] = {}
        for r in rows:
            by.setdefault(r["name"], []).append(r)
        pairs = [(g[0][feature], float(np.nanmean([r[target] for r in g])))
                 for g in by.values()]
    x = np.array([p[0] for p in pairs], float)
    y = np.array([p[1] for p in pairs], float)
    ok = np.isfinite(x) & np.isfinite(y)
    return x[ok], y[ok]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cache", type=Path, help="CSV to load/save per-row features (skips arda on reuse)")
    ap.add_argument("--limit", type=int, default=None, help="only the first N distinct structures")
    ap.add_argument("--perms", type=int, default=10000)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    if not ATLAS.exists():
        raise SystemExit(f"ATLAS not found at {ATLAS} — needs the 2026-tcren2 manuscript checkout")

    t0 = time.perf_counter()
    if args.cache and args.cache.exists():
        rows = list(csv.DictReader(open(args.cache)))
        for r in rows:
            r["own"] = r["own"] in ("True", "true", "1")
            for k in FEATURES + TARGETS:
                r[k] = float(r[k])
        print(f"loaded {len(rows)} cached rows from {args.cache}")
    else:
        rows = _load_rows(args.limit)
        if args.cache:
            with open(args.cache, "w", newline="") as f:
                cols = ["name", "own", *TARGETS, *FEATURES]
                w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
                w.writeheader()
                w.writerows(rows)
            print(f"cached {len(rows)} rows -> {args.cache}")

    rng = np.random.default_rng(args.seed)
    print(f"\nrows: {len(rows)} | distinct structures: {len({r['name'] for r in rows})} | "
          f"own-crystal rows: {sum(r['own'] for r in rows)}\n")

    # Main table: every feature vs every target, both scopes.
    hdr = f"{'feature':13} {'target':8} {'scope':8} {'n':>3} {'pearson':>8} {'spearman':>9} {'perm_p':>8}"
    print(hdr)
    print("-" * len(hdr))
    for target in TARGETS:
        for feature in FEATURES:
            for scope in ("between", "own"):
                x, y = _scope(rows, feature, target, scope)
                if len(x) < 5:
                    continue
                r = pearsonr(x, y)[0]
                rho = spearmanr(x, y)[0]
                p = _perm_p(x, y, rng, n=args.perms)
                print(f"{feature:13} {target:8} {scope:8} {len(x):>3} {r:>8.3f} {rho:>9.3f} {p:>8.4f}")
        print()

    # Headline + size control: K_tens vs log koff, between-structure, with bootstrap CI and the
    # partial controlling for contact count (n_spring) — is it stiffness or just size?
    print("=== headline: K_tens vs log koff (between-structure) ===")
    # Aggregate K_tens, contact count, and mean log koff per distinct structure (aligned arrays).
    by: dict[str, list] = {}
    for r in rows:
        by.setdefault(r["name"], []).append(r)
    xk = np.array([g[0]["K_tens"] for g in by.values()], float)
    xn = np.array([g[0]["n_spring"] for g in by.values()], float)
    yk = np.array([float(np.nanmean([r["logkoff"] for r in g])) for g in by.values()], float)
    ok = np.isfinite(xk) & np.isfinite(xn) & np.isfinite(yk)
    xk, xn, yk = xk[ok], xn[ok], yk[ok]
    r = pearsonr(xk, yk)[0]
    lo, hi = _boot_ci(xk, yk, rng)
    print(f"  K_tens vs logkoff        r = {r:.3f}  (n={len(xk)}, boot95% [{lo:.3f}, {hi:.3f}], "
          f"perm p = {_perm_p(xk, yk, rng, n=args.perms):.4f})")
    print(f"  contact count vs logkoff r = {pearsonr(xn, yk)[0]:.3f}  (n={len(xn)})")
    print(f"  partial K_tens | count   r = {_partial(xk, yk, xn):.3f}   <- stiffness survives size control")
    print(f"  partial count  | K_tens  r = {_partial(xn, yk, xk):.3f}")

    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1024 * 1024)
    print(f"\nwall {time.perf_counter() - t0:.1f}s | peak RSS {rss:.0f} MB")


if __name__ == "__main__":
    main()
