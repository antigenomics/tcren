"""CPL / TCRvdb binder-discrimination benchmark for candidate TCRen potentials.

Scores the TCR:peptide interface of modelled TCR-pMHC structures with each candidate
potential (reusing the production ``score_peptides`` scorer) and measures how well *low*
energy separates binders from non-binders. Two benchmarks, both TCR:peptide-only, both
"lower energy = better binder" (AUC of ``-energy``, positive class = binder):

* **CPL** — best-peptide vs worst-peptide structures per TCR, in
  ``<cpl_dir>/<tcr>_best/*.pdb`` and ``<tcr>_worst/*.pdb``. Reproduces
  ``notebooks_cleanup/<tcr>_contacts_clean.ipynb`` cell 41
  (``roc_auc_score(is_best, -tcren)``, pooled). Reported per-TCR and pooled.
* **TCRvdb** — ``padj < 1e-5`` = binder (``tcr_vdb_cleaned.ipynb`` cell 61). AUC of
  ``-energy`` over hash-named structures in ``<tcrvdb_structs>/*.pdb``. The notebook
  itself only runs Mann-Whitney; we add the matching AUC. TCRvdb structures are not yet
  on HF, so when the structure dir is empty only the label counts are reported and the
  scoring is flagged GATED (pull structures from aldan3/HF to complete it).

Contacts are computed **once per structure** and re-scored under every candidate (only
the 20x20 matrix lookup changes), so the same harness runs unchanged on the full
1,556-PDB CPL set. Run ``--verify`` to assert the fast path matches
``summarize_structure`` on the first structure.

    python3 scripts/binder_benchmark.py --verify
    python3 scripts/binder_benchmark.py --cpl-dir /path/to/data/cpl/pdb_cpl
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import mannwhitneyu, rankdata

from tcren.annotation import classify_chains
from tcren.contactmap import ContactMap
from tcren.mhc import annotate_mhc
from tcren.oracle import _native_peptide
from tcren.pipeline import _interface_energy
from tcren.potential import Potential
from tcren.potential import mj as mj_potential
from tcren.scoring import score_peptides
from tcren.scoring_rank import percentile_rank
from tcren.structure import import_structure

_MJ = mj_potential()  # fixed MJ for the two MHC interfaces (candidate-independent)

_TCREN_MS = Path(__file__).resolve().parents[1]
_CACHE = _TCREN_MS / "scratch" / "cache"
_MS = Path("/Users/mikesh/vcs/manuscripts/2026-tcren2")

# The decision set. legacy-paper reproduces the paper's published CPL AUC (the notebook's
# tcren_potentials.csv); legacy-shipped is the current production default (they differ,
# max |Δ|=1.62 — the known potential-provenance split). 2026-weighted/2026-off are the
# re-derivation candidates gated on this comparison.
CANDIDATES: dict[str, Path] = {
    "legacy-paper": _CACHE / "TCRen_legacy-paper.csv",
    "legacy-shipped": _TCREN_MS / "src" / "tcren" / "data" / "TCRen_potential.csv",
    "2026-weighted": _CACHE / "TCRen_2026-weighted-t6.csv",
    "2026-off": _CACHE / "TCRen_2026-off.csv",
}


def auc(scores: np.ndarray, labels: np.ndarray) -> float:
    """ROC-AUC via the rank-sum (Mann-Whitney) identity; higher score = positive."""
    labels = np.asarray(labels, dtype=bool)
    n1, n0 = int(labels.sum()), int((~labels).sum())
    if n1 == 0 or n0 == 0:
        return float("nan")
    r = rankdata(np.asarray(scores, dtype=float))
    return (r[labels].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)


def contacts_of(pdb: Path) -> tuple[ContactMap, str] | None:
    """Parse + classify + build the contact map once; None if it isn't a TCR-pMHC."""
    try:
        s = import_structure(str(pdb))
        classify_chains(s, organism="human", autodetect_species=True)
        annotate_mhc(s)  # types the MHC chain (needed for the two MHC interfaces)
        native = _native_peptide(s)
        if not native:
            return None
        return ContactMap.from_structure(s), native
    except Exception:
        return None


def energy(cm: ContactMap, native: str, pot: Potential) -> float | None:
    """TCR:peptide energy of the structure's native peptide under ``pot`` (sum of J)."""
    out = score_peptides(cm, [native], pot, interface="tcr_peptide")
    return float(out["score"][0]) if out.height else None


def iface_e(cm: ContactMap, pot: Potential, iface: str) -> float:
    """Native interface energy = sum of ``pot`` over the interface's observed contacts.

    The exact scorer ``run()``/``summarize_structure`` use for the ``scores`` frame (no
    peptide substitution). ``tcr_peptide`` takes the candidate TCRen; ``tcr_mhc`` and
    ``peptide_mhc`` take the fixed MJ.
    """
    return _interface_energy(cm.interface(iface, tcr_regions="all"), pot)


# --------------------------------------------------------------------------------------
# CPL: best vs worst per TCR
# --------------------------------------------------------------------------------------
def cpl_benchmark(cpl_dir: Path, pots: dict[str, Potential]) -> None:
    tcrs = sorted({d.name.rsplit("_", 1)[0] for d in cpl_dir.glob("*_*") if d.is_dir()})
    if not tcrs:
        print(f"  no <tcr>_best/<tcr>_worst folders under {cpl_dir} — nothing to score")
        return

    # one record per structure: tcr, is_best, {cand: tcr_peptide energy}, pep_mhc, tcr_mhc.
    # pep_mhc / tcr_mhc are MJ (candidate-independent) so they're computed once per structure.
    recs: list[dict] = []
    counts: dict[str, list[int]] = {t: [0, 0] for t in tcrs}
    skipped = 0
    for tcr in tcrs:
        for is_best, sub in ((1, f"{tcr}_best"), (0, f"{tcr}_worst")):
            folder = cpl_dir / sub
            if not folder.is_dir():
                continue
            for pdb in sorted(folder.glob("*.pdb")):
                got = contacts_of(pdb)
                if got is None:
                    skipped += 1
                    continue
                cm, native = got
                counts[tcr][1 - is_best] += 1
                rec = {
                    "tcr": tcr,
                    "is_best": is_best,
                    "tp": {c: iface_e(cm, pot, "tcr_peptide") for c, pot in pots.items()},
                    "pm": iface_e(cm, _MJ, "peptide_mhc"),
                    "tm": iface_e(cm, _MJ, "tcr_mhc"),
                }
                recs.append(rec)

    def auc_over(rows: list[dict], score) -> float:
        if not rows:
            return float("nan")
        s = np.array([score(r) for r in rows])
        lab = np.array([r["is_best"] for r in rows])
        return auc(-s, lab)  # lower energy -> better binder -> positive class

    cols = list(pots)

    # --- Table A: TCR:peptide only (the legacy CPL score), per-TCR + pooled ---
    print("\n=== CPL (A) — AUC of -TCRen(TCR:peptide) only, best(+) vs worst(-) ===")
    print(f"  {'TCR':<8}{'n_best':>7}{'n_worst':>8}  " + "".join(f"{c:>15}" for c in cols))
    for tcr in tcrs:
        nb, nw = counts[tcr][0], counts[tcr][1]
        sub = [r for r in recs if r["tcr"] == tcr]
        print(f"  {tcr:<8}{nb:>7}{nw:>8}  " + "".join(
            f"{auc_over(sub, lambda r, c=c: r['tp'][c]):>15.3f}" for c in cols))
    tb = sum(counts[t][0] for t in tcrs)
    tw = sum(counts[t][1] for t in tcrs)
    print(f"  {'POOLED':<8}{tb:>7}{tw:>8}  " + "".join(
        f"{auc_over(recs, lambda r, c=c: r['tp'][c]):>15.3f}" for c in cols))

    # --- Table B: does adding the MHC interfaces help? (pooled) ---
    print("\n=== CPL (B) — pooled AUC by score composition (lower energy = binder) ===")
    print(f"  {'score':<26}" + "".join(f"{c:>15}" for c in cols))
    variants = [
        ("TCR:pep (TCRen)", lambda r, c: r["tp"][c]),
        ("+ pep:MHC (MJ)", lambda r, c: r["tp"][c] + r["pm"]),
        ("+ pep:MHC + TCR:MHC", lambda r, c: r["tp"][c] + r["pm"] + r["tm"]),
    ]
    for name, fn in variants:
        print(f"  {name:<26}" + "".join(
            f"{auc_over(recs, lambda r, c=c, fn=fn: fn(r, c)):>15.3f}" for c in cols))
    # candidate-independent diagnostics (MHC interfaces alone)
    pm_auc = auc_over(recs, lambda r: r["pm"])
    tm_auc = auc_over(recs, lambda r: r["tm"])
    print(f"  {'pep:MHC only (MJ)':<26}{pm_auc:>15.3f}   [same for all candidates]")
    print(f"  {'TCR:MHC only (MJ)':<26}{tm_auc:>15.3f}   [same for all candidates]")
    if skipped:
        print(f"  ({skipped} structures skipped: not classifiable as TCR-pMHC / no peptide)")


# --------------------------------------------------------------------------------------
# TCRvdb: binder (padj<1e-5) vs non-binder
# --------------------------------------------------------------------------------------
def tcrvdb_benchmark(
    labels_csv: Path, structs_dir: Path | None, pots: dict[str, Potential], n_background: int = 1000
) -> None:
    """Authoritative TCRvdb binder-ID on the modelled structures (`isalgo/tcren_structures/tcrvdb/`).

    Labels = ``data/tcrvdb/tcrvdb_validation.csv`` (``TCR_hash``, ``is_binder`` = padj<1e-5,
    ``epitope_aa``); structures = ``<TCR_hash>.pdb`` (bootstrap ``--only tcrvdb``). For each
    structure computes, per candidate, both scores of the epitope on the TCR:peptide interface:
    the **raw energy** and the **%-rank** vs a random-peptide background (the manuscript's binder-ID
    metric — normalizes the per-structure baseline the raw energy can't). Reports AUC of
    binder-vs-non for both, pooled and per-epitope. Lower energy / lower %-rank = more binder-like.
    """
    print("\n=== TCRvdb binder-ID — padj<1e-5 = binder (authoritative set) ===")
    if not labels_csv.exists():
        print(f"  labels not found: {labels_csv}")
        return
    lab = pl.read_csv(labels_csv, infer_schema_length=5000)
    if "TCR_hash" not in lab.columns or "is_binder" not in lab.columns:
        print(f"  {labels_csv.name} lacks TCR_hash/is_binder — expected tcrvdb_validation.csv")
        return
    info = {
        r["TCR_hash"]: (bool(r["is_binder"]), r.get("epitope_aa", "?"))
        for r in lab.iter_rows(named=True)
    }
    pdbs = sorted(structs_dir.glob("*.pdb")) if structs_dir and structs_dir.is_dir() else []
    nb = sum(1 for v in info.values() if v[0])
    print(f"  labels: {len(info)} clonotypes (binders {nb}, non {len(info) - nb}); "
          f"structures present: {len(pdbs)}")
    if not pdbs:
        print("  no structures — run: python3 scripts/bootstrap_data.py --only tcrvdb --refresh")
        return

    recs: list[dict] = []  # {binder, epitope, tp:{cand:energy}, rank:{cand:rank_pct}}
    scored = skipped = 0
    for pdb in pdbs:
        h = pdb.stem
        if h not in info:
            skipped += 1
            continue
        got = contacts_of(pdb)
        if got is None:
            skipped += 1
            continue
        cm, native = got
        binder, epitope = info[h]
        rec = {"binder": binder, "epitope": epitope, "tp": {}, "rank": {}}
        for c, pot in pots.items():
            rec["tp"][c] = iface_e(cm, pot, "tcr_peptide")
            rec["rank"][c] = percentile_rank(
                cm, native, pot, interface="tcr_peptide", n_background=n_background, seed=0
            )["rank_pct"]
        recs.append(rec)
        scored += 1
    print(f"  scored {scored} ({skipped} skipped); background n={n_background}")

    def auc_over(rows: list[dict], score) -> float:
        if not rows:
            return float("nan")
        s = np.array([score(r) for r in rows])
        lab_ = np.array([r["binder"] for r in rows])
        return auc(s, lab_)  # score already signed so higher = more binder-like

    cols = list(pots)

    def report(rows: list[dict], title: str) -> None:
        b = sum(r["binder"] for r in rows)
        print(f"\n  --- {title} (n={len(rows)}, binders={b}, non={len(rows) - b}) ---")
        print(f"    {'metric':<18}" + "".join(f"{c:>15}" for c in cols))
        print(f"    {'raw energy AUC':<18}"
              + "".join(f"{auc_over(rows, lambda r, c=c: -r['tp'][c]):>15.3f}" for c in cols))
        print(f"    {'%-rank AUC':<18}"
              + "".join(f"{auc_over(rows, lambda r, c=c: -r['rank'][c]):>15.3f}" for c in cols))

    report(recs, "POOLED")
    for ep in sorted({r["epitope"] for r in recs}):
        report([r for r in recs if r["epitope"] == ep], f"epitope {ep}")


def verify_fast_path(cpl_dir: Path, pots: dict[str, Potential]) -> None:
    """Assert score_peptides(native) == summarize_structure scores[tcr_peptide]."""
    from tcren import summarize_structure

    pdb = next(cpl_dir.glob("*_best/*.pdb"), None) or next(cpl_dir.glob("*/*.pdb"), None)
    if pdb is None:
        print("verify: no CPL structure found to verify against")
        return
    for label, pot in pots.items():
        got = contacts_of(pdb)
        assert got is not None, f"could not build contacts for {pdb}"
        cm, native = got
        out = summarize_structure(
            str(pdb), superimpose=False, potentials={"tcr_peptide": str(CANDIDATES[label])}, background=10
        )
        s = out["scores"]
        checks = {
            "tcr_peptide": (iface_e(cm, pot, "tcr_peptide"), float(s["tcr_peptide"][0])),
            "peptide_mhc": (iface_e(cm, _MJ, "peptide_mhc"), float(s["peptide_mhc"][0])),
            "tcr_mhc": (iface_e(cm, _MJ, "tcr_mhc"), float(s["tcr_mhc"][0])),
        }
        for iface, (fast, slow) in checks.items():
            ok = abs(fast - slow) < 1e-9
            print(f"verify[{label}:{iface}] fast={fast:.6f} slow={slow:.6f} match={ok}")
            assert ok, f"fast path diverged from summarize_structure for {label}:{iface}"
    print(f"verify PASS on {pdb.name}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cpl-dir", type=Path, default=_MS / "data" / "cpl" / "pdb_cpl")
    ap.add_argument("--tcrvdb-labels", type=Path, default=_MS / "data" / "tcrvdb" / "tcrvdb_validation.csv")
    ap.add_argument("--tcrvdb-structs", type=Path, default=_MS / "data" / "tcrvdb")
    ap.add_argument("--tcrvdb-nbg", type=int, default=1000, help="random background size for %-rank")
    ap.add_argument("--skip-cpl", action="store_true", help="run only the TCRvdb benchmark")
    ap.add_argument("--verify", action="store_true", help="check the fast path then exit")
    args = ap.parse_args()

    missing = {k: v for k, v in CANDIDATES.items() if not v.exists()}
    if missing:
        print("WARNING: missing candidate tables:", {k: str(v) for k, v in missing.items()})
    pots = {k: Potential.from_csv(str(v)) for k, v in CANDIDATES.items() if v.exists()}
    print(f"candidates: {list(pots)}")

    t0 = time.perf_counter()
    if args.verify:
        verify_fast_path(args.cpl_dir, pots)
        return
    if not args.skip_cpl:
        cpl_benchmark(args.cpl_dir, pots)
    tcrvdb_benchmark(args.tcrvdb_labels, args.tcrvdb_structs, pots, n_background=args.tcrvdb_nbg)
    print(f"\nelapsed {time.perf_counter() - t0:.1f}s")


if __name__ == "__main__":
    main()
