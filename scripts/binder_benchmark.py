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
from tcren.oracle import _native_peptide
from tcren.potential import Potential
from tcren.scoring import score_peptides
from tcren.structure import import_structure

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


# --------------------------------------------------------------------------------------
# CPL: best vs worst per TCR
# --------------------------------------------------------------------------------------
def cpl_benchmark(cpl_dir: Path, pots: dict[str, Potential]) -> None:
    tcrs = sorted({d.name.rsplit("_", 1)[0] for d in cpl_dir.glob("*_*") if d.is_dir()})
    if not tcrs:
        print(f"  no <tcr>_best/<tcr>_worst folders under {cpl_dir} — nothing to score")
        return

    # rows: {tcr, label(1=best), {cand: energy}}
    pooled: dict[str, list[tuple[float, int]]] = {c: [] for c in pots}  # (energy, is_best)
    per_tcr: dict[str, dict[str, list[tuple[float, int]]]] = {}
    counts: dict[str, tuple[int, int]] = {}
    skipped = 0

    for tcr in tcrs:
        per_tcr[tcr] = {c: [] for c in pots}
        nb = nw = 0
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
                nb += is_best
                nw += 1 - is_best
                for c, pot in pots.items():
                    e = energy(cm, native, pot)
                    if e is None:
                        continue
                    per_tcr[tcr][c].append((e, is_best))
                    pooled[c].append((e, is_best))
        counts[tcr] = (nb, nw)

    def auc_of(pairs: list[tuple[float, int]]) -> float:
        if not pairs:
            return float("nan")
        e = np.array([p[0] for p in pairs])
        lab = np.array([p[1] for p in pairs])
        return auc(-e, lab)  # lower energy -> better binder -> positive

    cols = list(pots)
    print("\n=== CPL benchmark — AUC of -TCRen(TCR:peptide), best(+) vs worst(-) ===")
    header = f"  {'TCR':<10}{'n_best':>7}{'n_worst':>8}  " + "".join(f"{c:>16}" for c in cols)
    print(header)
    for tcr in tcrs:
        nb, nw = counts[tcr]
        line = f"  {tcr:<10}{nb:>7}{nw:>8}  " + "".join(
            f"{auc_of(per_tcr[tcr][c]):>16.3f}" for c in cols
        )
        print(line)
    tb = sum(c[0] for c in counts.values())
    tw = sum(c[1] for c in counts.values())
    print(f"  {'POOLED':<10}{tb:>7}{tw:>8}  " + "".join(f"{auc_of(pooled[c]):>16.3f}" for c in cols))
    if skipped:
        print(f"  ({skipped} structures skipped: not classifiable as TCR-pMHC / no peptide)")


# --------------------------------------------------------------------------------------
# TCRvdb: binder (padj<1e-5) vs non-binder
# --------------------------------------------------------------------------------------
def tcrvdb_benchmark(labels_csv: Path, structs_dir: Path | None, pots: dict[str, Potential]) -> None:
    print("\n=== TCRvdb benchmark — padj<1e-5 = binder ===")
    if not labels_csv.exists():
        print(f"  labels not found: {labels_csv}")
        return
    lab = pl.read_csv(labels_csv, infer_schema_length=5000)
    lab = lab.with_columns(pl.col("padj").cast(pl.Float64, strict=False))
    labelled = lab.filter(pl.col("padj").is_not_null())
    binders = labelled.filter(pl.col("padj") < 1e-5)
    print(
        f"  labels: {labelled.height} rows with padj  "
        f"(binders padj<1e-5: {binders.height}, non-binders: {labelled.height - binders.height})"
    )

    pdbs = sorted(structs_dir.glob("*.pdb")) if structs_dir and structs_dir.is_dir() else []
    if not pdbs:
        loc = structs_dir if structs_dir else "(no --tcrvdb-structs given)"
        print(f"  structures present: 0 under {loc}")
        print("  -> SCORING GATED: TCRvdb structures are not on HF yet (root tree has no")
        print("     tcrvdb folder); pull /projects/structures/TCRvdb from aldan3, name each")
        print("     <TCR_hash>.pdb, and re-run with --tcrvdb-structs to complete this row.")
        return

    # join structure hash -> is_validated. The hash join key is documented in
    # tcr_vdb_cleaned.ipynb (sha256 of cdr3a+va+ja+cdr3b+vb+jb+mhca+mhcb+epitope). If the
    # labels file already carries a matching key column we join on it; otherwise the caller
    # must supply a hash->label map. Here we expect a `TCR_hash` column alongside padj.
    if "TCR_hash" not in labelled.columns:
        print("  labels file has no TCR_hash column — cannot join structures to padj here.")
        print("  (the hash map lives in tcr_vdb_cleaned.ipynb; add TCR_hash to the padj CSV")
        print("   or pass a joined table to complete TCRvdb scoring.)")
        return

    is_val = {r["TCR_hash"]: (r["padj"] < 1e-5) for r in labelled.iter_rows(named=True)}
    rows: dict[str, list[tuple[float, bool]]] = {c: [] for c in pots}
    scored = skipped = 0
    for pdb in pdbs:
        h = pdb.stem
        if h not in is_val:
            skipped += 1
            continue
        got = contacts_of(pdb)
        if got is None:
            skipped += 1
            continue
        cm, native = got
        scored += 1
        for c, pot in pots.items():
            e = energy(cm, native, pot)
            if e is not None:
                rows[c].append((e, is_val[h]))

    print(f"  scored {scored} structures ({skipped} skipped)")
    for c in pots:
        pairs = rows[c]
        if not pairs:
            print(f"    {c:<16} n=0")
            continue
        e = np.array([p[0] for p in pairs])
        v = np.array([p[1] for p in pairs])
        a = auc(-e, v)
        if v.sum() and (~v).sum():
            mw = mannwhitneyu(e[v], e[~v], alternative="two-sided").pvalue
        else:
            mw = float("nan")
        print(f"    {c:<16} AUC={a:.3f}  MannWhitney_p={mw:.2e}  (binders={int(v.sum())}, non={int((~v).sum())})")


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
        fast = energy(cm, native, pot)
        out = summarize_structure(
            str(pdb), superimpose=False, potentials={"tcr_peptide": str(CANDIDATES[label])}, background=10
        )
        slow = float(out["scores"]["tcr_peptide"][0])
        ok = fast is not None and abs(fast - slow) < 1e-9
        print(f"verify[{label}] fast={fast!r} slow={slow!r} match={ok}")
        assert ok, f"fast path diverged from summarize_structure for {label}"
    print(f"verify PASS on {pdb.name}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cpl-dir", type=Path, default=_MS / "data" / "cpl" / "pdb_cpl")
    ap.add_argument("--tcrvdb-labels", type=Path, default=_MS / "data" / "tcrvdb" / "tcrvdb_padj.csv")
    ap.add_argument("--tcrvdb-structs", type=Path, default=_MS / "data" / "tcrvdb" / "structures")
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
    cpl_benchmark(args.cpl_dir, pots)
    tcrvdb_benchmark(args.tcrvdb_labels, args.tcrvdb_structs, pots)
    print(f"\nelapsed {time.perf_counter() - t0:.1f}s")


if __name__ == "__main__":
    main()
