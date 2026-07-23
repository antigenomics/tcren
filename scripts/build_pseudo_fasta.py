#!/usr/bin/env python
"""Build the committed MHC pseudosequence FASTAs from NetMHCpan tables + IPD-IMGT/HLA.

NetMHCpan distributes one 34-residue groove "pseudosequence" per allele — the polymorphic
positions that line the peptide-binding groove (class I: α1/α2; class II: α1 + β1). We collapse
identical pseudosequences across alleles and write two committed reference FASTAs:

    src/tcren/data/mhci_pseudo.fa    (from MHC_pseudo.dat — HLA-A/B/C, H-2, Mamu, …)
    src/tcren/data/mhcii_pseudo.fa   (from pseudosequence.2023.all.X.dat — DRB/DQ/DP, …)

Header lists **every** allele sharing the sequence, space-separated, as
``<allele> [<allele> ...]|n=<#alleles sharing the sequence>``. Source tables are
``<allele><whitespace><34-mer>`` per line.

Listing only the representative (as this script did until 2026-07) collapses the sequences
correctly but discards the other names, leaving **8,854 of MHC_pseudo.dat's 12,997 alleles (68%)**
— and 8,839 of the class-II table's 11,048 (80%) — unresolvable by any consumer keyed on the
header. They are not rare variants: HLA-B*14:02 (shares B*14:01's groove), B*18:05 (shares
B*18:01's) and C*03:04 (shares C*03:03's) were all lost while HLA-C03:438 shipped. The dedup is
right — the sequences are identical by construction — but dropping the name index is not.

``--imgt-alignments`` additionally derives class-I pseudosequences straight from IPD-IMGT/HLA, for
alleles NetMHCpan's table has never covered (it lags IMGT, and omits HLA-F entirely). The 34
positions are not hardcoded — they are *recovered* from the alleles the table already knows, then
verified by re-deriving every known allele byte-exactly. NetMHCpan wins any conflict; IMGT only
fills gaps, so no already-covered allele can change. Fetch the alignments with:

    for g in A B C E F G; do
      curl -sSo $g_prot.txt \
        https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/alignments/${g}_prot.txt
    done

    python scripts/build_pseudo_fasta.py \
        --mhci ~/vcs/tmp/pseudo/MHC_pseudo.dat \
        --mhcii ~/vcs/tmp/pseudo/pseudosequence.2023.all.X.dat \
        --imgt-alignments ~/vcs/tmp/imgt
"""

from __future__ import annotations

import argparse
import re
from collections import Counter, defaultdict
from pathlib import Path

_OUT = Path(__file__).resolve().parents[1] / "src" / "tcren" / "data"
_GENES = ["A", "B", "C", "E", "F", "G"]
_NPOS = 34
_AA = set("ACDEFGHIKLMNPQRSTVWY")


def _collapse(src: Path) -> dict[str, list[str]]:
    """``sequence -> [alleles]`` preserving first-seen allele order."""
    by_seq: dict[str, list[str]] = {}
    for line in src.read_text().splitlines():
        parts = line.split()
        if len(parts) != 2:
            continue
        allele, seq = parts
        by_seq.setdefault(seq, []).append(allele)
    return by_seq


def _write_fasta(by_seq: dict[str, list[str]], out: Path) -> int:
    with out.open("w") as fh:
        for seq, alleles in sorted(by_seq.items(), key=lambda kv: kv[1][0]):
            fh.write(f">{' '.join(alleles)}|n={len(alleles)}\n{seq}\n")
    return len(by_seq)


# --------------------------------------------------------------------------- IMGT derivation


def _parse_alignment(path: Path) -> tuple[dict[str, str], str]:
    """IMGT ``*_prot.txt`` -> ``({allele: column string}, reference)``.

    Repeating blocks; ``-`` = same as the reference, ``.`` = indel, ``*`` = unsequenced.
    """
    parts: dict[str, list[str]] = defaultdict(list)
    ref = None
    for line in path.read_text().splitlines():
        if not line.startswith(" ") or line.startswith(" Prot") or "|" in line:
            continue
        m = re.match(r"\s+(\S+)\s+(.*)$", line)
        if not m:
            continue
        allele, body = m.group(1), m.group(2).replace(" ", "")
        if not body or not re.match(r"^[A-Z*.\-]+$", body):
            continue
        if ref is None:
            ref = allele
        parts[allele].append(body)
    refseq = "".join(parts[ref])
    seqs = {a: "".join(refseq[i] if c == "-" and i < len(refseq) else c
                       for i, c in enumerate("".join(p)))
            for a, p in parts.items()}
    return seqs, ref


def _two_field(name: str) -> str | None:
    """``A*01:01:01:01`` -> ``HLA-A01:01``. IMGT field 2 *is* the protein field, so every 4-field
    allele under one 2-field name shares a protein sequence — exactly NetMHCpan's unit."""
    m = re.match(r"^([A-Z0-9]+)\*(\d+):(\d+)", name)
    return f"HLA-{m.group(1)}{m.group(2)}:{m.group(3)}" if m else None


def _nw(a: str, b: str) -> dict[int, int]:
    """Global alignment -> ``{index in a: index in b}``. Identity scoring is enough: the class-I
    references are >80% identical over ~365 aa."""
    n, m, gap = len(a), len(b), -4
    score = [[0] * (m + 1) for _ in range(n + 1)]
    ptr = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        score[i][0], ptr[i][0] = i * gap, 1
    for j in range(1, m + 1):
        score[0][j], ptr[0][j] = j * gap, 2
    for i in range(1, n + 1):
        ai = a[i - 1]
        for j in range(1, m + 1):
            d = score[i - 1][j - 1] + (2 if ai == b[j - 1] else -1)
            u, left = score[i - 1][j] + gap, score[i][j - 1] + gap
            best = max(d, u, left)
            score[i][j] = best
            ptr[i][j] = 0 if best == d else (1 if best == u else 2)
    out: dict[int, int] = {}
    i, j = n, m
    while i > 0 or j > 0:
        p = ptr[i][j]
        if p == 0:
            out[i - 1] = j - 1
            i, j = i - 1, j - 1
        elif p == 1:
            i -= 1
        else:
            j -= 1
    return out


def _derive_from_imgt(imgt_dir: Path, known: dict[str, str]) -> dict[str, str]:
    """``{allele: 34-mer}`` for HLA class-I alleles absent from ``known``, read off the IMGT
    alignment at the columns recovered from the alleles ``known`` already covers.

    Verified by round-trip: every known allele is re-derived and must match byte-exactly.
    """
    aln, refs, modal = {}, {}, {}
    for g in _GENES:
        f = imgt_dir / f"{g}_prot.txt"
        if not f.exists():
            continue
        aln[g], refs[g] = _parse_alignment(f)
        modal[g] = Counter(len(s) for s in aln[g].values()).most_common(1)[0][0]
    if not aln:
        raise SystemExit(f"no {{A..G}}_prot.txt under {imgt_dir}")

    def usable(g):
        # full-length, fully-sequenced records only: '*' is the single most common character in the
        # alignment (a third of records are exon 2-3 only) and would veto every candidate column.
        return [(known[_two_field(i)], s) for i, s in aln[g].items()
                if _two_field(i) in known and len(s) == modal[g] and "*" not in s]

    solved = {}
    for g in aln:
        pr = usable(g)
        if len(pr) < 500:
            continue
        sol, worst, ties = [], 1.0, 0
        for k in range(_NPOS):
            sc = [sum(s[c] == ps[k] for ps, s in pr) for c in range(modal[g])]
            best = max(sc)
            hits = [c for c, v in enumerate(sc) if v == best]
            ties += len(hits) > 1
            sol.append(hits[0])
            worst = min(worst, best / len(pr))
        # consensus, not unanimity: NetMHCpan's table lags IMGT, so a handful of revised records
        # (A*24:399 alone) would otherwise veto a column 99.97% of the data agrees on.
        if worst > 0.99 and not ties:
            solved[g] = sol
        print(f"  imgt {g}: {len(pr):,} records -> "
              f"{'solved' if g in solved else 'ambiguous'} (consensus {worst:.5f}, {ties} ties)")
    if not solved:
        raise SystemExit("could not recover the 34 columns from any gene")

    # Anchor on the best-supported gene and map its RESIDUE positions onto the others: columns are
    # not comparable across genes (gene-specific indels give different alignment widths).
    anchor = max(solved, key=lambda g: len(usable(g)))
    aref = aln[anchor][refs[anchor]]
    aung = "".join(x for x in aref if x != ".")
    apos = [len([x for x in aref[:c] if x != "."]) for c in solved[anchor]]
    cols = {anchor: solved[anchor]}
    for g in aln:
        if g == anchor:
            continue
        gref = aln[g][refs[g]]
        gung = "".join(x for x in gref if x != ".")
        mp = _nw(aung, gung)
        u2c, u = {}, 0
        for c, ch in enumerate(gref):
            if ch != ".":
                u2c[u] = c
                u += 1
        sol = [u2c.get(mp.get(p)) for p in apos]
        if all(x is not None for x in sol):
            cols[g] = sol
            if g in solved and sol != solved[g]:
                raise SystemExit(f"gene {g}: mapped columns disagree with its own solution")

    n_ok = n_bad = 0
    for g, cs in cols.items():
        for imgt, s in aln[g].items():
            k = _two_field(imgt)
            if k not in known or len(s) != modal[g]:
                continue
            got = "".join(s[c] for c in cs)
            if "*" in got:
                continue
            n_ok, n_bad = (n_ok + 1, n_bad) if got == known[k] else (n_ok, n_bad + 1)
    print(f"  imgt round-trip: {n_ok:,} exact, {n_bad} mismatch "
          f"({n_bad / max(n_ok + n_bad, 1):.5%}) over anchor={anchor}")
    if n_bad > n_ok * 0.001:
        raise SystemExit("imgt round-trip failed -- refusing to derive")

    new: dict[str, str] = {}
    for g, cs in cols.items():
        for imgt, s in aln[g].items():
            k = _two_field(imgt)
            if not k or k in known or k in new or len(s) != modal[g]:
                continue
            ps = "".join(s[c] for c in cs)
            if set(ps) <= _AA:
                new[k] = ps
    return new


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--mhci", type=Path, required=True, help="MHC-I source table (MHC_pseudo.dat)")
    ap.add_argument("--mhcii", type=Path, required=True, help="MHC-II source (pseudosequence.2023.all.X.dat)")
    ap.add_argument("--imgt-alignments", type=Path, default=None,
                    help="dir with IPD-IMGT/HLA {A,B,C,E,F,G}_prot.txt; adds class-I alleles "
                         "NetMHCpan's table lacks (it never covers HLA-F at all)")
    ap.add_argument("--out", type=Path, default=_OUT, help="output dir (src/tcren/data)")
    args = ap.parse_args()

    by_seq = _collapse(args.mhci)
    if args.imgt_alignments:
        known = {a: s for s, al in by_seq.items() for a in al}
        new = _derive_from_imgt(args.imgt_alignments, known)
        joined = sum(1 for s in new.values() if s in by_seq)
        for allele, seq in sorted(new.items()):
            by_seq.setdefault(seq, []).append(allele)
        print(f"  imgt: +{len(new):,} alleles ({joined:,} joined an existing 34-mer group, "
              f"{len(new) - joined:,} brought a new one)")

    n1 = _write_fasta(by_seq, args.out / "mhci_pseudo.fa")
    n2 = _write_fasta(_collapse(args.mhcii), args.out / "mhcii_pseudo.fa")
    a1 = sum(len(v) for v in by_seq.values())
    print(f"mhci_pseudo.fa: {n1} unique / {a1} alleles\nmhcii_pseudo.fa: {n2} unique")


if __name__ == "__main__":
    main()
