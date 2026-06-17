#!/usr/bin/env python
"""Build the committed MHC pseudosequence FASTAs from NetMHCpan source tables.

NetMHCpan distributes one 34-residue groove "pseudosequence" per allele — the polymorphic
positions that line the peptide-binding groove (class I: α1/α2; class II: α1 + β1). We collapse
identical pseudosequences across alleles and write two committed reference FASTAs:

    src/tcren/data/mhci_pseudo.fa    (from MHC_pseudo.dat — HLA-A/B/C, H-2, Mamu, …)
    src/tcren/data/mhcii_pseudo.fa   (from pseudosequence.2023.all.X.dat — DRB/DQ/DP, …)

Header is simplified to ``<representative-allele>|n=<#alleles sharing the sequence>``. Source
tables are ``<allele><whitespace><34-mer>`` per line.

    python scripts/build_pseudo_fasta.py \
        --mhci ~/vcs/tmp/pseudo/MHC_pseudo.dat \
        --mhcii ~/vcs/tmp/pseudo/pseudosequence.2023.all.X.dat
"""

from __future__ import annotations

import argparse
from pathlib import Path

_OUT = Path(__file__).resolve().parents[1] / "src" / "tcren" / "data"


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
            fh.write(f">{alleles[0]}|n={len(alleles)}\n{seq}\n")
    return len(by_seq)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--mhci", type=Path, required=True, help="MHC-I source table (MHC_pseudo.dat)")
    ap.add_argument("--mhcii", type=Path, required=True, help="MHC-II source (pseudosequence.2023.all.X.dat)")
    ap.add_argument("--out", type=Path, default=_OUT, help="output dir (src/tcren/data)")
    args = ap.parse_args()

    n1 = _write_fasta(_collapse(args.mhci), args.out / "mhci_pseudo.fa")
    n2 = _write_fasta(_collapse(args.mhcii), args.out / "mhcii_pseudo.fa")
    print(f"mhci_pseudo.fa: {n1} unique\nmhcii_pseudo.fa: {n2} unique")


if __name__ == "__main__":
    main()
