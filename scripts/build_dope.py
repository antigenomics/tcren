#!/usr/bin/env python
"""Build the bundled DOPE statistical-potential table from the MODELLER mean-force libraries.

DOPE (Discrete Optimized Protein Energy; Shen & Sali, Protein Science 2006) is an atom-level,
distance-dependent knowledge-based potential. tcren uses it as the energy for potential-guided
peptide refinement (`tcren.refine.refine_peptide`) — independent of the TCRen/MJ potentials used
for epitope scoring. The two source libraries are taken from the pymod/altmod distribution:

    atmcls-mf.lib  — (residue, atom name) -> atom class  (158 classes)
    dist-mf.lib    — per atom-class-pair spline of the potential vs distance
                     (upper triangle; 29 knots at 0.75..14.75 A, 0.5 A spacing)

This script fetches both (pinned altmod commit), parses them, and writes a compact
``src/tcren/data/dope_potential.npz`` (symmetric table + atom-class map + bin metadata).

    python scripts/build_dope.py
"""

from __future__ import annotations

import urllib.request
from pathlib import Path

import numpy as np

_OUT = Path(__file__).resolve().parents[1] / "src" / "tcren" / "data" / "dope_potential.npz"
_COMMIT = "9148cc7b9bf4ec0cdcd9e9ad4882a0a4aa0fd17b"
_BASE = f"https://raw.githubusercontent.com/pymodproject/altmod/{_COMMIT}/altmod/data/dope"
NBINS, X_START, DX = 29, 0.75, 0.5


def _fetch(name: str) -> str:
    with urllib.request.urlopen(f"{_BASE}/{name}", timeout=60) as fh:  # noqa: S310 - pinned host
        return fh.read().decode()


def _parse_atom_classes(text: str) -> tuple[list[str], dict[tuple[str, str], int]]:
    classes: list[str] = []
    index: dict[str, int] = {}
    atom_class: dict[tuple[str, str], int] = {}
    current = None
    for line in text.splitlines():
        s = line.strip()
        if s.startswith("ATMGRP"):
            current = s.split("'")[1]
            if current not in index:
                index[current] = len(classes)
                classes.append(current)
        elif s.startswith("ATOM") and current is not None:
            p = s.split("'")
            atom_class[(p[1], p[3])] = index[current]
    return classes, atom_class


def _parse_table(text: str, index: dict[str, int], n: int) -> np.ndarray:
    table = np.zeros((n, n, NBINS), dtype=np.float32)
    for line in text.splitlines():
        if not line.startswith("R"):
            continue
        f = line.split()
        ci, cj = index[f[8]], index[f[9]]
        y = np.asarray(f[10:][6 : 6 + NBINS], dtype=np.float32)  # 6 header floats, then NBINS knots
        table[ci, cj] = y
        table[cj, ci] = y  # symmetric (file stores the upper triangle)
    return table


def main() -> None:
    classes, atom_class = _parse_atom_classes(_fetch("atmcls-mf.lib"))
    index = {c: i for i, c in enumerate(classes)}
    table = _parse_table(_fetch("dist-mf.lib"), index, len(classes))

    keys = np.array([f"{r}:{a}" for (r, a) in atom_class], dtype="U8")
    vals = np.array(list(atom_class.values()), dtype=np.int16)
    _OUT.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(_OUT, table=table, classes=np.array(classes, dtype="U4"),
                        keys=keys, vals=vals,
                        x_start=np.float32(X_START), dx=np.float32(DX), nbins=np.int32(NBINS))
    print(f"wrote {_OUT}  ({len(classes)} classes, {len(atom_class)} atom maps, "
          f"table {table.shape}, {_OUT.stat().st_size/1e6:.2f} MB)")


if __name__ == "__main__":
    main()
