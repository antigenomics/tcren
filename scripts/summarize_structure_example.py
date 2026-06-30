"""Turn a TCR-pMHC PDB into the five tcren summary tables.

Runs the :func:`tcren.summarize_structure` facade on a structure and writes each of the
five returned frames to a CSV, demonstrating the one-call entry point the paper notebooks
use. The default structure is the bundled ``1ao7`` test fixture.

Usage::

    python scripts/summarize_structure_example.py [STRUCTURE.pdb] [OUT_DIR]
"""

from __future__ import annotations

import sys
from pathlib import Path

from tcren import summarize_structure

_DEFAULT_PDB = Path(__file__).resolve().parents[1] / "tests" / "assets" / "pdb" / "1ao7.pdb"


def main() -> None:
    pdb = Path(sys.argv[1]) if len(sys.argv) > 1 else _DEFAULT_PDB
    out_dir = Path(sys.argv[2]) if len(sys.argv) > 2 else Path("summary")
    out_dir.mkdir(parents=True, exist_ok=True)

    # One call composes S1-S4: pipeline scores + markup + contacts, percentile rank, and
    # the per-position alanine scan. superimpose=False skips canonical orientation here.
    tables = summarize_structure(pdb, superimpose=False, background=1000, alanine=True)

    for name, frame in tables.items():
        path = out_dir / f"{name}.csv"
        frame.write_csv(str(path))
        print(f"{name:9s} {frame.shape[0]:>4d} x {frame.shape[1]:<2d} -> {path}")


if __name__ == "__main__":
    main()
