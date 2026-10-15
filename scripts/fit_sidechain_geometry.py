"""Fit the side-chain internal-coordinate table from the Canonical2026 crystals. 2026-09-08.

`tcren.refine.sidechain.TOPOLOGY` declares the chemistry -- which atom hangs off which three, and
which chi drives its torsion -- and nothing numeric. This measures the numbers: per residue type and
atom, the bond length to its parent, the bond angle, and the torsion offset (the dihedral minus the
chi that drives it, or the bare dihedral where none does).

Circular quantities are averaged on the circle, never arithmetically. Writes
`src/tcren/data/sidechain_geometry.json`, which `build_sidechain` reads.

    python scripts/fit_sidechain_geometry.py --structures data/Canonical2026 -o src/tcren/data/sidechain_geometry.json
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from tcren.refine.sidechain import TOPOLOGY, bond_angle, dihedral  # noqa: E402
from tcren.structure import import_structure  # noqa: E402


def _coords(residue) -> dict[str, np.ndarray]:
    return {a.name: np.asarray(a.coord, float) for a in residue.atoms}


def chi_of(xyz: dict[str, np.ndarray], aa: str, k: int) -> float | None:
    """The k-th chi of a residue as built, read back through the same frames TOPOLOGY declares."""
    for atom, a, b, c, chi in TOPOLOGY[aa]:
        if chi == k and all(n in xyz for n in (a, b, c, atom)):
            return dihedral(xyz[a], xyz[b], xyz[c], xyz[atom])
    return None


def circmean(v: np.ndarray) -> float:
    return float(math.atan2(np.sin(v).mean(), np.cos(v).mean()))


def circsd(v: np.ndarray) -> float:
    r = math.hypot(np.sin(v).mean(), np.cos(v).mean())
    return float(math.sqrt(max(0.0, -2.0 * math.log(min(1.0, max(r, 1e-12))))))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--structures", default="data/Canonical2026")
    ap.add_argument("-o", "--out", default="src/tcren/data/sidechain_geometry.json")
    ap.add_argument("--limit", type=int, default=0)
    args = ap.parse_args()

    paths = sorted(p for p in Path(args.structures).iterdir() if p.suffix in (".pdb", ".gz", ".cif"))
    if args.limit:
        paths = paths[: args.limit]

    bonds: dict[tuple, list[float]] = defaultdict(list)
    angles: dict[tuple, list[float]] = defaultdict(list)
    offsets: dict[tuple, list[float]] = defaultdict(list)
    n_res = defaultdict(int)

    for i, p in enumerate(paths):
        try:
            s = import_structure(str(p))
        except Exception:
            continue
        for ch in s.chains:
            for r in ch.residues:
                rows = TOPOLOGY.get(r.aa)
                if not rows:
                    continue
                xyz = _coords(r)
                if not {"N", "CA", "C", "CB"} <= xyz.keys():
                    continue
                n_res[r.aa] += 1
                chis = {k: chi_of(xyz, r.aa, k) for k in range(1, 5)}
                for atom, a, b, c, chi in rows:
                    if not {atom, a, b, c} <= xyz.keys():
                        continue
                    bonds[(r.aa, atom)].append(float(np.linalg.norm(xyz[atom] - xyz[c])))
                    angles[(r.aa, atom)].append(bond_angle(xyz[b], xyz[c], xyz[atom]))
                    t = dihedral(xyz[a], xyz[b], xyz[c], xyz[atom])
                    if chi is not None:
                        if chis.get(chi) is None:
                            continue
                        t -= chis[chi]
                    offsets[(r.aa, atom)].append(math.atan2(math.sin(t), math.cos(t)))
        if (i + 1) % 100 == 0:
            print(f"  {i + 1}/{len(paths)} structures", file=sys.stderr)

    table: dict[str, dict] = {}
    print(f"{'aa':2s} {'atom':4s} {'n':>7s} {'bond':>7s} {'sd':>6s} "
          f"{'angle':>7s} {'sd':>6s} {'offset':>8s} {'sd':>6s}")
    for aa, rows in TOPOLOGY.items():
        if not rows:
            table[aa] = {}
            continue
        entry = {}
        for atom, a, b, c, chi in rows:
            k = (aa, atom)
            if not bonds[k]:
                raise SystemExit(f"no observations for {aa}:{atom} -- TOPOLOGY names a bad frame")
            bl = np.asarray(bonds[k]); an = np.asarray(angles[k]); of = np.asarray(offsets[k])
            entry[atom] = {"bond": float(bl.mean()), "bond_sd": float(bl.std()),
                           "angle": circmean(an), "angle_sd": circsd(an),
                           "offset": circmean(of), "offset_sd": circsd(of), "n": int(bl.size)}
            print(f"{aa:2s} {atom:4s} {bl.size:7d} {bl.mean():7.4f} {bl.std():6.4f} "
                  f"{math.degrees(circmean(an)):7.2f} {math.degrees(circsd(an)):6.2f} "
                  f"{math.degrees(circmean(of)):8.2f} {math.degrees(circsd(of)):6.2f}")
        table[aa] = entry

    out = Path(args.out)
    out.write_text(json.dumps({
        "_note": ("Side-chain internal coordinates fitted from Canonical2026 by "
                  "scripts/fit_sidechain_geometry.py. Angles and offsets in RADIANS, bonds in "
                  "angstrom. `offset` is the (a,b,c,atom) dihedral minus the chi that drives it, "
                  "or the bare dihedral where TOPOLOGY names no chi. Circular quantities are "
                  "circular means."),
        "n_structures": len(paths),
        "n_residues": {k: int(v) for k, v in sorted(n_res.items())},
        "geometry": table,
    }, indent=1, sort_keys=True) + "\n")
    print(f"\nwrote {out} over {len(paths)} structures", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
