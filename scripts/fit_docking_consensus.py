"""Fit the bundled crystal docking-geometry consensus → src/tcren/data/docking_consensus.json.

# 2026-07-24
Produces the multivariate-normal consensus of canonical class-I TCR:pMHC docking geometries used by
``tcren.orient.tcrdock_geometry.docking_zscore`` (the fit-free, AlphaFold-free pose-canonicality prior
of Bradley 2023 / Woods 2026). Source: the ``Canonical2026`` crystal reference set (genuine PDB
structures). Class II is skipped (``docking_geometry`` supports class I only).

Run (from the repo root, venv active):
    python scripts/fit_docking_consensus.py [CANONICAL2026_DIR]
"""
from __future__ import annotations

import glob
import json
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor
from importlib import resources
from pathlib import Path

warnings.filterwarnings("ignore")

DEFAULT_DIR = Path.home() / "hf" / "tcren_structures" / "Canonical2026"


def _geometry(path: str):
    """Annotate one crystal and return its 6-param docking geometry, or None (class II / failure)."""
    from tcren.annotation import classify_chains
    from tcren.mhc import annotate_mhc
    from tcren.orient.tcrdock_geometry import DockingGeometry, docking_geometry
    from tcren.structure import import_structure
    try:
        s = import_structure(path)
        classify_chains(s, autodetect_species=True)
        annotate_mhc(s)
        g: DockingGeometry = docking_geometry(s)
        return g.to_dict()
    except Exception:
        return None


def main() -> None:
    from tcren.orient.tcrdock_geometry import DockingGeometry, consensus_distribution

    src = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_DIR
    paths = sorted(glob.glob(str(src / "*.pdb.gz"))) + sorted(glob.glob(str(src / "*.pdb")))
    print(f"fitting from {len(paths)} structures in {src}")

    with ProcessPoolExecutor(max_workers=8) as ex:
        results = list(ex.map(_geometry, paths, chunksize=4))

    geoms = [DockingGeometry(**r) for r in results if r is not None]
    n_fail = sum(r is None for r in results)
    print(f"class-I geometries: {len(geoms)}   skipped (class II / fail): {n_fail}")

    consensus = consensus_distribution(geoms, mhc_class="I")

    # sanity: leave-none-out self Mahalanobis distances should look ~chi(7) (median ~sqrt(6.3)~2.5)
    import numpy as np
    from tcren.orient.tcrdock_geometry import _encode
    d = np.array([np.sqrt(consensus.mahalanobis2(_encode(g))) for g in geoms])
    print(f"self-distance: median={np.median(d):.2f} p90={np.quantile(d, 0.9):.2f} "
          f"max={d.max():.2f}  (mean geometry d={consensus.mean[0]:.1f} A)")

    out = Path(str(resources.files("tcren.data"))) / "docking_consensus.json"
    out.write_text(json.dumps({"I": consensus.to_dict()}, indent=1) + "\n")
    print(f"wrote {out}  (n={consensus.n})")


if __name__ == "__main__":
    main()
