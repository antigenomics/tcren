"""ProMod3 side-chain reconstruction engine (Apache-2.0, optional).

ProMod3 is the open (Apache-2.0, license-key-free) replacement for MODELLER. Its cleanest, most
directly useful capability for a groove peptide whose backbone is already placed (threaded or
CCD-closed) is **rotamer side-chain reconstruction** (``modelling.ReconstructSidechains``): it rebuilds
the side chains that :func:`substitute_peptide` strips, packed against the real receptor context with
ProMod3's rotamer library — the SCWRL/packer role. Full backbone *loop building* (``promod3.loop``
CCD/KIC + fragments, the MODELLER-loopmodel analog) is the documented upgrade; repack is the piece that
works via a stable API and complements the geometric (`ccd`) and physics (`openmm`) engines.

Optional dependency: ``conda install -c bioconda openstructure promod3``. Raises
:class:`EngineUnavailable` (never an ImportError at import) when absent. Reference oracle for the native
C++ rewrite (see ``CPP_REWRITE.md``): validates a future native rotamer packer.
"""

from __future__ import annotations

import tempfile
from pathlib import Path

import numpy as np

from ...structure.io import write_pdb
from ...structure.model import PEPTIDE_TYPE, Structure
from ..anchors import Decomposition
from ..oracle_flexpep import _map_peptide_back
from .base import EngineUnavailable, ModelResult


class ProMod3Engine:
    name = "promod3"

    def available(self) -> bool:
        try:
            import promod3  # noqa: F401
            import ost  # noqa: F401
        except ImportError:
            return False
        return True

    def run(self, structure: Structure, decomp: Decomposition, *, seed: int = 0,
            anchor_targets: np.ndarray | None = None) -> ModelResult:
        if not self.available():
            raise EngineUnavailable(
                "ProMod3/OpenStructure not installed; conda install -c bioconda openstructure promod3")
        import ost
        from promod3 import modelling

        pep = next((c for c in structure.chains if c.chain_type == PEPTIDE_TYPE), None)
        if pep is None:
            raise ValueError(f"no peptide chain in {structure.pdb_id!r}")

        with tempfile.TemporaryDirectory() as td:
            in_pdb = write_pdb(structure, Path(td) / "in.pdb", keep_hydrogens=False)
            ent = ost.io.LoadPDB(str(in_pdb))
            # Rebuild all side chains against the full-complex context (rotamer packing).
            modelling.ReconstructSidechains(ent, keep_sidechains=False)
            out_pdb = Path(td) / "packed.pdb"
            ost.io.SavePDB(ent, str(out_pdb))
            refined = _map_peptide_back(structure, out_pdb, pep.chain_id)

        return ModelResult(refined, float("nan"), self.name, decomp.anchors,
                           info={"mode": "ProMod3 rotamer side-chain reconstruction"})
