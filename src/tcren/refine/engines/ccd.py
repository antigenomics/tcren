"""CCD anchor-restrained closure engine (wraps the ``tcren._fold`` C++ kernel).

Cyclic Coordinate Descent drives the peptide's anchor Cα onto target positions (predicted MHC-groove
pocket centroids) while the rest of the backbone follows as a kinematic linkage. This is the
license-free geometric path: no Rosetta, no MODELLER, only the bundled stdlib-only ``_fold`` kernel.

Draft scope (ponytail: marked so the simplification reads as intent, not ignorance):
- Operates on the **Cα trace** with consecutive-Cα rotatable bonds, and writes the closed pose back by
  rigid per-residue translation (each residue's atoms shift by its Cα displacement). The kernel
  preserves Cα–Cα distances exactly, but because adjacent residues receive *different* translations
  the **inter-residue peptide-bond geometry (C(i)–N(i+1) ≈ 1.33 Å, φ/ψ) is only approximate** — intra-
  residue geometry is intact, the chain as a whole is a Cα-trace model, NOT a physically valid all-atom
  backbone. The output therefore MUST be followed by an energy refine (DOPE / OpenMM) to regularise the
  peptide bonds; the upgrade path is a full N–Cα–C kinematic chain + rotamer repack (OpenMM/ProMod3
  do this). Do not treat the ccd output, on its own, as a finished structure.
- Anchor *targets* must be supplied. In the self-reconstruction benchmark they are the native anchor
  Cα; predicting pocket centroids de novo (groove-pocket geometry) is the open piece flagged in
  STATUS.md and left to the scoring/orient layer.
"""

from __future__ import annotations

import numpy as np

from ...structure.model import PEPTIDE_TYPE, Atom, Chain, Residue, Structure
from ..anchors import Decomposition
from .base import EngineUnavailable, ModelResult


def _peptide_chain(structure: Structure) -> Chain:
    pep = next((c for c in structure.chains if c.chain_type == PEPTIDE_TYPE), None)
    if pep is None:
        raise ValueError(f"no peptide chain in {structure.pdb_id!r}")
    return pep


def _ca_trace(pep: Chain) -> np.ndarray:
    cas = [r.ca for r in pep.residues]
    if any(c is None for c in cas):
        raise ValueError("peptide residue missing a Cα atom")
    return np.asarray(cas, dtype=np.float64)


class CcdEngine:
    name = "ccd"

    def available(self) -> bool:
        try:
            from ... import _fold  # noqa: F401
        except ImportError:
            return False
        return True

    def run(self, structure: Structure, decomp: Decomposition, *, seed: int = 0,
            anchor_targets: np.ndarray | None = None, perturb: float = 0.0,
            max_iter: int = 1000, tol: float = 0.08) -> ModelResult:
        if not self.available():
            raise EngineUnavailable("tcren._fold extension not built (run pip install -e .)")
        from ... import _fold

        pep = _peptide_chain(structure)
        ca0 = _ca_trace(pep)  # (n, 3) native/threaded Cα
        n = len(ca0)
        anchors = [i for i in decomp.anchors if 0 <= i < n]
        if not anchors:
            raise ValueError(f"no in-range anchors for peptide of length {n}")

        # Targets: caller-supplied pocket centroids, else the current (native) anchor Cα.
        if anchor_targets is None:
            targets = ca0[anchors].copy()
        else:
            targets = np.asarray(anchor_targets, dtype=np.float64).reshape(-1, 3)
            if len(targets) != len(anchors):
                raise ValueError(f"got {len(targets)} targets for {len(anchors)} anchors")

        # Displaced start so CCD has work to do (a perturbed pose, deterministic in `seed`).
        start = ca0.copy()
        if perturb > 0.0:
            start = start + np.random.default_rng(seed).normal(0.0, perturb, size=start.shape)

        # Rotatable bonds = consecutive Cα pairs; rotating bond (i,i+1) moves Cα[i+2:].
        bonds = np.array([[i, i + 1] for i in range(n - 1)], dtype=np.int32)
        moving = np.asarray(anchors, dtype=np.int32)
        weights = np.ones(len(anchors), dtype=np.float64)

        closed, rmsd, iters = _fold.ccd_close(
            np.ascontiguousarray(start), bonds, moving,
            np.ascontiguousarray(targets), weights, max_iter, tol,
        )
        closed = np.asarray(closed)

        # Write back: rigid per-residue translation by the Cα displacement (keeps residues intact).
        delta = closed - ca0
        new_res = []
        for i, res in enumerate(pep.residues):
            atoms = tuple(Atom(a.name, a.element, a.coord + delta[i]) for a in res.atoms)
            new_res.append(Residue(res.seq_index, res.pdb_index, res.insertion_code,
                                   res.aa, res.resname, atoms))
        new_pep = Chain(pep.chain_id, new_res, chain_type=pep.chain_type,
                        chain_supertype=pep.chain_supertype)
        chains = [new_pep if c is pep else c for c in structure.chains]
        refined = Structure(structure.pdb_id, chains, complex_species=structure.complex_species,
                            cell_type=structure.cell_type)

        return ModelResult(refined, float(rmsd), self.name, tuple(anchors), iterations=int(iters),
                           info={"closure": "CCD Cα-trace", "perturb": perturb})
