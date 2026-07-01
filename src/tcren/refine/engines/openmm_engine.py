"""OpenMM anchor-restrained relaxation engine (physics, MIT-licensed, optional).

The open replacement for FlexPepDock's physics relaxation: build an OpenMM system for the complex,
freeze the receptor (TCR + MHC) by zeroing masses, pin the anchor Cα to their target pockets with
harmonic restraints (the physics analog of MODELLER's ``forms.gaussian`` distance restraints), and run
local energy minimisation. PDBFixer rebuilds the side chains that :func:`substitute_peptide` stripped,
so the output is a full-atom relaxed peptide (an improvement over the backbone-only threading).

Optional dependency: ``conda install -c conda-forge openmm pdbfixer``. Raises
:class:`EngineUnavailable` (never an ImportError at module import) when absent.

This engine is a **reference/accuracy oracle** for the native C++ rewrite (see ``CPP_REWRITE.md``): its
relaxed pose and relative energies validate the compact `_relax` minimiser that will replace it in the
dependency-free path — OpenMM's force field itself is not reimplemented.
"""

from __future__ import annotations

import tempfile
from pathlib import Path

import numpy as np

from ...structure.io import write_pdb
from ...structure.model import PEPTIDE_TYPE, Atom, Chain, Residue, Structure
from ..anchors import Decomposition
from .base import EngineUnavailable, ModelResult


class OpenMMEngine:
    name = "openmm"

    def available(self) -> bool:
        try:
            import openmm  # noqa: F401
            import pdbfixer  # noqa: F401
        except ImportError:
            return False
        return True

    def run(self, structure: Structure, decomp: Decomposition, *, seed: int = 0,
            anchor_targets: np.ndarray | None = None, restraint_k: float = 500.0,
            forcefield: str = "amber14-all.xml", max_iter: int = 0,
            tolerance: float = 5.0) -> ModelResult:
        if not self.available():
            raise EngineUnavailable(
                "OpenMM/pdbfixer not installed; conda install -c conda-forge openmm pdbfixer")
        import openmm
        from openmm import app, unit
        from pdbfixer import PDBFixer

        pep = next((c for c in structure.chains if c.chain_type == PEPTIDE_TYPE), None)
        if pep is None:
            raise ValueError(f"no peptide chain in {structure.pdb_id!r}")
        pep_id = pep.chain_id
        anchors = [i for i in decomp.anchors if 0 <= i < len(pep.residues)]

        with tempfile.TemporaryDirectory() as td:
            in_pdb = write_pdb(structure, Path(td) / "in.pdb", keep_hydrogens=False)

            # Rebuild the stripped side chains + add hydrogens.
            fixer = PDBFixer(filename=str(in_pdb))
            fixer.findMissingResidues()
            fixer.findMissingAtoms()
            fixer.addMissingAtoms()
            fixer.addMissingHydrogens(7.0)
            topology, positions = fixer.topology, fixer.positions

            ff = app.ForceField(forcefield)
            system = ff.createSystem(topology, nonbondedMethod=app.NoCutoff,
                                     constraints=app.HBonds)

            # Freeze the receptor: mass 0 => LocalEnergyMinimizer leaves it fixed; only the peptide moves.
            pep_atoms_by_res: dict[int, dict[str, int]] = {}
            for atom in topology.atoms():
                if atom.residue.chain.id != pep_id:
                    system.setParticleMass(atom.index, 0.0)
                else:
                    pep_atoms_by_res.setdefault(atom.residue.index, {})[atom.name] = atom.index

            # Anchor restraints: pull the anchor Cα toward their targets (native anchors if none given).
            pep_res_order = sorted(pep_atoms_by_res)  # topology residue index, in chain order
            if anchor_targets is not None:
                targets = np.asarray(anchor_targets, dtype=float).reshape(-1, 3)
            else:
                targets = np.array([[positions[pep_atoms_by_res[pep_res_order[i]]["CA"]].value_in_unit(unit.nanometer)]
                                    for i in anchors]).reshape(-1, 3)
                targets = targets * 10.0  # nm -> Å (targets handled in Å below)
            if len(targets) != len(anchors):
                raise ValueError(f"{len(targets)} targets for {len(anchors)} anchors")

            force = openmm.CustomExternalForce("0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
            force.addGlobalParameter("k", restraint_k * unit.kilojoule_per_mole / unit.nanometer**2)
            for p in ("x0", "y0", "z0"):
                force.addPerParticleParameter(p)
            for a_idx, tgt in zip(anchors, targets):
                ca = pep_atoms_by_res[pep_res_order[a_idx]].get("CA")
                if ca is None:
                    continue
                force.addParticle(ca, [float(tgt[0]) / 10.0, float(tgt[1]) / 10.0, float(tgt[2]) / 10.0])
            system.addForce(force)

            integrator = openmm.LangevinMiddleIntegrator(300 * unit.kelvin, 1 / unit.picosecond,
                                                         0.002 * unit.picosecond)
            integrator.setRandomNumberSeed(seed)
            context = openmm.Context(system, integrator)
            context.setPositions(positions)
            openmm.LocalEnergyMinimizer.minimize(
                context, tolerance * unit.kilojoule_per_mole / unit.nanometer, max_iter)
            state = context.getState(getPositions=True, getEnergy=True)
            energy = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
            final = state.getPositions(asNumpy=True).value_in_unit(unit.angstrom)

        refined = _rebuild_from_openmm(structure, pep, topology, pep_id, final, pep_res_order,
                                       pep_atoms_by_res)
        return ModelResult(refined, float(energy), self.name, tuple(anchors),
                           info={"forcefield": forcefield, "restraint_k": restraint_k,
                                 "mode": "receptor-frozen local minimisation"})


def _rebuild_from_openmm(structure, pep, topology, pep_id, coords_ang, pep_res_order,
                         pep_atoms_by_res):
    """Return a copy of ``structure`` with the peptide replaced by the OpenMM-relaxed heavy atoms."""
    idx_to_atom = {a.index: a for a in topology.atoms()}
    new_res = []
    for i, res in enumerate(pep.residues):
        top_res_idx = pep_res_order[i]
        atoms = []
        for name, aidx in pep_atoms_by_res[top_res_idx].items():
            elem = idx_to_atom[aidx].element
            if elem is not None and elem.symbol == "H":  # keep heavy atoms only (tcren convention)
                continue
            atoms.append(Atom(name, elem.symbol if elem else name[0],
                              np.asarray(coords_ang[aidx], dtype=np.float64)))
        new_res.append(Residue(res.seq_index, res.pdb_index, res.insertion_code,
                               res.aa, res.resname, tuple(atoms)))
    new_pep = Chain(pep.chain_id, new_res, chain_type=pep.chain_type,
                    chain_supertype=pep.chain_supertype)
    chains = [new_pep if c is pep else c for c in structure.chains]
    return Structure(structure.pdb_id, chains, complex_species=structure.complex_species,
                     cell_type=structure.cell_type)
