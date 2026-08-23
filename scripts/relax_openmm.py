#!/usr/bin/env python3
"""Full-complex OpenMM energy minimization of TCR-pMHC structures.

The physics relaxation :func:`tcren.refine_peptide` deliberately is not. ``refine_peptide`` moves
the peptide as a rigid body under DOPE, so crystal strain, added hydrogens and off-rotamer side
chains elsewhere in the complex stay exactly where the input file put them. This script minimizes
*every* atom of the complex in amber14 with GBn2 implicit solvent, which is the same energy family
an all-atom MD run evaluates, and it is roughly 10-30x faster than a Rosetta FastRelax.

Two things it is used for:

* **Forced poses.** AlphaFold/TCRmodel2 models of a non-binding pair are built to a plausible
  docking geometry regardless, and pay for it in interface strain. Minimizing relieves the strain
  without moving the model off its pose, so a geometry score trained on crystals can transfer.
* **Comparing a contact potential against MD.** MD never evaluates its energy on the deposited
  crystal: the first steps protonate the structure and let it settle. Everything that moves in that
  settling is geometry a residue-level contact map never sees, and it enters any crystal-vs-MD
  regression as noise. Scoring the minimized structure removes the difference from the tcren side.

This is minimization, not equilibration -- it reproduces the strain relief, not the nanosecond-scale
side-chain reorganisation of a production trajectory.

Needs ``openmm`` and ``pdbfixer``, which are not tcren dependencies::

    conda create -n tcren-fold -c conda-forge openmm pdbfixer
    conda run -n tcren-fold python scripts/relax_openmm.py in_dir/ out_dir/

Reads every ``*.pdb`` / ``*.pdb.gz`` in ``in_dir`` and writes ``<name>.pdb`` to ``out_dir``, skipping
any that is already there, so an interrupted run resumes. Minimization is single-structure and
independent, so shard it across cores -- one thread per shard beats one shard on all threads::

    for s in $(seq 0 11); do
      OPENMM_CPU_THREADS=1 python scripts/relax_openmm.py in/ out/ $s 12 &
    done
"""
import gzip
import shutil
import sys
import tempfile
import time
from pathlib import Path

import openmm
from openmm import app, unit
from pdbfixer import PDBFixer

src, out = Path(sys.argv[1]), Path(sys.argv[2])
out.mkdir(parents=True, exist_ok=True)
subset = sorted({p.name.split(".pdb")[0]: p for p in
                 sorted(src.glob("*.pdb")) + sorted(src.glob("*.pdb.gz"))}.items())
if len(sys.argv) > 4:
    shard, nsh = int(sys.argv[3]), int(sys.argv[4])
    subset = [x for i, x in enumerate(subset) if i % nsh == shard]

_FF = app.ForceField("amber14-all.xml", "implicit/gbn2.xml")


def minimize(in_pdb: str, out_pdb: str) -> None:
    fixer = PDBFixer(filename=in_pdb)
    fixer.removeHeterogens(keepWater=False)
    fixer.findMissingResidues(); fixer.findMissingAtoms(); fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)
    # a cutoff makes the nonbonded term O(N); NoCutoff is O(N^2) and hangs on a ~6.5k-atom complex
    system = _FF.createSystem(fixer.topology, nonbondedMethod=app.CutoffNonPeriodic,
                              nonbondedCutoff=1.0 * unit.nanometer, constraints=app.HBonds)
    integ = openmm.LangevinMiddleIntegrator(300 * unit.kelvin, 1 / unit.picosecond,
                                            0.002 * unit.picosecond)
    ctx = openmm.Context(system, integ, openmm.Platform.getPlatformByName("CPU"))
    ctx.setPositions(fixer.positions)
    openmm.LocalEnergyMinimizer.minimize(ctx, 10 * unit.kilojoule_per_mole / unit.nanometer, 300)
    pos = ctx.getState(getPositions=True).getPositions()
    with open(out_pdb, "w") as fh:
        app.PDBFile.writeFile(fixer.topology, pos, fh, keepIds=True)


for i, (name, path) in enumerate(subset):
    dst = out / f"{name}.pdb"
    if dst.exists():
        continue
    t0 = time.time()
    tmp = None
    try:
        if path.suffix == ".gz":
            with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as fh:
                with gzip.open(path, "rb") as gz:
                    shutil.copyfileobj(gz, fh)
                tmp = Path(fh.name)
        minimize(str(tmp or path), str(dst))
        print(f"  {name} minimized in {time.time() - t0:.1f}s ({i + 1}/{len(subset)})", flush=True)
    except Exception as e:
        print(f"  ERR {name}: {type(e).__name__}: {e}", flush=True)
    finally:
        if tmp is not None:
            tmp.unlink(missing_ok=True)
