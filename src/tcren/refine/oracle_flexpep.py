"""FlexPepDock oracle — the accuracy ceiling the open engines are measured against.

Rosetta FlexPepDock is the reference-standard peptide refiner. It is a protocol *inside* Rosetta, so
there are two ways to reach it:

* **PyRosetta API** (preferred here) — ``pyrosetta.rosetta.protocols.flexpep_docking``. Installed via
  ``pyrosetta-installer`` (academic license). No external binary needed.
* **External binary** — a licensed ``FlexPepDocking`` executable via ``$ROSETTA_BIN`` / ``rosetta_bin=``.

Either way this is an **oracle**, never a shipped tcren dependency: the RMSD FlexPepDock achieves from a
displaced start is the lower bound on error the license-free engines (ccd/dope/openmm) should approach
(see ``CPP_REWRITE.md``). Refinement is slow (tens of seconds to minutes per structure), so the oracle
column in ``fold_benchmark.py`` is meant for a subset, not the full sweep.
"""

from __future__ import annotations

import os
import shutil
import tempfile
from pathlib import Path

import numpy as np

from ..structure import parse_structure
from ..structure.io import write_pdb
from ..structure.model import PEPTIDE_TYPE, Atom, Chain, Residue, Structure

_PYROSETTA_INIT = False


def _peptide_last(structure: Structure) -> Structure:
    """Reorder chains so the peptide is last — FlexPepDock treats the final chain as the peptide."""
    others = [c for c in structure.chains if c.chain_type != PEPTIDE_TYPE]
    pep = [c for c in structure.chains if c.chain_type == PEPTIDE_TYPE]
    return Structure(structure.pdb_id, others + pep, complex_species=structure.complex_species,
                     cell_type=structure.cell_type)


def _map_peptide_back(annotated: Structure, refined_pdb: Path, pep_id: str) -> Structure:
    """Copy ``annotated`` but replace its peptide with the refined peptide coords (preserves typing).

    Parses ``refined_pdb`` (FlexPepDock output), pulls the peptide chain's heavy atoms by residue
    order, and threads them back onto ``annotated`` so the result stays chain-typed / MHC-annotated
    (what ``peptide_rmsd`` needs) without a second arda pass.
    """
    out = parse_structure(refined_pdb, pdb_id=annotated.pdb_id)
    ref_pep = next((c for c in out.chains if c.chain_id == pep_id), None)
    if ref_pep is None:  # chain id may be reassigned; fall back to the shortest chain
        ref_pep = min(out.chains, key=lambda c: len(c.residues))
    tgt_pep = next(c for c in annotated.chains if c.chain_type == PEPTIDE_TYPE)
    if len(ref_pep.residues) != len(tgt_pep.residues):
        raise RuntimeError(f"peptide length changed by FlexPepDock: "
                           f"{len(tgt_pep.residues)} -> {len(ref_pep.residues)}")
    new_res = []
    for tgt, ref in zip(tgt_pep.residues, ref_pep.residues):
        atoms = tuple(Atom(a.name, a.element, np.asarray(a.coord, dtype=np.float64))
                      for a in ref.atoms if a.element != "H")
        new_res.append(Residue(tgt.seq_index, tgt.pdb_index, tgt.insertion_code,
                               tgt.aa, tgt.resname, atoms))
    new_pep = Chain(tgt_pep.chain_id, new_res, chain_type=tgt_pep.chain_type,
                    chain_supertype=tgt_pep.chain_supertype)
    chains = [new_pep if c is tgt_pep else c for c in annotated.chains]
    return Structure(annotated.pdb_id, chains, complex_species=annotated.complex_species,
                     cell_type=annotated.cell_type)


def pyrosetta_available() -> bool:
    try:
        import pyrosetta  # noqa: F401
    except ImportError:
        return False
    return True


def _resolve_bin(rosetta_bin: str | None) -> str | None:
    cand = rosetta_bin or os.environ.get("ROSETTA_BIN")
    if cand and Path(cand).is_file():
        return cand
    for name in ("FlexPepDocking", "FlexPepDocking.default.macosclangrelease",
                 "FlexPepDocking.static.linuxgccrelease"):
        found = shutil.which(name)
        if found:
            return found
    return None


def flexpep_available(rosetta_bin: str | None = None) -> bool:
    """True if FlexPepDock is reachable (PyRosetta installed, or a binary resolvable)."""
    return pyrosetta_available() or _resolve_bin(rosetta_bin) is not None


def _init_pyrosetta(seed: int) -> None:
    global _PYROSETTA_INIT
    if _PYROSETTA_INIT:
        return
    import pyrosetta

    # -ex1/-ex2aro: extra rotamers; pep_refine set per-run via the protocol; keep it quiet + deterministic.
    pyrosetta.init(f"-mute all -ex1 -ex2aro -run:constant_seed -run:jran {int(seed)} "
                   "-flexPepDocking:pep_refine -flexPepDocking:flexpep_score_only false")
    _PYROSETTA_INIT = True


def _pyrosetta_refine(structure: Structure, pep_id: str, seed: int) -> Structure:
    import pyrosetta
    from pyrosetta.rosetta.protocols.flexpep_docking import FlexPepDockingProtocol

    _init_pyrosetta(seed)
    with tempfile.TemporaryDirectory() as td:
        in_pdb = write_pdb(_peptide_last(structure), Path(td) / "in.pdb", keep_hydrogens=False)
        pose = pyrosetta.pose_from_file(str(in_pdb))
        FlexPepDockingProtocol().apply(pose)  # pep_refine mode set via init flags
        out_pdb = Path(td) / "refined.pdb"
        pose.dump_pdb(str(out_pdb))
        pep_id = next(c.chain_id for c in structure.chains if c.chain_type == PEPTIDE_TYPE)
        return _map_peptide_back(structure, out_pdb, pep_id)


def _binary_refine(structure: Structure, binary: str, nstruct: int,
                   extra_flags: tuple[str, ...]) -> Structure:
    import subprocess

    with tempfile.TemporaryDirectory() as td:
        tdp = Path(td)
        in_pdb = write_pdb(_peptide_last(structure), tdp / "input.pdb")
        cmd = [binary, "-s", str(in_pdb), "-flexPepDocking:pep_refine",
               "-nstruct", str(nstruct), "-out:path:all", str(tdp), *extra_flags]
        subprocess.run(cmd, cwd=tdp, check=True, capture_output=True, text=True)
        out = sorted(tdp.glob("input_*.pdb"))
        if not out:
            raise RuntimeError(f"FlexPepDock produced no output in {tdp}")
        pep_id = next(c.chain_id for c in structure.chains if c.chain_type == PEPTIDE_TYPE)
        return _map_peptide_back(structure, out[0], pep_id)


def flexpep_refine(structure: Structure, *, rosetta_bin: str | None = None, seed: int = 0,
                   nstruct: int = 1, extra_flags: tuple[str, ...] = ()) -> Structure:
    """Refine the peptide pose with FlexPepDock; return the refined :class:`Structure`.

    Uses the PyRosetta API if installed, else a resolvable ``FlexPepDocking`` binary. The structure
    should carry full chains (imported with ``keep_c_gene=True``). Raises ``RuntimeError`` if neither
    route is available.
    """
    pep = next((c for c in structure.chains if c.chain_type == PEPTIDE_TYPE), None)
    if pep is None:
        raise ValueError(f"no peptide chain in {structure.pdb_id!r}")

    if pyrosetta_available():
        return _pyrosetta_refine(structure, pep.chain_id, seed)
    binary = _resolve_bin(rosetta_bin)
    if binary is not None:
        return _binary_refine(structure, binary, nstruct, extra_flags)
    raise RuntimeError("FlexPepDock unavailable: install pyrosetta-installer or set $ROSETTA_BIN")
