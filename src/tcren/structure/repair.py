"""Restore omitted modified-peptide atoms from a matching deposited coordinate copy."""
import gzip
import hashlib
import json
from pathlib import Path

import numpy as np
import typer

from ..docking import kabsch


def restore_modified_peptide(native_lines, deposited, *, peptide_chain="C", residue="CIR", max_rmsd=.002):
    """Rigidly match standard peptide Cα atoms and restore the specified modified residue.

    Preserve the established author-number/residue-name match and PDB-rounding RMSD gate.
    Nonpeptide coordinates stay byte-identical; this is source restoration, not refinement.
    """
    peptide = [line for line in native_lines if line.startswith("ATOM") and line[21] == peptide_chain]

    def ca_by_number(lines):
        return {(int(line[22:26]), line[26], line[17:20]): np.array([float(line[k:k+8]) for k in (30, 38, 46)])
                for line in lines if line.startswith("ATOM") and line[12:16].strip() == "CA"}

    native_ca = ca_by_number(peptide)
    if len(native_ca) < 3 or np.linalg.matrix_rank(np.array(list(native_ca.values())) - np.mean(list(native_ca.values()), axis=0)) < 2:
        raise ValueError("At least three noncollinear standard peptide C-alpha anchors required")
    copies = []
    for source_chain in sorted({line[21] for line in deposited if line.startswith("ATOM")}):
        source_ca = ca_by_number([line for line in deposited if line[21:22] == source_chain])
        if not native_ca.keys() <= source_ca.keys():
            continue
        moving = np.array([source_ca[k] for k in native_ca])
        fixed = np.array(list(native_ca.values()))
        if np.linalg.matrix_rank(moving - moving.mean(axis=0)) < 2:
            continue
        rot, tran, rmsd = kabsch(moving, fixed)
        copies.append((rmsd, source_chain, rot, tran))
    if not copies:
        raise ValueError("No matching deposited peptide")
    rmsd, chain, rot, tran = min(copies, key=lambda x: x[:2])
    if not rmsd < max_rmsd:
        raise ValueError(f"Coordinate identity exceeds PDB rounding precision: {rmsd}")
    modifications = []
    for line in deposited:
        if line.startswith("HETATM") and line[21] == chain and line[17:20] == residue:
            xyz = np.array([float(line[k:k+8]) for k in (30, 38, 46)]) @ rot + tran
            modifications.append(line[:21] + peptide_chain + line[22:30] + "".join(f"{x:8.3f}" for x in xyz) + line[54:])
    if not modifications:
        raise ValueError("No matching modified atoms in deposited peptide")
    others = [line for line in native_lines if line.startswith(("ATOM", "HETATM")) and line[21] != peptide_chain]
    peptide += [line for line in native_lines if line.startswith("HETATM") and line[21] == peptide_chain and line[17:20] != residue]
    peptide = sorted(peptide + modifications, key=lambda line: (int(line[22:26]), line[26]))
    details = {"deposited_chain": chain, "matched_standard_ca": len(native_ca),
               "restored_modified_atoms": len(modifications), "source_to_prepared_ca_RMSD_A": rmsd,
               "rotation": rot.tolist(), "translation": tran.tolist()}
    return "\n".join(others + peptide) + "\nEND\n", details


def command(
    structures: Path = typer.Option(..., exists=True, file_okay=False),
    deposits: Path = typer.Option(..., exists=True, file_okay=False),
    output: Path = typer.Option(..., "--output", "-o"),
    exclude: str = typer.Option("", help="Comma-separated source IDs intentionally not restored"),
    peptide_chain: str = typer.Option("C"),
    residue: str = typer.Option("CIR"),
):
    """Restore a batch from ID-deposited.pdb.gz sources; retain unmodified source paths."""
    if output.exists() or len(peptide_chain) != 1 or len(residue) != 3:
        raise ValueError("New output, one-character chain and three-character residue required")
    paths = sorted(q for q in structures.glob("*.pdb.gz") if not q.name.startswith("._"))
    if not paths:
        raise ValueError("No input PDB.gz files")
    output.mkdir(parents=True)
    selected, records, hashes = [], [], {}
    for path in paths:
        pid = path.name.split(".")[0]
        source = deposits / f"{pid}-deposited.pdb.gz"
        hashes[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
        if not source.exists() or pid in exclude.split(","):
            selected.append(path.resolve())
            continue
        hashes[str(source)] = hashlib.sha256(source.read_bytes()).hexdigest()
        with gzip.open(path, "rt") as f:
            native = f.read().splitlines()
        with gzip.open(source, "rt") as f:
            deposited = f.read().splitlines()
        text, record = restore_modified_peptide(native, deposited, peptide_chain=peptide_chain, residue=residue)
        target = output / path.name.removesuffix(".gz")
        target.write_text(text)
        selected.append(target.resolve())
        records.append({"pdb.id": pid, **record})
    (output / "structures.txt").write_text("".join(str(p)+"\n" for p in selected))
    (output / "corrections.json").write_text(json.dumps(records, indent=2)+"\n")
    (output / "SOURCES.json").write_text(json.dumps({"inputs": hashes, "excluded_ids": exclude.split(",") if exclude else [],
        "residue": residue, "peptide_chain": peptide_chain, "max_rmsd_A": .002}, indent=2)+"\n")
