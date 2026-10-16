"""Exact, fit-free table exports for benchmark and figure producers.

Input structures remain unchanged: receptor referencing either substitutes identities
on the same contact map or truncates a copied structure and rebuilds contacts.
Each export uses the existing numerical API, retaining its conventions.
"""
from __future__ import annotations

from dataclasses import asdict
import hashlib
import json
from importlib.metadata import PackageNotFoundError, version
import math
from pathlib import Path
import platform
import shutil
import subprocess

import numpy as np
import polars as pl

from .contactmap import ContactMap
from .potential import Potential

RECEPTOR_REGIONS = {"all": None, "cdr": ("CDR1", "CDR2", "CDR3"),
                    "cdr12": ("CDR1", "CDR2"), "cdr3": ("CDR3",)}
QUANTITIES = ("receptor-reference", "position", "residue", "contacts", "placement", "pmhc")


def receptor_reference(contact_map: ContactMap, peptide_potential: Potential,
                       mhc_potential: Potential, *, structure=None,
                       cutoff: float = 5.0) -> dict:
    """Φ(real receptor) − Φ(poly-Ala receptor), in each potential's units.

    Without ``structure``, all/CDR/CDR12/CDR3 substitutions retain the contact map.
    With ``structure``, truncate the selected receptor residues to alanine and rebuild
    contacts at ``cutoff``. Neither mode is a calibrated binding free energy.
    """
    from .pipeline import _interface_energy

    rebuilt = {}
    if structure is not None:
        from .refine.substitute import substitute_residues
        for name, regions in RECEPTOR_REGIONS.items():
            targets = {}
            for chain in structure.chains:
                if chain.chain_type not in ("TRA", "TRB"):
                    continue
                residues = (chain.residues if regions is None else
                            [r for region in chain.regions if region.region_type in regions
                             for r in region.residues])
                targets.update({(chain.chain_id, r.seq_index): "A" for r in residues})
            rebuilt[name] = ContactMap.from_structure(
                substitute_residues(structure, targets), cutoff=cutoff)
    row = {"pdb.id": contact_map.pdb_id}
    for interface, suffix, potential in (("tcr_peptide", "tcr_pep", peptide_potential),
                                          ("tcr_mhc", "tcr_mhc", mhc_potential)):
        contacts = contact_map.interface(interface)
        real = _interface_energy(contacts, potential)
        for name, regions in RECEPTOR_REGIONS.items():
            replacement = (pl.lit("A") if regions is None else
                           pl.when(pl.col("region.type.from").is_in(regions))
                           .then(pl.lit("A")).otherwise(pl.col("residue.aa.from")))
            reference = (rebuilt[name].interface(interface) if structure is not None else
                         contacts.with_columns(replacement.alias("residue.aa.from")))
            row[f"dFT_{suffix}_{name}"] = real - _interface_energy(reference, potential)
        if interface == "tcr_peptide":
            row["n_tcr_res_tp"] = contacts.select("chain.id.from", "residue.index.from").n_unique()
    for name in RECEPTOR_REGIONS:
        row[f"dFT_total_{name}"] = row[f"dFT_tcr_pep_{name}"] + row[f"dFT_tcr_mhc_{name}"]
    return row


def residue_identity(structure) -> pl.DataFrame:
    """Every physical residue, including noncontacting and unassigned residues."""
    rows = []
    for chain in structure.chains:
        regions = {r.seq_index: reg.region_type for reg in chain.regions for r in reg.residues}
        for pos, residue in enumerate(chain.residues, 1):
            rows.append({"pdb.id": structure.pdb_id, "chain.id": chain.chain_id,
                         "chain.type": chain.chain_type, "residue.index": residue.seq_index,
                         "residue.index.pdb": residue.pdb_index,
                         "insertion_code": residue.insertion_code, "residue.aa": residue.aa,
                         "region.type": regions.get(residue.seq_index),
                         "peptide.pos": pos if chain.chain_type == "PEPTIDE" else None})
    return pl.DataFrame(rows, infer_schema_length=None)


def placement_table(structure) -> pl.DataFrame:
    """Existing placement plus the figure's MHC-only signed normal, distances in Å.

    The normal is the MHCa-centroid→peptide-centroid vector orthogonalized against
    the peptide N→C PCA axis. It is independent of the receptor's side of the MHC.
    """
    from .docking.angles import _cdr_centroid, _chain_ca, _groove_frame, tcr_placement

    row = {"pdb.id": structure.pdb_id, **asdict(tcr_placement(structure))}
    u, _, _ = _groove_frame(structure)
    pc = _chain_ca(structure, ("PEPTIDE",)).mean(axis=0)
    mhc = _chain_ca(structure, ("MHCa",))
    if not len(mhc):
        raise ValueError("MHC-normal placement requires MHCa Cα atoms")
    v = pc - mhc.mean(axis=0)
    normal = v - np.dot(v, u) * u
    norm = np.linalg.norm(normal)
    if not np.isfinite(norm) or norm < 1e-12:
        raise ValueError("MHC-normal placement is degenerate")
    normal /= norm
    heights = (_chain_ca(structure, ("TRA", "TRB")) - pc) @ normal
    row.update(tcr_low=float(heights.min()), tcr_high=float(heights.max()),
               cdr_height_mhc=float(np.dot(_cdr_centroid(structure) - pc, normal)))
    for name, vector in (("peptide_center", pc), ("groove_axis", u), ("mhc_normal", normal)):
        row.update({f"{name}_{axis}": float(value) for axis, value in zip("xyz", vector)})
    return pl.DataFrame([row])


def export_structure(structure, quantity: str, peptide_potential: Potential,
                     mhc_potential: Potential, *, cutoff: float = 5.0,
                     interface: str = "tcr_peptide", reference_mode: str = "virtual") -> pl.DataFrame:
    """Export one annotated structure using the exact existing scientific APIs."""
    from .energetics.scoring import position_profile
    from .viz.pymol import residue_importance

    if reference_mode not in ("virtual", "structural"):
        raise ValueError("reference_mode must be virtual or structural")
    if reference_mode == "structural" and quantity != "receptor-reference":
        raise ValueError("structural reference requires receptor-reference quantity")
    if quantity not in QUANTITIES:
        raise ValueError(f"quantity must be one of {QUANTITIES}")
    if not math.isfinite(cutoff) or cutoff <= 0:
        raise ValueError("cutoff must be finite and positive")
    if quantity == "pmhc":
        from .topology.pmhc import pmhc_features
        return pl.DataFrame([{"pdb.id": structure.pdb_id, **pmhc_features(structure)}])
    types = {c.chain_type for c in structure.chains}
    if types & {"TRG", "TRD", "IGH", "IGK", "IGL"}:
        raise ValueError("export is restricted to alpha-beta TCR complexes")
    if not {"TRA", "TRB", "PEPTIDE", "MHCa"} <= types:
        raise ValueError("export requires annotated TRA/TRB, peptide and MHCa chains")
    peptide = structure.by_type("PEPTIDE")
    if len(peptide) != 1 or set(peptide[0].sequence()) - set("ACDEFGHIKLMNPQRSTVWY"):
        raise ValueError("export requires exactly one standard-amino-acid peptide")
    if quantity == "placement":
        return placement_table(structure)
    cm = ContactMap.from_structure(structure, cutoff=cutoff)
    potential = peptide_potential if interface == "tcr_peptide" else mhc_potential
    if quantity == "receptor-reference":
        return pl.DataFrame([receptor_reference(cm, peptide_potential, mhc_potential,
                            structure=structure if reference_mode == "structural" else None,
                            cutoff=cutoff)]).with_columns(pl.lit(reference_mode).alias("reference_mode"))
    if quantity == "position":
        if interface == "tcr_mhc":
            raise ValueError("position decomposition requires a peptide interface")
        table = position_profile(cm, potential, structure=structure, interface=interface)
        identity = residue_identity(structure).filter(pl.col("chain.type") == "PEPTIDE")
        return table.join(identity.select("peptide.pos", "chain.id", "residue.index",
                                         "residue.index.pdb", "insertion_code"), on="peptide.pos",
                          how="left", validate="1:1").with_columns(pl.lit(structure.pdb_id).alias("pdb.id"))
    if quantity == "residue":
        table = residue_importance(structure, interface=interface, cutoff=cutoff, potential=potential)
        return table.join(residue_identity(structure).select("chain.id", "residue.index",
                         "residue.index.pdb", "insertion_code", "chain.type"),
                         on=["chain.id", "residue.index"], how="left", validate="1:1").with_columns(
                             pl.lit(structure.pdb_id).alias("pdb.id"))
    table = cm.interface(interface)
    identity = residue_identity(structure)
    for side in ("from", "to"):
        keys = ["chain.id", "residue.index", "residue.index.pdb", "insertion_code"]
        table = table.join(identity.select(keys).rename({k: f"{k}.{side}" for k in keys}),
                           on=[f"chain.id.{side}", f"residue.index.{side}"],
                           how="left", validate="m:1")
    return table


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def annotation_lineage(organisms) -> dict:
    """Record actual installed versions and raw annotation-reference identities."""
    from arda.paths import vdj_dir
    from .mhc.reference import reference_fasta

    versions = {"python": platform.python_version()}
    for package in ("tcren", "arda-mapper", "numpy", "polars", "biopython"):
        try:
            versions[package] = version(package)
        except PackageNotFoundError:
            versions[package] = "unavailable"
    executable = shutil.which("mmseqs")
    mmseqs = (subprocess.run([executable, "version"], capture_output=True, text=True,
                            timeout=10, check=True).stdout.strip() if executable else "unavailable")
    refs = {}
    for organism in organisms:
        root = Path(vdj_dir(organism))
        refs[organism] = {p.name: _sha256(p) for p in sorted(root.iterdir())
                         if p.is_file() and p.suffix in {".fasta", ".tsv", ".json"}}
    refs["mhc"] = _sha256(Path(reference_fasta()))
    return {"versions": versions, "mmseqs_version": mmseqs, "reference_sha256": refs,
            "organisms": list(organisms), "annotation_policy": ("batched species autodetection" if len(organisms) > 1
                                  else "batched fixed species")}


from .runtime import resolve_threads, single_threaded


@single_threaded
def export_batch(sources, output, *, quantity, peptide_potential, mhc_potential,
                 cutoff=5.0, interface="tcr_peptide", organisms=("human", "mouse"), threads=1,
                 command="tcren export", reference_mode="virtual"):
    """Write exact tables plus complete per-input accounting and residue assignments.

    Inputs are files, directories, tar archives, globs or text manifests. Errors remain in the
    accounting sidecar; the caller must return nonzero when no structure succeeds.
    """
    from .annotation import classify_chains
    from .annotation.batch import annotate_batch
    from .mhc import annotate_mhc_batch
    from .provenance import annotation_indexes, stamp
    from .structure import iter_structures, parse_structure, structure_paths
    from .structure.io import _TAR_SUFFIXES, resolve_sources, structure_stem

    if reference_mode not in ("virtual", "structural"):
        raise ValueError("reference_mode must be virtual or structural")
    if reference_mode == "structural" and quantity != "receptor-reference":
        raise ValueError("structural reference requires receptor-reference quantity")
    threads = resolve_threads(threads)
    if not organisms or set(organisms) - {"human", "mouse"}:
        raise ValueError("positive threads and human/mouse organisms are required")
    if quantity not in QUANTITIES or interface not in {"tcr_peptide", "tcr_mhc", "peptide_mhc"}:
        raise ValueError("invalid export quantity or interface")
    if not math.isfinite(cutoff) or cutoff <= 0:
        raise ValueError("cutoff must be finite and positive")
    paths = [p for source in resolve_sources(sources) for p in structure_paths(source)]
    ids = {str(p): structure_stem(p) for p in paths}
    if len(set(ids.values())) != len(paths):
        raise ValueError("duplicate structure IDs in export inputs")
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    Path(str(output) + ".annotation.json").write_text("{}\n")
    accounting, parsed, tables, identities = [], [], [], []
    def parse_one(path, pdb_id, source, archive_hash=None):
        row = {"pdb.id": pdb_id, "source": source, "sha256": None,
               "archive_sha256": archive_hash, "status": "error", "error": "", "n_rows": 0}
        accounting.append(row)
        try:
            row["sha256"] = _sha256(path)
            structure = parse_structure(path, pdb_id=pdb_id)
            parsed.append((structure, row))
            return structure
        except Exception as exc:
            row["error"] = f"parse: {type(exc).__name__}: {exc}"
            return None

    for path in paths:
        if str(path).lower().endswith(_TAR_SUFFIXES) and path.is_file():
            try:
                archive_hash = _sha256(path)
                def importer(member_path, pdb_id):
                    return parse_one(member_path, pdb_id, f"{path}#{pdb_id}", archive_hash)
                for _ in iter_structures(path, importer=importer, full_ids=True):
                    pass
            except Exception as exc:
                accounting.append({"pdb.id": ids[str(path)], "source": str(path),
                                   "sha256": None, "archive_sha256": None, "status": "error",
                                   "error": f"archive: {type(exc).__name__}: {exc}", "n_rows": 0})
        else:
            parse_one(path, ids[str(path)], str(path))
    if len({row["pdb.id"] for row in accounting}) != len(accounting):
        raise ValueError("duplicate structure IDs in export inputs or archive members")
    lineage = {}
    if parsed:
        try:
            records = annotate_batch([s for s, _ in parsed], organisms=organisms, threads=threads)
            if len(records) != len(parsed):
                raise ValueError("annotation output count differs from parsed input count")
            Path(str(output) + ".annotation.json").write_text(json.dumps(
                {s.pdb_id: record for (s, _), record in zip(parsed, records)},
                indent=1, default=str) + "\n")
            typed = []
            for (structure, row), record in zip(parsed, records):
                try:
                    classify_chains(structure, organism=organisms[0],
                                    autodetect_species=len(organisms) > 1, precomputed_records=record)
                    typed.append((structure, row))
                except Exception as exc:
                    row["error"] = f"typing: {type(exc).__name__}: {exc}"
            annotate_mhc_batch([s for s, _ in typed], threads=threads)
            lineage = annotation_lineage(organisms)
            for structure, row in typed:
                try:
                    identity = residue_identity(structure)
                    identities.append(identity)
                    table = export_structure(structure, quantity, peptide_potential, mhc_potential,
                                             cutoff=cutoff, interface=interface,
                                             reference_mode=reference_mode)
                    if "pdb.id" not in table.columns:
                        table = table.with_columns(pl.lit(structure.pdb_id).alias("pdb.id"))
                    tables.append(table)
                    row.update(status="ok", n_rows=table.height)
                except Exception as exc:
                    row["error"] = f"export: {type(exc).__name__}: {exc}"
        except Exception as exc:
            for _, row in parsed:
                if not row["error"]:
                    row["error"] = f"annotation: {type(exc).__name__}: {exc}"
    table = pl.concat(tables, how="diagonal_relaxed") if tables else pl.DataFrame(schema={"pdb.id": pl.String})
    table.write_csv(output, separator="\t")
    ledger = pl.DataFrame(accounting, infer_schema_length=None) if accounting else pl.DataFrame(
        schema={"pdb.id": pl.String, "status": pl.String, "error": pl.String, "n_rows": pl.Int64})
    ledger.write_csv(str(output) + ".accounting.tsv", separator="\t")
    residue_table = pl.concat(identities, how="diagonal_relaxed") if identities else pl.DataFrame(
        schema={"pdb.id": pl.String})
    residue_table.write_csv(str(output) + ".residues.tsv", separator="\t")
    potentials = {}
    for label, potential in (() if quantity == "pmhc" else (("peptide", peptide_potential), ("mhc", mhc_potential))):
        path = Path(str(output) + f".{label}.potential.tsv")
        potential.matrix.write_csv(path, separator="\t")
        potentials[label] = {"name": potential.name, "sha256": _sha256(path),
                             "units": ("dimensionless contact log-odds" if potential.name in {"TCRen", "TCRen2", "TCRen2_em"}
                                       else "dimensionless MJ table scale" if potential.name == "MJ"
                                       else "input table units (not declared by Potential)"),
                             "thermodynamic_calibration": "none; not kcal/mol"}
    succeeded = sum(row["status"] == "ok" for row in accounting)
    from importlib.resources import files
    manifest_text = files("tcren.data").joinpath("potentials.json").read_text()
    stamp(output, command=command, columns=table.columns, extra={
        **lineage, "searched_indexes": annotation_indexes(s for s, _ in parsed),
        "organisms": list(organisms), "quantity": quantity, "interface": interface, "cutoff_A": cutoff,
        "threads": threads, "potentials": potentials, "bundled_potential_manifest": json.loads(manifest_text),
        "custom_potential_derivation": "not inferred; supplied matrix is exported and hashed",
        "n_input": len(accounting), "n_ok": succeeded,
        "n_error": len(accounting) - succeeded, "sources": [str(p) for p in paths],
        "reference_semantics": ("pMHC only; supplied conformation; no potential or receptor coordinates" if quantity == "pmhc"
                                else "structural alanine truncation; contacts rebuilt" if reference_mode == "structural"
                                else "virtual receptor identities on unchanged contact map"),
        "residue_phi_accounting": "each contact contributes to both endpoints; sum is 2*Phi"})
    return ledger
