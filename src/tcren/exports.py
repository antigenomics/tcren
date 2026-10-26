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
from .contacts.observations import CONFORMER_POLICY
from .potential import Potential

RECEPTOR_REGIONS = {"all": None, "cdr": ("CDR1", "CDR2", "CDR3"),
                    "cdr12": ("CDR1", "CDR2"), "cdr3": ("CDR3",)}
QUANTITIES = ("cdr3-neighborhood", "cdr3-contacts", "receptor-reference", "position", "residue", "contacts", "contact-context", "contact-atoms", "atom-inventory", "placement", "contact-topology", "recognition-topology", "recognition-betti", "weighted-recognition-topology", "persistence", "backbone", "pmhc", "pmhc-cells", "contact-exposure", "surface-contacts")


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


def backbone_table(structure) -> pl.DataFrame:
    """Intrinsic peptide/CDR C-alpha geometry, without fitting or a strain threshold.

    Windows retain physical residue identities. Missing C-alpha or discontinuous
    observed indices split windows instead of bridging missing coordinates.
    """
    import numpy as np
    from .geometry import _bond_angle
    from .docking.tcrdock_geometry import _dihedral

    rows = []
    for chain in structure.chains:
        regions = ([('PEPTIDE', chain.residues)] if chain.chain_type == 'PEPTIDE' else
                   [(r.region_type, r.residues) for r in chain.regions
                    if chain.chain_type in ('TRA', 'TRB') and r.region_type in ('CDR1', 'CDR2', 'CDR3')])
        for region, residues in regions:
            for size, metric, unit in ((2, 'ca_bond', 'angstrom'), (3, 'ca_angle', 'radian'),
                                       (4, 'ca_torsion', 'radian')):
                for start in range(len(residues)-size+1):
                    window = residues[start:start+size]
                    if any(r.ca is None for r in window) or any(
                            b.seq_index != a.seq_index+1 for a, b in zip(window, window[1:])):
                        continue
                    xyz = np.asarray([r.ca for r in window], dtype=float)
                    if np.any(np.linalg.norm(np.diff(xyz, axis=0), axis=1) <= 1e-8):
                        raise ValueError('Coincident C-alpha atoms in backbone window')
                    value = (float(np.linalg.norm(xyz[1]-xyz[0])) if size == 2 else
                             _bond_angle(*xyz) if size == 3 else _dihedral(*xyz))
                    row = {'chain.id': chain.chain_id, 'chain.type': chain.chain_type,
                           'region.type': region, 'metric': metric, 'units': unit, 'value': value}
                    for i, residue in enumerate(window, 1):
                        row.update({f'residue.index.{i}': residue.seq_index,
                                    f'residue.index.pdb.{i}': residue.pdb_index,
                                    f'insertion_code.{i}': residue.insertion_code,
                                    f'residue.aa.{i}': residue.aa})
                    rows.append(row)
    return (pl.DataFrame(rows, infer_schema_length=None) if rows else pl.DataFrame(
        schema={'chain.id': pl.String, 'chain.type': pl.String, 'region.type': pl.String,
                'metric': pl.String, 'units': pl.String, 'value': pl.Float64}))


def contact_topology(contact_map: ContactMap) -> pl.DataFrame:
    """Separate residue-contact graphs: CDR:peptide, CDR:MHC and peptide:MHC.

    Vertices are contacted residues (physical chain + observed index); edges are
    unique heavy-atom contacts. Betti numbers are graph components and cycles,
    not the Betti numbers of a spatial flag complex. Empty graphs have V=E=β0=β1=0.
    """
    import numpy as np
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components

    rows = []
    for interface in ("tcr_peptide", "tcr_mhc", "peptide_mhc"):
        d = contact_map.interface(interface)
        if interface.startswith("tcr_"):
            d = d.filter(pl.col("region.type.from").is_in(["CDR1", "CDR2", "CDR3"]))
        keys = [f"{key}.{side}" for side in ("from", "to")
                for key in ("chain.id", "residue.index")]
        edges = d.select(keys).unique().sort(keys)
        vertices = sorted({(c, int(r)) for c, r, _, _ in edges.iter_rows()} |
                          {(c, int(r)) for _, _, c, r in edges.iter_rows()})
        index = {key: i for i, key in enumerate(vertices)}
        n, e = len(vertices), edges.height
        if e:
            pairs = np.array([(index[(a, int(i))], index[(b, int(j))])
                              for a, i, b, j in edges.iter_rows()], dtype=int)
            graph = coo_matrix((np.ones(e), (pairs[:, 0], pairs[:, 1])), shape=(n, n))
            b0 = int(connected_components(graph, directed=False, return_labels=False))
        else:
            b0 = 0
        rows.append(dict(interface=interface, vertices=n, edges=e,
                         beta0=b0, beta1=e-n+b0, available=bool(e)))
    return pl.DataFrame(rows)


def recognition_betti(contact_map: ContactMap) -> pl.DataFrame:
    """Raw complex cycle fraction, keeping receptor-contact availability separate."""
    from .topology.graph import recognition_graph_features
    row = recognition_graph_features(contact_map)
    edges = row["g_recognition_edges"]
    row["complex_cyclo_frac"] = row["g_recognition_beta1"] / edges if edges else float("nan")
    row["receptor_contacts_available"] = bool(np.isfinite(row["g_recognition_cyclo_frac"]))
    return pl.DataFrame([row])


def _recognition_graphs(structure, contact_map, *, cutoff=5.0):
    """Shared physical edge ladder; preserve the released graph definitions."""
    from .contacts.geometry import all_atom_contacts
    from .structure.model import Structure
    from .topology.graph import _recognition_contact_edges

    # Reuse the native contact search on peptide atoms only. No sequence-distance
    # exclusion: adjacent residues qualify only through side-chain heavy atoms.
    internal = all_atom_contacts(Structure("peptide", structure.by_type("PEPTIDE")),
                                 cutoff=cutoff, scope="intra", sidechain=True)
    internal = internal.filter(pl.col("n_sc_pairs") > 0)
    peptide_contacts = {tuple(sorted(((a, int(i)), (b, int(j)))))
                        for a, i, b, j in internal.select(
                            "chain.id.from", "residue.index.from",
                            "chain.id.to", "residue.index.to").iter_rows()}

    parts = _recognition_contact_edges(contact_map)
    sequence = set()
    missing = 0
    for chain in structure.by_type("PEPTIDE"):
        for left, right in zip(chain.residues, chain.residues[1:]):
            c = next((a.coord for a in left.atoms if a.name == "C"), None)
            n = next((a.coord for a in right.atoms if a.name == "N"), None)
            if (right.seq_index != left.seq_index + 1 or c is None or n is None
                    or not 0 < np.linalg.norm(c-n) <= 1.8):
                missing += 1
                continue
            sequence.add(tuple(sorted(((chain.chain_id, left.seq_index),
                                       (chain.chain_id, right.seq_index)))))
    tp, tm, pm = (parts[n] for n in ("tcr_peptide", "tcr_mhc", "peptide_mhc"))
    graphs = (("released", tp | tm), ("tp_pm", tp | pm),
                        ("tp_pm_sequence", tp | pm | sequence),
                        ("tp_pm_sequence_tm", tp | pm | sequence | tm),
                        ("tp_pm_sidechain", tp | pm | peptide_contacts),
                        ("tp_pm_sidechain_tm", tp | pm | peptide_contacts | tm),
                        ("tp_pm_tm", tp | pm | tm))
    return graphs, parts, sequence, peptide_contacts, missing


def recognition_topology(structure, contact_map: ContactMap, *, cutoff: float = 5.0) -> pl.DataFrame:
    """Fixed exploratory graph ladder; beta1/E is redundancy, not rupture energy.

    Peptide sequence edges connect consecutive observed residues only when their
    backbone C--N distance is <=1.8 Angstrom. Missing atoms/breaks remain explicit.
    All nodes have physical chain/index keys. No TCR or MHC backbone edges.
    """
    from .topology.graph import _clique_homology, _contact_graph_counts

    graphs, parts, sequence, peptide_contacts, missing = _recognition_graphs(
        structure, contact_map, cutoff=cutoff)
    tp, tm, pm = (parts[n] for n in ("tcr_peptide", "tcr_mhc", "peptide_mhc"))
    rows = []
    for name, edges in graphs:
        e, v, components, beta1 = _contact_graph_counts(edges)
        homology = _clique_homology(edges, beta1)
        rows.append(dict(graph=name, vertices=v, edges=e, beta0=components,
                         beta1=beta1, cyclo_frac=beta1/e if e else None,
                         tcr_peptide_edges=len(tp), tcr_mhc_edges=len(tm),
                         peptide_mhc_edges=len(pm), sequence_edges=len(sequence),
                         peptide_sidechain_edges=len(peptide_contacts),
                         sequence_breaks=missing, available=bool(tp or tm), **homology))
    return pl.DataFrame(rows)


def export_structure(structure, quantity: str, peptide_potential: Potential,
                     mhc_potential: Potential, *, cutoff: float = 5.0,
                     interface: str = "tcr_peptide", reference_mode: str = "virtual",
                     pmhc_potential: Potential | None = None,
                     graph_weighting="positive", graph_beta=1.0) -> pl.DataFrame:
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
    if quantity == "surface-contacts":
        from .contacts.exposure import contact_surface
        return contact_surface(structure)[1]
    if quantity == "contact-exposure":
        from .contacts.exposure import context_exposure
        return context_exposure(structure)
    if quantity == "pmhc-cells":
        from .topology.pmhc import pmhc_cells
        return pmhc_cells(structure)
    if quantity == "pmhc":
        from .topology.pmhc import pmhc_features
        return pl.DataFrame([{"pdb.id": structure.pdb_id, **pmhc_features(structure)}])
    types = {c.chain_type for c in structure.chains}
    if types & {"TRG", "TRD", "IGH", "IGK", "IGL"}:
        raise ValueError("export is restricted to alpha-beta TCR complexes")
    if quantity == "backbone":
        if not ({"TRA", "TRB"} <= types or "PEPTIDE" in types):
            raise ValueError("backbone export requires annotated alpha-beta TCR or peptide")
        return backbone_table(structure)
    if not {"TRA", "TRB", "PEPTIDE", "MHCa"} <= types:
        raise ValueError("export requires annotated TRA/TRB, peptide and MHCa chains")
    peptide = structure.by_type("PEPTIDE")
    if len(peptide) != 1 or set(peptide[0].sequence()) - set("ACDEFGHIKLMNPQRSTVWY"):
        raise ValueError("export requires exactly one standard-amino-acid peptide")
    if quantity == "placement":
        return placement_table(structure)
    if quantity == "atom-inventory":
        from .contacts.observations import atom_inventory
        return atom_inventory(structure)
    if quantity == "contact-atoms":
        from .contacts.observations import contact_observations
        return contact_observations(structure, cutoff=cutoff, interface=interface)
    if quantity == "cdr3-contacts":
        from .topology.neighborhood import neighborhood_contacts
        return neighborhood_contacts(structure, cutoff=cutoff)
    if quantity == "cdr3-neighborhood":
        from .topology.neighborhood import neighborhood_features
        return pl.DataFrame([neighborhood_features(structure, cutoff=cutoff)])
    cm = ContactMap.from_structure(structure, cutoff=cutoff)
    if quantity == "contact-topology":
        return contact_topology(cm)
    if quantity == "recognition-betti":
        return recognition_betti(cm)
    if quantity == "recognition-topology":
        return recognition_topology(structure, cm, cutoff=cutoff)
    if quantity == "weighted-recognition-topology":
        from .topology.weighted import weighted_graph_features
        return weighted_graph_features(structure, cm, cutoff=cutoff,
            peptide_potential=peptide_potential, mhc_potential=mhc_potential,
            pmhc_potential=pmhc_potential, weighting=graph_weighting, beta=graph_beta)
    if quantity == "persistence":
        from .topology.persistence import persistence_features
        return pl.DataFrame([persistence_features(structure, contact_map=cm)])
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
    if quantity == "contact-context" and interface != "tcr_peptide":
        raise ValueError("contact-context requires tcr_peptide")
    table = cm.interface(interface)
    if quantity == "contact-context":
        from .topology.contact_context import contact_context
        table = contact_context(structure, table, potential)
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
                 command="tcren export", reference_mode="virtual", register_reference_id=None,
                 register_p1=None, peptide_selection=None, pmhc_potential=None,
                 graph_weighting="positive", graph_beta=1.0):
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
    if (register_reference_id is None) != (register_p1 is None):
        raise ValueError("register reference ID and one-based reference P1 are required together")
    selections = {}
    if peptide_selection is not None:
        selection = pl.read_csv(peptide_selection, separator="\t",
                                schema_overrides={"first_pdb": pl.Int64, "last_pdb": pl.Int64})
        required = ["pdb.id", "chain.id", "first_pdb", "last_pdb"]
        if not set(required) <= set(selection.columns) or selection["pdb.id"].n_unique() != selection.height:
            raise ValueError("peptide selection requires unique pdb.id, chain.id, first_pdb, last_pdb")
        selections = {r["pdb.id"]: r for r in selection.iter_rows(named=True)}
    threads = resolve_threads(threads)
    if not organisms or set(organisms) - {"human", "mouse"}:
        raise ValueError("positive threads and human/mouse organisms are required")
    if (quantity not in QUANTITIES or interface not in {"tcr_peptide", "tcr_mhc", "peptide_mhc", "all"}
            or (interface == "all" and quantity not in {"contact-atoms", "atom-inventory", "cdr3-neighborhood", "cdr3-contacts"})):
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
    accounting, parsed, tables, identities, surfaces = [], [], [], [], []
    def parse_one(path, pdb_id, source, archive_hash=None):
        row = {"pdb.id": pdb_id, "source": source, "sha256": None,
               "archive_sha256": archive_hash, "status": "error", "error": "", "n_rows": 0}
        accounting.append(row)
        try:
            row["sha256"] = _sha256(path)
            structure = (parse_structure(path, pdb_id=pdb_id, keep_modified=True)
                         if register_reference_id is not None or quantity in {"contact-atoms", "atom-inventory"}
                         else parse_structure(path, pdb_id=pdb_id))
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
    if quantity == "weighted-recognition-topology":
        from .potential import mj, tcren2
        peptide_potential = peptide_potential if peptide_potential is not None else tcren2()
        mhc_potential = mhc_potential if mhc_potential is not None else mj()
    lineage, registers = {}, []
    if register_reference_id is not None and register_reference_id not in {s.pdb_id for s, _ in parsed}:
        raise ValueError("register reference must be included in the same input batch")
    if selections.keys() - {s.pdb_id for s, _ in parsed}:
        raise ValueError("peptide selection contains IDs absent from input batch")
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
            from .refine.structural_register import select_peptide, structural_register
            for structure, row in typed:
                if structure.pdb_id in selections:
                    r = selections[structure.pdb_id]
                    select_peptide(structure, r["chain.id"], int(r["first_pdb"]), int(r["last_pdb"]))
            reference = next((s for s, _ in typed if s.pdb_id == register_reference_id), None)
            if register_reference_id is not None and reference is None:
                raise ValueError("register reference annotation failed")
            lineage = annotation_lineage(organisms)
            for structure, row in typed:
                try:
                    if reference is not None and any(c.chain_supertype == "MHCII" for c in structure.chains):
                        registers.append(structural_register(structure, reference, register_p1))
                    identity = residue_identity(structure)
                    identities.append(identity)
                    if quantity == "surface-contacts":
                        from .contacts.exposure import contact_surface
                        atom_surface, table = contact_surface(structure)
                        surfaces.append(atom_surface)
                    else:
                        table = export_structure(structure, quantity, peptide_potential, mhc_potential,
                                                 cutoff=cutoff, interface=interface,
                                                 reference_mode=reference_mode, pmhc_potential=pmhc_potential,
                                                 graph_weighting=graph_weighting, graph_beta=graph_beta)
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
    if quantity in {"contact-atoms", "atom-inventory", "contact-exposure", "surface-contacts"} and output.suffix == ".parquet":
        table.write_parquet(output)
    else:
        table.write_csv(output, separator="\t")
    if quantity == "surface-contacts":
        (pl.concat(surfaces, how="diagonal_relaxed") if surfaces else pl.DataFrame(schema={"pdb.id": pl.String})).write_parquet(str(output) + ".atoms.parquet")
    if quantity == "pmhc-cells" and table.height:
        from .topology.pmhc import spatial_features
        spatial_features(table).write_csv(str(output) + ".summary.tsv", separator="\t")
    if quantity == "contact-atoms" and tables:
        from .contacts.observations import summarize_contact_observations
        summarize_contact_observations(table).write_csv(str(output) + ".summary.tsv", separator="\t")
    ledger = pl.DataFrame(accounting, infer_schema_length=None) if accounting else pl.DataFrame(
        schema={"pdb.id": pl.String, "status": pl.String, "error": pl.String, "n_rows": pl.Int64})
    ledger.write_csv(str(output) + ".accounting.tsv", separator="\t")
    residue_table = pl.concat(identities, how="diagonal_relaxed") if identities else pl.DataFrame(
        schema={"pdb.id": pl.String})
    residue_table.write_csv(str(output) + ".residues.tsv", separator="\t")
    if register_reference_id is not None:
        (pl.concat(registers) if registers else pl.DataFrame(schema={"pdb.id": pl.String})).write_csv(
            str(output) + ".registers.tsv", separator="\t")
    potentials = {}
    for label, potential in (() if quantity in {"pmhc", "pmhc-cells", "contact-topology", "recognition-topology", "recognition-betti", "cdr3-neighborhood", "cdr3-contacts", "persistence", "backbone", "contact-atoms", "atom-inventory", "contact-exposure", "surface-contacts"} else (("peptide", peptide_potential), ("mhc", mhc_potential))):
        path = Path(str(output) + f".{label}.potential.tsv")
        potential.matrix.write_csv(path, separator="\t")
        potentials[label] = {"name": potential.name, "sha256": _sha256(path),
                             "units": ("dimensionless contact log-odds" if potential.name in {"TCRen", "TCRen2", "TCRen2_em"}
                                       else "dimensionless MJ table scale" if potential.name == "MJ"
                                       else "input table units (not declared by Potential)"),
                             "thermodynamic_calibration": "none; not kcal/mol"}
    if pmhc_potential is not None and quantity == "weighted-recognition-topology":
        path = Path(str(output) + ".pmhc.potential.tsv")
        pmhc_potential.matrix.write_csv(path, separator="\t")
        potentials["pmhc"] = {"name": pmhc_potential.name, "sha256": _sha256(path),
                              "units": "input table units; see supplied derivation",
                              "thermodynamic_calibration": "none; not kcal/mol"}
    succeeded = sum(row["status"] == "ok" for row in accounting)
    from importlib.resources import files
    manifest_text = files("tcren.data").joinpath("potentials.json").read_text()
    stamp(output, command=command, columns=table.columns, extra={
        **lineage, "searched_indexes": annotation_indexes(s for s, _ in parsed),
        "organisms": list(organisms), "quantity": quantity, "interface": "all" if quantity == "atom-inventory" else interface, "cutoff_A": None if quantity == "atom-inventory" else cutoff,
        "threads": threads, "potentials": potentials, "bundled_potential_manifest": json.loads(manifest_text),
        "custom_potential_derivation": "not inferred; supplied matrix is exported and hashed",
        "n_input": len(accounting), "n_ok": succeeded,
        "n_error": len(accounting) - succeeded, "sources": [str(p) for p in paths],
        "weighted_topology": ({"goodness": "exp(-beta*energy)" if graph_weighting == "boltzmann" else "max(-energy,0)",
                              "weighting": graph_weighting, "beta": graph_beta,
                              "orientation": ("TCR to peptide; TCR to MHC; peptide to MHC" if pmhc_potential is not None
                                              else "TCR/MHC to peptide; TCR to MHC"),
                              "per_vertex": "E/V and W/V; bounded E/(E+V), W/(W+V); V incident residues",
                              "missing": "edge retained, weighted descriptors withheld",
                              "internal_peptide": "unweighted; missing_weight incremented"}
                             if quantity == "weighted-recognition-topology" else None),
        "persistence": ({"backend": "gudhi", "field": "GF(2)", "max_edge_A": 12.,
                         "skeleton": 2, "atoms": "CA and CB; Gly CA fallback",
                         "censoring": "finite bars only; infinite counts retained"}
                        if quantity == "persistence" else None),
        "reference_semantics": ("directed conserved surface allocation; intact target; separate endpoint areas; no score" if quantity == "surface-contacts"
                                else "atom SASA in isolated binder/intact pMHC/bound complex; unchanged coordinates; no energy" if quantity == "contact-exposure"
                                else "selected observed heavy atoms across full protein components; zero-contact atoms retained; no energy" if quantity == "atom-inventory"
                                else "observed heavy-atom pairs; no energy or fitted model" if quantity == "contact-atoms"
                                else "all heavy-atom CDR3 peptide and complete MHC contacts; fixed v2 chemistry, Bondi overlaps, separate interfaces; no fitted potential" if quantity in {"cdr3-neighborhood", "cdr3-contacts"}
                                else "intrinsic C-alpha geometry; not a strain score" if quantity == "backbone"
                                else "separate residue-contact graphs; CDR-only receptor; no fitted score" if quantity in {"contact-topology", "recognition-topology", "recognition-betti"}
                                else "pMHC only; supplied conformation; no potential or receptor coordinates" if quantity in {"pmhc", "pmhc-cells"}
                                else "structural alanine truncation; contacts rebuilt" if reference_mode == "structural"
                                else "virtual receptor identities on unchanged contact map"),
        "surface_summary_sha256": (_sha256(Path(str(output) + ".summary.tsv")) if quantity == "pmhc-cells" and table.height else None),
        "contact_summary_sha256": (_sha256(Path(str(output) + ".summary.tsv")) if quantity == "contact-atoms" and tables else None),
        "contact_observations": (dict(conformer_policy=CONFORMER_POLICY,
            scope="CDR1/2/3, peptide, annotated MHC helices; interchain only",
            chemistry="existing v2 distance/type heuristics; no angular hydrogen-bond assignment",
            symmetry="one physical atom pair; unordered residue+atom chemical endpoint tuples",
            angles="trRosetta N/CA/observed CB definitions; radians; missing/degenerate null",
            fitted=False) if quantity == "contact-atoms" else None),
        "atom_inventory": ({"conformer_policy": CONFORMER_POLICY,
            "source_index": "zero-based heavy-atom enumeration across all parsed chains before selection; not PDB serial",
            "scope": "TRA/TRB, peptide, full MHCa/MHCb/B2M; parsed observed protein atoms only",
            "noncontact_atoms_retained": True, "units": "angstrom", "model": "parser-selected model",
            "missing_atoms": "not reconstructed", "hydrogens": "H/D excluded",
            "chemical_state": "no waters, ions, bonded topology, charges or protonation inferred",
            "omission_count": "heavy atoms removed by conformer/name selection only; repeated per retained atom"}
            if quantity == "atom-inventory" else None),
        "surface_allocation": ({"n_points": 100, "probe_A": 1.4, "sphere": "native Fibonacci-y", "owner": "minimum signed distance to opposing expanded sphere; exact ties split equally", "atoms_sha256": _sha256(Path(str(output) + ".atoms.parquet")), "energy": "none"} if quantity == "surface-contacts" else None),
        "exposure_parameters": ({"n_points": 100, "probe_A": 1.4, "radii": "Bio.PDB.SASA", "units": "A^2", "hydrogens": "excluded", "frame": "supplied coordinates", "context": "intact peptide/MHC/B2M target; peptide alone only for MHC shielding"} if quantity == "contact-exposure" else None),
        "surface_parameters": ({"grid_yx": [64,32], "extent_A": [-20,20,-25,25], "probe_A": 1.4, "smooth": "8-neighbour, padded", "scale": "Kyte-Doolittle", "ridge_percentile": 90, "frame": "pMHC only; receptor removed"} if quantity == "pmhc-cells" else None),
        "structural_register": ({"reference_id": register_reference_id, "reference_p1": register_p1,
            "method": "MHC-groove-only superposition; forward/reverse nine-residue C-alpha RMSD",
            "sha256": _sha256(Path(str(output) + ".registers.tsv"))} if register_reference_id else None),
        "peptide_selection": ({"sha256": _sha256(Path(peptide_selection)), "rows": selections}
                              if peptide_selection is not None else None),
        "residue_phi_accounting": ("not applicable: atom inventory, no energy" if quantity == "atom-inventory"
                                   else "not applicable: atom surface areas, no energy" if quantity in {"contact-exposure", "surface-contacts"}
                                   else "each contact contributes to both endpoints; sum is 2*Phi")})
    return ledger
