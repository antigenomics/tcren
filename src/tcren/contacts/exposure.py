"""Atom-resolved exposure in unchanged TCR:pMHC coordinates, in square angstroms.

This supplies observed surface opportunity, not a normalized contact potential.
The existing Bio.PDB Shrake-Rupley kernel owns sampling and atomic radii.
"""
from __future__ import annotations

import numpy as np
import polars as pl


def _exposure_atoms(structure):
    """Expose every selected heavy atom, including noncontact atoms.

    Isolated partners are the full TRA/TRB assembly and full peptide/MHC/B2M
    assembly. Peptide-alone exposure additionally measures shielding by its MHC
    context. No relaxation, water energy, fitted weights or residue pooling.
    """
    types = {c.chain_type for c in structure.chains}
    peptide = structure.by_type("PEPTIDE")
    mhc = structure.by_type("MHCa", "MHCb")
    classes = {c.chain_supertype for c in mhc if c.chain_supertype in {"MHCI", "MHCII"}}
    if (not {"TRA", "TRB", "MHCa"} <= types or len(peptide) != 1 or len(classes) != 1
            or types & {"TRG", "TRD", "IGH", "IGK", "IGL"}
            or set(peptide[0].sequence()) - set("ACDEFGHIKLMNPQRSTVWY")
            or (classes == {"MHCII"} and "MHCb" not in types)):
        raise ValueError("exposure requires alpha-beta TCR and one standard peptide in an unambiguous MHC groove")
    selected = structure.by_type("TRA", "TRB", "PEPTIDE", "MHCa", "MHCb", "B2M")
    if len({c.chain_id for c in selected}) != len(selected):
        raise ValueError("duplicate physical chain IDs")
    rows, atoms = [], []
    for chain in selected:
        regions = {r.seq_index: reg.region_type for reg in chain.regions for r in reg.residues}
        for residue in chain.residues:
            seen = set()
            for atom in residue.atoms:
                if atom.element in {"H", "D"} or atom.name in seen:
                    continue
                seen.add(atom.name)  # Same first-conformer atom-name policy as the writer.
                if np.asarray(atom.coord).shape != (3,) or not np.isfinite(atom.coord).all():
                    raise ValueError("finite three-dimensional atom coordinates required")
                atoms.append(atom)
                rows.append({"pdb.id": structure.pdb_id, "chain.id": chain.chain_id,
                    "chain.type": chain.chain_type, "residue.index": residue.seq_index,
                    "residue.index.pdb": residue.pdb_index, "insertion_code": residue.insertion_code,
                    "residue.aa": residue.aa, "resname": residue.resname,
                    "region.type": regions.get(residue.seq_index), "atom.name": atom.name,
                    "element": atom.element, "role": "binder" if chain.chain_type in {"TRA", "TRB"} else "target",
                    "mhc_class": "II" if classes == {"MHCII"} else "I"})
    if not atoms:
        raise ValueError("no heavy atoms in selected assembly")
    table = pl.DataFrame(rows, infer_schema_length=None)
    if table.select("chain.id", "residue.index", "atom.name").is_duplicated().any():
        raise ValueError("duplicate physical atom keys")
    return table, atoms


def context_exposure(structure):
    """Bio.PDB reference exposure; see module documentation for physical states."""
    from Bio.PDB.Atom import Atom
    from Bio.PDB.Chain import Chain
    from Bio.PDB.Model import Model
    from Bio.PDB.Residue import Residue
    from Bio.PDB.SASA import ShrakeRupley

    table, atoms = _exposure_atoms(structure)
    rows = table.to_dicts()
    sampler = ShrakeRupley(n_points=100, probe_radius=1.4)
    if any(a.element not in sampler.radii_dict for a in atoms):
        raise ValueError("unknown atomic radius; no default radius substitution")

    def areas(mask):
        # One call per physical state into the established SASA implementation.
        model = Model(0)
        indices = np.flatnonzero(mask)
        refs = []
        for i in indices:
            row, atom = rows[i], atoms[i]
            cid, rid = row["chain.id"], (" ", row["residue.index"], " ")
            if cid not in model:
                model.add(Chain(cid))
            if rid not in model[cid]:
                model[cid].add(Residue(rid, row["resname"], " "))
            bioatom = Atom(atom.name, np.asarray(atom.coord, dtype=float), 0., 1., " ",
                           atom.name, int(i), element=atom.element)
            model[cid][rid].add(bioatom)
            refs.append(bioatom)
        sampler.compute(model, level="A")
        result = np.full(len(atoms), np.nan)
        result[indices] = [a.sasa for a in refs]
        return result

    binder = table["role"].to_numpy() == "binder"
    is_peptide = table["chain.type"].to_numpy() == "PEPTIDE"
    bound = areas(np.ones(len(atoms), dtype=bool))
    isolated = np.where(binder, areas(binder), areas(~binder))
    peptide_alone = areas(is_peptide)
    burial = isolated - bound
    shielding = peptide_alone - isolated
    if np.min(burial) < -1e-10 or np.nanmin(shielding) < -1e-10:
        raise RuntimeError("surface exposure increased on adding occluding atoms")
    return table.drop("resname").with_columns(
        pl.Series("isolated_sasa_A2", isolated), pl.Series("bound_sasa_A2", bound),
        pl.Series("partner_burial_A2", burial),
        pl.Series("peptide_alone_sasa_A2", peptide_alone, nan_to_null=True),
        pl.Series("mhc_shielding_A2", shielding, nan_to_null=True))


def _boundary_allocation(result, roles):
    """Assign each lost patch to the nearest sampled opposite-partner boundary.

    One batched native spatial query per endpoint role; only exact distance ties
    require ragged processing. Equal nearest owners split the patch once each.
    """
    from scipy.spatial import cKDTree

    points = np.asarray(result["sample_xyz"], dtype=np.float64)
    owners = np.asarray(result["sample_atom"], dtype=np.int64)
    lost = np.asarray(result["sample_lost"], dtype=bool)
    areas = np.asarray(result["sample_area_A2"], dtype=np.float64)
    source_parts, target_parts, area_parts = [], [], []
    for role in (0, 1):
        query = np.flatnonzero(lost & (roles[owners] == role))
        if not len(query):
            continue
        boundary = np.flatnonzero(roles[owners] != role)
        if not len(boundary):
            raise ValueError("lost patches require a sampled opposite boundary")
        tree = cKDTree(points[boundary])
        distance, nearest = tree.query(points[query], k=2, workers=1)
        source = owners[query]
        target = owners[boundary[nearest[:, 0]]]
        weight = areas[query].copy()
        # The inflated radius retrieves floating-point boundary candidates only;
        # equality of recomputed squared distances, not a tolerance, defines ties.
        tied = np.flatnonzero(distance[:, 0] == distance[:, 1])
        candidates = tree.query_ball_point(points[query[tied]],
            np.nextafter(distance[tied, 0], np.inf), workers=1)
        for row, neighbours in zip(tied, candidates):
            indices = boundary[np.asarray(neighbours, dtype=np.int64)]
            squared = np.sum((points[indices]-points[query[row]])**2, axis=1)
            targets = np.unique(owners[indices[squared == squared.min()]])
            target[row] = targets[0]
            weight[row] /= len(targets)
            if len(targets) > 1:
                source_parts.append(np.full(len(targets)-1, source[row], dtype=np.int64))
                target_parts.append(targets[1:])
                area_parts.append(np.full(len(targets)-1, weight[row]))
        source_parts.append(source); target_parts.append(target); area_parts.append(weight)
    if not source_parts:
        return np.empty(0, dtype=np.int64), np.empty(0, dtype=np.int64), np.empty(0)
    pairs = np.column_stack([np.concatenate(source_parts), np.concatenate(target_parts)])
    unique, inverse = np.unique(pairs, axis=0, return_inverse=True)
    area = np.bincount(inverse, weights=np.concatenate(area_parts))
    if np.any(np.asarray(result["isolated_A2"])[unique[:, 1]] <= 0):
        raise RuntimeError("boundary owner has no sampled isolated exposure")
    return unique[:, 0], unique[:, 1], area


def contact_surface(structure, *, ownership="deepest_sphere"):
    """Return native atom exposures and directed area-conserving atom edges.

    Sphere sampling is the existing native Fibonacci-y convention (100 points).
    Each exposed endpoint patch is assigned to the deepest opposing expanded
    sphere by default, with exact ties split equally. ``accessible_boundary``
    instead uses the nearest sampled isolated opposite boundary, dividing exact
    ties equally among distinct atom owners. Binder/target endpoint areas are
    retained separately and never averaged or folded at residue level.
    """
    from Bio.PDB.SASA import ATOMIC_RADII
    from .. import _geom

    if ownership not in {"deepest_sphere", "accessible_boundary"}:
        raise ValueError("unknown surface ownership method")
    table, atoms = _exposure_atoms(structure)
    if any(a.element not in ATOMIC_RADII for a in atoms):
        raise ValueError("unknown atomic radius; no default radius substitution")
    xyz = np.array([a.coord for a in atoms], dtype=np.float64)
    radii = np.array([ATOMIC_RADII[a.element] for a in atoms], dtype=np.float64)
    role = np.asarray(table["role"] == "target", dtype=np.int32)
    result = _geom.surface_partition(xyz, radii, role, 1.4, 100, samples=ownership=="accessible_boundary")
    isolated, bound = result["isolated_A2"], result["bound_A2"]
    peptide = np.asarray(table["chain.type"] == "PEPTIDE")
    alone = np.full(len(atoms), np.nan)
    alone[peptide] = _geom.shrake_rupley(xyz[peptide], radii[peptide], 1.4, 100)
    atom_table = table.drop("resname").with_columns(
        pl.Series("isolated_sasa_A2", isolated), pl.Series("bound_sasa_A2", bound),
        pl.Series("partner_burial_A2", isolated-bound),
        pl.Series("peptide_alone_sasa_A2", alone, nan_to_null=True),
        pl.Series("mhc_shielding_A2", alone-isolated, nan_to_null=True))
    source, target = np.asarray(result["from"], dtype=np.int64), np.asarray(result["to"], dtype=np.int64)
    area = np.asarray(result["area_A2"], dtype=np.float64)
    if ownership == "accessible_boundary":
        source, target, area = _boundary_allocation(result, role)
    binder_source = role[source] == 0
    edges = pl.DataFrame({"binder_atom": np.where(binder_source,source,target),
        "target_atom": np.where(binder_source,target,source),
        "binder_area_A2": np.where(binder_source,area,0.),
        "target_area_A2": np.where(binder_source,0.,area)}).group_by("binder_atom","target_atom").agg(
            pl.col("binder_area_A2","target_area_A2").sum()).sort("binder_atom","target_atom")
    identity = atom_table.with_row_index("atom_index").with_columns(pl.col("atom_index").cast(pl.Int64))
    for key, side in (("binder_atom","from"),("target_atom","to")):
        mapping = {c: c+"."+side for c in identity.columns if c not in {"pdb.id","mhc_class","atom_index"}}
        meta = identity.drop("pdb.id","mhc_class").rename({"atom_index":key, **mapping})
        edges = edges.join(meta,on=key,how="left",validate="m:1")
    edges = edges.with_columns(pl.lit(structure.pdb_id).alias("pdb.id"),
        pl.lit(atom_table["mhc_class"][0]).alias("mhc_class"))
    # No unallocated partner burial and no silent endpoint reversal.
    for side, field in (("binder","binder_area_A2"),("target","target_area_A2")):
        expected = atom_table.filter(pl.col("role")==side)["partner_burial_A2"].sum()
        if not np.isclose(edges[field].sum(), expected, rtol=1e-12, atol=1e-10):
            raise RuntimeError("surface allocation violates endpoint area conservation")
    return atom_table, edges
