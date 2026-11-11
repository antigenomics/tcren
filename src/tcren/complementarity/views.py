"""Canonical-frame topography and faceted distances: typed tables own all values."""
from __future__ import annotations

from collections import defaultdict

import numpy as np
import polars as pl
from scipy.spatial.distance import cdist

from ..topology.surface import (
    KYTE_DOOLITTLE, _BONDI, _DEFAULT_RADIUS, _V_REGIONS, _height_candidates,
)

INT, F, S, B = pl.Int64, pl.Float64, pl.String, pl.Boolean
LOOPS = ("alpha:CDR1", "alpha:CDR2", "alpha:CDR3", "beta:CDR3", "beta:CDR2", "beta:CDR1")
PARTNERS = ("peptide", "helix_1", "helix_2")
COMMON = dict(structure_id=S, reference=S, task=S)
SCHEMAS = {
    "topography_residues": {**COMMON, "residue_id": INT, "chain": S, "role": S, "region": S,
                            "seq_index": INT, "number": INT, "insertion": S, "aa": S,
                            "ca_atom": INT, "ca_x": F, "ca_y": F, "ca_z": F,
                            "atom_z_min": F, "atom_z_max": F, "face_apex_z": F,
                            "hydropathy": F, "heavy_atoms": INT},
    "topography_cells": {**COMMON, "ix": INT, "iy": INT, "x": F, "y": F, "occupied": B,
                         "surface_z": F, "hydropathy": F, "owner_atom": INT, "owner_residue": INT,
                         "chain": S, "role": S, "region": S, "aa": S},
    "contact_heatmap": {**COMMON, "loop_facet": S, "loop_order": INT, "partner_facet": S,
                        "partner_order": INT, "tcr_index": INT, "partner_index": INT,
                        "tcr_residue_id": INT, "tcr_chain": S, "tcr_number": INT,
                        "tcr_insertion": S, "tcr_aa": S, "tcr_seq_index": INT,
                        "partner_residue_id": INT, "partner_chain": S, "partner_number": INT,
                        "partner_insertion": S, "partner_aa": S, "partner_region": S,
                        "partner_seq_index": INT, "ca_distance": F, "has_both_ca": B,
                        "atom_contact_count": INT, "contacting": B},
}


def validate_grid(grid, extent, probe):
    if len(grid) != 2 or any(isinstance(n, bool) or not isinstance(n, int) or not 1 <= n <= 256 for n in grid):
        raise ValueError("grid must have two integer dimensions in [1,256]")
    if len(extent) != 4 or not np.isfinite(extent).all() or extent[0] >= extent[1] or extent[2] >= extent[3]:
        raise ValueError("extent must be finite x0<x1,y0<y1")
    if not np.isfinite((extent[1]-extent[0], extent[3]-extent[2])).all():
        raise ValueError("extent widths must be finite")
    if not np.isfinite(probe) or probe < 0:
        raise ValueError("probe must be finite and nonnegative")


def _frame(records, name):
    return pl.DataFrame(records, schema=SCHEMAS[name])


def topography_tables(scene, *, grid=(64, 64), extent=(-40., 40., -40., 40.), probe=1.4):
    """Unsmoothed sphere-envelope surface, in the owning scene's exact frame.

    Epitope ranking maps the peptide upper face; TCR ranking maps the variable
    domains' lower face. The existing surface ray caster and KD constants are
    reused; the legacy surface command's independent frame is never called.
    Every grid cell is retained, with null chemistry/height for empty cells.
    """
    validate_grid(grid, extent, probe)
    task = scene.get("task")
    if task not in {"tcr", "epitope"}:
        raise ValueError("scene task must be tcr or epitope")
    common = dict(structure_id=scene["id"], reference=scene["frame"]["reference"], task=task)
    selected = [r for r in scene["residues"] if
                (task == "epitope" and r["role"] == "peptide") or
                (task == "tcr" and r["role"] in {"alpha", "beta"} and r["region"] in _V_REGIONS)]
    if not selected:
        raise ValueError("no residues on the requested varying partner")
    selected_ids = {r["id"] for r in selected}
    atoms = [a for a in scene["atoms"] if a["residue"] in selected_ids]
    if not atoms:
        raise ValueError("no heavy atoms on the requested varying partner")
    by_residue = defaultdict(list)
    for atom in atoms:
        by_residue[atom["residue"]].append(atom)
    residues = []
    for r in selected:
        members = by_residue[r["id"]]
        z = [a["xyz"][2] for a in members]
        ca = scene["atoms"][r["ca_atom"]]["xyz"] if r["ca_atom"] is not None else [None]*3
        residues.append({**common, "residue_id": r["id"], **{k: r[k] for k in
                         ("chain", "role", "region", "seq_index", "number", "insertion", "aa", "ca_atom")},
                         **dict(zip(("ca_x", "ca_y", "ca_z"), ca)),
                         "atom_z_min": min(z) if z else None, "atom_z_max": max(z) if z else None,
                         "face_apex_z": (min(z) if task == "tcr" else max(z)) if z else None,
                         "hydropathy": KYTE_DOOLITTLE.get(r["aa"]), "heavy_atoms": len(z)})
    xyz = np.array([a["xyz"] for a in atoms], dtype=float)
    if not np.isfinite(xyz).all():
        raise ValueError("topography requires finite coordinates")
    local = xyz.copy()
    direction = -1 if task == "tcr" else 1
    local[:, 2] *= direction
    radii = np.array([_BONDI.get(a["element"], _DEFAULT_RADIUS) for a in atoms]) + probe
    # The shared sampler enumerates lateral offsets for every atom. Grid size
    # alone does not bound that work when the probe is huge or the extent tiny.
    ny, nx = grid
    dx, dy = (extent[1]-extent[0])/nx, (extent[3]-extent[2])/ny
    radius = float(radii.max())
    offsets = (radius/dx if dx else float("inf"), radius/dy if dy else float("inf"))
    if (not np.isfinite(offsets).all() or max(offsets) > 20_000_000
            or (2*np.ceil(offsets[0])+1)*(2*np.ceil(offsets[1])+1)*len(atoms) > 20_000_000):
        raise ValueError("topography exceeds 20,000,000 atom/offset evaluations; "
                         "use a smaller grid/probe or a wider extent")
    cells, owners, heights = _height_candidates(local, radii, grid, extent)
    order = np.lexsort((heights, cells))
    cells, owners, heights = cells[order], owners[order], heights[order]
    last = np.r_[cells[1:] != cells[:-1], True] if len(cells) else np.zeros(0, bool)
    exposed = {int(cell): (int(owner), float(height)*direction)
               for cell, owner, height in zip(cells[last], owners[last], heights[last])}
    rows = []
    for iy in range(ny):
        for ix in range(nx):
            exposed_atom = exposed.get(iy*nx+ix)
            atom = atoms[exposed_atom[0]] if exposed_atom else None
            residue = scene["residues"][atom["residue"]] if atom else None
            rows.append({**common, "ix": ix, "iy": iy, "x": extent[0]+(ix+.5)*dx,
                         "y": extent[2]+(iy+.5)*dy, "occupied": atom is not None,
                         "surface_z": exposed_atom[1] if exposed_atom else None,
                         "hydropathy": KYTE_DOOLITTLE.get(residue["aa"]) if residue else None,
                         "owner_atom": atom["id"] if atom else None,
                         "owner_residue": residue["id"] if residue else None,
                         **{k: residue[k] if residue else None for k in ("chain", "role", "region", "aa")}})
    return _frame(residues, "topography_residues"), _frame(rows, "topography_cells")


def heatmap_table(scene):
    """Dense C-alpha distances and existing contact multiplicities, in 18 facets.

    Beta facet order is reversed; residues within each loop retain their natural
    sequence order. Missing C-alpha distances remain null, never zero. This does
    not detect contacts: witnesses are joined from the existing scene graph.
    """
    common = dict(structure_id=scene["id"], reference=scene["frame"]["reference"], task=scene["task"])
    residues = scene["residues"]
    helices = {r["region"] for r in residues if r["role"] == "mhc"}
    if {"HELIX_A2", "HELIX_B1"} <= helices:
        raise ValueError("ambiguous second groove helix: both class-I A2 and class-II B1 present")
    second = "HELIX_B1" if "HELIX_B1" in helices else "HELIX_A2"
    partners = [[r for r in residues if r["role"] == "peptide"],
                [r for r in residues if r["role"] == "mhc" and r["region"] == "HELIX_A1"],
                [r for r in residues if r["role"] == "mhc" and r["region"] == second]]
    counts = {tuple(sorted((e["source"], e["target"]))): e["contact_count"] for e in scene["residue_contacts"]}
    rows, facets = [], []
    for column, loop in enumerate(LOOPS):
        role, region = loop.split(":")
        tcr = sorted((r for r in residues if r["role"] == role and r["region"] == region),
                     key=lambda r: (r["chain"], r["seq_index"]))
        for row, partner in enumerate(partners):
            partner = sorted(partner, key=lambda r: (r["chain"], r["seq_index"]))
            facets.append(dict(loop_facet=loop, partner_facet=PARTNERS[row],
                               tcr_residues=len(tcr), partner_residues=len(partner),
                               partner_regions=sorted({r["region"] for r in partner})))
            if len(rows) + len(tcr)*len(partner) > 200_000:
                raise ValueError("heatmap exceeds 200,000 residue pairs for one structure")
            # Only resolved CA coordinates enter cdist; missing entries are scattered as null.
            ti = [i for i, r in enumerate(tcr) if r["ca_atom"] is not None]
            pi = [i for i, r in enumerate(partner) if r["ca_atom"] is not None]
            distances = np.full((len(tcr), len(partner)), np.nan)
            if ti and pi:
                coordinates = [scene["atoms"][r["ca_atom"]]["xyz"]
                               for r in [*(tcr[i] for i in ti), *(partner[i] for i in pi)]]
                if not np.isfinite(coordinates).all():
                    raise ValueError("heatmaps require finite C-alpha coordinates")
                distances[np.ix_(ti, pi)] = cdist(
                    [scene["atoms"][tcr[i]["ca_atom"]]["xyz"] for i in ti],
                    [scene["atoms"][partner[i]["ca_atom"]]["xyz"] for i in pi])
            for i, t in enumerate(tcr):
                for j, p in enumerate(partner):
                    count = counts.get(tuple(sorted((t["id"], p["id"]))), 0)
                    row_values = {**common, "loop_facet": loop, "loop_order": column,
                                  "partner_facet": PARTNERS[row], "partner_order": row,
                                  "tcr_index": i, "partner_index": j,
                                  "ca_distance": float(distances[i, j]) if np.isfinite(distances[i, j]) else None,
                                  "has_both_ca": t["ca_atom"] is not None and p["ca_atom"] is not None,
                                  "atom_contact_count": count, "contacting": count > 0}
                    for prefix, residue in [("tcr", t), ("partner", p)]:
                        row_values.update({prefix+"_"+key: residue[key] for key in
                                           ("chain", "number", "insertion", "aa", "seq_index")})
                        row_values[prefix+"_residue_id"] = residue["id"]
                    row_values["partner_region"] = p["region"]
                    rows.append(row_values)
    return _frame(rows, "contact_heatmap"), facets


def view_tables(scene, **options):
    residues, cells = topography_tables(scene, **options)
    heatmap, facets = heatmap_table(scene)
    return {"topography_residues": residues, "topography_cells": cells, "contact_heatmap": heatmap}, facets
