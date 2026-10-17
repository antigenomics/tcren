"""Local contact graph support and pMHC topography; no label-dependent filters."""
import numpy as np
import polars as pl

from ..structure.model import Structure
from .surface import KYTE_DOOLITTLE, RIDGE_PERCENTILE, SOURCE_NAMES, _groove_frame, surface_map


def contact_context(structure, contacts, potential):
    """Annotate TCR:peptide residue-pair contacts in their receptor-free pMHC frame.

    Four-cycles measure local bipartite redundancy, not physical interaction truth.
    Surface cells use the closest peptide atom's x/y projection; ownership is by
    peptide versus MHC, not individual peptide residue. Heights are in angstroms,
    hydropathy is the KD index and energies retain the supplied potential's units.
    """
    selected = Structure(structure.pdb_id, structure.by_type("PEPTIDE", "MHCa", "MHCb", "B2M"))
    groove = selected.by_type("MHCa", "MHCb")
    classes = {c.chain_supertype for c in groove if c.chain_supertype in {"MHCI", "MHCII"}}
    if len(classes) != 1 or (classes == {"MHCII"} and {c.chain_type for c in groove} != {"MHCa", "MHCb"}):
        raise ValueError("contact context requires an unambiguous complete MHC groove")
    origin, basis = _groove_frame(selected)
    smap = surface_map(selected, scale="kd")
    cells = smap.to_frame()
    helix = cells.filter(pl.col("source").is_in(["mhc_helix_a1", "mhc_helix_a2", "mhc_helix_b1"]))
    ridge = float(np.quantile(helix["h"].to_numpy(), RIDGE_PERCENTILE/100)) if helix.height else float("nan")
    left = list(zip(contacts["chain.id.from"], contacts["residue.index.from"]))
    right = list(zip(contacts["chain.id.to"], contacts["residue.index.to"]))
    li = {k: i for i, k in enumerate(sorted(set(left)))}
    ri = {k: i for i, k in enumerate(sorted(set(right)))}
    u = np.array([li[k] for k in left], dtype=int)
    v = np.array([ri[k] for k in right], dtype=int)
    B = np.zeros((len(li), len(ri)), dtype=np.int64)
    B[u, v] = 1
    dl, dr = B.sum(1), B.sum(0)
    cycles = (B @ B.T @ B)[u, v] - dl[u] - dr[v] + 1
    atoms = {(c.chain_id, r.seq_index, a.name): a.coord
             for c in selected.by_type("PEPTIDE") for r in c.residues for a in r.atoms}
    xyz = np.array([atoms[(k[0], k[1], a)] for k, a in zip(right, contacts["atom.to"])])
    xyz = xyz.reshape((-1, 3))
    local = (xyz - origin) @ basis.T
    ny, nx = smap.grid
    x0, x1, y0, y1 = smap.extent
    ix = np.floor((local[:, 0]-x0)*nx/(x1-x0)).astype(int)
    iy = np.floor((local[:, 1]-y0)*ny/(y1-y0)).astype(int)
    inside = (ix >= 0) & (ix < nx) & (iy >= 0) & (iy < ny)
    h = np.full(len(u), np.nan); kd = h.copy()
    source = np.full(len(u), "outside", dtype=object)
    h[inside] = smap.channels["h"][iy[inside], ix[inside]]
    kd[inside] = smap.channels["phobic"][iy[inside], ix[inside]]
    source[inside] = [SOURCE_NAMES[i] for i in smap.source[iy[inside], ix[inside]]]
    matrix, alphabet = potential.as_matrix()
    energy = np.array([matrix[alphabet[a], alphabet[b]] if a in alphabet and b in alphabet
                       else np.nan for a, b in zip(contacts["residue.aa.from"], contacts["residue.aa.to"])])
    return contacts.with_columns(
        pl.Series("tcr_degree", dl[u]), pl.Series("peptide_degree", dr[v]),
        pl.Series("four_cycles", cycles), pl.Series("surface_source", source.tolist(), dtype=pl.String),
        pl.Series("surface_h_A", h), pl.Series("surface_hydropathy", kd),
        pl.Series("atom_height_A", local[:, 2]), pl.lit(ridge).alias("ridge_height_A"),
        pl.Series("surface_excess_A", h - ridge),
        pl.Series("peptide_hydropathy", [KYTE_DOOLITTLE[a] for a in contacts["residue.aa.to"]]),
        pl.Series("pair_energy", energy))
