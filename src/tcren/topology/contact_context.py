"""Local contact graph support and pMHC topography; no label-dependent filters."""
import numpy as np
import polars as pl

from ..structure.model import Structure
from .surface import KYTE_DOOLITTLE, RIDGE_PERCENTILE, SOURCE_NAMES, _groove_frame, surface_map


def summarize_contacts(contacts: pl.DataFrame, accounting: pl.DataFrame) -> pl.DataFrame:
    """Canonical pair-count reduction; identical compositions remain tied.

    Fixed masks: all contacts, edges in four-cycles, peptide surface above the
    existing ridge, and their intersection. No fitting or label inputs. Complete
    accounting retains zero-contact models. Missing topography stays undefined.
    """
    keys = ["pdb.id", "chain.id.from", "residue.index.from", "chain.id.to", "residue.index.to"]
    required = set(keys + ["four_cycles", "surface_source", "surface_excess_A",
                           "pair_energy", "residue.aa.from", "residue.aa.to"])
    if required - set(contacts.columns) or accounting["pdb.id"].is_duplicated().any():
        raise ValueError("complete contact columns and unique accounting IDs required")
    if not accounting.height or set(accounting["status"]) != {"ok"}:
        raise ValueError("complete successful accounting required")
    if contacts.is_empty():
        # A header-only CSV has string columns even for these numeric quantities.
        contacts = contacts.with_columns(
            pl.col("four_cycles", "pair_energy", "surface_excess_A").cast(pl.Float64))
    if contacts.select(keys).is_duplicated().any() or contacts.filter(pl.col("four_cycles") < 0).height:
        raise ValueError("unique residue-pair edges with nonnegative cycle counts required")
    ids = sorted(accounting["pdb.id"].to_list())
    index = {s: i for i, s in enumerate(ids)}
    if set(contacts["pdb.id"]) - set(ids):
        raise ValueError("unaccounted contact IDs")
    counts = contacts.group_by("pdb.id").len()
    ledger = accounting.select("pdb.id", pl.col("n_rows").cast(pl.UInt32).alias("expected"))
    if ledger.join(counts, on="pdb.id", how="left").filter(
            pl.col("expected") != pl.col("len").fill_null(0)).height:
        raise ValueError("contact rows disagree with accounting")
    pairs = contacts.select("residue.aa.from", "residue.aa.to", "pair_energy").unique().sort(
        "residue.aa.from", "residue.aa.to")
    if pairs.select("residue.aa.from", "residue.aa.to").is_duplicated().any():
        raise ValueError("inconsistent potential values for a residue pair")
    codes = {(a, b): i for i, (a, b) in enumerate(pairs.select(
        "residue.aa.from", "residue.aa.to").iter_rows())}
    u = np.array([index[s] for s in contacts["pdb.id"]], dtype=int)
    v = np.array([codes[k] for k in contacts.select("residue.aa.from", "residue.aa.to").iter_rows()], dtype=int)
    if np.isinf(pairs["pair_energy"].to_numpy()).any():
        raise ValueError("infinite potential entries")
    energy = np.nan_to_num(pairs["pair_energy"].to_numpy(), nan=0)
    known = np.isfinite(contacts["surface_excess_A"].to_numpy())
    cycle = contacts["four_cycles"].to_numpy() > 0
    bulge = (contacts["surface_source"].to_numpy() == "peptide") & (
        contacts["surface_excess_A"].to_numpy() > 0)
    result = {"pdb.id": ids}
    missing = np.bincount(u[~known], minlength=len(ids))
    for name, mask in (("all", np.ones(len(u), dtype=bool)), ("cycle", cycle),
                       ("bulge", bulge), ("joint", cycle & bulge)):
        matrix = np.bincount(u[mask]*len(pairs)+v[mask], minlength=len(ids)*len(pairs)).reshape(len(ids), len(pairs))
        values = matrix @ energy
        number = matrix.sum(1).astype(float)
        if name in {"bulge", "joint"}:
            values[missing > 0] = np.nan; number[missing > 0] = np.nan
        result["n_"+name] = number
        result["context_energy" if name == "all" else name+"_energy"] = values
    result["n_missing_topography"] = missing
    result["n_missing_energy"] = np.bincount(u[~np.isfinite(contacts["pair_energy"].to_numpy())], minlength=len(ids))
    return pl.DataFrame(result)


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
