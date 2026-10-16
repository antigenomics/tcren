"""Receptor-independent pMHC topography and peptide solvent exposure.

The supplied coordinates are unchanged; removing a receptor does not undo induced
fit. SASA is Shrake–Rupley (100 sphere points, 1.4 Å probe) in the groove frame.
"""
from __future__ import annotations

import numpy as np

from .surface import KYTE_DOOLITTLE, _groove_frame, surface_map, surface_stats


def pmhc_cells(structure):
    """Exact KD top-face cells in a receptor-independent groove frame.

    Supplied conformations are unchanged. Numeric channels retain the existing
    8-neighbour smoothing; source ownership is the original topmost atom's region.
    Allele calls are derived annotations, never experimental allele assignments.
    """
    import polars as pl
    from ..structure.model import Structure

    peptide = structure.by_type("PEPTIDE")
    mhc = structure.by_type("MHCa", "MHCb")
    if len(peptide) != 1 or not mhc or set(peptide[0].sequence()) - set(KYTE_DOOLITTLE):
        raise ValueError("pMHC cells require one standard-amino-acid peptide and annotated MHC")
    classes = {c.chain_supertype for c in mhc if c.chain_supertype in {"MHCI", "MHCII"}}
    if len(classes) != 1 or (classes == {"MHCII"} and {c.chain_type for c in mhc} != {"MHCa", "MHCb"}):
        raise ValueError("pMHC cells require unambiguous class and both class-II groove chains")
    selected = Structure(structure.pdb_id, structure.by_type("PEPTIDE", "MHCa", "MHCb", "B2M"))
    return surface_map(selected, scale="kd").to_frame().rename({"structure.id": "pdb.id"}).with_columns(
        pl.lit(peptide[0].sequence()).alias("peptide"),
        pl.lit("II" if classes == {"MHCII"} else "I").alias("mhc_class"),
        pl.lit(";".join(c.allele_info or "unknown" for c in mhc)).alias("mhc_call_derived"))


def spatial_features(cells):
    """Reduce authoritative occupied-cell exports with population covariances.

    For common finite peptide cells, e=max(h-r,0), where r is the existing
    90th percentile helix ridge. Cov(e,k) is the spatial excess-height/KD
    co-localization; Cov(h,k) retains height variation below the ridge too.
    Both have units Å times KD index. No fit, direction selection or new cutoff.
    E(e*k)=E(e)*E(k)+Cov(e,k); neither term is a binding or pose-quality score.
    """
    import polars as pl
    from .surface import RIDGE_PERCENTILE

    required = {"pdb.id", "ix", "iy", "source", "h", "phobic", "peptide", "mhc_class"}
    if not cells.height or required - set(cells.columns):
        raise ValueError("nonempty complete pMHC cell export required")
    if cells.select("pdb.id", "ix", "iy").is_duplicated().any():
        raise ValueError("duplicate physical grid cells")
    if cells.select(pl.any_horizontal(pl.col("pdb.id", "ix", "iy", "source", "peptide", "mhc_class").is_null()).any()).item():
        raise ValueError("null pMHC cell keys")
    if not set(cells["mhc_class"]) <= {"I", "II"}:
        raise ValueError("MHC class must be I or II")
    meta = [c for c in ("pdb.id", "peptide", "mhc_class", "mhc_call_derived") if c in cells.columns]
    base = cells.select(meta).unique()
    if base["pdb.id"].is_duplicated().any():
        raise ValueError("inconsistent per-structure metadata")
    helix = cells.filter(pl.col("source").is_in(["mhc_helix_a1", "mhc_helix_a2", "mhc_helix_b1"]) & pl.col("h").is_finite())
    ridge = helix.group_by("pdb.id").agg(pl.col("h").quantile(RIDGE_PERCENTILE/100, interpolation="linear").alias("ridge_height_A"))
    peptide = cells.filter((pl.col("source") == "peptide") & pl.col("h").is_finite())
    counts = peptide.group_by("pdb.id").agg(pl.len().alias("n_peptide_cells"))
    q = peptide.filter(pl.col("phobic").is_finite()).join(ridge, on="pdb.id", how="left")
    q = q.with_columns((pl.col("h")-pl.col("ridge_height_A")).clip(lower_bound=0).alias("excess"))
    q = q.with_columns(*[(pl.col(c)-pl.col(c).mean().over("pdb.id")).alias(c+"_centred") for c in ("h", "phobic", "excess")])
    values = q.group_by("pdb.id").agg(
        pl.len().alias("n_joint_cells"), pl.col("h").std(ddof=0).alias("relief"),
        pl.col("phobic").mean().alias("phobic_mean"),
        pl.col("excess").mean().alias("ridge_excess_mean_A"),
        (pl.col("h_centred")*pl.col("phobic_centred")).mean().alias("height_hydropathy_cov"),
        (pl.col("excess_centred")*pl.col("phobic_centred")).mean().alias("ridge_excess_hydropathy_cov"),
        (pl.col("excess")*pl.col("phobic")).mean().alias("ridge_excess_hydropathy_mean"),
        (pl.col("excess") > 0).sum().alias("n_above_ridge_cells"),
        (pl.col("excess") > 0).mean().alias("frac_above_ridge"),
        pl.col("phobic").filter(pl.col("excess") > 0).mean().alias("above_ridge_hydropathy"))
    # A missing helix ridge is undefined, rather than a zero above-ridge count.
    out = base.join(ridge, on="pdb.id", how="left").join(counts, on="pdb.id", how="left").join(values, on="pdb.id", how="left")
    return out.with_columns(pl.when(pl.col("ridge_height_A").is_finite()).then(pl.col("n_above_ridge_cells")).otherwise(None).alias("n_above_ridge_cells")).sort("pdb.id")


def pmhc_features(structure) -> dict:
    """Return pMHC-only surface scalars and peptide/side-chain SASA in Å².

    Receptor coordinates enter neither the surface normal nor the occluding atoms.
    ``sasa_hydropathy`` weights Kyte–Doolittle values by exposed side-chain area;
    glycine has zero side-chain area. No class-II core/register is guessed.
    All MHC chains retain their physical identities. Nonstandard peptides fail.
    """
    from Bio.PDB.Atom import Atom
    from Bio.PDB.Chain import Chain
    from Bio.PDB.Model import Model
    from Bio.PDB.Residue import Residue
    from Bio.PDB.SASA import ShrakeRupley
    from scipy.spatial import cKDTree

    from ..structure.model import Structure

    peptide = structure.by_type("PEPTIDE")
    mhc = structure.by_type("MHCa", "MHCb")
    if len(peptide) != 1 or not mhc or set(peptide[0].sequence()) - set(KYTE_DOOLITTLE):
        raise ValueError("pMHC features require one standard-amino-acid peptide and annotated MHC")
    classes = {c.chain_supertype for c in mhc if c.chain_supertype in {"MHCI", "MHCII"}}
    if len(classes) != 1:
        raise ValueError("pMHC requires one unambiguous annotated MHC class")
    if classes == {"MHCII"} and {c.chain_type for c in mhc} != {"MHCa", "MHCb"}:
        raise ValueError("class II requires both MHC groove chains")
    selected = Structure(structure.pdb_id, structure.by_type("PEPTIDE", "MHCa", "MHCb", "B2M"))
    origin, basis = _groove_frame(selected)
    values = surface_stats(surface_map(selected, scale="kd"))
    values.update(peptide=peptide[0].sequence(), peptide_length=len(peptide[0].residues),
                  mhc_class="II" if classes == {"MHCII"} else "I")
    sr = ShrakeRupley(n_points=100, probe_radius=1.4)
    # Atoms farther than the sum of expanded radii cannot occlude peptide SASA.
    # Filter before invoking the existing native-neighbour SASA implementation.
    peptide_xyz = np.array([a.coord for r in peptide[0].residues for a in r.atoms
                            if a.element not in {"H", "D"}])
    if not len(peptide_xyz):
        raise ValueError("peptide has no atoms")
    tree = cKDTree(peptide_xyz)
    atoms = [a for c in selected.chains for r in c.residues for a in r.atoms
             if a.element not in {"H", "D"}]
    radius = max(sr.radii_dict[a.element] for a in atoms) + sr.probe_radius
    model = Model(0)
    neighbours = tree.query(np.array([a.coord for a in atoms]))[0] <= 2 * radius
    nearby = {id(a) for a, keep in zip(atoms, neighbours) if keep}
    for chain in selected.chains:
        biochain = Chain(chain.chain_id)
        model.add(biochain)
        for residue in chain.residues:
            biores = Residue((" ", residue.seq_index, " "), residue.resname, " ")
            biochain.add(biores)
            for atom in residue.atoms:
                if atom.element in {"H", "D"}:
                    continue
                if chain.chain_type != "PEPTIDE" and id(atom) not in nearby:
                    continue
                if atom.name in biores:  # Match the writer's single-conformer atom convention.
                    continue
                biores.add(Atom(atom.name, (atom.coord - origin) @ basis.T, 0., 1., " ",
                                atom.name, len(biores), element=atom.element))
    sr.compute(model, level="A")
    backbone = {"N", "CA", "C", "O", "OXT", "OT1", "OT2", "OC1", "OC2"}
    total, side, weighted = 0., 0., 0.
    for residue in peptide[0].residues:
        atoms = model[peptide[0].chain_id][(" ", residue.seq_index, " ")]
        area = sum(a.sasa for a in atoms if a.name not in backbone)
        total += sum(a.sasa for a in atoms)
        side += area
        weighted += area * KYTE_DOOLITTLE[residue.aa]
    values.update(peptide_sasa_A2=float(total), sidechain_sasa_A2=float(side),
                  sasa_hydropathy=float(weighted / side) if side else float("nan"))
    return values
