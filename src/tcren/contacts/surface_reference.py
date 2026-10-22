"""Directed chemistry against exposed-surface opportunity, without assay fitting.

Scores are binder-area-weighted composite log odds (A^2 times table scale),
not binding free energies or independent-contact log likelihoods.
"""
from __future__ import annotations

import numpy as np
import polars as pl
from scipy.special import logsumexp

CDRS = ["CDR1", "CDR2", "CDR3"]
PAIR = ["residue.aa.from", "residue.aa.to"]


def surface_reference_score(atoms, edges, potential):
    """One structure: exposed noncontacts define Q; binder patch area defines U.

    Q is the product of role-specific exposed AA area fractions, conditional on
    finite ordered matrix support. C=U+A*logZ, Z=sum Q*exp(-V). Missing cells and
    zero-opportunity observed area are reported and omitted, without imputation.
    Empty supported area gives empty sums of zero; its mean is undefined.
    """
    if atoms["pdb.id"].n_unique() != 1 or (edges.height and set(edges["pdb.id"]) != set(atoms["pdb.id"])):
        raise ValueError("surface reference requires one matching structure")
    if atoms.select("chain.id", "residue.index", "atom.name").is_duplicated().any():
        raise ValueError("unique physical atoms required")
    if edges.select("binder_atom", "target_atom").is_duplicated().any():
        raise ValueError("unique directed surface edges required")
    if (atoms.filter((~pl.col("isolated_sasa_A2").is_finite() | (pl.col("isolated_sasa_A2") < 0)).fill_null(True)).height
            or edges.filter((~pl.col("binder_area_A2").is_finite() | (pl.col("binder_area_A2") < 0)).fill_null(True)).height):
        raise ValueError("finite nonnegative exposed and allocated areas required")
    if edges.filter((pl.col("role.from") != "binder") | (pl.col("role.to") != "target")).height:
        raise ValueError("surface pairs must preserve binder-to-target roles")
    binder = atoms.filter((pl.col("role") == "binder") & pl.col("region.type").is_in(CDRS))
    target = atoms.filter((pl.col("role") == "target") & (pl.col("chain.type") == "PEPTIDE"))
    selected = edges.filter(pl.col("region.type.from").is_in(CDRS) & (pl.col("chain.type.to") == "PEPTIDE") & (pl.col("binder_area_A2") > 0))
    matrix = potential.matrix.filter(pl.col("value").is_finite())
    if matrix.select(PAIR).is_duplicated().any():
        raise ValueError("unique ordered potential cells required")
    def composition(table, side):
        return table.group_by("residue.aa").agg(pl.col("isolated_sasa_A2").sum().alias("area_"+side)).rename({"residue.aa":"residue.aa."+side})
    q = composition(binder,"from").join(composition(target,"to"),how="cross").with_columns(
        (pl.col("area_from")*pl.col("area_to")).alias("mass"))
    mass = q["mass"].sum()
    q = q.join(matrix,on=PAIR,how="inner",validate="1:1").filter(pl.col("mass") > 0)
    supported_mass = q["mass"].sum()
    logz = None
    if supported_mass > 0:
        q = q.with_columns((pl.col("mass")/supported_mass).alias("q"))
        logz = float(logsumexp(np.log(q["q"].to_numpy())-q["value"].to_numpy()))
    observed = selected.join(matrix,on=PAIR,how="left",validate="m:1").join(
        q.select(*PAIR,"mass"),on=PAIR,how="left",validate="m:1")
    supported = observed.filter(pl.col("value").is_finite() & (pl.col("mass") > 0))
    area = supported["binder_area_A2"].sum()
    raw = supported.select((pl.col("binder_area_A2")*pl.col("value")).sum()).item()
    corrected = raw + area*logz if logz is not None else 0.
    return dict(surface_raw_sum=float(raw), surface_context_sum=float(corrected),
        surface_context_mean=float(corrected/area) if area else None,
        surface_logZ=logz, surface_supported_area_A2=float(area),
        surface_cdr_peptide_area_A2=float(selected["binder_area_A2"].sum()),
        surface_missing_potential_area_A2=float(observed.filter(pl.col("value").is_null())["binder_area_A2"].sum()),
        surface_out_of_opportunity_area_A2=float(observed.filter(pl.col("value").is_finite() & pl.col("mass").is_null())["binder_area_A2"].sum()),
        surface_opportunity_support_fraction=float(supported_mass/mass) if mass else None,
        surface_binder_exposed_A2=float(binder["isolated_sasa_A2"].sum()),
        surface_target_exposed_A2=float(target["isolated_sasa_A2"].sum()),
        surface_reference_factor_mean=float(np.exp(np.log(q["q"].to_numpy())-q["value"].to_numpy()-logz).sum()) if logz is not None else None,
        surface_n_atom_edges=selected.height, surface_n_supported_atom_edges=supported.height,
        surface_status="supported" if area else "no_supported_area")
