"""Observed atom contacts for sparse reference statistics; no fitted scoring."""
from dataclasses import replace

import numpy as np
import polars as pl

from .geometry import BACKBONE_ATOMS, all_atom_contacts
from ..contact_types import classify_contacts, stacked_pairs

INTERFACES = ("tcr_peptide", "tcr_mhc", "peptide_mhc")
CONFORMER_POLICY = "per residue: lexicographically first nonblank altloc, plus blank atoms; selected label wins name ties"
_CDRS = {"CDR1", "CDR2", "CDR3"}
_HELICES = {"HELIX_A1", "HELIX_A2", "HELIX_B1"}


def _validate_complex(structure):
    """Shared alpha-beta receptor, standard peptide and MHC scope validation."""
    types = {c.chain_type for c in structure.chains}
    if not {"TRA", "TRB", "PEPTIDE", "MHCa"} <= types or types & {"TRG", "TRD", "IGH", "IGK", "IGL"}:
        raise ValueError("contact observations require an alpha-beta TCR:pMHC complex")
    peptides = structure.by_type("PEPTIDE")
    if len(peptides) != 1 or set(peptides[0].sequence()) - set("ACDEFGHIKLMNPQRSTVWY"):
        raise ValueError("exactly one standard-amino-acid peptide required")
    classes = {c.chain_supertype for c in structure.chains if c.chain_type in {"MHCa", "MHCb"}}
    if classes not in ({"MHCI"}, {"MHCII"}) or (classes == {"MHCII"} and "MHCb" not in types):
        raise ValueError("consistent annotated MHC-I or two-chain MHC-II groove required")
    return classes


_INVENTORY_SCHEMA = {
    "pdb.id": pl.String, "chain.id": pl.String, "chain.type": pl.String,
    "residue.index": pl.Int64, "residue.index.pdb": pl.Int64,
    "insertion_code": pl.String, "residue.aa": pl.String, "residue.resname": pl.String,
    "region.type": pl.String, "peptide.pos": pl.Int64, "atom.name": pl.String,
    "atom.source_index": pl.Int64, "element": pl.String,
    "altloc": pl.String, "altloc.selected": pl.String,
    "x": pl.Float64, "y": pl.Float64, "z": pl.Float64,
    "role": pl.String, "contact_role": pl.String, "contact_eligible": pl.Boolean,
    "n_altloc_omitted": pl.Int64,
}


def atom_inventory(structure) -> pl.DataFrame:
    """Selected observed heavy atoms, including uncontacted full-component context.

    Coordinates are in Angstroms. ``peptide.pos`` is the one-based ordinal in the
    complete observed peptide, independent of author numbering and seq_index.
    ``atom.source_index`` enumerates heavy atoms across all input chains before
    selection, exactly as for contact observations. ``n_altloc_omitted`` counts
    heavy atoms removed by residue conformer/name selection only, repeated per
    retained atom; no Gly or region filtering is applied to inventory membership.
    ``contact_role`` uses the existing CDR/helix role; ``contact_eligible`` adds
    the existing Gly backbone restriction. No geometry or missing atoms are built.
    """
    _validate_complex(structure)
    chain_keys = set()
    for chain in structure.chains:
        if chain.chain_id in chain_keys:
            raise ValueError("duplicate chain identity in atom inventory")
        chain_keys.add(chain.chain_id)
        sequence_keys, author_keys = set(), set()
        for residue in chain.residues:
            author = (residue.pdb_index, residue.insertion_code)
            if residue.seq_index in sequence_keys or author in author_keys:
                raise ValueError("duplicate residue identity in atom inventory")
            sequence_keys.add(residue.seq_index)
            author_keys.add(author)
    roles = {"TRA": "tcr", "TRB": "tcr", "PEPTIDE": "peptide",
             "MHCa": "mhc", "MHCb": "mhc", "B2M": "b2m"}
    records = []
    for chain, residue, position, region, contact_role, selected, heavy, kept in _residue_atoms(structure):
        if chain.chain_type not in roles:
            continue
        for source_index, atom in kept:
            try:
                xyz = np.asarray(atom.coord, dtype=float)
            except (TypeError, ValueError) as exc:
                raise ValueError("atom inventory requires finite three-dimensional atom coordinates") from exc
            if xyz.shape != (3,) or not np.isfinite(xyz).all():
                raise ValueError("atom inventory requires finite three-dimensional atom coordinates")
            records.append({
                "pdb.id": structure.pdb_id, "chain.id": chain.chain_id,
                "chain.type": chain.chain_type, "residue.index": residue.seq_index,
                "residue.index.pdb": residue.pdb_index, "insertion_code": residue.insertion_code,
                "residue.aa": residue.aa, "residue.resname": residue.resname,
                "region.type": region, "peptide.pos": position if chain.chain_type == "PEPTIDE" else None,
                "atom.name": atom.name, "atom.source_index": source_index, "element": atom.element,
                "altloc": atom.altloc, "altloc.selected": selected,
                **{axis: float(xyz[i]) for i, axis in enumerate("xyz")},
                "role": roles[chain.chain_type], "contact_role": contact_role,
                "contact_eligible": contact_role is not None and (residue.aa != "G" or atom.name in BACKBONE_ATOMS),
                "n_altloc_omitted": len(heavy) - len(kept),
            })
    return pl.DataFrame(records, schema=_INVENTORY_SCHEMA)


def _residue_atoms(structure):
    """Yield conformer-selected atoms with indices over *all* source heavy atoms."""
    source_index = 0
    for chain in structure.chains:
        regions = {r.seq_index: reg.region_type for reg in chain.regions for r in reg.residues}
        for position, residue in enumerate(chain.residues, 1):
            heavy = []
            for atom in residue.atoms:
                if atom.element.strip().upper() in {"H", "D"}:
                    continue
                heavy.append((source_index, atom))
                source_index += 1
            region = regions.get(residue.seq_index)
            role = ("tcr" if chain.chain_type in {"TRA", "TRB"} and region in _CDRS
                    else "peptide" if chain.chain_type == "PEPTIDE"
                    else "mhc" if chain.chain_type in {"MHCa", "MHCb"} and region in _HELICES
                    else None)
            labels = sorted({a.altloc for _, a in heavy if a.altloc})
            selected = labels[0] if labels else ""
            by_name = {}
            for index, atom in sorted(heavy, key=lambda x: bool(x[1].altloc)):
                if atom.altloc in {"", selected}:
                    by_name[atom.name] = (index, atom)
            yield chain, residue, position, region, role, selected, heavy, list(by_name.values())


def select_conformer(structure):
    """Copy full components with the existing per-residue atom/conformer policy."""
    selected = {c.chain_id: [] for c in structure.chains}
    for chain, residue, _, _, _, _, _, kept in _residue_atoms(structure):
        selected[chain.chain_id].append(replace(residue, atoms=tuple(a for _, a in kept)))
    return replace(structure, chains=[replace(c, residues=selected[c.chain_id])
                                      for c in structure.chains])


def _selected_atoms(structure):
    """Copy only supported component atoms; retain original heavy-atom indices."""
    residues_by_chain = {id(chain): [] for chain in structure.chains}
    records = []
    for chain, residue, _, region, role, selected, heavy, conformer in _residue_atoms(structure):
        kept = [(i, a) for i, a in conformer if role is not None
                and (residue.aa != "G" or a.name in BACKBONE_ATOMS)]
        residues_by_chain[id(chain)].append(replace(residue, atoms=tuple(a for _, a in kept)))
        coords = {a.name: a.coord for _, a in kept}
        for original_index, atom in kept:
            xyz = np.asarray(atom.coord, dtype=float)
            if xyz.shape != (3,) or not np.isfinite(xyz).all():
                raise ValueError("contact observations require finite three-dimensional atom coordinates")
            record = {"atom.index": len(records), "atom.source_index": original_index,
                      "residue.index.pdb": residue.pdb_index, "insertion_code": residue.insertion_code,
                      "chain.type": chain.chain_type, "region.type": region, "role": role,
                      "altloc": atom.altloc, "altloc.selected": selected,
                      "n_atoms_omitted": len(heavy) - len(kept), "element": atom.element,
                      "chem": residue.aa + ":" + atom.name}
            for prefix, value in (("", xyz), ("n.", coords.get("N")),
                                  ("ca.", coords.get("CA")), ("cb.", coords.get("CB"))):
                record.update({prefix + axis: None if value is None else float(value[k])
                               for k, axis in enumerate("xyz")})
            records.append(record)
    chains = [replace(chain, residues=residues_by_chain[id(chain)]) for chain in structure.chains]
    schema = {"atom.index": pl.Int64, "atom.source_index": pl.Int64, "residue.index.pdb": pl.Int64,
              "insertion_code": pl.String, "chain.type": pl.String, "region.type": pl.String,
              "role": pl.String, "altloc": pl.String, "altloc.selected": pl.String,
              "n_atoms_omitted": pl.Int64, "element": pl.String, "chem": pl.String,
              **{prefix + axis: pl.Float64 for prefix in ("", "n.", "ca.", "cb.") for axis in "xyz"}}
    return replace(structure, chains=chains), pl.DataFrame(records, schema=schema)


def _angle(a, b, c):
    u, v = a - b, c - b
    with np.errstate(divide="ignore", invalid="ignore"):
        cosine = (u*v).sum(1) / (np.linalg.norm(u, axis=1)*np.linalg.norm(v, axis=1))
    return np.arccos(np.clip(cosine, -1, 1))


def _dihedral(a, b, c, d):
    b0, axis, b2 = a-b, c-b, d-c
    with np.errstate(divide="ignore", invalid="ignore"):
        axis = axis / np.linalg.norm(axis, axis=1)[:, None]
    v = b0 - (b0*axis).sum(1)[:, None]*axis
    w = b2 - (b2*axis).sum(1)[:, None]*axis
    value = np.arctan2((np.cross(axis, v)*w).sum(1), (v*w).sum(1))
    return np.where((np.linalg.norm(v, axis=1) > 0) & (np.linalg.norm(w, axis=1) > 0), value, np.nan)


def contact_observations(structure, *, cutoff=5.0, interface="all") -> pl.DataFrame:
    """Unique heavy-atom contacts between CDRs, peptide and annotated MHC helices.

    Select one deterministic alternate conformer per residue (see
    ``CONFORMER_POLICY``); do not alter the input or legacy contact calculation.
    Glycine contributes observed backbone atoms only. Biological roles remain
    separate from unordered ``chem.a``/``chem.b`` endpoint identities. Orientation
    angles use actual N/CA/CB atoms; missing or degenerate geometry stays null.
    Chemical flags are existing v2 distance/type heuristics, not angular validation.
    """
    if interface not in ("all", *INTERFACES):
        raise ValueError("invalid contact observation interface")
    classes = _validate_complex(structure)
    for chain in structure.by_type("MHCa", "MHCb"):
        if not any(r.region_type in _HELICES for r in chain.regions):
            raise ValueError("MHC helix region annotation required for contact observations")
    selected, metadata = _selected_atoms(structure)
    table = all_atom_contacts(selected, cutoff=cutoff, atom_pairs=True, atom_indices=True)
    for side in ("from", "to"):
        table = table.join(metadata.rename({c: f"{c}.{side}" for c in metadata.columns}),
                           on=f"atom.index.{side}", how="left", validate="m:1")
    a, b = pl.col("role.from"), pl.col("role.to")
    table = table.with_columns(
        pl.when(((a == "tcr") & (b == "peptide")) | ((b == "tcr") & (a == "peptide")))
        .then(pl.lit("tcr_peptide"))
        .when(((a == "tcr") & (b == "mhc")) | ((b == "tcr") & (a == "mhc")))
        .then(pl.lit("tcr_mhc"))
        .when(((a == "peptide") & (b == "mhc")) | ((b == "peptide") & (a == "mhc")))
        .then(pl.lit("peptide_mhc")).otherwise(None).alias("interface"))
    table = table.filter(pl.col("interface").is_not_null())
    if interface != "all":
        table = table.filter(pl.col("interface") == interface)
    table = classify_contacts(table, scheme="v2", stacking=stacked_pairs(selected))
    table = table.with_columns(
        pl.lit(structure.pdb_id).alias("pdb.id"),
        pl.lit("I" if classes == {"MHCI"} else "II").alias("mhc_class"),
        pl.min_horizontal("chem.from", "chem.to").alias("chem.a"),
        pl.max_horizontal("chem.from", "chem.to").alias("chem.b"),
        pl.when(pl.col("role.from") == "tcr").then(pl.col("region.type.from"))
        .when(pl.col("role.to") == "tcr").then(pl.col("region.type.to"))
        .otherwise(pl.lit("none")).alias("receptor_region"))
    def xyz(prefix, side):
        return table.select([f"{prefix}{axis}.{side}" for axis in "xyz"]).to_numpy()
    ca, cb, n = [{side: xyz(prefix, side) for side in ("from", "to")}
                 for prefix in ("ca.", "cb.", "n.")]
    angles = {"omega_rad": _dihedral(ca["from"], cb["from"], cb["to"], ca["to"]),
              "cb_distance_A": np.linalg.norm(cb["to"] - cb["from"], axis=1)}
    for side, other in (("from", "to"), ("to", "from")):
        angles[f"theta.{side}_rad"] = _dihedral(n[side], ca[side], cb[side], cb[other])
        angles[f"phi.{side}_rad"] = _angle(ca[side], cb[side], cb[other])
    return table.with_columns([pl.Series(k, v, dtype=pl.Float64).fill_nan(None)
                               for k, v in angles.items()]).sort("atom.index.from", "atom.index.to")


def summarize_contact_observations(table: pl.DataFrame) -> pl.DataFrame:
    """Summarize unordered atom pairs and directed binder→target residue pairs.

    Residue ``chem.a``/``chem.b`` mean binder/target respectively: TCR→peptide,
    TCR→MHC and peptide→MHC. Atom identities alone use canonical min/max order.

    Distances are atom distances in atom rows and closest retained atom distances
    in residue rows. Counts describe contacts, not independent observations.
    """
    context = ["mhc_class", "interface", "receptor_region"]
    keys = ["pdb.id", "chain.id.from", "residue.index.from", "chain.id.to", "residue.index.to"]
    frames = []
    flags = [c for c in table.columns if c.startswith("is_")]
    for level in ("atom", "residue"):
        frame = table
        if level == "residue":
            binder = pl.when(pl.col("interface") == "peptide_mhc").then(pl.lit("peptide")).otherwise(pl.lit("tcr"))
            swap = pl.col("role.from") != binder
            frame = frame.with_columns(*[
                pl.when(swap).then(pl.col(field + "." + other)).otherwise(pl.col(field + "." + side)).alias(field + "." + side)
                for field in ("chain.id", "residue.index", "residue.aa")
                for side, other in (("from", "to"), ("to", "from"))])
            frame = frame.group_by(*keys, *context, "residue.aa.from", "residue.aa.to").agg(
                pl.col("dist").min(), *[pl.col(c).any() for c in flags]).with_columns(
                pl.col("residue.aa.from").alias("chem.a"),
                pl.col("residue.aa.to").alias("chem.b"))
        frames.append(frame.group_by(*context, "chem.a", "chem.b").agg(
            pl.len().alias("n_contacts"), pl.col("pdb.id").n_unique().alias("n_structures"),
            pl.col("dist").mean().alias("mean_distance_A"),
            pl.col("dist").min().alias("min_distance_A"),
            pl.col("dist").max().alias("max_distance_A"),
            *[pl.col(c).sum().alias("n_" + c[3:]) for c in flags])
            .with_columns(pl.lit(level).alias("level")))
    return pl.concat(frames).sort("level", *context, "chem.a", "chem.b")
