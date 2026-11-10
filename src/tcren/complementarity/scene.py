"""Versioned residue/atom scenes for annotated structures; no scoring side effects."""
from __future__ import annotations

from collections import Counter

import numpy as np
from scipy.spatial import cKDTree

from .frame import CanonicalFrame

_REGIONS = {"CDR1", "CDR2", "CDR3", "HELIX_A1", "HELIX_A2", "HELIX_B1", "GROOVE_FLOOR"}
_BACKBONE = {"N", "CA", "C", "O", "OXT"}
_ROLES = {"TRA": "alpha", "TRB": "beta", "PEPTIDE": "peptide", "MHCa": "mhc", "MHCb": "mhc"}
_RADII = {"C": .76, "N": .71, "O": .66, "S": 1.05, "P": 1.07, "SE": 1.20}


def build_scene(structure, frame: CanonicalFrame, *, cutoff: float = 5.0) -> dict:
    """Export all interface atoms and exact heavy-atom contact witnesses.

    Nodes have scene-local stable integer IDs. Residue keys retain author insertion
    codes, and atom IDs distinguish repeated names/alternate locations. Contacts
    use inclusive Euclidean distance in Angstroms. Hydrogen/deuterium are omitted.
    Trace links require both consecutive sequence indices and CA distance <=4.5 A.
    Atom backbone links require a named backbone bond and distance <=1.9 A;
    these links are a visual trace, not an inferred complete chemical bond graph.
    """
    if not np.isfinite(cutoff) or cutoff <= 0:
        raise ValueError("cutoff must be finite and positive")
    residues, atoms, trace, backbone, bonds = [], [], [], [], []
    warnings = []
    source_coords = []
    for chain in structure.chains:
        role = _ROLES.get(chain.chain_type)
        if role is None:
            continue
        regions = {r.seq_index: reg.region_type for reg in chain.regions for r in reg.residues}
        if role != "peptide" and not regions:
            warnings.append(f"{chain.chain_id}: missing region annotation")
        previous = None
        previous_c = None
        for res in chain.residues:
            region = "PEPTIDE" if role == "peptide" else regions.get(res.seq_index)
            if region not in _REGIONS | {"PEPTIDE"}:
                previous = previous_c = None
                continue
            rid = len(residues)
            ca = res.ca
            node = dict(id=rid, chain=chain.chain_id, role=role, region=region,
                        seq_index=res.seq_index, number=res.pdb_index,
                        insertion=res.insertion_code, aa=res.aa,
                        xyz=None if ca is None else frame.transform(ca).tolist())
            residues.append(node)
            if ca is None:
                warnings.append(f"{chain.chain_id}:{res.pdb_index}{res.insertion_code}: missing CA")
            if previous is not None and node["xyz"] is not None:
                if (res.seq_index == previous["seq_index"] + 1
                        and np.linalg.norm(np.array(node["xyz"]) - previous["xyz"]) <= 4.5):
                    trace.append([previous["id"], rid])
            previous = node if ca is not None else None
            named = {}
            residue_atoms = []
            for atom in res.atoms:
                if atom.element.upper() in {"H", "D"}:
                    continue
                aid = len(atoms)
                atoms.append(dict(id=aid, residue=rid, name=atom.name, element=atom.element,
                                  altloc=atom.altloc, backbone=atom.name in _BACKBONE,
                                  xyz=frame.transform(atom.coord).tolist()))
                source_coords.append(atom.coord)
                named.setdefault(atom.name, []).append(aid)
                residue_atoms.append(aid)

            for offset, i in enumerate(residue_atoms):
                for j in residue_atoms[offset + 1:]:
                    left, right = atoms[i], atoms[j]
                    if left["name"] == right["name"]:
                        continue
                    if left["altloc"] and right["altloc"] and left["altloc"] != right["altloc"]:
                        continue
                    radii = [_RADII.get(a["element"].upper()) for a in (left, right)]
                    if None in radii:
                        continue
                    distance = np.linalg.norm(np.array(left["xyz"]) - right["xyz"])
                    if .4 < distance <= 1.2 * sum(radii):
                        bonds.append([i, j])

            def link(a, b):
                for i in a:
                    for j in b:
                        left, right = atoms[i], atoms[j]
                        if left["altloc"] and right["altloc"] and left["altloc"] != right["altloc"]:
                            continue
                        distance = np.linalg.norm(np.array(left["xyz"]) - right["xyz"])
                        if 0.4 < distance <= 1.9:
                            backbone.append([i, j])

            for a, b in [("N", "CA"), ("CA", "C"), ("C", "O"), ("C", "OXT")]:
                link(named.get(a, []), named.get(b, []))
            if previous_c is not None and previous_c[0] + 1 == res.seq_index:
                link(previous_c[1], named.get("N", []))
            previous_c = (res.seq_index, named.get("C", []))
    contacts = []
    counts = Counter()
    if atoms:
        coords = np.array(source_coords)
        for i, j in sorted(cKDTree(coords).query_pairs(cutoff)):
            a, b = atoms[i], atoms[j]
            ra, rb = residues[a["residue"]], residues[b["residue"]]
            if ra["chain"] == rb["chain"]:
                continue
            # Alternate locations are alternative conformers, not simultaneous atoms.
            if a["altloc"] and b["altloc"] and a["altloc"] != b["altloc"]:
                continue
            contacts.append([i, j, float(np.linalg.norm(coords[i] - coords[j]))])
            counts[tuple(sorted((ra["id"], rb["id"])))] += 1
    return dict(schema="tcren.complementarity/1", id=structure.pdb_id,
                provenance=structure.annotation_provenance,
                components=[dict(chain=c.chain_id, role=_ROLES[c.chain_type],
                                 sequence=c.sequence(), allele=c.allele_info,
                                 regions={r.region_type: r.sequence for r in c.regions})
                            for c in structure.chains if c.chain_type in _ROLES],
                frame=dict(reference=frame.reference_id, origin=frame.origin,
                           rotation=frame.rotation.tolist(), translation=frame.translation.tolist(),
                           rmsd=frame.rmsd, anchors=frame.anchor_count),
                cutoff=cutoff, residues=residues, atoms=atoms, trace=trace,
                backbone=backbone, bonds=sorted({tuple(b) for b in bonds + backbone}), atom_contacts=contacts,
                residue_contacts=[[a, b, n] for (a, b), n in sorted(counts.items())],
                warnings=warnings)
