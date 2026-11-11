"""One normalized molecular graph owns contacts, residue aggregation and counts."""
from __future__ import annotations

from collections import Counter, defaultdict
import json

import numpy as np
from scipy.spatial import cKDTree

from .frame import CanonicalFrame
from .topology import BACKBONE, compatible, conformer_rank, polymer_bonds

_ROLES = {"TRA": "alpha", "TRB": "beta", "PEPTIDE": "peptide", "MHCa": "mhc", "MHCb": "mhc", "B2M": "b2m"}


def _key(*values):
    return json.dumps(values, separators=(",", ":"))


def summarize_graph(atoms, residues, bonds, contacts):
    """Populate exact global degrees and aggregate the single atom-contact table.

    Contact count = atom-pair witnesses incident on a node. Contact degree =
    distinct partner atoms for atoms, distinct partner residues for residues.
    Covalent degree counts distinct bonded atom/residue neighbours respectively.
    These are global counts; the viewer computes visibly filtered counts separately.
    """
    atom_bonds, atom_contacts = defaultdict(set), defaultdict(set)
    residue_bonds, residue_partners = defaultdict(set), defaultdict(set)
    residue_atoms, residue_contact_atoms = defaultdict(set), defaultdict(set)
    counts = Counter()
    for atom in atoms:
        residue_atoms[atom["residue"]].add(atom["id"])
    for edge in bonds:
        a, b = edge["source"], edge["target"]
        atom_bonds[a].add(b)
        atom_bonds[b].add(a)
        ra, rb = atoms[a]["residue"], atoms[b]["residue"]
        if ra != rb:
            residue_bonds[ra].add(rb)
            residue_bonds[rb].add(ra)
    for edge in contacts:
        a, b = edge["source"], edge["target"]
        atom_contacts[a].add(b)
        atom_contacts[b].add(a)
        ra, rb = atoms[a]["residue"], atoms[b]["residue"]
        residue_partners[ra].add(rb)
        residue_partners[rb].add(ra)
        residue_contact_atoms[ra].add(a)
        residue_contact_atoms[rb].add(b)
        counts[tuple(sorted((ra, rb)))] += 1
    residue_counts = Counter()
    for (a, b), count in counts.items():
        residue_counts[a] += count
        residue_counts[b] += count
    for atom in atoms:
        partners = atom_contacts[atom["id"]]
        atom.update(covalent_degree=len(atom_bonds[atom["id"]]),
                    contact_count=len(partners), contact_degree=len(partners),
                    partner_residue_count=len({atoms[p]["residue"] for p in partners}),
                    partner_chain_count=len({residues[atoms[p]["residue"]]["chain"] for p in partners}))
    for residue in residues:
        rid = residue["id"]
        residue.update(atom_count=len(residue_atoms[rid]),
                       contacting_atom_count=len(residue_contact_atoms[rid]),
                       contact_count=residue_counts[rid], contact_degree=len(residue_partners[rid]),
                       covalent_degree=len(residue_bonds[rid]),
                       partner_chain_count=len({residues[p]["chain"] for p in residue_partners[rid]}))
    return [dict(source=a, target=b, contact_count=n) for (a, b), n in sorted(counts.items())]


def build_scene(structure, frame: CanonicalFrame, *, cutoff: float = 5.0) -> dict:
    """Build all typed polymer atoms once; residue views resolve owning CA atoms.

    Heavy-atom contacts use inclusive raw-coordinate distances, compatible alternate
    locations and distinct chains; directly covalently bonded atoms are excluded.
    Frameworks and other resolved regions are retained for chain/region selection.
    No independent residue-coordinate or residue-contact calculation exists.
    """
    if not np.isfinite(cutoff) or cutoff <= 0:
        raise ValueError("cutoff must be finite and positive")
    residues, atoms, source_coords, warnings = [], [], [], []
    chain_atoms = defaultdict(list)
    for chain in structure.chains:
        role = _ROLES.get(chain.chain_type)
        if role is None:
            warnings.append(f"untyped_chain:{chain.chain_id}")
            continue
        regions = {r.seq_index: reg.region_type for reg in chain.regions for r in reg.residues}
        if role not in {"peptide", "b2m"} and not regions:
            warnings.append(f"missing_region_annotation:{chain.chain_id}")
        for res in chain.residues:
            rid = len(residues)
            region = "PEPTIDE" if role == "peptide" else regions.get(res.seq_index, "UNANNOTATED")
            residue = dict(id=rid, key=_key(structure.pdb_id, chain.chain_id, res.seq_index),
                           chain=chain.chain_id, role=role, region=region, seq_index=res.seq_index,
                           number=res.pdb_index, insertion=res.insertion_code, aa=res.aa,
                           resname=res.resname, ca_atom=None)
            residues.append(residue)
            for ordinal, atom in enumerate(res.atoms):
                if atom.element.upper() in {"H", "D"}:
                    continue
                aid = len(atoms)
                atoms.append(dict(id=aid, key=_key(structure.pdb_id, chain.chain_id, res.seq_index, ordinal),
                                  ordinal=ordinal, residue=rid, name=atom.name, element=atom.element,
                                  altloc=atom.altloc, backbone=atom.name in BACKBONE,
                                  xyz=None))
                source_coords.append(atom.coord)
                chain_atoms[chain.chain_id].append(aid)
                if atom.name == "CA":
                    previous_ca = residue["ca_atom"]
                    if previous_ca is None or conformer_rank(atom.altloc) < conformer_rank(atoms[previous_ca]["altloc"]):
                        residue["ca_atom"] = aid
            if residue["ca_atom"] is None:
                warnings.append(f"missing_CA:{residue['key']}")
    coordinates = np.asarray(source_coords, dtype=float).reshape((-1, 3))
    for atom, xyz in zip(atoms, frame.transform(coordinates)):
        atom["xyz"] = xyz.tolist()
    bonds, bond_warnings = polymer_bonds(atoms, residues, coordinates)
    warnings.extend(bond_warnings)
    covalent_pairs = {(b["source"], b["target"]) for b in bonds}
    contacts = []
    chains = sorted(chain_atoms)
    trees = {c: cKDTree(coordinates[chain_atoms[c]]) for c in chains}
    for index, left in enumerate(chains):
        for right in chains[index + 1:]:
            for i, neighbours in enumerate(trees[left].query_ball_tree(trees[right], cutoff)):
                a = chain_atoms[left][i]
                for j in neighbours:
                    b = chain_atoms[right][j]
                    if not compatible(atoms[a], atoms[b]) or tuple(sorted((a, b))) in covalent_pairs:
                        continue
                    contacts.append(dict(source=min(a, b), target=max(a, b),
                                         distance=float(np.linalg.norm(coordinates[a] - coordinates[b]))))
    contacts.sort(key=lambda edge: (edge["source"], edge["target"]))
    residue_contacts = summarize_graph(atoms, residues, bonds, contacts)
    contact_groups = Counter()
    for edge in residue_contacts:
        roles = sorted((residues[edge["source"]]["role"], residues[edge["target"]]["role"]))
        contact_groups[":".join(roles)] += edge["contact_count"]
    tcr_peptide_contacts = contact_groups["alpha:peptide"] + contact_groups["beta:peptide"]
    trace = []
    for edge in bonds:
        a, b = atoms[edge["source"]], atoms[edge["target"]]
        ra, rb = a["residue"], b["residue"]
        if edge["kind"] == "backbone" and ra != rb and residues[ra]["ca_atom"] is not None and residues[rb]["ca_atom"] is not None:
            trace.append([min(ra, rb), max(ra, rb)])
    return dict(schema="tcren.complementarity/2", id=structure.pdb_id,
                provenance=structure.annotation_provenance,
                components=[dict(chain=c.chain_id, role=_ROLES[c.chain_type], sequence=c.sequence(),
                                 allele=c.allele_info, regions={r.region_type: r.sequence for r in c.regions})
                            for c in structure.chains if c.chain_type in _ROLES],
                frame=dict(reference=frame.reference_id, origin=frame.origin,
                           origin_component=frame.origin_component,
                           anchor_role_positions=[list(key) for key in frame.anchor_role_positions],
                           rotation=frame.rotation.tolist(), translation=frame.translation.tolist(),
                           rmsd=frame.rmsd, anchors=frame.anchor_count),
                cutoff=cutoff, residues=residues, atoms=atoms, bonds=bonds,
                trace=[list(pair) for pair in sorted({tuple(pair) for pair in trace})],
                atom_contacts=contacts, residue_contacts=residue_contacts,
                summary=dict(atoms=len(atoms), residues=len(residues), covalent_bonds=len(bonds),
                             atom_contacts=len(contacts), residue_contacts=len(residue_contacts),
                             tcr_peptide_atom_contacts=tcr_peptide_contacts,
                             no_tcr_peptide_interface=tcr_peptide_contacts == 0,
                             contacts_by_role_pair=dict(contact_groups)),
                warnings=warnings)
