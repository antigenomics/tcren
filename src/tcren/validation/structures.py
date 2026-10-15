"""Sequence-verified coordinate inspection and explicit component comparisons."""

from __future__ import annotations

import hashlib
from itertools import combinations
import json
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

from ..contacts.geometry import _atom_arrays, all_atom_contacts
from ..structure import parse_structure


def read_chains(path):
    """Read exact polymer sequences and Cα arrays through tcren's structure parser."""
    s = parse_structure(path, keep_hydrogens=False)
    rows = []
    for chain in s.chains:
        if any(r.ca is None for r in chain.residues):
            raise ValueError(f"Missing CA: {path}:{chain.chain_id}")
        rows.append((chain.chain_id, chain.sequence(), np.asarray([r.ca for r in chain.residues])))
    return rows


def map_sequence_roles(path, roles, sequences, *, max_terminal_missing=4, allow_extended=False):
    """Match components uniquely, with explicit opt-in for longer reference controls.

    Extended output chains are compared only over the exact requested sequence
    span. Both reference and output starts are reported; default matching stays
    strict apart from the existing terminal-missing allowance.
    """
    mapped, identity = {}, []
    chains = read_chains(path)
    for role in roles:
        sequence = sequences[role]
        matches = []
        for cid, seq, coords in chains:
            if seq == sequence:
                matches.append((cid, coords, 0, len(seq), 0))
            elif seq in sequence and len(seq) >= len(sequence) - max_terminal_missing:
                start = sequence.index(seq)
                if sequence.find(seq, start + 1) >= 0:
                    raise ValueError(f"Ambiguous repeated sequence span: {path}:{role}")
                matches.append((cid, coords, start, len(seq), 0))
            elif allow_extended and sequence in seq:
                start = seq.index(sequence)
                if seq.find(sequence, start + 1) >= 0:
                    raise ValueError(f"Ambiguous repeated output span: {path}:{role}")
                matches.append(
                    (cid, coords[start : start + len(sequence)], 0, len(sequence), start)
                )
        if len(matches) != 1:
            raise ValueError(
                f"Absent/ambiguous sequence match {path}: {role}, {len(matches)} matches"
            )
        cid, coords, start, length, output_start = matches[0]
        mapped[role] = (coords, start, length)
        identity.append(
            dict(
                role=role,
                output_chain=cid,
                reference_start=start,
                output_start=output_start,
                matched_residues=length,
                input_residues=len(sequence),
            )
        )
    if len({row["output_chain"] for row in identity}) != len(identity):
        raise ValueError("multiple roles matched the same coordinate chain")
    return mapped, identity


def parse_pdb_text(text, source="structure"):
    """Parse streamed PDB coordinates with the same parser used by tcren file inputs."""
    from ..structure.io import _parse_pdb_fast
    from tempfile import TemporaryDirectory

    s = _parse_pdb_fast(text.encode(), source, False)
    if s is not None:
        return s
    with TemporaryDirectory() as temporary:
        path = Path(temporary) / "structure.pdb"
        path.write_text(text)
        return parse_structure(path, pdb_id=source, keep_hydrogens=False)


def coordinate_diagnostics(path) -> dict:
    """Finite coordinates, missing backbone atoms and close-pair counts; not MolProbity.

    First model, every alternate location, no hydrogens, matching tcren's contact convention.
    Bond chemistry is not inferred: same-chain nonadjacent close pairs include disulfides.
    The checksum covers the actual source file bytes, including compression when present.
    """
    return _diagnostics(
        parse_structure(path, keep_hydrogens=False),
        str(path),
        hashlib.sha256(Path(path).read_bytes()).hexdigest(),
        "source file bytes",
    )


def close_atom_pairs(path, *, cutoff=2.0):
    """Identify nonlocal heavy-atom proximity pairs; these are not MolProbity clashes.

    Same-residue and same-chain neighboring residues are excluded. SG–SG CYS
    pairs are labelled possible disulfides, without asserting a covalent bond.
    Source alternate-location labels are retained; differing labels alone do not
    establish incompatible conformers across residues. Bond classification is
    explicitly unknown: adjacency exclusions are a sequence heuristic, not a bond graph.
    """
    if not np.isfinite(cutoff) or cutoff <= 0:
        raise ValueError("Expected positive finite cutoff in Å")
    structure = parse_structure(path, keep_hydrogens=False)
    atoms = [
        (c.chain_id, r, a)
        for c in structure.chains
        for r in c.residues
        for a in r.atoms
        if a.element != "H"
    ]
    xyz = np.asarray([a.coord for _, _, a in atoms])
    if not len(xyz) or not np.isfinite(xyz).all():
        raise ValueError("Empty/nonfinite atom coordinates")
    pairs = cKDTree(xyz).query_pairs(cutoff, output_type="ndarray")
    rows = []
    for i, j in pairs:
        c1, r1, a1 = atoms[i]
        c2, r2, a2 = atoms[j]
        if c1 == c2 and abs(r1.seq_index - r2.seq_index) <= 1:
            continue
        distance = float(np.linalg.norm(xyz[i] - xyz[j]))
        if distance >= cutoff:
            continue
        row = {
            "distance_A": distance,
            "interchain": c1 != c2,
            "possible_disulfide": a1.name == a2.name == "SG" and r1.aa == r2.aa == "C",
            "alternate_location_involved": bool(a1.altloc or a2.altloc),
            "bond_classification": "unknown",
        }
        for suffix, cid, residue, atom in [("from", c1, r1, a1), ("to", c2, r2, a2)]:
            row.update(
                {
                    f"chain.{suffix}": cid,
                    f"residue_index.{suffix}": residue.seq_index,
                    f"author_residue.{suffix}": residue.pdb_index,
                    f"insertion_code.{suffix}": residue.insertion_code,
                    f"residue.{suffix}": residue.aa,
                    f"atom.{suffix}": atom.name,
                    f"altloc.{suffix}": atom.altloc,
                }
            )
        rows.append(row)
    return sorted(
        rows,
        key=lambda row: (
            row["distance_A"],
            row["chain.from"],
            row["residue_index.from"],
            row["atom.from"],
        ),
    )


def diagnose_pdb_text(text, source):
    """Coordinate diagnostics on a streamed PDB, with a checksum of the actual text."""
    return _diagnostics(
        parse_pdb_text(text, source),
        source,
        hashlib.sha256(text.encode()).hexdigest(),
        "source PDB text",
    )


def _diagnostics(s, source, sha256, representation):
    xyz, chains, residues, _, _ = _atom_arrays(s)
    if not len(xyz) or not np.isfinite(xyz).all():
        raise ValueError(f"Empty/nonfinite coordinates: {source}")
    pairs = cKDTree(xyz).query_pairs(2.0, output_type="ndarray")
    inter = chains[pairs[:, 0]] != chains[pairs[:, 1]]
    distant = (~inter) & (np.abs(residues[pairs[:, 0]] - residues[pairs[:, 1]]) > 1)
    distances = np.linalg.norm(xyz[pairs[:, 0]] - xyz[pairs[:, 1]], axis=1)
    missing = 0
    for c in s.chains:
        for r in c.residues:
            names = {a.name for a in r.atoms}
            if names & {"OC1", "OT1"}:
                names.add("O")
            missing += len({"N", "CA", "C", "O"} - names)
    return dict(
        source=source,
        coordinate_sha256=sha256,
        coordinate_sha256_representation=representation,
        chain_lengths=json.dumps({c.chain_id: len(c.residues) for c in s.chains}, sort_keys=True),
        residues=sum(len(c.residues) for c in s.chains),
        heavy_atoms=len(xyz),
        missing_backbone_atoms=missing,
        interchain_heavy_pairs_lt2A=int(inter.sum()),
        samechain_nonadjacent_heavy_pairs_lt2A=int(distant.sum()),
        shortest_interchain_pair_lt2A_A=float(distances[inter].min()) if inter.any() else None,
    )


def geometry_and_contacts(path, identity, *, cutoff=4.5):
    """Explicit chain-pair coverage including zero-contact pairs, using tcren contacts."""
    s = parse_structure(path, keep_hydrogens=False)
    ct = all_atom_contacts(s, cutoff=cutoff, count_atoms=True)
    atom_sets = {
        c.chain_id: np.asarray([a.coord for r in c.residues for a in r.atoms]) for c in s.chains
    }
    by_role = {row["role"]: row["output_chain"] for row in identity}
    rows = []
    for left, right in combinations(sorted(by_role), 2):
        a, b = by_role[left], by_role[right]
        select = ct.filter(
            ((ct["chain.id.from"] == a) & (ct["chain.id.to"] == b))
            | ((ct["chain.id.from"] == b) & (ct["chain.id.to"] == a))
        )
        row = dict(
            role_from=left,
            role_to=right,
            cutoff_A=float(cutoff),
            heavy_atom_pairs=int(select["n_atom_contacts"].sum()),
            residue_pairs=select.height,
            shortest_heavy_atom_distance_A=float(
                cKDTree(atom_sets[b]).query(atom_sets[a])[0].min()
            ),
        )
        if cutoff == 4.5:
            row.update(
                heavy_atom_pairs_le4_5A=row["heavy_atom_pairs"],
                residue_pairs_le4_5A=row["residue_pairs"],
            )
        rows.append(row)
    return coordinate_diagnostics(path), rows


def compatible_reference_models(models, reference_roles, construct_roles):
    """Select by covered physical roles, preserving complete-query/partial-reference matches.

    ``construct_roles`` explicitly defines each caller's construct labels. Subsequent
    ``all`` fits must use only ``reference_roles``; unmatched components are not scored.
    """
    required = set(reference_roles)
    if not required:
        raise ValueError("Reference must identify at least one component role")
    return [model for model in models if required <= set(construct_roles[model["construct"]])]
