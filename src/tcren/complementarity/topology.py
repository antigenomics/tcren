"""Named polymer connectivity; distance validation never substitutes for atom identity."""
from __future__ import annotations

import numpy as np

_SIDECHAINS = {
    "A": "CA-CB", "G": "", "V": "CA-CB CB-CG1 CB-CG2",
    "L": "CA-CB CB-CG CG-CD1 CG-CD2", "I": "CA-CB CB-CG1 CB-CG2 CG1-CD1",
    "S": "CA-CB CB-OG", "T": "CA-CB CB-OG1 CB-CG2", "C": "CA-CB CB-SG",
    "M": "CA-CB CB-CG CG-SD SD-CE", "P": "CA-CB CB-CG CG-CD CD-N",
    "D": "CA-CB CB-CG CG-OD1 CG-OD2", "N": "CA-CB CB-CG CG-OD1 CG-ND2",
    "E": "CA-CB CB-CG CG-CD CD-OE1 CD-OE2", "Q": "CA-CB CB-CG CG-CD CD-OE1 CD-NE2",
    "K": "CA-CB CB-CG CG-CD CD-CE CE-NZ", "R": "CA-CB CB-CG CG-CD CD-NE NE-CZ CZ-NH1 CZ-NH2",
    "H": "CA-CB CB-CG CG-ND1 CG-CD2 ND1-CE1 CD2-NE2 CE1-NE2",
    "F": "CA-CB CB-CG CG-CD1 CG-CD2 CD1-CE1 CD2-CE2 CE1-CZ CE2-CZ",
    "Y": "CA-CB CB-CG CG-CD1 CG-CD2 CD1-CE1 CD2-CE2 CE1-CZ CE2-CZ CZ-OH",
    "W": "CA-CB CB-CG CG-CD1 CG-CD2 CD1-NE1 NE1-CE2 CD2-CE2 CD2-CE3 CE2-CZ2 CE3-CZ3 CZ2-CH2 CZ3-CH2",
}
BACKBONE = frozenset({"N", "CA", "C", "O", "OXT"})


def conformer_rank(altloc):
    """Deterministic representative: blank, A, then remaining labels."""
    return altloc not in {"", "A"}, altloc


def compatible(left, right):
    """Blank atoms are shared; different nonblank conformers are alternatives."""
    return not (left["altloc"] and right["altloc"] and left["altloc"] != right["altloc"])


def polymer_bonds(atoms, residues, coordinates):
    """Return unique distance-validated named bonds and explicit missing/long flags.

    Standard amino-acid templates specify connectivity, without assigning bond
    orders. MSE uses the methionine template with SE in place of SD. Peptide C–N
    bonds require consecutive resolved sequence positions and <=1.9 A. Unknown
    residues retain a backbone template and an explicit sidechain warning.
    """
    by_residue = [[] for _ in residues]
    for atom in atoms:
        by_residue[atom["residue"]].append(atom)
    bonds, warnings = {}, []

    def add(left, right, kind, maximum):
        if not compatible(left, right):
            return
        distance = float(np.linalg.norm(coordinates[left["id"]] - coordinates[right["id"]]))
        if not .4 < distance <= maximum:
            warnings.append(f"bond_unavailable:{left['id']}:{right['id']}:{distance:.3f}")
            return
        key = tuple(sorted((left["id"], right["id"])))
        bonds[key] = dict(source=key[0], target=key[1], kind=kind, distance=distance,
                          method="named-polymer-template")

    previous = {}
    for residue, members in zip(residues, by_residue):
        named = {}
        for atom in members:
            named.setdefault(atom["name"], []).append(atom)
        missing_backbone = {"N", "CA", "C", "O"} - named.keys()
        if missing_backbone:
            warnings.append(f"missing_backbone_atoms:{residue['id']}:" + ",".join(sorted(missing_backbone)))
        sidechain = _SIDECHAINS.get(residue["aa"])
        if sidechain is None:
            warnings.append(f"unknown_sidechain_template:{residue['id']}:{residue['resname']}")
            sidechain = ""
        if residue["resname"] == "MSE":
            sidechain = sidechain.replace("SD", "SE")
        supported = BACKBONE | {name for pair in sidechain.split() for name in pair.split("-")}
        for atom in members:
            if atom["name"] not in supported:
                warnings.append(f"untemplated_atom:{atom['id']}:{atom['name']}:{residue['resname']}")
        for pair in ("N-CA CA-C C-O C-OXT " + sidechain).split():
            a, b = pair.split("-")
            for left in named.get(a, []):
                for right in named.get(b, []):
                    backbone = a in BACKBONE and b in BACKBONE
                    add(left, right, "backbone" if backbone else "sidechain", 1.9 if backbone else 2.3)
        last = previous.get(residue["chain"])
        if last and last[0] + 1 == residue["seq_index"]:
            for left in last[1]:
                for right in named.get("N", []):
                    add(left, right, "backbone", 1.9)
        previous[residue["chain"]] = residue["seq_index"], named.get("C", [])
    # Disulfides are covalent even when they bridge distinct chains.
    sulfurs = [a for a in atoms if a["name"] == "SG" and residues[a["residue"]]["aa"] == "C"]
    for i, left in enumerate(sulfurs):
        for right in sulfurs[i + 1:]:
            if compatible(left, right):
                distance = np.linalg.norm(coordinates[left["id"]] - coordinates[right["id"]])
                if 1.7 <= distance <= 2.3:
                    add(left, right, "disulfide", 2.3)
    return [bonds[key] for key in sorted(bonds)], warnings
