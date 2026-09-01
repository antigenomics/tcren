"""Backbone-preserving peptide substitution.

``score_peptides`` scores a candidate peptide *virtually* (it re-indexes the potential matrix over
the native contact map — no atoms move). When you want to actually re-dock / refine a candidate you
first need its coordinates: :func:`substitute_peptide` threads an equal-length sequence onto the
existing peptide backbone, keeping N/Cα/C/O(+Cβ) and dropping the old side-chain atoms beyond Cβ
(a refiner / rotamer repack rebuilds them). Pure data-model manipulation; returns a new structure.

**Alanine is the case this gets exactly right, and it matters.** Alanine's heavy atoms are N, Cα,
C, O and Cβ and no others, so a substitution to alanine needs no rotamer, no relaxation and no
choice: truncating at Cβ *is* the alanine. Every other target is left as a Cβ stub whose reach is
therefore under-stated, and needs a side-chain builder
(:func:`tcren.rotamers.repack` rotates existing atoms but cannot create them; ProMod3's
``ReconstructSidechains``, wired up in :mod:`tcren.refine.engines.promod3_engine`, can).

Glycine has no Cβ to keep, so mutating a glycine to anything else needs one built.
:func:`virtual_cb` places it from N, Cα and C by ideal L-amino-acid geometry; measured against
1,679 real crystallographic Cβ atoms the construction lands a median 0.09 Å away (99th percentile
0.31 Å), which is far inside the 5 Å contact definition it feeds.
"""

from __future__ import annotations

from collections.abc import Mapping

import numpy as np

from ..structure.model import PEPTIDE_TYPE, Atom, Chain, RegionMarkup, Residue, Structure

# Atoms retained on substitution: backbone + Cβ (the rest of the side chain is identity-specific).
_KEEP = {"N", "CA", "C", "O", "CB"}

_ONE_TO_THREE = {
    "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS", "Q": "GLN", "E": "GLU",
    "G": "GLY", "H": "HIS", "I": "ILE", "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE",
    "P": "PRO", "S": "SER", "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL",
}


def virtual_cb(n: np.ndarray, ca: np.ndarray, c: np.ndarray) -> np.ndarray:
    """Cβ position from the backbone, by ideal L-amino-acid geometry.

    The standard construction: with ``b = Cα − N`` and ``c = C − Cα``, the Cβ sits along a fixed
    combination of ``b``, ``c`` and their cross product, which fixes both the ~1.53 Å bond length
    and the tetrahedral chirality. Used only where there is no Cβ to keep — i.e. mutating away
    from glycine.
    """
    b, cc = ca - n, c - ca
    a = np.cross(b, cc)
    return -0.58273431 * a + 0.56802827 * b - 0.54067466 * cc + ca


def _rewrite_residue(res: Residue, aa: str) -> Residue:
    """One residue re-typed to ``aa``, truncated to backbone + Cβ.

    Alanine's heavy atoms are exactly N, Cα, C, O and Cβ, so for an alanine target this *is* the
    mutant: no rotamer, no relaxation, no choice. Mutating to glycine drops the Cβ; mutating away
    from glycine builds one by ideal geometry, because a residue with no Cβ reaches no further
    than its backbone and the contact map would under-count it.
    """
    resname = _ONE_TO_THREE.get(aa)
    if resname is None:
        raise ValueError(f"non-standard amino acid {aa!r}")
    keep = _KEEP - ({"CB"} if aa == "G" else set())
    atoms = tuple(a for a in res.atoms if a.name in keep)
    if aa != "G" and not any(a.name == "CB" for a in atoms):
        bb = {a.name: a.coord for a in res.atoms}
        if {"N", "CA", "C"} <= bb.keys():
            atoms += (Atom(name="CB", element="C",
                           coord=virtual_cb(bb["N"], bb["CA"], bb["C"])),)
    return Residue(res.seq_index, res.pdb_index, res.insertion_code, aa, resname, atoms)


def _rebuild_chain(chain: Chain, new_residues: list[Residue]) -> Chain:
    """A chain carrying ``new_residues``, with its region markup re-pointed at them.

    Without the re-point the contact map has null ``pos.from``/``pos.to`` and every downstream
    scorer fails.
    """
    by_index = {r.seq_index: r for r in new_residues}
    return Chain(chain_id=chain.chain_id, residues=new_residues,
                 chain_type=chain.chain_type, chain_supertype=chain.chain_supertype,
                 allele_info=chain.allele_info,
                 regions=[RegionMarkup(
                     region_type=reg.region_type, start_seq_index=reg.start_seq_index,
                     end_seq_index=reg.end_seq_index,
                     sequence="".join(by_index[r.seq_index].aa for r in reg.residues),
                     residues=[by_index[r.seq_index] for r in reg.residues],
                 ) for reg in chain.regions])


def substitute_residues(
    structure: Structure, mutations: "Mapping[tuple[str, int], str]"
) -> Structure:
    """Return a copy of ``structure`` with the named residues re-typed, in 3D.

    The general primitive behind every structural substitution: it moves no backbone, drops
    side-chain atoms past Cβ on the residues it touches, and leaves every other residue and every
    other chain byte-identical. **Any chain may be targeted** — this is what makes a receptor-side
    alanine scan possible, where :func:`substitute_peptide` threads the peptide chain alone.

    Args:
        structure: The structure to mutate.
        mutations: ``(chain_id, residue.seq_index) -> one-letter target``. An empty mapping
            returns the structure unchanged.

    Returns:
        A new :class:`~tcren.structure.model.Structure`. Region markup on every touched chain is
        re-pointed at the rewritten residues, so the result goes straight into a contact map.

    Raises:
        ValueError: if a chain id or residue index is absent, or a code is non-standard.
    """
    if not mutations:
        return structure
    by_chain: dict[str, dict[int, str]] = {}
    for (chain_id, seq_index), aa in mutations.items():
        by_chain.setdefault(chain_id, {})[seq_index] = aa.upper()

    chains = []
    for chain in structure.chains:
        want = by_chain.pop(chain.chain_id, None)
        if not want:
            chains.append(chain)
            continue
        present = {r.seq_index for r in chain.residues}
        missing = sorted(set(want) - present)
        if missing:
            raise ValueError(
                f"chain {chain.chain_id!r} has no residue at seq_index {missing} "
                f"(it runs {min(present)}..{max(present)})"
            )
        chains.append(_rebuild_chain(
            chain,
            [_rewrite_residue(r, want[r.seq_index]) if r.seq_index in want else r
             for r in chain.residues],
        ))
    if by_chain:
        raise ValueError(f"no such chain(s) in {structure.pdb_id!r}: {sorted(by_chain)}")
    return Structure(pdb_id=structure.pdb_id, chains=chains,
                     complex_species=structure.complex_species, cell_type=structure.cell_type)


def substitute_peptide(structure: Structure, new_peptide: str,
                       chain_type: str = PEPTIDE_TYPE) -> Structure:
    """Return a copy of ``structure`` with the peptide chain threaded to ``new_peptide``.

    The peptide backbone (and Cβ) is preserved; side-chain atoms beyond Cβ are dropped (and Cβ
    too for any position mutated to glycine, while a position mutated *from* glycine has one built
    by :func:`virtual_cb`). ``new_peptide`` must equal the peptide length and use the 20 standard
    one-letter amino acids. Region markup is carried over onto the new residues, so the result can
    go straight into a contact map and be scored.

    For ``new_peptide`` all-alanine this is exact, which is what makes a *structural* poly-alanine
    reference possible (:func:`tcren.ddg.reference_delta` with a structure). For other targets the
    result is a Cβ stub: correct in position, short in reach.

    Raises:
        ValueError: if there is no peptide chain, the length differs, or a code is non-standard.
    """
    new_peptide = new_peptide.upper()
    pep = next((c for c in structure.chains if c.chain_type == chain_type), None)
    if pep is None:
        raise ValueError(f"no {chain_type} chain in structure {structure.pdb_id!r}")
    if len(new_peptide) != len(pep.residues):
        raise ValueError(
            f"length mismatch: peptide has {len(pep.residues)} residues, got {len(new_peptide)}"
        )

    new_residues = [_rewrite_residue(res, aa) for res, aa in zip(pep.residues, new_peptide)]

    new_pep = _rebuild_chain(pep, new_residues)
    chains = [new_pep if c is pep else c for c in structure.chains]
    return Structure(pdb_id=structure.pdb_id, chains=chains,
                     complex_species=structure.complex_species, cell_type=structure.cell_type)
