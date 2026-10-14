"""Select a single TCR-pMHC complex and rename its chains to the canonical A–E scheme."""

from __future__ import annotations

import copy

from ..contacts.geometry import all_atom_contacts
from ..structure.model import Structure

# Canonical chain id by chain role.
CHAIN_RENAME = {
    "TRA": "A", "TRG": "A",   # VJ chain
    "TRB": "B", "TRD": "B",   # VDJ chain
    "PEPTIDE": "C",
    "MHCa": "D",
    "MHCb": "E", "B2M": "E",
}
_ROLE_ORDER = ("A", "B", "C", "D", "E")


def _has_multiple_copies(structure: Structure) -> bool:
    seen: dict[str, int] = {}
    for c in structure.chains:
        target = CHAIN_RENAME.get(c.chain_type)
        if target:
            seen[target] = seen.get(target, 0) + 1
    return any(n > 1 for n in seen.values())


def select_primary_complex(structure: Structure) -> Structure:
    """Keep one mutually-contacting chain per canonical role (one TCR-pMHC complex).

    Chosen as a connected unit so chains that do not touch the peptide (notably β2m) are not
    grabbed from another copy: the primary peptide is the one most embedded in an MHC-α groove
    (then most TCR contacts, then shortest); the TCR α/β and MHC-α are taken by contacts to
    that peptide, retaining chains connected through the selected core; β2m / MHC-β prefer
    contacts to the chosen MHC-α. Disconnected copies are omitted. No-op for single complexes.
    """
    if not _has_multiple_copies(structure):
        return structure
    contacts = all_atom_contacts(structure, cutoff=8.0)
    pair_counts: dict[tuple[str, str], int] = {}
    for row in contacts.iter_rows(named=True):
        a, b = row["chain.id.from"], row["chain.id.to"]
        pair_counts[(a, b)] = pair_counts.get((a, b), 0) + 1
        pair_counts[(b, a)] = pair_counts.get((b, a), 0) + 1

    by_role: dict[str, list] = {}
    for c in structure.chains:
        target = CHAIN_RENAME.get(c.chain_type)
        if target:
            by_role.setdefault(target, []).append(c)
    peptides = by_role.get("C", [])
    if not peptides:
        return structure

    def contacts_to(chain, others) -> int:
        return sum(pair_counts.get((chain.chain_id, o.chain_id), 0) for o in others)

    mhca = by_role.get("D", [])
    tcr_all = by_role.get("A", []) + by_role.get("B", [])
    # primary peptide = the groove peptide (most MHC-α contacts), then TCR contacts, then shortest
    primary_pep = max(peptides, key=lambda p: (contacts_to(p, mhca), contacts_to(p, tcr_all),
                                               -len(p.residues)))
    kept = [primary_pep]
    pending = [r for r in ("D", "A", "B", "E") if by_role.get(r)]
    while pending:
        progressed = False
        for role in pending[:]:
            anchor = (next((c for c in kept if CHAIN_RENAME.get(c.chain_type) == "D"), primary_pep)
                      if role == "E" else primary_pep)
            selected = max(by_role[role], key=lambda c: (contacts_to(c, [anchor]), contacts_to(c, kept)))
            if contacts_to(selected, kept) > 0:
                kept.append(selected)
                pending.remove(role)
                progressed = True
        if not progressed:
            break
    kept_ids = {c.chain_id for c in kept}
    new = copy.copy(structure)
    new.chains = [c for c in structure.chains if c.chain_id in kept_ids]
    return new


def rename_chains(structure: Structure) -> tuple[Structure, dict[str, str]]:
    """Return a copy with **only** the canonical complex, chain ids remapped per
    :data:`CHAIN_RENAME`, plus the old→new map.

    Chains with no canonical role (tags, additives, unrelated proteins) are dropped so the
    output is exactly the A–E TCR-pMHC complex. Raises ``ValueError`` if two source chains map
    to the same canonical id (unresolved multi-copy — run :func:`select_primary_complex` first).
    """
    chain_map: dict[str, str] = {}
    used: dict[str, str] = {}
    kept = []
    for c in structure.chains:
        target = CHAIN_RENAME.get(c.chain_type)
        if target is None:
            continue
        if target in used:
            raise ValueError(
                f"chain collision on {target!r}: {used[target]} and {c.chain_id} "
                f"(run select_primary_complex first)"
            )
        used[target] = c.chain_id
        chain_map[c.chain_id] = target
        nc = copy.copy(c)
        nc.chain_id = target
        kept.append(nc)
    new = copy.copy(structure)
    order = {cid: i for i, cid in enumerate(_ROLE_ORDER)}
    new.chains = sorted(kept, key=lambda c: order.get(c.chain_id, 99))
    return new, chain_map
