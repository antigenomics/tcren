"""Chemical typing of interface contacts — a DSSP-style annotation layer for TCR:pMHC contact maps.

Types each heavy-atom contact in a :class:`~tcren.contactmap.ContactMap` interface from geometry and
atom identity alone: no hydrogens (models and many crystals lack them), no external DSSP binary.

Two schemes ship, and the difference matters:

``"v2"`` (default)
    ``salt_bridge``, ``hydrogen_bond``, ``cation_pi``, ``stacking``, ``aromatic``, ``hydrophobic``,
    ``polar``, ``vdw``, ``other`` — where ``other`` now means only "too far to be anything",
    never "unrecognised". Apolarity is decided **per atom** (a carbon with no bonded N/O) rather than
    per residue, donors and acceptors are typed, and a contact may carry more than one type — the
    ``is_<type>`` booleans are independent and ``contact.type`` is only the highest-priority label.
``"v1"``
    The original five-type, residue-level, winner-takes-all scheme, kept byte-for-byte because the
    frozen recognition models in :mod:`tcren.recognition` were trained on its ``ct_*`` counts.

**Why v2 exists.** Measured on ``tests/assets/pdb/{1ao7,1bd2,2ckb,5m01,6bj3}``, v1 typed **72.3%** of
TCR:peptide contacts as ``other`` — not because those contacts are featureless but because of four
specific gaps. 59% of the ``other`` rows were mixed C–O or C–N pairs, for which v1 had no class at
all. Every C–C ``other`` row failed a residue-level apolarity test that excludes Tyr — the most
common TCR interface residue — along with the aliphatic CB/CG/CD of Arg, Lys and Gln. A further 29
N–O and O–O pairs sat between 3.5 Å and 5 Å, outside a heavy-atom H-bond cutoff that is strict when
hydrogens are absent. And ring stacking was measured in :mod:`tcren.stacking` but never joined here.

**Seeing more than one atom pair.** :func:`~tcren.contacts.geometry.all_atom_contacts` collapses each
residue pair to its closest atom pair, which hides a salt bridge whose nearest contact happens to be
between two carbons. Pass ``atom_pairs=True`` there, or use :func:`residue_pair_types`, to type from
every atom pair a residue pair makes.
"""
from __future__ import annotations

import polars as pl

# --- residue-level sets (v1 and, where still correct, v2) -----------------------------------------
_CATIONIC_RES = {"K", "R", "H"}
_ANIONIC_RES = {"D", "E"}
_AROMATIC_RES = {"F", "Y", "W", "H"}
_APOLAR_RES = {"A", "V", "L", "I", "M", "F", "W", "P", "C"}     # v1 only

# side-chain atom names carrying formal charge
_CATIONIC_ATOMS = {"NZ", "NE", "NH1", "NH2", "ND1", "NE2"}
_ANIONIC_ATOMS = {"OD1", "OD2", "OE1", "OE2"}
# aromatic ring atoms (Phe/Tyr/Trp/His)
_RING_ATOMS = {"CG", "CD1", "CD2", "CE1", "CE2", "CZ", "NE1", "CE3", "CZ2", "CZ3", "CH2", "ND1", "NE2"}

# Carbons bonded to a backbone or side-chain N/O, i.e. the ones that are NOT apolar. Everything else
# with element C or S is. Deriving apolarity per atom rather than per residue is what recovers Tyr's
# ring carbons, Arg/Lys/Gln's aliphatic stems, and Thr's CG2.
_BACKBONE_POLAR_C = {"C", "CA"}
_POLAR_SIDECHAIN_C = {
    "D": {"CG"}, "E": {"CD"}, "N": {"CG"}, "Q": {"CD"}, "S": {"CB"}, "T": {"CB"}, "Y": {"CZ"},
    "R": {"CD", "CZ"}, "K": {"CE"}, "H": {"CG", "CD2", "CE1"}, "W": {"CD1", "CE2"}, "P": {"CD"},
}
# Hydroxyl oxygens donate as well as accept; His ring nitrogens do both (tautomer-dependent). Every
# other N is a donor and every other O an acceptor, which is what makes O-O contacts an H-bond only
# when a hydroxyl is involved and N-N only when a histidine is.
_HYDROXYL = {("S", "OG"), ("T", "OG1"), ("Y", "OH")}

_HBOND_MAX = 3.9        # heavy-atom D...A; the permissive end, because there are no hydrogens to
_SALT_MAX = 4.0         # place. v1 used 3.5, which dropped 29 N-O/O-A pairs into `other`.
_AROMATIC_MAX = 5.0
_CATION_PI_MAX = 6.0
_HYDROPHOBIC_MAX = 4.5
_POLAR_MAX = 4.5
_VDW_MAX = 4.5          # a close contact that carries no other chemistry is still van der Waals

#: Ring-pair geometry accepted as a stack (see :func:`tcren.stacking.ring_stacking`): centroids close
#: enough to interact, and either near-parallel (a face-to-face or parallel-displaced stack) or
#: near-perpendicular (edge-to-face / T-shaped).
STACK_CENTROID_MAX = 5.5
STACK_PARALLEL_MAX = 30.0
STACK_PERPENDICULAR_MIN = 60.0

#: v2 types, in the priority order :func:`classify_contacts` uses for the ``contact.type`` label.
TYPES_V2 = ("salt_bridge", "hydrogen_bond", "cation_pi", "stacking", "aromatic", "hydrophobic",
            "polar", "vdw", "other")
#: v1 types, frozen — the recognition models' ``ct_*`` features are counts over these.
TYPES_V1 = ("salt_bridge", "hydrogen_bond", "aromatic", "hydrophobic", "other")
TYPES = TYPES_V2


def _elem(atom: str) -> str:
    """Element from a PDB atom name (first alphabetic character)."""
    for ch in atom:
        if ch.isalpha():
            return ch
    return "?"


def _is_apolar(aa: str, atom: str) -> bool:
    """True for a carbon with no bonded N/O, or a sulfur (Met SD / Cys SG pack like carbon)."""
    e = _elem(atom)
    if e == "S":
        return True
    if e != "C":
        return False
    return atom not in _BACKBONE_POLAR_C and atom not in _POLAR_SIDECHAIN_C.get(aa, ())


def _donor_acceptor(aa: str, atom: str) -> tuple[bool, bool]:
    e = _elem(atom)
    if e == "O":
        return (aa, atom) in _HYDROXYL, True
    if e == "N":
        if aa == "P" and atom == "N":
            return False, False       # proline's ring nitrogen carries no hydrogen to donate
        return True, aa == "H" and atom in ("ND1", "NE2")
    return False, False


def _charged(aa: str, atom: str) -> tuple[bool, bool]:
    return (aa in _CATIONIC_RES and atom in _CATIONIC_ATOMS,
            aa in _ANIONIC_RES and atom in _ANIONIC_ATOMS)


def _types_v2(aa_a: str, aa_b: str, atom_a: str, atom_b: str, dist: float) -> set[str]:
    """Every type a single atom-atom contact satisfies (may be more than one, may be empty)."""
    out: set[str] = set()
    cat_a, ani_a = _charged(aa_a, atom_a)
    cat_b, ani_b = _charged(aa_b, atom_b)
    don_a, acc_a = _donor_acceptor(aa_a, atom_a)
    don_b, acc_b = _donor_acceptor(aa_b, atom_b)
    ring_a = aa_a in _AROMATIC_RES and atom_a in _RING_ATOMS
    ring_b = aa_b in _AROMATIC_RES and atom_b in _RING_ATOMS

    if dist <= _SALT_MAX and ((cat_a and ani_b) or (cat_b and ani_a)):
        out.add("salt_bridge")
    if dist <= _HBOND_MAX and ((don_a and acc_b) or (don_b and acc_a)):
        out.add("hydrogen_bond")
    if dist <= _CATION_PI_MAX and ((cat_a and ring_b) or (cat_b and ring_a)):
        out.add("cation_pi")
    if dist <= _AROMATIC_MAX and ring_a and ring_b:
        out.add("aromatic")
    if dist <= _HYDROPHOBIC_MAX and _is_apolar(aa_a, atom_a) and _is_apolar(aa_b, atom_b):
        out.add("hydrophobic")
    if dist <= _POLAR_MAX and not out and (_elem(atom_a) in "NO" or _elem(atom_b) in "NO"):
        out.add("polar")
    if dist <= _VDW_MAX and not out:
        # e.g. Ser CB (bonded to OG, so not apolar) against Leu CD1 (apolar): no polar atom, so not
        # `polar`, and not a hydrophobic pair either -- but at 3.7 A it is plainly a contact.
        out.add("vdw")
    return out


def _classify(aa_a: str, aa_b: str, atom_a: str, atom_b: str, dist: float) -> str:
    """The frozen v1 winner-takes-all classifier. Do not change: trained models depend on it."""
    ea, eb = _elem(atom_a), _elem(atom_b)
    cat = (aa_a in _CATIONIC_RES and atom_a in _CATIONIC_ATOMS,
           aa_b in _CATIONIC_RES and atom_b in _CATIONIC_ATOMS)
    ani = (aa_a in _ANIONIC_RES and atom_a in _ANIONIC_ATOMS,
           aa_b in _ANIONIC_RES and atom_b in _ANIONIC_ATOMS)
    if dist <= _SALT_MAX and ((cat[0] and ani[1]) or (cat[1] and ani[0])):
        return "salt_bridge"
    if dist <= 3.5 and ea in ("N", "O") and eb in ("N", "O"):
        return "hydrogen_bond"
    if (dist <= _AROMATIC_MAX and aa_a in _AROMATIC_RES and aa_b in _AROMATIC_RES
            and atom_a in _RING_ATOMS and atom_b in _RING_ATOMS):
        return "aromatic"
    if dist <= _HYDROPHOBIC_MAX and ea == "C" and eb == "C" and aa_a in _APOLAR_RES and aa_b in _APOLAR_RES:
        return "hydrophobic"
    return "other"


_PAIR_KEYS = ["chain.id.from", "residue.index.from", "chain.id.to", "residue.index.to"]


def stacked_pairs(structure) -> set[tuple]:
    """Residue pairs whose aromatic rings are arranged as a stack, keyed like a contact table.

    A contact potential scores a residue pair by identity, so it cannot tell two rings face to face
    at 3.5 Å from the same two residues brushing past edge-on. :func:`tcren.stacking.ring_stacking`
    measures that geometry; this turns it into the predicate :func:`classify_contacts` needs.

    Returns:
        A set of ``(chain.id.from, residue.index.from, chain.id.to, residue.index.to)`` tuples,
        canonically ordered to match :func:`~tcren.contacts.geometry.all_atom_contacts`.
    """
    from .stacking import ring_stacking

    df = ring_stacking(structure, cutoff=STACK_CENTROID_MAX)
    if df.height == 0:
        return set()
    df = df.filter(
        (pl.col("centroid_distance") <= STACK_CENTROID_MAX)
        & ((pl.col("interplanar_angle") <= STACK_PARALLEL_MAX)
           | (pl.col("interplanar_angle") >= STACK_PERPENDICULAR_MIN))
    )
    out = set()
    for row in df.iter_rows(named=True):
        a = (row["chain.id.from"], row["residue.index.from"])
        b = (row["chain.id.to"], row["residue.index.to"])
        lo, hi = (a, b) if a <= b else (b, a)          # ring_stacking orders by discovery, not by key
        out.add((lo[0], lo[1], hi[0], hi[1]))
    return out


def classify_contacts(interface_df: pl.DataFrame, scheme: str = "v2",
                      stacking: set[tuple] | None = None) -> pl.DataFrame:
    """Add chemical typing to an interface frame.

    Args:
        interface_df: an interface frame from :meth:`ContactMap.interface`, carrying
            ``residue.aa.from/to``, ``atom.from/to`` and ``dist``.
        scheme: ``"v2"`` (default) or the frozen ``"v1"`` — see the module docstring.
        stacking: residue pairs to mark as ``stacking``, from :func:`stacked_pairs`. Ring geometry
            needs coordinates, which a contact frame does not carry, so it is passed in.

    Returns:
        The frame with a ``contact.type`` column (the highest-priority label) and, for ``"v2"``, an
        ``is_<type>`` boolean per type — a contact can be both a salt bridge and a hydrogen bond,
        and collapsing that to one label loses a real interaction.

    Raises:
        ValueError: for an unknown ``scheme``.
    """
    if scheme not in ("v1", "v2"):
        raise ValueError(f"scheme must be 'v1' or 'v2', got {scheme!r}")
    types = TYPES_V1 if scheme == "v1" else TYPES_V2
    if interface_df.height == 0:
        out = interface_df.with_columns(pl.lit(None, dtype=pl.Utf8).alias("contact.type"))
        if scheme == "v2":
            out = out.with_columns([pl.lit(None, dtype=pl.Boolean).alias(f"is_{t}") for t in types])
        return out

    cols = [interface_df[c].to_list() for c in
            ("residue.aa.from", "residue.aa.to", "atom.from", "atom.to", "dist")]
    if scheme == "v1":
        return interface_df.with_columns(
            pl.Series("contact.type", [_classify(*row) for row in zip(*cols)]))

    stacking = stacking or set()
    keys = list(zip(*[interface_df[c].to_list() for c in _PAIR_KEYS])) if stacking else None
    per_row = []
    for n, row in enumerate(zip(*cols)):
        hit = _types_v2(*row)
        if keys is not None and keys[n] in stacking:
            hit.add("stacking")
        per_row.append(hit)

    label = [next((t for t in TYPES_V2 if t in hit), "other") for hit in per_row]
    return interface_df.with_columns(
        [pl.Series("contact.type", label)]
        + [pl.Series(f"is_{t}", [t in hit for hit in per_row]) for t in TYPES_V2]
    )


def _typed_atom_pairs(structure, interface: str, tcr_regions: str, cutoff: float) -> pl.DataFrame:
    """Every heavy-atom pair of one interface, typed under v2 with ring stacking joined."""
    from .contactmap import ContactMap

    cm = ContactMap.from_structure(structure, cutoff=cutoff, atom_pairs=True)
    return classify_contacts(cm.interface(interface, tcr_regions=tcr_regions), "v2",
                             stacked_pairs(structure))


def residue_pair_types(structure, interface: str = "tcr_peptide", tcr_regions: str = "all",
                       cutoff: float = 5.0) -> pl.DataFrame:
    """Per residue pair, the union of the types its atom pairs make.

    This is the form the collapsed contact table cannot give: a residue pair is credited with a salt
    bridge if *any* of its atom pairs makes one, not only if its closest pair happens to.

    Args:
        structure: a chain-typed, annotated :class:`~tcren.structure.model.Structure`.
        interface: interface name (``"tcr_peptide"``, ``"tcr_mhc"``, ``"peptide_mhc"``).
        tcr_regions: passed through to :meth:`ContactMap.interface`.
        cutoff: contact distance threshold (Å).

    Returns:
        One row per residue pair with ``dist`` (the closest atom pair's), ``contact.type`` and the
        ``is_<type>`` booleans.
    """
    typed = _typed_atom_pairs(structure, interface, tcr_regions, cutoff)
    if typed.height == 0:
        return typed
    agg = (typed.group_by(_PAIR_KEYS, maintain_order=True)
           .agg([pl.col("dist").min(),
                 pl.col("residue.aa.from").first(), pl.col("residue.aa.to").first()]
                + [pl.col(f"is_{t}").any() for t in TYPES_V2]))
    label = [next((t for t in TYPES_V2 if row[f"is_{t}"]), "other")
             for row in agg.iter_rows(named=True)]
    return agg.with_columns(pl.Series("contact.type", label))


def contact_type_counts(cm, interface: str = "tcr_peptide", tcr_regions: str = "all",
                        scheme: str = "v1", structure=None, cutoff: float = 5.0) -> dict[str, int]:
    """Per-type contact counts + distinct residue-pair counts for one interface.

    Args:
        cm: a :class:`~tcren.contactmap.ContactMap` (used by ``"v1"``; ``"v2"`` needs ``structure``).
        interface: interface name (``"tcr_peptide"``, ``"tcr_mhc"``, ``"peptide_mhc"``).
        tcr_regions: passed through to :meth:`ContactMap.interface`.
        scheme: defaults to the frozen ``"v1"`` here, unlike :func:`classify_contacts`, because
            these counts feed the trained recognition models. Pass ``"v2"`` for the current typing.
        structure: the source structure — required for ``"v2"``, which needs every atom pair and
            the ring geometry, neither of which a collapsed contact map carries.
        cutoff: contact distance threshold for the ``"v2"`` rebuild (Å).

    Returns:
        Mapping with ``n_<type>`` (contacts of each type) and ``pairs_<type>`` (distinct
        residue-residue pairs with >=1 contact of that type), e.g. ``pairs_hydrogen_bond`` is the
        documented ``n_hbond`` feature. Under ``"v2"`` a contact is counted under **every** type it
        satisfies, so the ``n_*`` values do not sum to the contact count.

    Raises:
        ValueError: for ``scheme="v2"`` without a ``structure``.
    """
    types = TYPES_V1 if scheme == "v1" else TYPES_V2
    out = {f"n_{t}": 0 for t in types} | {f"pairs_{t}": 0 for t in types}
    if scheme == "v1":
        df = classify_contacts(cm.interface(interface, tcr_regions=tcr_regions), "v1")
        if df.height == 0:
            return out
        for t in types:
            sub = df.filter(pl.col("contact.type") == t)
            out[f"n_{t}"] = sub.height
            if sub.height:
                out[f"pairs_{t}"] = sub.select(_PAIR_KEYS).unique().height
        return out

    if structure is None:
        raise ValueError("scheme='v2' needs the source structure (every atom pair + ring geometry)")
    typed = _typed_atom_pairs(structure, interface, tcr_regions, cutoff)
    if typed.height == 0:
        return out
    for t in types:
        sub = typed.filter(pl.col(f"is_{t}"))
        out[f"n_{t}"] = sub.height
        if sub.height:
            out[f"pairs_{t}"] = sub.select(_PAIR_KEYS).unique().height
    return out


#: Types that carry no chemistry beyond proximity — what :func:`type_weights` drops.
UNTYPED = ("vdw", "other")


def type_weights(typed: pl.DataFrame, drop: "tuple[str, ...]" = UNTYPED) -> "np.ndarray":
    """0/1 per-contact weights that keep only chemically-typed contacts.

    The review's fallback, and the cheap half of a type-aware potential: rather than re-derive the
    matrix conditioned on the contact type, use the type to *discard* pairs that are within 5 Å but
    make no interaction — a contact map built on proximity alone counts them the same as a salt
    bridge. Feed the result to ``score_peptides(..., weights=...)``.

    Args:
        typed: a frame carrying the ``is_<type>`` booleans, from :func:`classify_contacts` with
            ``scheme="v2"`` or from :func:`residue_pair_types`.
        drop: types to zero out. The default drops only the two that mean "nothing but proximity".

    Returns:
        One float (0.0 or 1.0) per row of ``typed``, in its row order.

    Raises:
        ValueError: if the frame carries no ``is_<type>`` columns (it was typed under ``"v1"``).
    """
    import numpy as np

    keep_types = [t for t in TYPES_V2 if t not in drop]
    have = [f"is_{t}" for t in keep_types if f"is_{t}" in typed.columns]
    if not have:
        raise ValueError("no is_<type> columns; classify with scheme='v2' or use residue_pair_types")
    if typed.height == 0:
        return np.zeros(0, dtype=np.float64)
    keep = np.zeros(typed.height, dtype=bool)
    for col in have:
        keep |= np.asarray(typed[col].to_list(), dtype=bool)
    return keep.astype(np.float64)
