"""Compare distinct engaged atoms with the observed local atom inventory."""
from dataclasses import dataclass
from math import sqrt

import polars as pl

from ..clashes import _radius
from ..structure.model import Structure
from .geometry import representative_atom_contacts
from .observations import atom_inventory, contact_observations, select_conformer
from .table import directed_views

_RESIDUE = ["chain.id", "residue.index"]
_PAIR = [f"{key}.{side}" for side in ("from", "to") for key in _RESIDUE]
_INTERFACES = {
    "tcr_peptide": ("tcr", "peptide"),
    "tcr_mhc": ("tcr", "mhc"),
    "peptide_mhc": ("peptide", "mhc"),
}


@dataclass(frozen=True, slots=True)
class CoverageCounts:
    """Physical atom counts in biological interface order."""

    from_available: int
    to_available: int
    from_engaged: int
    to_engaged: int
    contact_pairs: int
    in_shell_pairs: int
    clashing_pairs: int
    from_entropy_bits: float | None = None
    to_entropy_bits: float | None = None
    from_information_bits: float | None = None
    to_information_bits: float | None = None

    @property
    def information_bits(self) -> float | None:
        """Mean directional divergence from shell partner availability."""
        if self.from_information_bits is None or self.to_information_bits is None:
            return None
        return (self.from_information_bits + self.to_information_bits) / 2

    @property
    def coverage(self) -> float | None:
        """Geometric mean of the two fractions, or None if opportunity is absent."""
        if not self.from_available or not self.to_available:
            return None
        return sqrt(self.from_engaged / self.from_available * self.to_engaged / self.to_available)


def _partner_information(opportunities: pl.DataFrame, edges: pl.DataFrame) -> tuple:
    if edges.is_empty():
        return (None, None, None, None)
    witnesses = edges.group_by(_PAIR).len()
    entropy, divergence = [], []
    for side, other in (("from", "to"), ("to", "from")):
        key = [f"{name}.{side}" for name in _RESIDUE]
        capacity = pl.col(f"n_atoms.{other}")
        reference = opportunities.with_columns((capacity / capacity.sum().over(key)).alias("q"))
        contact = witnesses.join(reference.select(*_PAIR, "q"), on=_PAIR, validate="1:1")
        p = pl.col("len") / pl.col("len").sum().over(key)
        w = pl.col("len") / edges.height
        h, kl = contact.select((-w * p.log(2)).sum().alias("entropy"),
            (w * (p / pl.col("q")).log(2)).sum().alias("divergence")).row(0)
        entropy.append(h)
        divergence.append(max(0.0, kl))
    return (*entropy, *divergence)


def reciprocal_coverage(
    structure: Structure, *, partner_information: bool = False
) -> dict[str, CoverageCounts]:
    """Measure engagement on the 12 Å Cβ/Cα shell with 5 Å heavy-atom contacts.

    Atom inventory and contacts use the same conformer and CDR/helix scope.
    Contacts outside the shell and overlaps above 0.4 Å are counted for audit.
    Only nonclashing atom pairs whose residue pair is in the shell engage atoms.
    The optional entropy and divergence use directed partner-residue contact
    counts, with a reference proportional to observed shell partner atom count.
    """
    inventory = atom_inventory(structure).filter(pl.col("contact_eligible"))
    residues = inventory.group_by(_RESIDUE).agg(
        pl.col("contact_role").first(), pl.len().alias("n_atoms")
    )
    shell = directed_views(representative_atom_contacts(
        select_conformer(structure), kind="cb", cutoff=12.0
    )).select(_PAIR)
    for side in ("from", "to"):
        metadata = residues.rename({key: f"{key}.{side}" for key in residues.columns})
        shell = shell.join(metadata, on=[f"{key}.{side}" for key in _RESIDUE], validate="m:1")

    contacts = contact_observations(structure)
    swapped = contacts.rename({
        name: name[:-4] + "to" if name.endswith(".from") else name[:-2] + "from"
        for name in contacts.columns if name.endswith((".from", ".to"))
    })
    contacts = pl.concat([contacts, swapped.select(contacts.columns)])
    radii = {element: _radius(element) for element in inventory["element"].unique()}
    overlap = (
        pl.col("element.from").replace_strict(radii, return_dtype=pl.Float64)
        + pl.col("element.to").replace_strict(radii, return_dtype=pl.Float64)
        - pl.col("dist")
    )
    result = {}
    for interface, (left, right) in _INTERFACES.items():
        opportunities = shell.filter(
            (pl.col("contact_role.from") == left) & (pl.col("contact_role.to") == right)
        )
        edges = contacts.filter((pl.col("interface") == interface) & (pl.col("role.from") == left))
        in_shell = edges.join(opportunities.select(_PAIR), on=_PAIR, how="semi")
        valid = in_shell.filter(overlap <= 0.4)
        available = []
        engaged = []
        for side in ("from", "to"):
            available.append(opportunities.unique(subset=[f"{key}.{side}" for key in _RESIDUE])[
                f"n_atoms.{side}"
            ].sum())
            engaged.append(valid[f"atom.source_index.{side}"].n_unique())
        information = _partner_information(opportunities, valid) if partner_information else ()
        result[interface] = CoverageCounts(
            *available, *engaged, edges.height, in_shell.height,
            in_shell.height - valid.height, *information
        )
    return result
