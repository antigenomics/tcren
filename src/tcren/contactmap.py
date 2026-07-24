"""Residue-level contact map and interface partitioning.

A :class:`ContactMap` wraps the annotated, symmetrised contact table and exposes the
three biological interfaces (TCR↔peptide, TCR↔MHC, peptide↔MHC). The TCR↔peptide
interface is the central object for scoring and reproduces the schema of
``data/contact_maps_PDB.csv`` once chains and regions are annotated.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal

import polars as pl

from .contacts.table import tidy_contacts
from .structure.model import MHC_TYPES, PEPTIDE_TYPE, RECEPTOR_TYPES, Structure

Interface = Literal["tcr_peptide", "tcr_mhc", "peptide_mhc"]

#: TCR region sets selectable on the ``from`` (TCR) side of an interface. ``"all"`` (no
#: filter) is the default and reproduces the legacy behaviour byte-for-byte; ``"cdr"`` keeps
#: only the three CDRs; ``"cdr+fr"`` adds the FR1–FR3 framework regions (FR4 excluded).
TCR_REGIONS: dict[str, set[str] | None] = {
    "cdr": {"CDR1", "CDR2", "CDR3"},
    "cdr+fr": {"CDR1", "CDR2", "CDR3", "FR1", "FR2", "FR3"},
    "all": None,
}


@dataclass(slots=True)
class ContactMap:
    """Annotated, symmetrised residue contacts for one structure."""

    pdb_id: str
    contacts: pl.DataFrame
    peptide_length: int | None = None
    # Per-(interface, tcr_regions) result cache; the table is immutable and the recognition
    # path re-requests the same interface many times per structure.
    _iface_cache: dict = field(default_factory=dict, init=False, repr=False, compare=False)

    @classmethod
    def from_structure(
        cls, structure: Structure, cutoff: float = 5.0, count_atoms: bool = False
    ) -> "ContactMap":
        """Build a contact map from an (annotated) structure.

        When ``count_atoms`` is set, the annotated table carries an ``n_atom_contacts``
        per-residue-pair heavy-atom count column (needed for atomic-weighted scoring).
        Default ``False`` keeps the contacts table byte-identical to the legacy output.
        """
        df = tidy_contacts(
            structure, cutoff=cutoff, count_atoms=count_atoms
        ).with_columns(pl.lit(structure.pdb_id).alias("pdb.id"))
        peptide_length = next(
            (len(c.residues) for c in structure.chains if c.chain_type == PEPTIDE_TYPE),
            None,
        )
        return cls(pdb_id=structure.pdb_id, contacts=df, peptide_length=peptide_length)

    def _interface(self, from_types: tuple[str, ...], to_types: tuple[str, ...]) -> pl.DataFrame:
        sel = self.contacts.filter(
            pl.col("chain.type.from").is_in(list(from_types))
            & pl.col("chain.type.to").is_in(list(to_types))
        )
        # pos = residue.index - region.start (0-based position within a region).
        return sel.with_columns(
            (pl.col("residue.index.from") - pl.col("region.start.from")).alias("pos.from"),
            (pl.col("residue.index.to") - pl.col("region.start.to")).alias("pos.to"),
        )

    def interface(self, which: Interface, tcr_regions: str = "all") -> pl.DataFrame:
        """Return the contacts of one interface with within-region positions.

        Args:
            which: ``"tcr_peptide"``, ``"tcr_mhc"`` or ``"peptide_mhc"``.
            tcr_regions: which TCR regions to keep on the ``from`` (TCR) side —
                ``"all"`` (default, no filter; legacy behaviour), ``"cdr"`` (CDR1–CDR3
                only), or ``"cdr+fr"`` (CDR1–CDR3 plus FR1–FR3). Has no effect on
                ``"peptide_mhc"`` (no TCR side).

        Returns:
            Filtered contacts with added ``pos.from``/``pos.to`` columns.
        """
        if tcr_regions not in TCR_REGIONS:
            raise ValueError(f"unknown tcr_regions {tcr_regions!r}")
        key = (which, tcr_regions)
        cached = self._iface_cache.get(key)
        if cached is not None:
            return cached
        if which == "tcr_peptide":
            sel = self._interface(RECEPTOR_TYPES, (PEPTIDE_TYPE,))
        elif which == "tcr_mhc":
            sel = self._interface(RECEPTOR_TYPES, MHC_TYPES)
        elif which == "peptide_mhc":
            sel = self._interface((PEPTIDE_TYPE,), MHC_TYPES)
            self._iface_cache[key] = sel
            return sel
        else:
            raise ValueError(f"unknown interface {which!r}")

        keep = TCR_REGIONS[tcr_regions]
        if keep is not None:  # TCR is on the 'from' side for these interfaces
            sel = sel.filter(pl.col("region.type.from").is_in(list(keep)))
        self._iface_cache[key] = sel
        return sel

    def tcr_peptide(self) -> pl.DataFrame:
        """Convenience accessor for the TCR↔peptide interface."""
        return self.interface("tcr_peptide")

    def to_csv(self, path: str | Path) -> None:
        """Write the full annotated contact table to CSV."""
        self.contacts.write_csv(str(path))


# --- Graphon featurisation: CDR3xpeptide maps in normalised (x, y) coordinates -----------------
# These reduce a variable-size CDR3xpeptide Ca map to a fixed-shape descriptor so loops of different
# length live in one vector space. They are structure->geometry featurisation, NOT binder scores:
# the discriminative signal collapses to chance under epitope matching (see the algem monograph E13).

from dataclasses import dataclass as _dataclass  # noqa: E402


@_dataclass(slots=True)
class ModeCentroid:
    """Length-invariant binding-mode centroid in graphon coordinates (see :func:`binding_mode`)."""

    apex_x: float           #: contact-weighted mean loop position, pooled over both loops (~0.47)
    y_alpha: float          #: contact-weighted mean peptide position read by CDR3alpha (~0.41)
    y_beta: float           #: contact-weighted mean peptide position read by CDR3beta (~0.62)
    sigma_sum: float        #: ``y_alpha + y_beta`` — the sigma involution reads ~1.0
    footprint_width_alpha: float
    footprint_width_beta: float
    n_contacts_alpha: int
    n_contacts_beta: int


def _cdr3_ca(structure, ctype: str):
    """Resolved CDR3 Cα array for one TCR chain type (``"TRA"``/``"TRB"``), or ``None``."""
    import numpy as np
    for c in structure.chains:
        if getattr(c, "chain_type", None) != ctype:
            continue
        for reg in getattr(c, "regions", []) or []:
            if reg.region_type == "CDR3":
                pts = [r.ca for r in reg.residues if r.ca is not None]
                return np.asarray(pts, float) if len(pts) >= 3 else None
    return None


def _cdr3_peptide_dmaps(structure):
    """``{"cdr3a": (L_a x M) dist, "cdr3b": (L_b x M) dist}`` Cα maps against the peptide, via markup.

    Uses the region markup (``classify_chains``) for CDR3 spans — cleaner than substring matching —
    and :func:`tcren.orient.docking._chain_ca` for the peptide. Missing loops are omitted.
    """
    import numpy as np

    from .orient.docking import _chain_ca
    pep = _chain_ca(structure, (PEPTIDE_TYPE,))
    out: dict[str, np.ndarray] = {}
    if pep is None or len(pep) < 2:
        return out
    for loop, ctype in (("cdr3a", "TRA"), ("cdr3b", "TRB")):
        cas = _cdr3_ca(structure, ctype)
        if cas is not None:
            out[loop] = np.linalg.norm(cas[:, None, :] - pep[None, :, :], axis=2)
    return out


def _register(dmap, grid: int):
    """Resample an (L x M) map onto a (grid x grid) grid in graphon coords ``(i/(L+1), j/(M+1))``."""
    import numpy as np
    from scipy.ndimage import map_coordinates
    u, v = np.meshgrid(np.linspace(0, dmap.shape[0] - 1, grid),
                       np.linspace(0, dmap.shape[1] - 1, grid), indexing="ij")
    return map_coordinates(dmap, [u.ravel(), v.ravel()], order=1, mode="nearest").reshape(grid, grid)


def registered_map(structure, *, grid: int = 8, target: str = "peptide",
                   metric: str = "distance", cutoff: float = 5.0):
    """CDR3×``target`` Cα map resampled onto a fixed ``grid``×``grid`` graphon grid.

    Bilinearly resamples the variable-length CDR3α/β × peptide Cα distance map onto normalised
    coordinates ``(i/(L+1), j/(M+1))`` so loops of different length are directly comparable, averageable
    and linear-model-ready. **Featurisation, not a binder score** — the discriminative signal is
    epitope-identity provenance and collapses to chance under epitope matching (algem monograph E13).

    Args:
        structure: a chain-typed structure (``classify_chains`` run).
        grid: output grid size ``G`` per axis.
        target: ``"peptide"`` (only target currently supported).
        metric: ``"distance"`` (Cα distances) or ``"contact"`` (binary at ``cutoff`` Å).
        cutoff: contact threshold in Å when ``metric="contact"``.

    Returns:
        ``(2, G, G)`` stacking the α then β blocks when both loops resolve; ``(G, G)`` when only one
        does; ``None`` when neither loop or the peptide is available.
    """
    import numpy as np
    if target != "peptide":
        raise ValueError("registered_map currently supports target='peptide' only")
    if metric not in ("distance", "contact"):
        raise ValueError("metric must be 'distance' or 'contact'")
    dmaps = _cdr3_peptide_dmaps(structure)
    blocks = []
    for loop in ("cdr3a", "cdr3b"):
        if loop in dmaps:
            reg = _register(dmaps[loop], grid)                             # resample the continuous distance
            blocks.append((reg <= cutoff).astype(float) if metric == "contact" else reg)  # then threshold
    if not blocks:
        return None
    return blocks[0] if len(blocks) == 1 else np.stack(blocks)


def binding_mode(structure, *, contact: float = 5.0) -> ModeCentroid | None:
    """Length-invariant binding-mode centroid of the CDR3α/β × peptide contact map.

    Reduces each loop's contact map to graphon-coordinate centroids: ``apex_x = i/(L+1)`` (where on the
    loop the contacts sit, ~0.47 = the apex) and ``y = j/(M+1)`` (where on the peptide each loop reads).
    The σ involution shows up as ``y_alpha + y_beta ≈ 1`` (α reads the N-terminal half, β the
    C-terminal). The two loops are **never pooled** for ``y`` — the split is the point — but ``apex_x``
    is contact-weighted over both.

    **Characterises the mode, does not discriminate specificity** — the gross mode is universal across
    HLA-A*02:01 9-mers, so do not use it as a same-epitope classifier (algem monograph E19).

    Args:
        structure: a chain-typed structure.
        contact: closest-Cα contact threshold in Å (the reference used an 8 Å Cα proxy; 5 Å is tighter).

    Returns:
        A :class:`ModeCentroid`, or ``None`` if neither loop makes ``>= 3`` contacts with the peptide.
    """
    import numpy as np
    dmaps = _cdr3_peptide_dmaps(structure)
    stats: dict[str, tuple[float, float, float, int]] = {}
    for loop in ("cdr3a", "cdr3b"):
        if loop not in dmaps:
            continue
        d = dmaps[loop]
        L, M = d.shape
        C = (d <= contact).astype(float)
        n = int(C.sum())
        if n < 3:
            continue
        x = (np.arange(1, L + 1) / (L + 1))[:, None]
        y = (np.arange(1, M + 1) / (M + 1))[None, :]
        xc = float((C * x).sum() / n)
        yc = float((C * y).sum() / n)
        yvals = np.repeat(y.ravel(), C.sum(axis=0).astype(int))
        width = float(np.sqrt(((yvals - yc) ** 2).mean())) if len(yvals) else float("nan")
        stats[loop] = (xc, yc, width, n)
    if not stats:
        return None
    a, b = stats.get("cdr3a"), stats.get("cdr3b")
    n_a, n_b = (a[3] if a else 0), (b[3] if b else 0)
    apex_x = float(np.average([s[0] for s in (a, b) if s], weights=[n for n in (n_a, n_b) if n]))
    return ModeCentroid(
        apex_x=apex_x,
        y_alpha=a[1] if a else float("nan"), y_beta=b[1] if b else float("nan"),
        sigma_sum=(a[1] + b[1]) if (a and b) else float("nan"),
        footprint_width_alpha=a[2] if a else float("nan"),
        footprint_width_beta=b[2] if b else float("nan"),
        n_contacts_alpha=n_a, n_contacts_beta=n_b)
