"""pMHC surface topology — the height and chemistry of the face a TCR actually sees.

A contact potential scores an interface that already exists. This module describes the pMHC
*before* a TCR arrives: how the presented surface is shaped and what it is made of. The peptide
sits in a groove between two helices, and a TCR approaching from above meets one surface — so the
descriptor is a **height field** ``h(x, y)`` over the groove plane, with per-cell chemistry painted
on. Two epitopes are then comparable as two rasters, which is what
:func:`surface_distance` exploits.

The method follows SURFMAP (Schweke et al., *J Chem Inf Model* 2022, 62:1595) — surface shell,
per-point feature, grid, 8-neighbour smoothing, Manhattan map distance — with one deliberate
departure. SURFMAP projects a globular protein onto an equal-area *spherical* chart because a
closed surface has no undistorted plane. The TCR-facing pMHC surface is an open, near-planar patch
sitting in a groove frame we can define from the coordinates, so a flat raster is both simpler and
undistorting. Protein Surface Topography (Berkut et al., *JBC* 2019) supplies the other idea taken
here: centre the chart on the functional site, so maps of different molecules are registered.

**The frame is refit from every structure, not inherited.** ``_groove_frame`` takes the SVD of the
MHC groove-floor Cα and signs the axes from the peptide and the TCR, so ``x`` = groove width,
``y`` = peptide N→C, ``z`` = toward the TCR, always. Maps are therefore comparable without
prealigning the inputs — SURFMAP's standing caveat — and without depending on whether the caller
ran :func:`tcren.docking.canonicalize_structure` first.

**What "featureless" means numerically.** :func:`surface_stats` reports ``relief`` (the height
spread over the peptide's own footprint), ``peak_to_valley`` and ``frac_above_ridge`` (how much
peptide surface clears the helix rims). A flat, MHC-dominated landscape — the "featureless" epitope
of Tynan et al. (*Nat Immunol* 2007) and Motozono et al. (*J Immunol* 2014) — scores low on all
three; a bulged epitope scores high.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import polars as pl

# Kyte & Doolittle (1982) J Mol Biol 157:105-132, Table 1 — the standard hydropathy scale.
KYTE_DOOLITTLE: dict[str, float] = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5, "G": -0.4,
    "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6, "S": -0.8,
    "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}

# Formal side-chain charge at pH 7. His is given a fraction, not a whole unit: its pKa sits near
# physiological pH, so it is neither reliably charged nor reliably neutral.
SIDE_CHAIN_CHARGE: dict[str, float] = {"D": -1.0, "E": -1.0, "K": 1.0, "R": 1.0, "H": 0.1}

# Bondi vdW radii, matching tcren.binder.features and tcren.clashes.
_BONDI = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "P": 1.80}
_DEFAULT_RADIUS = 1.70

_MHC_TYPES = ("MHCa", "MHCb", "MHC")
_TCR_TYPES = ("TRA", "TRB", "TRD", "TRG")
_HELIX_REGIONS = ("HELIX_A1", "HELIX_A2", "HELIX_B1")
_GROOVE_REGIONS = _HELIX_REGIONS + ("GROOVE_FLOOR",)
_V_REGIONS = ("FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "FR4")

#: Default map window in Å, (x0, x1, y0, y1), centred on the groove-floor centroid. Fixed rather
#: than fitted per structure, so every map shares one grid and cells correspond across epitopes.
#: Wide enough for a class-II 15-mer with its flanking overhangs.
DEFAULT_EXTENT = (-20.0, 20.0, -25.0, 25.0)
DEFAULT_GRID = (64, 32)          # (n_y, n_x) cells — SURFMAP's 72x36, sized to this window

CHANNELS = ("h", "phobic", "charge")

#: Percentile of MHC-helix cell heights taken as the groove rim in :func:`surface_stats`.
RIDGE_PERCENTILE = 90.0

#: Z cutoff for :func:`surface_complementarity` — the largest ``h_tcr − h_pmhc`` clearance, in Å,
#: at which a cell still counts as surface facing surface. Calibrated over 60 Native2026 crystals:
#: inside :data:`COMPARE_WINDOW` the cutoff reaches 0.895 of occupied pMHC cells at 4 Å, 0.951 at
#: 10 Å and 0.962 with no cutoff at all, so 10 Å sits where the curve has gone flat. One-sided on
#: purpose: the median gap is **−1.7 Å** and 71% of cells are interdigitated (the receptor's lowest
#: point in a cell lies below the groove's highest point in the same cell), because the two faces
#: interlock rather than stack.
MAX_GAP = 10.0

#: Half-widths ``(x, y)`` in Å of the window :func:`surface_complementarity` compares over.
#: :data:`DEFAULT_EXTENT` is sized for a class-II 15-mer with overhangs, which is much wider than
#: any receptor's footprint: over the full extent a TCR projection reaches only 0.741 of occupied
#: pMHC cells however large the Z cutoff, and the shortfall is nearly all at the far groove end
#: (coverage 0.348 beyond y = +15 Å against 0.987 near y = 0). Cropping to ±12 Å lifts coverage to
#: 0.951 without a Z cutoff doing the work. The peptide's own cells are covered at 0.917 even over
#: the full extent, so no part of the epitope surface is being discarded here.
COMPARE_WINDOW = (12.0, 12.0)


@dataclass(slots=True)
class SurfaceMap:
    """A gridded height + chemistry map of one pMHC's TCR-facing surface."""

    structure_id: str
    grid: tuple[int, int]                                   # (n_y, n_x)
    extent: tuple[float, float, float, float]               # x0, x1, y0, y1 in Å
    channels: dict[str, np.ndarray]                         # each (n_y, n_x), NaN where unoccupied
    source: np.ndarray                                      # (n_y, n_x) int codes, see SOURCE_CODES
    scale: str = "kd"
    n_atoms: int = 0
    peptide: str = ""
    side: str = "pmhc"                                      # which face was mapped

    def occupancy(self) -> float:
        """Fraction of grid cells that any surface point reached."""
        return float(np.isfinite(self.channels["h"]).mean())

    def to_frame(self) -> pl.DataFrame:
        """Long form: one row per occupied cell, with the cell centre in Å."""
        n_y, n_x = self.grid
        yy, xx = np.meshgrid(_centres(self.extent[2], self.extent[3], n_y),
                             _centres(self.extent[0], self.extent[1], n_x), indexing="ij")
        keep = np.isfinite(self.channels["h"]).ravel()
        data = {
            "structure.id": np.full(keep.sum(), self.structure_id),
            "iy": np.repeat(np.arange(n_y), n_x)[keep],
            "ix": np.tile(np.arange(n_x), n_y)[keep],
            "x": xx.ravel()[keep], "y": yy.ravel()[keep],
            "source": np.array([SOURCE_NAMES[c] for c in self.source.ravel()[keep]]),
        }
        data.update({ch: self.channels[ch].ravel()[keep] for ch in self.channels})
        return pl.DataFrame(data)


SOURCE_NAMES = ("none", "peptide", "mhc_helix_a1", "mhc_helix_a2", "mhc_helix_b1", "mhc_floor",
                "cdr1a", "cdr2a", "cdr3a", "cdr1b", "cdr2b", "cdr3b", "tcr_fr")
SOURCE_CODES = {name: i for i, name in enumerate(SOURCE_NAMES)}


def _centres(lo: float, hi: float, n: int) -> np.ndarray:
    edges = np.linspace(lo, hi, n + 1)
    return 0.5 * (edges[:-1] + edges[1:])


# =========================================================================================
# frame
# =========================================================================================
def _floor_ca(structure) -> np.ndarray:
    """Groove-floor Cα; falls back to all MHC Cα when the groove is not region-annotated."""
    pts = [r.ca for c in structure.chains if c.chain_type in _MHC_TYPES
           for reg in c.regions if reg.region_type == "GROOVE_FLOOR"
           for r in reg.residues if r.ca is not None]
    if len(pts) < 3:
        pts = [r.ca for c in structure.chains if c.chain_type in _MHC_TYPES
               for r in c.residues if r.ca is not None]
    return np.asarray(pts, float)


def _centroid_of(structure, chain_types) -> np.ndarray | None:
    pts = [r.ca for c in structure.chains if c.chain_type in chain_types
           for r in c.residues if r.ca is not None]
    return np.mean(pts, axis=0) if pts else None


def _peptide_ca(structure) -> np.ndarray:
    pep = next((c for c in structure.chains if c.chain_type == "PEPTIDE"), None)
    return np.asarray([r.ca for r in pep.residues if r.ca is not None], float) if pep else np.zeros((0, 3))


def _groove_frame(structure) -> tuple[np.ndarray, np.ndarray]:
    """Return ``(origin, basis)`` with basis rows x = groove width, y = peptide N→C, z = up.

    ``y`` comes from the **peptide**, not from the floor's principal axis. The class-I groove floor
    is a β-sheet whose strands run across the groove and whose GROOVE_FLOOR markup drags in the
    loops descending to α3, so its longest principal axis is not the groove axis — using it puts the
    two helices diagonally across the map instead of on either side of it. The peptide's N→C vector
    *is* the groove axis by construction. ``z`` is the floor plane's normal, which is well
    determined even when its in-plane axes are not, orthogonalised against ``y`` and signed toward
    the TCR (or, absent a TCR, toward the peptide).

    The origin puts the **peptide centroid at (0, 0)** while heights stay measured from the floor —
    Protein Surface Topography's "centre the chart on the functional site", which is what registers
    maps of different alleles and peptide lengths against each other.
    """
    floor = _floor_ca(structure)
    if len(floor) < 3:
        raise ValueError("need >=3 MHC groove Cα to fit the groove plane")
    floor_centre = floor.mean(axis=0)
    pep = _peptide_ca(structure)

    up_ref = _centroid_of(structure, _TCR_TYPES)
    if up_ref is None and len(pep):
        up_ref = pep.mean(axis=0)
    if up_ref is None:
        raise ValueError("need a TCR or peptide chain to orient the surface normal")

    _u, _s, vt = np.linalg.svd(floor - floor_centre, full_matrices=True)
    normal = vt[2]
    if np.dot(up_ref - floor_centre, normal) < 0:
        normal = -normal

    axis = pep[-1] - pep[0] if len(pep) >= 2 else vt[0]      # fall back to the floor's long axis
    y = axis - normal * np.dot(axis, normal)
    if np.linalg.norm(y) < 1e-6:                              # degenerate: peptide ~parallel to z
        y = vt[0] - normal * np.dot(vt[0], normal)
    y /= np.linalg.norm(y)
    x = np.cross(y, normal)

    basis = np.vstack([x, y, normal])
    origin = floor_centre
    if len(pep):                                              # slide x/y onto the peptide centre
        off = (pep.mean(axis=0) - floor_centre) @ basis.T
        origin = floor_centre + np.array([off[0], off[1], 0.0]) @ basis
    return origin, basis


# =========================================================================================
# surface points
# =========================================================================================
def _surface_atoms(structure) -> tuple[np.ndarray, np.ndarray, list, list]:
    """Heavy atoms of the pMHC groove face: coords, radii, per-atom residue, per-atom source code.

    The TCR is excluded on purpose — this is the surface presented *to* it. Buried MHC domains
    (α3, β2m, the class-II Ig domains) are excluded too when region markup is available, because
    they add atoms that can never be reached from above and only slow the SASA down.
    """
    coords, radii, residues, sources = [], [], [], []

    def _take(residue, code):
        for a in residue.atoms:
            if a.element == "H":
                continue
            coords.append(a.coord)
            radii.append(_BONDI.get(a.element, _DEFAULT_RADIUS))
            residues.append(residue)
            sources.append(code)

    for chain in structure.chains:
        if chain.chain_type == "PEPTIDE":
            for r in chain.residues:
                _take(r, SOURCE_CODES["peptide"])
        elif chain.chain_type in _MHC_TYPES:
            groove = [reg for reg in chain.regions if reg.region_type in _GROOVE_REGIONS]
            if groove:
                for reg in groove:
                    code = SOURCE_CODES.get(f"mhc_{reg.region_type.lower()}", SOURCE_CODES["mhc_floor"])
                    for r in reg.residues:
                        _take(r, code)
            else:
                for r in chain.residues:
                    _take(r, SOURCE_CODES["mhc_floor"])

    if not coords:
        raise ValueError("no peptide or MHC atoms found; is the structure chain-typed?")
    return np.asarray(coords, float), np.asarray(radii, float), residues, sources


def _tcr_atoms(structure) -> tuple[np.ndarray, np.ndarray, list, list]:
    """Heavy atoms of the TCR's **underside** — the face that descends onto the groove.

    The complement of :func:`_surface_atoms`: pMHC excluded, TCR kept. Only the V domain is taken
    when region markup is available (FR1-4 + CDR1-3), because the constant domains sit far above
    the groove and can never own a cell's lowest point; dropping them is a speed-up, not a change.
    """
    coords, radii, residues, sources = [], [], [], []
    loop_code = {("TRA", "CDR1"): "cdr1a", ("TRA", "CDR2"): "cdr2a", ("TRA", "CDR3"): "cdr3a",
                 ("TRB", "CDR1"): "cdr1b", ("TRB", "CDR2"): "cdr2b", ("TRB", "CDR3"): "cdr3b"}

    def _take(residue, code):
        for a in residue.atoms:
            if a.element == "H":
                continue
            coords.append(a.coord)
            radii.append(_BONDI.get(a.element, _DEFAULT_RADIUS))
            residues.append(residue)
            sources.append(code)

    for chain in structure.chains:
        if chain.chain_type not in _TCR_TYPES:
            continue
        vdom = [reg for reg in chain.regions if reg.region_type in _V_REGIONS]
        if vdom:
            for reg in vdom:
                code = SOURCE_CODES[loop_code.get((chain.chain_type, reg.region_type), "tcr_fr")]
                for r in reg.residues:
                    _take(r, code)
        else:
            for r in chain.residues:
                _take(r, SOURCE_CODES["tcr_fr"])

    if not coords:
        raise ValueError("no TCR atoms found; is the structure chain-typed?")
    return np.asarray(coords, float), np.asarray(radii, float), residues, sources


def _height_candidates(local: np.ndarray, radii: np.ndarray, grid, extent):
    """Every (cell, height, atom) a top-down ray cast produces.

    The height of the solvent-accessible surface directly above a point ``(x, y)`` is
    ``max_i z_i + sqrt(R_i² − r_lat²)`` over the atoms whose expanded sphere that column passes
    through. Casting rays in the groove frame — rather than sampling each atom's sphere and keeping
    the points that survive, as Shrake-Rupley does — makes the map **exactly equivariant**: the
    same structure rotated gives the same map, because nothing depends on how the molecule happens
    to sit relative to a fixed sampling sphere. (Measured on 1ao7, sphere sampling moved the median
    cell by 1.35 Å under a rigid rotation and shifted ``relief`` by 19%.)

    Accessibility needs no separate probe test: the highest surface in a column is, by definition,
    the one nothing else is above.
    """
    n_y, n_x = grid
    x0, x1, y0, y1 = extent
    dx, dy = (x1 - x0) / n_x, (y1 - y0) / n_y
    if not len(radii):
        return np.zeros(0, int), np.zeros(0, int), np.zeros(0)
    x, y, z = local[:, 0], local[:, 1], local[:, 2]
    r2 = radii ** 2

    ax = np.floor((x - x0) / dx).astype(int)
    ay = np.floor((y - y0) / dy).astype(int)
    half_x = int(np.ceil(radii.max() / dx))
    half_y = int(np.ceil(radii.max() / dy))

    cells, heights, owners = [], [], []
    atom = np.arange(len(x))
    for oy in range(-half_y, half_y + 1):
        ty = ay + oy
        cy = y0 + (ty + 0.5) * dy
        for ox in range(-half_x, half_x + 1):
            tx = ax + ox
            cx = x0 + (tx + 0.5) * dx
            lat2 = (cx - x) ** 2 + (cy - y) ** 2
            keep = (tx >= 0) & (tx < n_x) & (ty >= 0) & (ty < n_y) & (lat2 < r2)
            if not keep.any():
                continue
            cells.append(ty[keep] * n_x + tx[keep])
            heights.append(z[keep] + np.sqrt(r2[keep] - lat2[keep]))
            owners.append(atom[keep])
    if not cells:
        return (np.zeros(0, int),) * 2 + (np.zeros(0),)
    return np.concatenate(cells), np.concatenate(owners), np.concatenate(heights)


# =========================================================================================
# the map
# =========================================================================================
def _smooth(a: np.ndarray) -> np.ndarray:
    """Average each cell with its 8 neighbours, ignoring empty ones (SURFMAP step 5).

    Padded, not wrapped. The map is a flat window over one groove, not a torus, and both helices run
    the full length of it — ``np.roll`` here averaged the far end of the groove into the near one
    (90 of 1820 occupied cells on 1ao7, by up to 4.3 Å).
    """
    n_y, n_x = a.shape
    filled = np.pad(np.nan_to_num(a, nan=0.0), 1)
    mask = np.pad(np.isfinite(a).astype(float), 1)
    acc = sum(filled[dy:dy + n_y, dx:dx + n_x] for dy in (0, 1, 2) for dx in (0, 1, 2))
    cnt = sum(mask[dy:dy + n_y, dx:dx + n_x] for dy in (0, 1, 2) for dx in (0, 1, 2))
    out = np.where(cnt > 0, acc / np.maximum(cnt, 1), np.nan)
    return np.where(np.isfinite(a), out, np.nan)      # never invent a cell the surface never reached


def surface_map(
    structure,
    *,
    grid: tuple[int, int] = DEFAULT_GRID,
    extent: tuple[float, float, float, float] = DEFAULT_EXTENT,
    scale: str = "kd",
    probe: float = 1.4,
    smooth: bool = True,
    side: str = "pmhc",
) -> SurfaceMap:
    """Build the TCR-facing height + chemistry map of one pMHC, or the TCR face that meets it.

    Each cell keeps the **highest** point of the solvent-accessible surface above it — the first
    thing a TCR descending onto the groove would touch — and takes its chemistry from the atom that
    height belongs to. A mean over the column would average the exposed tip together with the flank
    behind it and blur exactly the relief the map exists to measure.

    Args:
        structure: chain-typed (and ideally MHC-annotated) TCR-pMHC or pMHC structure.
        grid: ``(n_y, n_x)`` cell counts.
        extent: ``(x0, x1, y0, y1)`` window in Å, registered on the peptide centroid.
        scale: hydropathy scale for the ``phobic`` channel — ``"kd"`` (Kyte-Doolittle) or
            ``"mj"`` (the hydrophobicity axis recovered from the MJ 1996 contact matrix by
            :meth:`tcren.potential.Potential.hydrophobicity_fit`).
        probe: solvent probe radius in Å, added to each atom's vdW radius.
        smooth: apply the 8-neighbour average to the numeric channels.
        side: ``"pmhc"`` (default) keeps the **highest** surface point per cell — the groove face a
            TCR descends onto. ``"tcr"`` keeps the **lowest** point of the TCR V domains in the same
            groove frame, i.e. the receptor's underside. Both maps carry the same grid, extent and
            frame, so they register cell-for-cell and can be compared by
            :func:`surface_complementarity`.

    Returns:
        A :class:`SurfaceMap`.

    Raises:
        ValueError: if the groove plane or the pMHC atoms cannot be located.
    """
    if side not in ("pmhc", "tcr"):
        raise ValueError(f"side must be 'pmhc' or 'tcr', got {side!r}")
    xyz, radii, residues, sources = (_surface_atoms if side == "pmhc" else _tcr_atoms)(structure)
    origin, basis = _groove_frame(structure)
    local = (xyz - origin) @ basis.T
    # The ray cast always keeps the topmost surface. Flipping z turns it into a cast from below,
    # which is what the TCR underside is; the heights are flipped back into the groove frame after.
    if side == "tcr":
        local = local * np.array([1.0, 1.0, -1.0])

    cell, owner, height = _height_candidates(local, radii + probe, grid, extent)
    n_y, n_x = grid

    # Highest surface per cell, without a Python loop: sort by (cell, height), take each run's last.
    order = np.lexsort((height, cell))
    cell_s, owner_s, height_s = cell[order], owner[order], height[order]
    last = np.r_[cell_s[1:] != cell_s[:-1], True] if cell_s.size else np.zeros(0, bool)
    top_cell, top_owner, top_z = cell_s[last], owner_s[last], height_s[last]

    phobic_by_aa = KYTE_DOOLITTLE if scale == "kd" else _mj_hydropathy()
    flat = {ch: np.full(n_y * n_x, np.nan) for ch in CHANNELS}
    src = np.zeros(n_y * n_x, dtype=np.int8)
    flat["h"][top_cell] = -top_z if side == "tcr" else top_z
    flat["phobic"][top_cell] = [phobic_by_aa.get(residues[o].aa, np.nan) for o in top_owner]
    flat["charge"][top_cell] = [SIDE_CHAIN_CHARGE.get(residues[o].aa, 0.0) for o in top_owner]
    src[top_cell] = [sources[o] for o in top_owner]

    channels = {ch: flat[ch].reshape(n_y, n_x) for ch in CHANNELS}
    if smooth:
        channels = {ch: _smooth(a) for ch, a in channels.items()}

    pep = next((c for c in structure.chains if c.chain_type == "PEPTIDE"), None)
    return SurfaceMap(
        structure_id=getattr(structure, "pdb_id", "") or "",
        grid=grid, extent=extent, channels=channels, source=src.reshape(n_y, n_x),
        scale=scale, n_atoms=len(radii), side=side,
        peptide="".join(r.aa for r in pep.residues) if pep else "",
    )


def _mj_hydropathy() -> dict[str, float]:
    """Per-residue hydrophobicity recovered from the MJ 1996 contact matrix itself.

    :meth:`Potential.hydrophobicity_fit` reads the hydrophobicity axis straight off the matrix
    (R² = 0.98 on ``mj1996``), so this is a structure-derived alternative to a tabulated scale.
    TCRen itself cannot be used: it is directed (TCR→peptide) and the fit refuses an asymmetric
    matrix.
    """
    from ..potential import mj1996

    fit = mj1996().hydrophobicity_fit()
    return {aa: float(fit.q[i]) for aa, i in fit.index.items()}


# =========================================================================================
# scalars
# =========================================================================================
def surface_stats(smap: SurfaceMap) -> dict[str, float]:
    """Reduce a map to the scalars that say how featured the epitope's surface is.

    Returns a dict with:

    ``relief``
        Standard deviation of height over the cells the peptide owns — the spread of the
        peptide's own topography. A flat epitope is small here.
    ``peak_to_valley``
        Max minus min height over the same cells.
    ``frac_above_ridge``
        Fraction of peptide cells that clear the MHC helix crest (:data:`RIDGE_PERCENTILE` of the
        helix cell heights). This is the one that separates a bulged epitope (much of it above the
        rims) from a featureless one (buried between them), and it needs no reference structure.
    ``phobic_mean`` / ``phobic_centre``
        Mean hydropathy over all peptide cells, and over the central third along the groove axis
        (the TCR-facing bulge, where the Chowell et al. 2015 immunogenicity signal sits).
    ``charge_mean``
        Mean formal charge over peptide cells.
    ``area_frac_peptide``
        Share of occupied cells the peptide owns rather than the MHC.
    """
    h = smap.channels["h"]
    pep = smap.source == SOURCE_CODES["peptide"]
    helix = np.isin(smap.source, [SOURCE_CODES[n] for n in ("mhc_helix_a1", "mhc_helix_a2", "mhc_helix_b1")])
    occupied = np.isfinite(h)
    pep_h = h[pep & occupied]
    out = {k: float("nan") for k in ("relief", "peak_to_valley", "frac_above_ridge", "phobic_mean",
                                     "phobic_centre", "charge_mean", "area_frac_peptide")}
    if pep_h.size == 0:
        return out

    # The ridge is the helix *crest*, taken as a high percentile rather than a mean. A helix's outer
    # flank slopes away from the groove and is still the top surface in its own cells, so its mean
    # height sits well below the rims the peptide actually has to clear (on 1ao7: mean ~12 Å against
    # a crest near 18 Å, which would put 85% of a flat 9-mer "above the ridge").
    helix_h = h[helix & occupied]
    ridge = float(np.percentile(helix_h, RIDGE_PERCENTILE)) if helix_h.size else float("nan")
    n_y = smap.grid[0]
    band = np.zeros_like(pep)
    band[n_y // 3: 2 * n_y // 3, :] = True

    out["relief"] = float(np.std(pep_h))
    out["peak_to_valley"] = float(pep_h.max() - pep_h.min())
    out["frac_above_ridge"] = float(np.mean(pep_h > ridge)) if np.isfinite(ridge) else float("nan")
    out["phobic_mean"] = float(np.nanmean(smap.channels["phobic"][pep & occupied]))
    centre = smap.channels["phobic"][pep & occupied & band]
    out["phobic_centre"] = float(np.nanmean(centre)) if centre.size else float("nan")
    out["charge_mean"] = float(np.nanmean(smap.channels["charge"][pep & occupied]))
    out["area_frac_peptide"] = float(pep_h.size / max(int(occupied.sum()), 1))
    return out


# =========================================================================================
# comparison
# =========================================================================================
def surface_distance(maps: list[SurfaceMap], channel: str = "h",
                     region: str | None = None) -> tuple[list[str], np.ndarray]:
    """Pairwise Manhattan distance between maps (SURFMAP eq. 1), normalised per shared cell.

    Cells occupied in only one of the two maps carry no comparison, so the sum runs over the
    intersection and is divided by its size. Without that, a map with fewer occupied cells would
    look closer to everything.

    Args:
        maps: maps sharing one ``grid`` and ``extent``.
        channel: which channel to compare (``"h"``, ``"phobic"``, ``"charge"``).
        region: restrict to cells of one source, e.g. ``"peptide"``; ``None`` uses every cell.

    Returns:
        ``(ids, distances)`` — the structure ids in row order and a square ``(n, n)`` matrix.

    Raises:
        ValueError: if the maps do not share a grid and extent.
    """
    if not maps:
        return [], np.zeros((0, 0))
    grid, extent, side = maps[0].grid, maps[0].extent, maps[0].side
    if any(m.grid != grid or m.extent != extent for m in maps):
        raise ValueError("all maps must share the same grid and extent to be comparable")
    if any(m.side != side for m in maps):
        raise ValueError("all maps must map the same face; use surface_complementarity across faces")

    stack = np.stack([m.channels[channel] for m in maps])
    if region is not None:
        code = SOURCE_CODES[region]
        stack = np.where(np.stack([m.source for m in maps]) == code, stack, np.nan)

    n = len(maps)
    d = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            both = np.isfinite(stack[i]) & np.isfinite(stack[j])
            d[i, j] = d[j, i] = (float(np.abs(stack[i][both] - stack[j][both]).sum() / both.sum())
                                 if both.any() else np.nan)
    return [m.structure_id for m in maps], d


def surface_complementarity(pmhc: SurfaceMap, tcr: SurfaceMap, *,
                            max_gap: float = MAX_GAP,
                            window: tuple[float, float] | None = COMPARE_WINDOW,
                            region: str | tuple[str, ...] | None = None,
                            tcr_region: str | tuple[str, ...] | None = None) -> dict[str, float]:
    """Cell-for-cell agreement between a pMHC face and the TCR underside that meets it.

    Both maps must come from the same structure (same groove frame, grid and extent), one built
    with ``side="pmhc"`` and one with ``side="tcr"``. A cell enters the comparison when both maps
    reach it **and** the vertical clearance between them is at most ``max_gap`` — the Z cutoff that
    keeps the footprint to surface actually facing surface, rather than the map corners where the
    receptor overhangs nothing. See :data:`MAX_GAP` for how the default was chosen.

    Returned keys, all over the retained cells:

    ``n_cells`` / ``coverage``
        Cells retained, and their share of the occupied pMHC cells inside ``window``.
    ``gap_mean`` / ``gap_sd``
        Mean and spread of ``h_tcr − h_pmhc`` in Å. A tight, even gap is a well-packed interface;
        a wide or ragged one is a receptor resting on a few high points.
    ``shape_r``
        Pearson *r* between the two height fields. **Positive** is complementary: where the groove
        rises the receptor must ride up over it.
    ``charge_r`` / ``charge_product``
        Pearson *r* between the two charge fields, and the mean of their product.
        **Negative** is complementary — plus meeting minus.
    ``phobic_r`` / ``phobic_product``
        The same for hydropathy. **Positive** is complementary — apolar meeting apolar.
    ``d_h`` / ``d_charge`` / ``d_phobic``
        Mean absolute difference per cell in each channel, the SURFMAP map distance of
        :func:`surface_distance` applied across the two faces instead of across two structures.

    Raises:
        ValueError: if the maps disagree on grid/extent, or are not one of each side.
    """
    if pmhc.grid != tcr.grid or pmhc.extent != tcr.extent:
        raise ValueError("maps must share a grid and extent to be compared cell-for-cell")
    if (pmhc.side, tcr.side) != ("pmhc", "tcr"):
        raise ValueError(f"expected (pmhc, tcr) sides, got ({pmhc.side}, {tcr.side})")

    hp, ht = pmhc.channels["h"], tcr.channels["h"]
    occupied = np.isfinite(hp)
    if window is not None:
        n_y, n_x = pmhc.grid
        x0, x1, y0, y1 = pmhc.extent
        yy = np.abs(_centres(y0, y1, n_y))[:, None] <= window[1]
        xx = np.abs(_centres(x0, x1, n_x))[None, :] <= window[0]
        occupied = occupied & yy & xx
    for smap, sel in ((pmhc, region), (tcr, tcr_region)):
        if sel is not None:
            names = (sel,) if isinstance(sel, str) else tuple(sel)
            occupied = occupied & np.isin(smap.source, [SOURCE_CODES[n] for n in names])
    keep = occupied & np.isfinite(ht) & ((ht - hp) <= max_gap)
    out = {k: float("nan") for k in ("gap_mean", "gap_sd", "shape_r", "charge_r", "charge_product",
                                     "phobic_r", "phobic_product", "d_h", "d_charge", "d_phobic")}
    out["n_cells"] = float(keep.sum())
    out["coverage"] = float(keep.sum() / max(int(occupied.sum()), 1))
    if keep.sum() < 3:
        return out

    gap = ht[keep] - hp[keep]
    out["gap_mean"], out["gap_sd"] = float(gap.mean()), float(gap.std())
    for name, a, b in (("shape", hp[keep], ht[keep]),
                       ("charge", pmhc.channels["charge"][keep], tcr.channels["charge"][keep]),
                       ("phobic", pmhc.channels["phobic"][keep], tcr.channels["phobic"][keep])):
        m = np.isfinite(a) & np.isfinite(b)
        # A constant field (an all-neutral charge patch) has no correlation to report, not a zero.
        if m.sum() >= 3 and a[m].std() > 0 and b[m].std() > 0:
            out[f"{name}_r"] = float(np.corrcoef(a[m], b[m])[0, 1])
        if name != "shape" and m.any():
            out[f"{name}_product"] = float(np.mean(a[m] * b[m]))
        out["d_h" if name == "shape" else f"d_{name}"] = (
            float(np.abs(a[m] - b[m]).mean()) if m.any() else float("nan"))
    return out


def surface_tree(maps: list[SurfaceMap], channel: str = "h", region: str | None = None,
                 method: str = "complete"):
    """Hierarchical clustering of maps by :func:`surface_distance` (SURFMAP's distance tree).

    Returns:
        ``(ids, linkage)`` — the ids in row order and a :func:`scipy.cluster.hierarchy.linkage`
        matrix, ready for ``dendrogram`` or ``fcluster``.
    """
    from scipy.cluster.hierarchy import linkage
    from scipy.spatial.distance import squareform

    ids, d = surface_distance(maps, channel=channel, region=region)
    if len(ids) < 2:
        return ids, np.zeros((0, 4))
    return ids, linkage(squareform(np.nan_to_num(d, nan=float(np.nanmax(d))), checks=False), method=method)


def surface_table(maps: list[SurfaceMap]) -> pl.DataFrame:
    """One row of :func:`surface_stats` per map, with the peptide and grid occupancy."""
    return pl.DataFrame([
        {"structure.id": m.structure_id, "peptide": m.peptide, "scale": m.scale,
         "n.atoms": m.n_atoms, "occupancy": round(m.occupancy(), 4), **surface_stats(m)}
        for m in maps
    ])
