"""Interactive 3D peptide-pocket view with CDR1-3 overlay (py3Dmol / 3Dmol.js).

Renders the MHC groove (cartoon + optional translucent surface, histo.fyi style), the
peptide as sticks, and the TCR CDR1-3 loops as Cα traces coloured by the shared palette,
all in the canonical (oriented) frame. ``py3Dmol`` is imported lazily so the rest of the
package has no hard 3D dependency.
"""

from __future__ import annotations


from ..orient.align import align_to_native
from ..structure.model import Structure
from .palette import REGION_COLOR

_CDR_REGIONS = ("CDR1", "CDR2", "CDR3")
_TCR_TYPES = ("TRA", "TRB", "TRD", "TRG", "IGH", "IGK", "IGL")
_GROOVE = ("HELIX_A1", "HELIX_A2", "HELIX_B1", "GROOVE_FLOOR")


def _oriented_coords(structure: Structure, reference_id):
    """Return a function mapping a Cα coord to the canonical frame (identity on failure)."""
    try:
        result = align_to_native(structure, reference_id=reference_id)
        return lambda c: c @ result.rotation + result.translation
    except Exception:
        return lambda c: c


def _pdb_block(structure: Structure, transform) -> str:
    """Well-formed PDB text of the heavy atoms in the (oriented) frame for py3Dmol.

    One conformer per atom name per residue (drops duplicate altlocs) and hydrogens are
    omitted, so 3Dmol.js parses clean backbone connectivity for cartoon/surface drawing.
    """
    from ..structure.io import pdb_lines

    return "\n".join(pdb_lines(structure, transform=transform, keep_hydrogens=False))


def view_pocket_cdr(
    structure: Structure,
    reference_id: str | None = None,
    surface: bool = True,
    width: int = 700,
    height: int = 500,
):
    """Build a py3Dmol view of the groove + peptide + CDR loops (oriented frame).

    Args:
        structure: a chain-typed, MHC-annotated structure.
        reference_id: native reference for orientation (identity frame if unavailable).
        surface: draw a translucent groove surface.
        width, height: viewer size.

    Returns:
        A ``py3Dmol.view`` (call ``.show()`` in a notebook).
    """
    import py3Dmol

    transform = _oriented_coords(structure, reference_id)
    view = py3Dmol.view(width=width, height=height)
    view.addModel(_pdb_block(structure, transform), "pdb")
    view.setStyle({}, {"cartoon": {"color": "lightgray"}})

    mhc_chains = [c.chain_id for c in structure.chains if c.chain_type in ("MHCa", "MHCb")]
    pep_chains = [c.chain_id for c in structure.chains if c.chain_type == "PEPTIDE"]

    for cid in pep_chains:
        view.setStyle({"chain": cid}, {"stick": {"colorscheme": "yellowCarbon"}})
    if surface and mhc_chains:
        view.addSurface(
            py3Dmol.VDW, {"opacity": 0.45, "color": "lightgray"},
            {"chain": mhc_chains},
        )

    # CDR1-3 loops as colored Cα spheres+lines.
    for chain in structure.chains:
        if chain.chain_type not in _TCR_TYPES:
            continue
        for region in chain.regions:
            if region.region_type not in _CDR_REGIONS:
                continue
            color = REGION_COLOR.get(region.region_type.lower(), "magenta")
            resi = [r.pdb_index for r in region.residues]
            view.addStyle(
                {"chain": chain.chain_id, "resi": resi},
                {"cartoon": {"color": color}, "stick": {"color": color, "radius": 0.15}},
            )
    view.zoomTo()
    return view
