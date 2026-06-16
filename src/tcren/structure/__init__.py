"""Structure data model and parsing."""

from .io import import_structure, parse_structure
from .model import Atom, Chain, RegionMarkup, Residue, Structure

__all__ = [
    "parse_structure", "import_structure",
    "Structure", "Chain", "Residue", "Atom", "RegionMarkup",
]
