"""Structure data model and parsing."""

from .io import (
    import_structure,
    is_structure_file,
    iter_structures,
    parse_structure,
    structure_id_from_path,
    structure_paths,
)
from .model import Atom, Chain, RegionMarkup, Residue, Structure

__all__ = [
    "parse_structure", "import_structure",
    "iter_structures", "structure_paths", "structure_id_from_path", "is_structure_file",
    "Structure", "Chain", "Residue", "Atom", "RegionMarkup",
]
