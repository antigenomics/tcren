"""2D complementarity-map projection, canonical tables, and pocket markers."""

from .frame import ProjectionResult, project_structure
from .pockets import pocket_markers
from .tables import (
    ca_contacts_table,
    classify_contact,
    contacts_table,
    region_pair_contacts,
    region_pair_summary,
    residue_markup_table,
)

__all__ = [
    "project_structure", "ProjectionResult",
    "residue_markup_table", "contacts_table", "ca_contacts_table", "classify_contact",
    "region_pair_contacts", "region_pair_summary",
    "pocket_markers",
]
