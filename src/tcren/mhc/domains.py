"""Canonical MHC groove region definitions.

Loads the bundled ``mhc_canonical.json``: for each ``"<class>|<role>"`` key it holds a
canonical mature chain sequence and the 0-based positions of each groove region
(``HELIX_A1``/``HELIX_A2`` for class I, ``HELIX_A1``/``HELIX_B1`` for class II, and
``GROOVE_FLOOR``). Region boundaries follow established mature-numbering ranges for the
α1/α2 (class I) and α1/β1 (class II) groove domains. Regions are projected onto query
chains in :mod:`tcren.mhc.regions`.
"""

from __future__ import annotations

import json
from functools import lru_cache
from importlib import resources


@lru_cache(maxsize=1)
def canonical_groove() -> dict:
    """Return the bundled canonical groove definitions."""
    text = resources.files("tcren.data").joinpath("mhc_canonical.json").read_text()
    return json.loads(text)


def groove_for(mhc_class: str, chain_role: str) -> dict | None:
    """Canonical groove definition for a ``(class, role)``, or ``None`` if none exists."""
    return canonical_groove().get(f"{mhc_class}|{chain_role}")
