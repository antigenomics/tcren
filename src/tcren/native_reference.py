"""The versioned Native387 reference used by cohort ranking and contact statistics."""

from __future__ import annotations

from hashlib import sha256
from importlib.resources import files
import json
from pathlib import Path

FORMAT = "tcren-native-reference-1"
VERSION = "native387-v1"
SHA256 = "da552cff09e85ffd66cf524be4ff6b90f4524611537917c4dcd48395dad6f8c6"


def load_native_reference(
    path: str | Path | None = None, *, expected_sha256: str = SHA256,
) -> dict:
    """Load the frozen reference, checking its content identity before decoding.

    An explicit path can be used for a local copy. By default it must match the
    shipped bytes. Include the returned version and digest in derived tables.
    """
    source = files("tcren.data").joinpath("native387_v1.json") if path is None else Path(path)
    raw = source.read_bytes()
    digest = sha256(raw).hexdigest()
    if digest != expected_sha256:
        raise ValueError("native reference SHA-256 does not match the expected digest")
    reference = json.loads(raw)
    if reference.get("format") != FORMAT or reference.get("version") != VERSION:
        raise ValueError("unsupported native reference format or version")
    return {**reference, "reference_sha256": digest}


def contact_probability(reference: dict, alphabet: str, scope: str, cell: str) -> float | None:
    """Look up a cell on the frozen support. Unknown cells return ``None``."""
    return reference["contact_profiles"][f"{alphabet}/{scope}"]["probability"].get(cell)
