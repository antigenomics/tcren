"""Provenance stamps for generated tables, and the check that a table is not stale.

The failure this exists to prevent: a feature or score table written by an older tcren is read back
months later, a number is quoted from it, and nothing anywhere says the descriptor that produced it
has since been renamed, redefined or removed. The number is then wrong and irreproducible, and the
only symptom is that it does not match a fresh run.

:func:`stamp` writes a sidecar JSON beside every table a command emits. :func:`check` reads it back
and **raises** unless the installed package would produce the same columns from the same registry.
The registry digest is the load-bearing part: it changes whenever any descriptor is added, dropped,
or has its units or definition edited, so a table computed under a different catalogue cannot pass
silently.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

__all__ = ["registry_digest", "sidecar_path", "stamp", "read", "check", "StaleTableError"]


class StaleTableError(RuntimeError):
    """A generated table does not match the installed descriptor catalogue."""


def registry_digest() -> str:
    """SHA-256 over the descriptor catalogue: names, families, invariance classes, units.

    Definitions are included, so editing what a descriptor *means* invalidates every table that
    carries it even when the number would not change. That is deliberate -- a redefined column is a
    different quantity under the same name, which is the case a value comparison cannot catch.
    """
    from .recognition import DESCRIPTORS, DETAIL, INVARIANCE

    payload = [
        [name, fam, bool(tcr), INVARIANCE.get(name, ""), *DETAIL.get(name, ("", ""))]
        for name, (fam, tcr) in sorted(DESCRIPTORS.items())
    ]
    blob = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(blob).hexdigest()


def sidecar_path(table: str | Path) -> Path:
    """Where the stamp for ``table`` lives: the table's path with ``.provenance.json`` appended."""
    return Path(str(table) + ".provenance.json")


def stamp(table: str | Path, *, command: str, columns=None, extra: dict | None = None) -> Path:
    """Write the sidecar for ``table`` and return its path.

    Args:
        table: the file just written.
        command: the invocation that produced it, verbatim enough to re-run.
        columns: the table's column names, or ``None`` to omit them.
        extra: anything else worth recording (structure count, potentials, cutoff).
    """
    from . import __version__

    out = sidecar_path(table)
    out.write_text(json.dumps({
        "tcren": __version__,
        "registry_digest": registry_digest(),
        "command": command,
        "columns": list(columns) if columns is not None else None,
        **(extra or {}),
    }, indent=1) + "\n")
    return out


def read(table: str | Path) -> dict | None:
    """The stamp for ``table``, or ``None`` if it has none."""
    p = sidecar_path(table)
    return json.loads(p.read_text()) if p.exists() else None


def check(table: str | Path, *, require: bool = True) -> dict:
    """Raise :class:`StaleTableError` unless ``table`` was written by this catalogue.

    Args:
        table: the table about to be read.
        require: raise when the stamp is missing entirely. ``True`` (default) is right for anything
            whose numbers will be reported -- an unstamped table is one written before stamping
            existed, which is exactly the era this guard is aimed at.

    Returns:
        The stamp.
    """
    from . import __version__

    meta = read(table)
    if meta is None:
        if require:
            raise StaleTableError(
                f"{table} has no provenance stamp ({sidecar_path(table).name} is missing). It was "
                f"written before this tcren, or by something that does not stamp. Recompute it "
                f"with the installed {__version__} rather than reading it.")
        return {}
    want = registry_digest()
    if meta.get("registry_digest") != want:
        raise StaleTableError(
            f"{table} was written under a different descriptor catalogue "
            f"(tcren {meta.get('tcren')}, digest {str(meta.get('registry_digest'))[:12]}; "
            f"installed {__version__}, digest {want[:12]}). Descriptors have been renamed, "
            f"redefined or removed since. Recompute it: {meta.get('command')!r}")
    return meta
