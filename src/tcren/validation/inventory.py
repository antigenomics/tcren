"""Stream private/local coordinate collections and retain sequence-match provenance."""

import gzip
from pathlib import Path
import tarfile

from ..structure.io import is_structure_file


def coordinate_stream(root):
    """Yield source identity and PDB text from folders or streaming tar/gzip archives."""
    root = Path(root)
    for path in [root] if root.is_file() else sorted(root.rglob("*")):
        if not path.is_file() or path.name.startswith("._"):
            continue
        if path.name.endswith(".pdb"):
            yield str(path), path.read_text()
        elif path.name.endswith(".pdb.gz"):
            with gzip.open(path, "rt") as handle:
                yield str(path), handle.read()
        elif path.name.endswith((".tgz", ".tar.gz")):
            with tarfile.open(path, "r|gz") as archive:
                for member in archive:
                    if (
                        member.isfile()
                        and member.name.endswith(".pdb")
                        and is_structure_file(member.name)
                    ):
                        yield f"{path}::{member.name}", archive.extractfile(member).read().decode()


def compatible_sequence_span(sequence, query, terminal_budget=0):
    """Unique exact shared span, allowing only a stated number of missing query termini.

    No substitutions or internal deletions; repeated/ambiguous spans return None.
    Ambiguity is rejected at the first matching truncation level, including repeated
    alternatives beside a unique span. Extra reference termini are allowed.
    """
    if not query or terminal_budget < 0 or terminal_budget >= len(query):
        raise ValueError("nonempty query and terminal budget smaller than query length required")
    for missing in range(terminal_budget + 1):
        matches = []
        for left in range(missing + 1):
            right = missing - left
            fragment = query[left : len(query) - right]
            start = sequence.find(fragment)
            if start >= 0:
                if sequence.find(fragment, start + 1) >= 0:
                    return None
                matches.append(dict(reference_start=start, input_start=left, length=len(fragment)))
        if matches:
            return matches[0] if len({tuple(sorted(row.items())) for row in matches}) == 1 else None
    return None
