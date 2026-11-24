#!/usr/bin/env python3
"""Pack the recorded model coefficients and development effects for local inference."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from tempfile import NamedTemporaryFile

from tcren.contacts.language import FORMAT, VERSION, load_contact_language_reference
import polars as pl


def _source(directory: Path, name: str) -> Path:
    paths = [p for suffix in (".tsv", ".tsv.gz") if (p := directory / (name + suffix)).is_file()]
    if len(paths) != 1:
        raise ValueError(f"expected one {name} TSV in {directory}")
    return paths[0]


def _mapping(path: Path, key: str, value: str) -> dict[str, float]:
    table = pl.read_csv(path, separator="\t").select(key, value)
    if table[key].n_unique() != table.height or table[key].null_count():
        raise ValueError(f"duplicate or missing {key} in {path}")
    return dict(table.iter_rows())


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ensemble", type=Path, required=True, help="verified effects and coefficient TSVs")
    parser.add_argument("--atom-effects", type=Path, required=True, help="directional cell/effect TSV")
    parser.add_argument("--residue-effects", type=Path, required=True, help="directional cell/effect TSV")
    parser.add_argument("-o", "--output", type=Path, required=True)
    args = parser.parse_args()

    names = sorted(f"_sv{pose}_b{bins}_{scope}" for pose in range(6)
                   for bins in range(1, 4) for scope in ("pooled", "split"))
    models = []
    for name in names:
        coefficients = _source(args.ensemble, f"l3_coefficients{name}")
        effect_path = _source(args.ensemble, f"l2_cell_effects{name}")
        table = pl.read_csv(effect_path, separator="\t")
        real = (pl.col("numerator_real") + .5) / (pl.col("denominator_real") + 1)
        mock = (pl.col("numerator_mock") + .5) / (pl.col("denominator_mock") + 1)
        complement = (1 - real).clip(lower_bound=1e-6).log() - (1 - mock).clip(lower_bound=1e-6).log()
        effects = table.select("interface", "block", "cell", "beta_shrunk", complement)
        models.append({
            "name": name,
            "coefficients": _mapping(coefficients, "feature", "coefficient"),
            "effects": effects.rows(),
            "coefficient_sha256": sha256(coefficients.read_bytes()).hexdigest(),
            "effects_sha256": sha256(effect_path.read_bytes()).hexdigest(),
        })
    packed = {
        "format": FORMAT, "version": VERSION, "models": models,
        "atom_effects": _mapping(args.atom_effects, "cell", "effect"),
        "residue_effects": _mapping(args.residue_effects, "cell", "effect"),
        "atom_effects_sha256": sha256(args.atom_effects.read_bytes()).hexdigest(),
        "residue_effects_sha256": sha256(args.residue_effects.read_bytes()).hexdigest(),
    }
    raw = (json.dumps(packed, separators=(",", ":"), allow_nan=False) + "\n").encode()
    with NamedTemporaryFile(dir=args.output.parent, delete=False) as stream:
        temp = Path(stream.name)
        stream.write(raw)
    try:
        reference = load_contact_language_reference(temp)
        temp.replace(args.output)
        print(reference.sha256)
    finally:
        temp.unlink(missing_ok=True)


if __name__ == "__main__":
    main()
