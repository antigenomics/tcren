"""Preserve declared MHC context separately from sequence-search representatives."""
from __future__ import annotations

import warnings
import re
from pathlib import Path

from ..metadata import read_metadata
from . import reference
from .imgt import _HLA_LOCUS_ROLE, _read_fasta


def _normalise(value):
    return str(value).strip().removeprefix("HLA-")


def _compatible(declared, candidate):
    declared, candidate = _normalise(declared), _normalise(candidate)
    return candidate == declared or candidate.startswith(declared + ":")


def _role(allele):
    return _HLA_LOCUS_ROLE.get(_normalise(allele).split("*")[0])


def _validate(allele, reference_names):
    normal = _normalise(allele)
    human = re.fullmatch(r"[A-Z]+[0-9]*\*\d{2,3}(?::\d{2,3}){0,3}[NLSCAQ]?", normal)
    if not (human and _role(allele)) and allele not in reference_names:
        raise ValueError(f"Unsupported or malformed MHC declaration: {allele}")


def set_supplied_allele(call, allele, sequence, sequences=None):
    """Retain supplied resolution and report support, without changing chain typing.

    Exact observed-sequence support cannot establish distinctions outside a cropped
    structure. Absence from co-best hits is not proof of a genotype conflict.
    """
    allele = str(allele).strip()
    role = _role(allele)
    if role is not None and role != (call.mhc_class, call.chain_role):
        raise ValueError(f"Declared MHC {allele} conflicts with mapped class/chain role")
    if sequences is None:
        try:
            sequences = list(_read_fasta(reference.reference_fasta()))
        except FileNotFoundError:
            sequences = []
    if isinstance(sequences, dict):
        matching = sequences[_normalise(allele)]
    else:
        names = [(reference.parse_header(header)["allele"], seq) for header, seq in sequences]
        _validate(allele, {name for name, _ in names})
        matching = [seq for name, seq in names if _compatible(allele, name)]
    support = "not_supported_by_returned_cobest"
    if any(_compatible(allele, item) for item in call.allele_candidates):
        support = "compatible_returned_cobest"
    if matching:
        support = ("exact_observed_sequence" if sequence and any(sequence in seq for seq in matching)
                   else "observed_sequence_not_exact")
    call.representative_allele = call.representative_allele or call.allele
    call.allele = allele
    call.allele_source = "supplied"
    call.supplied_support = support
    if support in {"not_supported_by_returned_cobest", "observed_sequence_not_exact"}:
        warnings.warn(f"Declared MHC {allele}: {support}; inspect sequence evidence",
                      UserWarning, stacklevel=2)


def apply_metadata(structures, specs):
    """Apply unique id/complex.id keyed mhc or mhc.a/mhc.b metadata after mapping."""
    rows = {}
    for spec in ([specs] if isinstance(specs, (str, bytes, Path)) else specs):
        table = read_metadata(spec)
        if table is None:
            continue
        columns = [c for c in ("mhc", "mhc.a", "mhc.b") if c in table.columns]
        for row in table.select("id", *columns).iter_rows(named=True):
            values = {k: str(v).strip() for k, v in row.items() if k != "id" and v is not None and str(v).strip()}
            previous = rows.setdefault(row["id"], {})
            for key, value in values.items():
                if key in previous and previous[key] != value:
                    raise ValueError(f"Conflicting {key} metadata for structure {row['id']}")
                previous[key] = value
    if not any(rows.values()):
        return
    reference_rows = [(reference.parse_header(header)["allele"], seq)
                      for header, seq in _read_fasta(reference.reference_fasta())]
    declarations = {allele.strip() for values in rows.values()
                    for value in values.values() for allele in value.split('/')}
    reference_names = {name for name, _ in reference_rows}
    sequences = {}
    for allele in declarations:
        _validate(allele, reference_names)
        sequences[_normalise(allele)] = [seq for name, seq in reference_rows
                                       if _compatible(allele, name)]
    for structure in structures:
        declared = {}
        for column, value in rows.get(structure.pdb_id, {}).items():
            for allele in value.split("/"):
                allele = allele.strip()
                role_info = _role(allele)
                role = {"mhc.a": "MHCa", "mhc.b": "MHCb"}.get(column)
                if role_info and role and role_info[1] != role:
                    raise ValueError(f"Declared MHC {allele} conflicts with {column}")
                role = role or (role_info[1] if role_info else "MHCa")
                if role in declared and _normalise(declared[role]) != _normalise(allele):
                    raise ValueError(f"Conflicting declared MHC alleles for {role}")
                declared[role] = allele
        for role, allele in declared.items():
            calls = [c for c in structure.mhc_calls or [] if c.chain_role == role]
            if len(calls) != 1:
                raise ValueError(f"Declared MHC {allele} requires one mapped {role} chain")
            call = calls[0]
            chain = structure.chain(call.chain_id)
            set_supplied_allele(call, allele, chain.sequence(), sequences)
            chain.allele_info = call.allele
