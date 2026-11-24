"""Prediction-to-prediction residue contacts; the comparator is not experimental truth."""

import math

import polars as pl


def contact_cutoff_from_provenance(metadata):
    """Read current/legacy inspection cutoff names, rejecting absent or conflicting evidence."""
    names = ('cutoff_A', 'contact_cutoff_A')
    values = []
    for name in names:
        if name not in (metadata or {}):
            continue
        raw = metadata[name]
        if isinstance(raw, bool):
            raise ValueError('Contact cutoff must be finite and positive')
        try:
            value = float(raw)
        except (TypeError, ValueError) as error:
            raise ValueError('Contact cutoff must be finite and positive') from error
        if not math.isfinite(value) or value <= 0:
            raise ValueError('Contact cutoff must be finite and positive')
        values.append(value)
    if not values:
        raise ValueError('Contact cutoff is absent from inspection provenance')
    if len(set(values)) != 1:
        raise ValueError('Conflicting contact cutoff aliases in inspection provenance')
    return values[0]


def compare_peptide_contact_sets(contacts, chains, *, reference_id):
    """Compare tcren inspection tables over exact TRA/TRB/peptide sequences.

    Tables must come from the same inclusive heavy-atom cutoff. Contacts are the
    directed residue-pair export, not an atom-pair table. Models missing a required
    chain or carrying a different sequence stay explicitly unavailable. An empty
    contact set on a complete matching complex is valid and receives zero overlap
    against a nonempty reference. No confidence, native accuracy or affinity is
    inferred from Jaccard agreement with another predicted structure.
    """
    required = {"TRA", "TRB", "PEPTIDE"}
    chain_map = {}
    for row in chains.iter_rows(named=True):
        model = chain_map.setdefault(row["complex.id"], {})
        if row["chain.type"] in required:
            if row["chain.type"] in model:
                raise ValueError("Ambiguous repeated chain role: " + row["complex.id"])
            model[row["chain.type"]] = row["sequence"]
    reference = chain_map.get(reference_id)
    if reference is None or set(reference) != required:
        raise ValueError("Reference requires one TRA, TRB and peptide")
    selected = contacts.filter(
        pl.col("chain.type.from").is_in(["TRA", "TRB"]) & (pl.col("chain.type.to") == "PEPTIDE")
    )
    sets = {key: set() for key in chain_map}
    for row in selected.iter_rows(named=True):
        model = row["complex.id"]
        role = row["chain.type.from"]
        i, j = row["residue.index.from"], row["residue.index.to"]
        sequences = chain_map[model]
        if not (0 <= i < len(sequences[role]) and 0 <= j < len(sequences["PEPTIDE"])):
            raise ValueError("Contact index outside chain sequence")
        if (
            row["residue.aa.from"] != sequences[role][i]
            or row["residue.aa.to"] != sequences["PEPTIDE"][j]
        ):
            raise ValueError("Contact amino acids differ from chain sequences")
        key = (role, i, j)
        if key in sets[model]:
            raise ValueError(
                "Repeated directed residue contact; require collapsed residue-pair table"
            )
        sets[model].add(key)
    baseline = sets[reference_id]
    rows = []
    for model, sequences in chain_map.items():
        status = "comparable" if sequences == reference else "missing_or_different_sequence"
        candidate = sets[model]
        common = candidate & baseline
        union = candidate | baseline
        row = {
            "complex.id": model,
            "contact_reference": reference_id,
            "contact_comparison_status": status,
            "reference_residue_pairs": len(baseline),
            "candidate_residue_pairs": None,
            "shared_residue_pairs": None,
            "reference_only_pairs": None,
            "candidate_only_pairs": None,
            "contact_jaccard": None,
        }
        if status == "comparable":
            row.update(
                candidate_residue_pairs=len(candidate),
                shared_residue_pairs=len(common),
                reference_only_pairs=len(baseline - candidate),
                candidate_only_pairs=len(candidate - baseline),
                contact_jaccard=len(common) / len(union) if union else None,
            )
        rows.append(row)
    return pl.DataFrame(rows)
