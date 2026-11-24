"""Sequence-verified model replay comparisons with explicit component fit scopes."""

from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np

from .md import component_relative_rmsd
from .structures import map_sequence_roles, read_chains


def compare_replay(
    reference,
    model,
    sequences,
    *,
    reference_terminal_missing=4,
    model_terminal_missing=4,
    reference_extended=True,
    model_extended=False,
    reference_spans=None,
):
    """Return Cα RMSDs and literal source/span identity, without biological claims.

    Requested roles A/B/C/D mean TRA/TRB/peptide/MHC. Atom correspondence requires
    unique exact sequence spans. Terminal omissions/extensions are declared; no
    internal substitution alignment, pose-based chain assignment or independent
    refit of the measured component is performed after selecting the fit frame.
    """
    roles = tuple("ABCD")
    if set(sequences) != set(roles) or any(not sequences[k] for k in roles):
        raise ValueError("Replay requires explicit nonempty TRA/TRB/peptide/MHC sequences")
    overrides = reference_spans or {}
    if not set(overrides) <= set(roles):
        raise ValueError("Unknown explicit reference role")
    ref, ref_identity = map_sequence_roles(
        reference,
        [r for r in roles if r not in overrides],
        sequences,
        max_terminal_missing=reference_terminal_missing,
        allow_extended=reference_extended,
    )
    candidate, model_identity = map_sequence_roles(
        model,
        roles,
        sequences,
        max_terminal_missing=model_terminal_missing,
        allow_extended=model_extended,
    )
    explicit = {}
    chains = (
        {cid: (seq, coords) for cid, seq, coords in read_chains(reference)} if overrides else {}
    )
    for role, segments in overrides.items():
        points = {}
        used = set()
        for segment in segments:
            cid = segment["output_chain"]
            start, output_start, length = (
                segment[k] for k in ("input_start", "output_start", "length")
            )
            if (
                any(type(x) is not int for x in (start, output_start, length))
                or min(start, output_start) < 0
                or length <= 0
            ):
                raise ValueError("Invalid explicit reference span")
            seq, coords = chains[cid]
            query = sequences[role][start : start + length]
            if len(query) != length or query != seq[output_start : output_start + length]:
                raise ValueError("Explicit reference span sequence differs")
            for i in range(length):
                if start + i in points:
                    raise ValueError("Overlapping explicit reference spans")
                points[start + i] = coords[output_start + i]
            used.add(cid)
            ref_identity.append(
                dict(
                    role=role,
                    output_chain=cid,
                    reference_start=start,
                    output_start=output_start,
                    matched_residues=length,
                    input_residues=len(sequences[role]),
                    explicit_span=True,
                )
            )
        if len(used) != 1 or not points:
            raise ValueError("A reference role requires one chain and nonempty spans")
        explicit[role] = points
    chain_roles = {}
    for item in ref_identity:
        previous = chain_roles.setdefault(item["output_chain"], item["role"])
        if previous != item["role"]:
            raise ValueError("Multiple roles use the same reference chain")
    arrays = [[], []]
    selections = {}
    spans = []
    offset = 0
    names = {"A": "TRA", "B": "TRB", "C": "peptide", "D": "MHC"}
    for role in roles:
        if role in explicit:
            reference_points = explicit[role]
        else:
            rc, rs, rn = ref[role]
            reference_points = {rs + i: rc[i] for i in range(rn)}
        mc, ms, mn = candidate[role]
        positions = sorted(set(reference_points) & set(range(ms, ms + mn)))
        if not positions:
            raise ValueError("No shared exact coordinate span: " + role)
        arrays[0].append(np.asarray([reference_points[i] for i in positions]))
        arrays[1].append(np.asarray([mc[i - ms] for i in positions]))
        selections[names[role]] = np.arange(offset, offset + len(positions), dtype=int)
        spans.append(
            dict(
                role=role,
                input_indices=positions,
                matched_ca=len(positions),
                input_residues=len(sequences[role]),
            )
        )
        offset += len(positions)
    selections["TCR"] = np.concatenate([selections["TRA"], selections["TRB"]])
    selections["all"] = np.arange(offset, dtype=int)
    rows = component_relative_rmsd(np.concatenate(arrays[0]), np.concatenate(arrays[1]), selections)
    rows = [r for r in rows if r["fit_component"] in {"all", "TCR", "MHC"}]
    return dict(
        rows=rows,
        spans=spans,
        reference_identity=ref_identity,
        model_identity=model_identity,
        atom_selection="CA",
        units="angstrom",
        reference_sha256=hashlib.sha256(Path(reference).read_bytes()).hexdigest(),
        model_sha256=hashlib.sha256(Path(model).read_bytes()).hexdigest(),
        scope="Coordinate agreement on exact shared spans; reference provenance determines accuracy scope",
    )


def inspect_replay(
    reference,
    model,
    sequences,
    *,
    alerts,
    contact_cutoff=4.5,
    minimum_interchain_distance_A=1.0,
    **mapping,
):
    """Apply declared diagnostic alerts without selecting or repairing a model.

    RMSD alerts are engineering triggers, not universal accuracy thresholds.
    Contact counts describe whole mapped output chains, including any declared
    extensions; RMSD coverage is separately reported at exact shared positions.
    """
    from .structures import geometry_and_contacts

    if not np.isfinite(contact_cutoff) or contact_cutoff <= 0:
        raise ValueError("Contact cutoff must be finite and positive")
    if (
        isinstance(minimum_interchain_distance_A, bool)
        or not np.isfinite(minimum_interchain_distance_A)
        or not 0 < minimum_interchain_distance_A <= 2
    ):
        raise ValueError("Extreme-clash distance must be positive and at most 2 A")
    result = compare_replay(reference, model, sequences, **mapping)
    observed = {
        r["fit_component"] + "/" + r["measured_component"]: r["rmsd_A"] for r in result["rows"]
    }
    if not alerts or any(
        k not in observed or isinstance(v, bool) or not np.isfinite(v) or v <= 0
        for k, v in alerts.items()
    ):
        raise ValueError(
            "Declare finite positive RMSD alert thresholds for known fit/measurement pairs"
        )
    diagnostics, contacts = geometry_and_contacts(
        model, result["model_identity"], cutoff=contact_cutoff
    )
    peptide_contacts = sum(
        r["residue_pairs"] for r in contacts if r["role_from"] in {"A", "B"} and r["role_to"] == "C"
    )
    findings = [
        dict(kind="rmsd_alert", scope=k, observed_A=observed[k], threshold_A=v)
        for k, v in alerts.items()
        if observed[k] > v
    ]
    shortest = diagnostics["shortest_interchain_pair_lt2A_A"]
    if shortest is not None and shortest < minimum_interchain_distance_A:
        findings.append(
            dict(
                kind="extreme_interchain_clash",
                observed_A=shortest,
                threshold_A=minimum_interchain_distance_A,
            )
        )
    if peptide_contacts == 0:
        findings.append(dict(kind="zero_tcr_peptide_contacts"))
    if diagnostics["missing_backbone_atoms"]:
        findings.append(
            dict(kind="missing_backbone_atoms", count=diagnostics["missing_backbone_atoms"])
        )
    return dict(
        comparison=result,
        coordinate_diagnostics=diagnostics,
        chain_pair_contacts=contacts,
        contact_cutoff_A=contact_cutoff,
        minimum_interchain_distance_A=minimum_interchain_distance_A,
        contact_scope="whole mapped output chains",
        tcr_peptide_residue_pairs=peptide_contacts,
        alert_thresholds_A=alerts,
        findings=findings,
        status="requires_review" if findings else "diagnostic_checks_passed",
    )


def compare_replay_contacts(contacts, chains, sequences, mappings, *, reference_id):
    """Compare peptide contacts on the exact query positions shared by every model.

    Mappings are the sequence-verified identity rows from compare_replay or
    map_sequence_roles, keyed by complex.id. No unresolved residues are imputed.
    Native substitutions must be excluded through declared exact mapping spans.
    """
    import polars as pl
    from .prediction_comparison import compare_peptide_contact_sets

    names = {"A": "TRA", "B": "TRB", "C": "PEPTIDE"}
    common = {role: set(range(len(sequences[role]))) for role in names}
    maps = {}
    expected_ids = {}
    source_sequences = {}
    chain_rows = []
    for identifier, identity in mappings.items():
        maps[identifier] = {}
        expected_ids[identifier] = {}
        source_sequences[identifier] = {}
        for role, kind in names.items():
            fragments = [r for r in identity if r["role"] == role]
            if not fragments or len({r["output_chain"] for r in fragments}) != 1:
                raise ValueError("Require one mapped chain per contact role")
            cid = fragments[0]["output_chain"]
            source = chains.filter(
                (pl.col("complex.id") == identifier) & (pl.col("chain.id") == cid)
            )
            if source.height != 1 or source["chain.type"][0] != kind:
                raise ValueError("Contact annotation differs from declared chain role")
            sequence = source["sequence"][0]
            source_sequences[identifier][kind] = sequence
            positions = {}
            query_positions = set()
            for fragment in fragments:
                i, j, n = (
                    fragment[k] for k in ("reference_start", "output_start", "matched_residues")
                )
                if (
                    any(type(v) is not int for v in (i, j, n))
                    or min(i, j) < 0
                    or n <= 0
                    or len(sequences[role][i : i + n]) != n
                    or sequences[role][i : i + n] != sequence[j : j + n]
                ):
                    raise ValueError("Contact mapping is not an exact sequence span")
                for offset in range(n):
                    if j + offset in positions or i + offset in query_positions:
                        raise ValueError("Overlapping contact mapping spans")
                    positions[j + offset] = i + offset
                    query_positions.add(i + offset)
            expected_ids[identifier][kind] = cid
            maps[identifier][kind] = positions
            common[role] &= set(positions.values())
            chain_rows.append(
                {"complex.id": identifier, "chain.type": kind, "sequence": sequences[role]}
            )
    if not mappings or reference_id not in mappings or any(not v for v in common.values()):
        raise ValueError("Require nonempty common TRA/TRB/peptide contact spans and reference")
    accepted = {names[k]: v for k, v in common.items()}
    projected = []
    excluded = {key: 0 for key in mappings}
    selected = contacts.filter(
        pl.col("chain.type.from").is_in(["TRA", "TRB"]) & (pl.col("chain.type.to") == "PEPTIDE")
    )
    for row in selected.iter_rows(named=True):
        identifier = row["complex.id"]
        if identifier not in maps:
            continue
        kind = row["chain.type.from"]
        if (
            row["chain.id.from"] != expected_ids[identifier][kind]
            or row["chain.id.to"] != expected_ids[identifier]["PEPTIDE"]
        ):
            raise ValueError("Contact row chain differs from declared role")
        for direction, role in [("from", kind), ("to", "PEPTIDE")]:
            index = row["residue.index." + direction]
            sequence = source_sequences[identifier][role]
            if (
                type(index) is not int
                or not 0 <= index < len(sequence)
                or row["residue.aa." + direction] != sequence[index]
            ):
                raise ValueError("Contact row differs from source sequence/index")
        i = maps[identifier][kind].get(row["residue.index.from"])
        j = maps[identifier]["PEPTIDE"].get(row["residue.index.to"])
        if i not in accepted[kind] or j not in accepted["PEPTIDE"]:
            excluded[identifier] += 1
            continue
        row.update({"residue.index.from": i, "residue.index.to": j})
        projected.append(row)
    table = pl.DataFrame(projected, schema=contacts.schema) if projected else contacts.head(0)
    result = compare_peptide_contact_sets(
        table, pl.DataFrame(chain_rows), reference_id=reference_id
    )
    return result, dict(
        common_query_indices={k: sorted(v) for k, v in common.items()},
        excluded_residue_pairs=excluded,
        scope="Exact common query positions across all supplied complexes",
    )


def replay_ensemble_consistency(models, sequences, *, threshold_A, confidence_selected_id=None):
    """Flag confidence/pose disagreement in a small producer ensemble, without reranking.

    Complete-linkage clusters bound every within-cluster MHC-fitted TCR Cα RMSD
    by the declared threshold. A dominant cluster is only an ensemble-agreement
    diagnostic: every model can agree on an inaccurate pose. No native structure
    or previous prediction is used to construct the clusters.
    """
    from itertools import combinations
    from scipy.cluster.hierarchy import fcluster, linkage
    from scipy.spatial.distance import squareform

    if isinstance(threshold_A, bool) or not np.isfinite(threshold_A) or threshold_A <= 0:
        raise ValueError("Declare a finite positive docking-consistency threshold")
    if not 2 <= len(models) <= 64 or any(not isinstance(k, str) or not k for k in models):
        raise ValueError("Require 2 to 64 uniquely identified producer models")
    if confidence_selected_id is not None and confidence_selected_id not in models:
        raise ValueError("Confidence-selected model is absent from the ensemble")
    identifiers = sorted(models)
    distances = np.zeros((len(identifiers), len(identifiers)))
    pairs = []
    spans = None
    hashes = {}
    for i, j in combinations(range(len(identifiers)), 2):
        left, right = identifiers[i], identifiers[j]
        result = compare_replay(models[left], models[right], sequences, reference_extended=False)
        if spans is None:
            spans = result["spans"]
        elif spans != result["spans"]:
            raise ValueError("Ensemble pairs require identical shared query-position coverage")
        value = next(
            r["rmsd_A"]
            for r in result["rows"]
            if r["fit_component"] == "MHC" and r["measured_component"] == "TCR"
        )
        distances[i, j] = distances[j, i] = value
        pairs.append(dict(model_from=left, model_to=right, rmsd_A=value))
        hashes[left] = result["reference_sha256"]
        hashes[right] = result["model_sha256"]
    labels = fcluster(
        linkage(squareform(distances), method="complete"), t=threshold_A, criterion="distance"
    )
    groups = sorted(
        (
            [identifiers[i] for i in range(len(identifiers)) if labels[i] == label]
            for label in set(labels)
        ),
        key=lambda group: (-len(group), group),
    )
    clusters = []
    for number, members in enumerate(groups, 1):
        indices = [identifiers.index(name) for name in members]
        clusters.append(
            dict(
                cluster=number,
                members=members,
                maximum_pair_rmsd_A=float(distances[np.ix_(indices, indices)].max()),
            )
        )
    dominant = clusters[0] if len(clusters) == 1 or len(groups[0]) > len(groups[1]) else None
    selected_in = (
        None
        if dominant is None or confidence_selected_id is None
        else confidence_selected_id in dominant["members"]
    )
    return dict(
        clusters=clusters,
        dominant_cluster=dominant["cluster"] if dominant else None,
        dominant_fraction=len(dominant["members"]) / len(models) if dominant else None,
        confidence_selected_id=confidence_selected_id,
        confidence_selected_in_dominant_cluster=selected_in,
        review_required=dominant is None or selected_in is False,
        pairs=pairs,
        spans=spans,
        source_sha256=hashes,
        threshold_A=float(threshold_A),
        atom_selection="CA",
        fit_component="MHC",
        measured_component="TCR",
        method="complete linkage",
        scope="Ensemble consistency only; no native accuracy, binding claim or automatic rank substitution",
    )
