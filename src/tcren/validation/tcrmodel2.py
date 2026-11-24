"""Retain residue confidence from locally generated TCRmodel2 outputs before cleanup."""

import hashlib
import json
from pathlib import Path
import pickle

import numpy as np

from .structures import parse_pdb_text


def export_generated_confidence(result_pickle, coordinate_file, output):
    """Export a *trusted, locally generated* result pickle to portable confidence JSON.

    Pickle is a producer format, not an interchange input for arbitrary downloaded files.
    Token order is verified against the producer's unrelaxed PDB. Chain sequences and 0-based
    residue positions permit exact remapping after TCRmodel2 changes chain IDs/author numbering.
    """
    result_pickle, coordinate_file, output = map(Path, (result_pickle, coordinate_file, output))
    # Copy both producer payloads into memory before decoding, so cleanup cannot
    # invalidate later path reads during confidence validation/serialisation.
    coordinate_bytes = coordinate_file.read_bytes()
    raw = result_pickle.read_bytes()
    result = pickle.loads(raw)
    coordinate_text = coordinate_bytes.decode()
    if not coordinate_text.rstrip().endswith("END"):
        raise ValueError("TCRmodel2 producer PDB is not complete yet")
    s = parse_pdb_text(coordinate_text, str(coordinate_file))
    chain_ids = [c.chain_id for c in s.chains for _ in c.residues]
    positions = [i for c in s.chains for i, _ in enumerate(c.residues)]
    n = len(positions)
    pae = np.asarray(result["predicted_aligned_error"], dtype=float)
    plddt = np.asarray(result["plddt"], dtype=float)
    if (
        pae.shape != (n, n)
        or plddt.shape != (n,)
        or not np.isfinite(pae).all()
        or not np.isfinite(plddt).all()
    ):
        raise ValueError("TCRmodel2 confidence dimensions do not match resolved protein residues")
    if (pae < 0).any() or (plddt < 0).any() or (plddt > 100).any():
        raise ValueError("TCRmodel2 PAE/pLDDT outside allowed range")
    ca_confidence = np.asarray(
        [
            float(line[60:66])
            for line in coordinate_text.splitlines()
            if line.startswith("ATOM") and line[12:16].strip() == "CA"
        ]
    )
    # The producer writes pLDDT into the PDB B column. Check every token against the
    # rounded Cα value to detect the wrong model or an unexpected residue order.
    if ca_confidence.shape != plddt.shape or not np.allclose(
        ca_confidence, plddt, rtol=0, atol=0.0051
    ):
        raise ValueError("TCRmodel2 coordinate token order/pLDDT differs from prediction result")
    payload = dict(
        generator="tcrmodel2",
        token_order_verified_by_ca_plddt=True,
        residue_indexing="sequence_index",
        token_chain_ids=chain_ids,
        token_res_ids=positions,
        chain_sequences={c.chain_id: c.sequence() for c in s.chains},
        pae=pae.tolist(),
        residue_plddts=plddt.tolist(),
        producer_result_sha256=hashlib.sha256(raw).hexdigest(),
        producer_coordinate_sha256=hashlib.sha256(coordinate_bytes).hexdigest(),
        producer_result=str(result_pickle),
        producer_coordinates=str(coordinate_file),
        scalars={
            key: float(result[key])
            for key in ("ranking_confidence", "iptm", "ptm", "num_recycles")
            if key in result
        },
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, separators=(",", ":")))
    temporary.replace(output)
    return payload


def native_ranking_confidence(row):
    """Read either literal TCRmodel2 score spelling and reject conflicting aliases."""
    values = [float(row[key]) for key in ("ranking_confidence", "model_confidence") if key in row]
    if not values:
        raise ValueError("Missing TCRmodel2 ranking confidence")
    if not np.isfinite(values).all() or len(set(values)) != 1:
        raise ValueError("Conflicting or nonfinite TCRmodel2 ranking confidence")
    return values[0]


def native_ranking_metric(row):
    """The upstream pMHC wrapper replaces its reported ranking score with pTM."""
    score = native_ranking_confidence(row)
    if "iptm" not in row and "ptm" in row and float(row["ptm"]) == score:
        return "ptm"
    return "ranking_confidence"


def ranked_detail_file(directory, statistics, rank):
    """Select retained confidence by the producer's exact ranking scalar, refusing ambiguity."""
    target = native_ranking_confidence(statistics[rank])
    metric = native_ranking_metric(statistics[rank])
    matches = [
        path
        for path in Path(directory).glob("*.json")
        if json.loads(path.read_text()).get("scalars", {}).get(metric) == target
    ]
    if len(matches) != 1:
        raise ValueError(
            f"{rank}: {len(matches)} retained confidence matches; require unique producer identity"
        )
    return matches[0]


def retain_generated_inputs(features, msa, output):
    """Retain trusted complete producer pickles before upstream cleanup.

    Only locally generated inputs are accepted by this internal workflow. Parsing
    verifies completed writes and matching token dimensions; the bytes are copied
    unchanged so later backend/recycle controls use the same producer features.
    """
    features, msa, output = map(Path, (features, msa, output))
    feature_bytes, msa_bytes = features.read_bytes(), msa.read_bytes()
    feature_data, msa_data = pickle.loads(feature_bytes), pickle.loads(msa_bytes)
    tokens = np.asarray(feature_data["aatype"]).shape[0]
    if (
        tokens <= 0
        or np.asarray(feature_data["template_aatype"]).shape != (4, tokens)
        or np.asarray(msa_data["msa"]).ndim != 2
        or np.asarray(msa_data["msa"]).shape[1] != tokens
    ):
        raise ValueError("Producer feature/MSA token dimensions differ")
    output.mkdir(parents=True, exist_ok=True)
    record = dict(
        tokens=tokens,
        features_source=str(features),
        msa_source=str(msa),
        features_sha256=hashlib.sha256(feature_bytes).hexdigest(),
        msa_sha256=hashlib.sha256(msa_bytes).hexdigest(),
        trusted_generated_inputs=True,
    )
    for name, data in [("producer-features.pkl", feature_bytes), ("producer-msa.pkl", msa_bytes)]:
        path = output / name
        if path.exists() and path.read_bytes() != data:
            raise ValueError("Refuse to replace different retained producer inputs")
        if not path.exists():
            temporary = path.with_suffix(path.suffix + ".tmp")
            temporary.write_bytes(data)
            temporary.replace(path)
    receipt = output / "producer-inputs.json"
    temporary = receipt.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(record, indent=2) + "\n")
    temporary.replace(receipt)
    return record


def compare_generated_feature_arrays(left, right):
    """Audit already-decoded trusted producer inputs without exposing sequence values.

    Pickle-byte hashes alone do not establish array identity across Python/numpy
    versions. Report missing keys, shape/dtype and exact element equality; numeric
    differences retain their native feature units. This does not infer why a
    prediction changed or deserialize arbitrary external pickles.
    """
    if any(not isinstance(key, str) for key in set(left) | set(right)):
        raise ValueError("Feature names must be strings")
    rows = []
    for key in sorted(set(left) | set(right)):
        if key not in left or key not in right:
            rows.append(
                dict(name=key, missing_from="left" if key not in left else "right", identical=False)
            )
            continue
        a, b = np.asarray(left[key]), np.asarray(right[key])
        same_shape = a.shape == b.shape
        equal = bool(same_shape and np.array_equal(a, b))
        maximum = None
        if same_shape and a.size and a.dtype.kind in "biuf" and b.dtype.kind in "biuf":
            af, bf = a.astype(float), b.astype(float)
            if np.isfinite(af).all() and np.isfinite(bf).all():
                maximum = float(np.max(np.abs(af - bf)))
        rows.append(
            dict(
                name=key,
                left_shape=list(a.shape),
                right_shape=list(b.shape),
                left_dtype=str(a.dtype),
                right_dtype=str(b.dtype),
                values_equal=equal,
                maximum_absolute_difference=maximum,
                identical=equal and a.dtype == b.dtype,
            )
        )
    return dict(
        identical=all(row["identical"] for row in rows),
        arrays=rows,
        scope="Exact producer input arrays; hashes/metadata and causal inference remain separate",
    )


def compare_msa_row_inventory(left, right, *, chain_ids):
    """Distinguish row ordering, sequence inventory and cross-chain pairing changes.

    Encoded MSAs are compared as multisets, including duplicate/padding rows.
    Equal per-chain inventories but different full rows isolate changed pairing
    or cross-chain ordering; this does not establish its upstream cause.
    """
    from collections import Counter

    arrays = [np.asarray(value) for value in (left, right)]
    ids = np.asarray(chain_ids)
    if any(a.ndim != 2 or a.dtype.kind not in "iu" for a in arrays):
        raise ValueError("Require two integer-encoded MSA matrices")
    if arrays[0].shape[1] != arrays[1].shape[1] or ids.shape != (arrays[0].shape[1],):
        raise ValueError("MSA columns and chain IDs must match")
    if (
        ids.dtype.kind not in "iuf"
        or not np.isfinite(ids).all()
        or not np.equal(ids, np.floor(ids)).all()
        or (ids < -(2**63)).any()
        or (ids >= 2**63).any()
    ):
        raise ValueError("Require finite integer-valued producer asymmetry IDs")
    ids = ids.astype(np.int64)

    def inventory(columns):
        tables = [a[:, columns].astype(np.int64) for a in arrays]
        counts = [Counter(row.tobytes() for row in table) for table in tables]
        return dict(
            order_equal=np.array_equal(*tables),
            multiset_equal=counts[0] == counts[1],
            left_rows=len(tables[0]),
            right_rows=len(tables[1]),
            left_unique_rows=len(counts[0]),
            right_unique_rows=len(counts[1]),
            common_rows_with_multiplicity=sum((counts[0] & counts[1]).values()),
        )

    return dict(
        full=inventory(np.ones(len(ids), dtype=bool)),
        chains={str(chain): inventory(ids == chain) for chain in np.unique(ids)},
        scope="Literal encoded rows including duplicates/padding; no sequence or pairing repair",
    )
