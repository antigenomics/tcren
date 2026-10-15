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


def ranked_detail_file(directory, statistics, rank):
    """Select retained confidence by the producer's exact ranking scalar, refusing ambiguity."""
    target = float(statistics[rank]["ranking_confidence"])
    matches = [
        path
        for path in Path(directory).glob("*.json")
        if json.loads(path.read_text()).get("scalars", {}).get("ranking_confidence") == target
    ]
    if len(matches) != 1:
        raise ValueError(
            f"{rank}: {len(matches)} retained confidence matches; require unique producer identity"
        )
    return matches[0]
