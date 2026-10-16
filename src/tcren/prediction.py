"""Inspect generated structures using tcren annotation, descriptors and contacts.

Confidence is generator evidence, attached to coordinate observations without changing their
definition. A chain-pair ipTM is repeated on that chain pair's contacts, never interpreted as a
residue-contact probability. AF3 token-pair contact probability and directional PAE are separate.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
import hashlib
import json
from pathlib import Path

import numpy as np
import polars as pl


@dataclass
class PredictionConfidence:
    """Literal confidence values and their source; absent quantities stay absent."""

    source: str
    generator: str
    sha256: str = ""
    token_mapping_status: str = "unavailable"
    detail_sha256: str | None = None
    producer_structure_sha256: str | None = None
    global_values: dict = field(default_factory=dict)
    chain_pairs: dict = field(default_factory=dict)
    token_index: dict = field(default_factory=dict)
    token_sequences: dict = field(default_factory=dict)
    residue_plddts: dict = field(default_factory=dict)
    pae: np.ndarray | None = None
    contact_probs: np.ndarray | None = None


def confidence_chain_pairs(confidence: PredictionConfidence) -> pl.DataFrame:
    """Directional interface confidence, including pairs with no observed contacts.

    AF2/TCRmodel2 exports lacking pair matrices produce an empty typed table.
    Monomer diagonal entries are excluded; no pair confidence is inferred.
    """
    schema = {
        "chain.id.from": pl.String,
        "chain.id.to": pl.String,
        "chain_pair_iptm": pl.Float64,
        "chain_pair_pae_min": pl.Float64,
        "iptm": pl.Float64,
        "confidence_sha256": pl.String,
    }
    rows = [
        {
            "chain.id.from": a,
            "chain.id.to": b,
            "chain_pair_iptm": values.get("chain_pair_iptm"),
            "chain_pair_pae_min": values.get("chain_pair_pae_min"),
            "iptm": confidence.global_values.get("iptm"),
            "confidence_sha256": confidence.sha256,
        }
        for (a, b), values in sorted(confidence.chain_pairs.items())
        if a != b
    ]
    return pl.DataFrame(rows, schema=schema)


def _matrix(values, size: int, name: str, *, probability: bool = False):
    array = np.asarray(values, dtype=float)
    if array.shape != (size, size) or not np.isfinite(array).all():
        raise ValueError(f"{name}: expected finite {size} x {size} matrix")
    if (array < 0).any() or (probability and (array > 1).any()):
        raise ValueError(f"{name}: values outside allowed range")
    return array


def read_prediction_confidence(
    path: str | Path,
    *,
    generator: str,
    model_key: str | None = None,
    detail_path: str | Path | None = None,
    producer_structure_path: str | Path | None = None,
) -> PredictionConfidence:
    """Read native AF3, TCRmodel2, Boltz-2 or retained tFold-TCR/AF2-Multimer confidence.

    TCRmodel2 requires its literal ranked key (e.g. ``ranked_0``). Only scalar quantities are
    available from statistics after upstream cleanup; no chain-pair confidence is inferred.
    AF3 matrices retain their order and direction. Token identities are chain + author residue
    number, with a unique token required per standard protein residue.
    Boltz-2 requires the processed structure NPZ and uses paired writer PAE/pLDDT
    NPZ files for detail; modified/masked entities are refused.
    """
    path = Path(path)
    raw = json.loads(path.read_text())
    result = PredictionConfidence(
        str(path.resolve()), generator, sha256=hashlib.sha256(path.read_bytes()).hexdigest()
    )
    if generator in {"tfold-tcr", "af2-multimer"}:
        if detail_path is not None:
            raise ValueError(f"{generator} confidence is one retained JSON bundle")
        return _read_retained_confidence(result, raw)
    if generator == "tcrmodel2":
        if model_key is None or model_key not in raw:
            raise ValueError(f"{path}: missing TCRmodel2 ranked key {model_key!r}")
        raw = raw[model_key]
    elif generator not in {"af3", "boltz2", "protenix", "rf3"}:
        raise ValueError("generator must be af3, tcrmodel2, boltz2, tfold-tcr, af2-multimer, protenix or rf3")
    # TCRmodel2 statistics serializes some scalar values as strings; preserve the metric
    # names while accepting those numeric strings. Arrays / descriptive fields stay separate.
    for key, value in raw.items():
        if isinstance(value, (int, float, str)):
            try:
                result.global_values[key] = float(value)
            except ValueError:
                if key in (
                    "iptm",
                    "ptm",
                    "plddt",
                    "model_confidence",
                    "ranking_confidence",
                    "IpLDDT",
                ):
                    raise ValueError(f"{path}: nonnumeric {key}") from None
    if not all(np.isfinite(v) for v in result.global_values.values()):
        raise ValueError(f"{path}: nonfinite confidence scalar")
    for name in ("iptm", "ptm"):
        if name in result.global_values and not 0 <= result.global_values[name] <= 1:
            raise ValueError(f"{path}: {name} outside [0, 1]")
    if generator in {"protenix", "rf3"}:
        # Native scalar names verified against Protenix sample_confidence.py and
        # RF3 RF3Output/compile_af3_style_confidence_outputs. Atom/token identity
        # needs a version-pinned producer fixture before residue joins are safe.
        if detail_path is not None or producer_structure_path is not None:
            raise ValueError(f"{generator}: token mapping is unverified; only native summary scalars are supported")
        result.token_mapping_status = "unavailable_unverified_token_mapping"
        return result
    if generator == "boltz2":
        if producer_structure_path is None:
            raise ValueError("Boltz-2 requires explicit processed producer_structure NPZ")
        return _read_boltz2_confidence(result, raw, producer_structure_path, detail_path)
    if generator == "af3":
        chains = raw.get("chain_ids", [])
        if len(set(chains)) != len(chains):
            raise ValueError("duplicate AF3 confidence chain IDs")
        for name in ("chain_pair_iptm", "chain_pair_pae_min"):
            if name not in raw:
                continue
            matrix = _matrix(raw[name], len(chains), name, probability=name == "chain_pair_iptm")
            for i, a in enumerate(chains):
                for j, b in enumerate(chains):
                    result.chain_pairs.setdefault((a, b), {})[name] = float(matrix[i, j])
    if detail_path:
        result.detail_sha256 = hashlib.sha256(Path(detail_path).read_bytes()).hexdigest()
        detail = json.loads(Path(detail_path).read_text())
        ids = list(zip(detail["token_chain_ids"], detail["token_res_ids"]))
        if len(detail["token_chain_ids"]) != len(detail["token_res_ids"]):
            raise ValueError("token chain/residue arrays differ in length")
        if len(set(ids)) != len(ids):
            raise ValueError("multiple tokens per residue: standard protein tokens required")
        result.token_index = {key: i for i, key in enumerate(ids)}
        result.token_mapping_status = "producer_identity_available"
        if generator == "tcrmodel2":
            if (
                detail.get("generator") != generator
                or detail.get("residue_indexing") != "sequence_index"
            ):
                raise ValueError(
                    "TCRmodel2 detail requires retained sequence-index confidence schema"
                )
            result.token_sequences = detail["chain_sequences"]
            expected = {
                (cid, i) for cid, seq in result.token_sequences.items() for i in range(len(seq))
            }
            if set(ids) != expected:
                raise ValueError("TCRmodel2 token identities differ from retained chain sequences")
            if (
                "ranking_confidence" in result.global_values
                and detail["scalars"]["ranking_confidence"]
                != result.global_values["ranking_confidence"]
            ):
                raise ValueError("TCRmodel2 detail ranking confidence differs from requested rank")
            plddts = np.asarray(detail["residue_plddts"], dtype=float)
            if (
                plddts.shape != (len(ids),)
                or not np.isfinite(plddts).all()
                or (plddts < 0).any()
                or (plddts > 100).any()
            ):
                raise ValueError("TCRmodel2 residue pLDDT dimensions/range invalid")
            result.residue_plddts = dict(zip(ids, plddts.tolist()))
        if "pae" in detail:
            result.pae = _matrix(detail["pae"], len(ids), "pae")
        if "contact_probs" in detail:
            result.contact_probs = _matrix(
                detail["contact_probs"], len(ids), "contact_probs", probability=True
            )
    return result


def _read_retained_confidence(result, raw):
    """Retained AF2/tFold scalar/token evidence with explicit producer units."""
    if (
        raw.get("generator") != result.generator
        or raw.get("residue_indexing") != "sequence_index"
        or raw.get("plddt_units")
        != ("score_0_100" if result.generator == "af2-multimer" else "probability_0_1")
    ):
        raise ValueError(f"Expected retained {result.generator} sequence-index confidence schema")
    values = {k: float(v) for k, v in raw["scalars"].items()}
    if not all(np.isfinite(v) for v in values.values()) or any(
        not 0 <= values[k] <= 1 for k in ("iptm", "ptm") if k in values
    ):
        raise ValueError(f"Invalid {result.generator} confidence scalar")
    sequences = raw["chain_sequences"]
    if not sequences or any(
        not seq or any(c not in "ACDEFGHIKLMNPQRSTVWY" for c in seq) for seq in sequences.values()
    ):
        raise ValueError(f"{result.generator} requires explicit canonical protein sequences")
    ids = list(zip(raw["token_chain_ids"], raw["token_res_ids"]))
    expected = {(cid, i) for cid, seq in sequences.items() for i in range(len(seq))}
    maximum = 100 if result.generator == "af2-multimer" else 1
    if (
        len(raw["token_chain_ids"]) != len(raw["token_res_ids"])
        or len(set(ids)) != len(ids)
        or set(ids) != expected
    ):
        raise ValueError(f"{result.generator} token identities differ from retained sequences")
    plddt = np.asarray(raw["residue_plddts"], dtype=float)
    if (
        plddt.shape != (len(ids),)
        or not np.isfinite(plddt).all()
        or ((plddt < 0) | (plddt > maximum)).any()
    ):
        raise ValueError(f"{result.generator} pLDDT dimensions/range invalid")
    result.global_values = values
    result.token_sequences = sequences
    result.token_index = {key: i for i, key in enumerate(ids)}
    result.token_mapping_status = "producer_identity_available"
    result.residue_plddts = dict(zip(ids, (100 / maximum * plddt).tolist()))
    if result.generator == "af2-multimer" and "pae" in raw:
        result.pae = _matrix(raw["pae"], len(ids), "pae")
    result.detail_sha256 = result.sha256
    return result


def _read_boltz2_confidence(result, raw, producer_structure_path, detail_path):
    """Canonical all-protein Boltz outputs, mapped through the producer NPZ.

    Reject masked/modified entities instead of guessing token order. The Boltz
    writer orders tokens by chain table then residue span; one token per standard
    protein residue. NPZ inputs are loaded without Python pickle support.
    """
    from Bio.Data.IUPACData import protein_letters_3to1

    producer = Path(producer_structure_path)
    result.producer_structure_sha256 = hashlib.sha256(producer.read_bytes()).hexdigest()
    with np.load(producer, allow_pickle=False) as archive:
        chains, residues, mask = archive["chains"], archive["residues"], archive["mask"]
    if len(mask) != len(chains) or not mask.all():
        raise ValueError("Boltz-2 masked chains are unsupported for exact token mapping")
    chain_map, ids = {}, []
    covered = []
    for i, chain in enumerate(chains):
        cid = str(chain["name"])
        if int(chain["asym_id"]) != i or cid in chain_map.values() or int(chain["mol_type"]) != 0:
            raise ValueError("Boltz-2 requires unique, ordered protein chains")
        chain_map[str(i)] = cid
        start, length = int(chain["res_idx"]), int(chain["res_num"])
        span = residues[start : start + length]
        if start < 0 or length < 1 or len(span) != length or not span["is_standard"].all():
            raise ValueError("Boltz-2 requires canonical protein residue spans")
        if list(span["res_idx"]) != list(range(length)):
            raise ValueError("Boltz-2 residue indices do not match token order")
        try:
            sequence = "".join(protein_letters_3to1[str(r["name"]).title()] for r in span)
        except KeyError:
            raise ValueError("Boltz-2 noncanonical residue name") from None
        result.token_sequences[cid] = sequence
        ids.extend((cid, j) for j in range(length))
        covered.extend(range(start, start + length))
    if sorted(covered) != list(range(len(residues))):
        raise ValueError("Boltz-2 residue spans overlap or omit producer residues")
    pairs = raw.get("pair_chains_iptm", {})
    if set(pairs) != set(chain_map) or any(set(v) != set(chain_map) for v in pairs.values()):
        raise ValueError("Boltz-2 chain-pair matrix differs from producer chains")
    for a, values in pairs.items():
        for b, value in values.items():
            value = float(value)
            if not np.isfinite(value) or not 0 <= value <= 1:
                raise ValueError("Boltz-2 chain-pair ipTM outside [0, 1]")
            result.chain_pairs[(chain_map[a], chain_map[b])] = {"chain_pair_iptm": value}
    if detail_path is not None:
        detail = Path(detail_path)
        if not detail.name.startswith("pae_") or detail.suffix != ".npz":
            raise ValueError("Boltz-2 detail requires writer's pae_MODEL.npz")
        plddt = detail.with_name("plddt_" + detail.name[4:])
        result.detail_sha256 = hashlib.sha256(detail.read_bytes() + plddt.read_bytes()).hexdigest()
        with np.load(detail, allow_pickle=False) as archive:
            result.pae = _matrix(archive["pae"], len(ids), "pae")
        with np.load(plddt, allow_pickle=False) as archive:
            values = np.asarray(archive["plddt"], dtype=float)
        if (
            values.shape != (len(ids),)
            or not np.isfinite(values).all()
            or ((values < 0) | (values > 1)).any()
        ):
            raise ValueError("Boltz-2 token pLDDT dimensions/range invalid")
        result.token_index = {key: i for i, key in enumerate(ids)}
        result.token_mapping_status = "producer_identity_available"
        # Boltz stores probabilities; expose the same 0..100 pLDDT units as AF2.
        result.residue_plddts = dict(zip(ids, (100 * values).tolist()))
    return result


def confidence_on_contacts(
    contacts: pl.DataFrame, structure, confidence: PredictionConfidence | None
) -> pl.DataFrame:
    """Attach global/chain-pair/token-pair evidence to directed tcren residue contacts.

    Author IDs are looked up from each residue, never equated to a sequential index. This is
    important for AHO-renumbered TCRmodel2 output. Unmatched AF3 token IDs raise instead of
    silently producing a plausible wrong PAE. Native structures pass ``None``.
    """
    author = {
        (c.chain_id, r.seq_index): (c.chain_id, r.pdb_index)
        for c in structure.chains
        for r in c.residues
    }
    if confidence and confidence.generator == "af2-multimer":
        sequences = {c.chain_id: c.sequence() for c in structure.chains}
        mapping = {cid: cid for cid in sequences}
        if sequences != confidence.token_sequences:
            # AF2's PDB writer can rename one-based asym_id values B, C, ... .
            # Preserve the retained input token IDs and require an exact bijection.
            mapping = {}
            for cid, sequence in sequences.items():
                matches = [
                    key for key, seq in confidence.token_sequences.items() if seq == sequence
                ]
                if len(matches) != 1:
                    raise ValueError("Absent/ambiguous exact AF2 producer chain sequence")
                mapping[cid] = matches[0]
            if len(mapping) != len(confidence.token_sequences) or set(mapping.values()) != set(
                confidence.token_sequences
            ):
                raise ValueError("AF2 producer/coordinate chain mapping is not bijective")
        author = {
            (c.chain_id, r.seq_index): (mapping[c.chain_id], i)
            for c in structure.chains
            for i, r in enumerate(c.residues)
        }
    elif confidence and confidence.generator in {"boltz2", "tfold-tcr"}:
        if {c.chain_id: c.sequence() for c in structure.chains} != confidence.token_sequences:
            raise ValueError(
                f"{confidence.generator} coordinates differ from producer chain sequences"
            )
        author = {
            (c.chain_id, r.seq_index): (c.chain_id, i)
            for c in structure.chains
            for i, r in enumerate(c.residues)
        }
    elif confidence and confidence.token_sequences:
        # TCRmodel2 rewrites author numbering and chain IDs after inference. Match exact
        # sequence spans to the retained producer chain; never guess from final IDs.
        author = {}
        for chain in structure.chains:
            sequence = chain.sequence()
            matches = [
                (cid, seq.index(sequence))
                for cid, seq in confidence.token_sequences.items()
                if sequence in seq and seq.find(sequence, seq.index(sequence) + 1) < 0
            ]
            if len(matches) != 1:
                raise ValueError(f"absent/ambiguous retained TCRmodel2 sequence: {chain.chain_id}")
            cid, start = matches[0]
            author.update(
                {
                    (chain.chain_id, r.seq_index): (cid, start + i)
                    for i, r in enumerate(chain.residues)
                }
            )
    extra = []
    for row in contacts.iter_rows(named=True):
        item = {
            "confidence.generator": confidence.generator if confidence else None,
            "confidence.source": confidence.source if confidence else None,
            "confidence.sha256": confidence.sha256 if confidence else None,
            "confidence.detail_sha256": confidence.detail_sha256 if confidence else None,
            "iptm": confidence.global_values.get("iptm") if confidence else None,
            "tcr-pmhc_iptm": confidence.global_values.get("tcr-pmhc_iptm") if confidence else None,
            "chain_pair_iptm": None,
            "chain_pair_pae_min": None,
            "token_pair_pae": None,
            "token_pair_contact_probability": None,
            "residue_plddt.from": None,
            "residue_plddt.to": None,
        }
        if confidence:
            a, b = row["chain.id.from"], row["chain.id.to"]
            item.update(confidence.chain_pairs.get((a, b), {}))
            if confidence.token_index:
                keys = [
                    author[(row[f"chain.id.{side}"], row[f"residue.index.{side}"])]
                    for side in ("from", "to")
                ]
                if any(key not in confidence.token_index for key in keys):
                    raise ValueError(
                        f"coordinate residue missing from confidence token identities: {keys}"
                    )
                i, j = (confidence.token_index[key] for key in keys)
                if confidence.pae is not None:
                    item["token_pair_pae"] = float(confidence.pae[i, j])
                if confidence.contact_probs is not None:
                    item["token_pair_contact_probability"] = float(confidence.contact_probs[i, j])
                for side, key in zip(("from", "to"), keys):
                    item[f"residue_plddt.{side}"] = confidence.residue_plddts.get(key)
        extra.append(item)
    schema = {
        "confidence.generator": pl.String,
        "confidence.source": pl.String,
        "confidence.sha256": pl.String,
        "confidence.detail_sha256": pl.String,
        "iptm": pl.Float64,
        "tcr-pmhc_iptm": pl.Float64,
        "chain_pair_iptm": pl.Float64,
        "chain_pair_pae_min": pl.Float64,
        "token_pair_pae": pl.Float64,
        "token_pair_contact_probability": pl.Float64,
        "residue_plddt.from": pl.Float64,
        "residue_plddt.to": pl.Float64,
    }
    return contacts.hstack(pl.DataFrame(extra, schema=schema))


from .runtime import resolve_threads, single_threaded


@single_threaded
def inspect_predictions(
    manifest: pl.DataFrame, *, cutoff: float = 5.0, organism: str = "human", threads: int = 1,
    historical_assessment: bool = False,
) -> dict[str, pl.DataFrame]:
    """One batch -> model rows, residue markup, descriptors and confidence-bearing contacts.

    Required manifest columns: ``complex.id``, ``structure``. Optional: ``generator`` (af3,
    tcrmodel2, boltz2, tfold-tcr or af2-multimer), ``confidence``, ``detail_confidence``, ``producer_structure``,
    ``model_key``. Paths should be absolute;
    the CLI resolves relative paths against the manifest's directory. IDs must be globally unique.
    TCR-only / pMHC-only structures get annotation and contacts, but ternary descriptors are
    explicitly unavailable. Fitted readouts are retired and only calculated with
    ``historical_assessment=True`` for historical reproduction. ``threads`` budgets both batched annotation and descriptor workers.
    """
    from .annotation import classify_chains
    from .annotation.batch import annotate_batch
    from .contacts.table import residue_annotation
    from .contactmap import ContactMap
    from .descriptors.table import recognition_table
    from .docking import docking_angles, tcr_placement
    from .mhc import annotate_mhc_batch
    from .structure import import_structure

    threads = resolve_threads(threads)
    if not np.isfinite(cutoff) or cutoff <= 0:
        raise ValueError("threads and cutoff must be positive")
    if not {"complex.id", "structure"} <= set(manifest.columns) or not manifest.height:
        raise ValueError("nonempty manifest with complex.id and structure required")
    ids = manifest["complex.id"].to_list()
    if any(not isinstance(id_, str) or not id_ for id_ in ids) or len(set(ids)) != len(ids):
        raise ValueError("complex.id must contain unique nonempty strings")
    entries = manifest.to_dicts()
    structures = [
        import_structure(e["structure"], pdb_id=e["complex.id"], keep_hydrogens=False)
        for e in entries
    ]
    records = annotate_batch(structures, organisms=(organism,), threads=threads)
    for s, rec in zip(structures, records):
        classify_chains(s, organism=organism, autodetect_species=False, precomputed_records=rec)
    annotate_mhc_batch(structures, threads=threads)
    contacts, residues, models, complete, angles, chain_rows, interfaces = (
        [],
        [],
        [],
        [],
        [],
        [],
        [],
    )
    for entry, s in zip(entries, structures):
        chain_rows.extend(
            {
                "complex.id": s.pdb_id,
                "chain.id": chain.chain_id,
                "chain.type": chain.chain_type,
                "chain.supertype": chain.chain_supertype,
                "allele_info": chain.allele_info,
                "sequence": chain.sequence(),
                "residues": len(chain.residues),
                "cdr3": next((r.sequence for r in chain.regions if r.region_type == "CDR3"), None),
                "author_residue_ids": json.dumps(
                    [[r.pdb_index, r.insertion_code] for r in chain.residues]
                ),
            }
            for chain in s.chains
        )
        types = {c.chain_type for c in s.chains}
        ternary = {"TRA", "TRB", "PEPTIDE", "MHCa"} <= types
        confidence = None
        if entry.get("confidence"):
            confidence = read_prediction_confidence(
                entry["confidence"],
                generator=entry["generator"],
                model_key=entry.get("model_key"),
                detail_path=entry.get("detail_confidence"),
                producer_structure_path=entry.get("producer_structure"),
            )
            if (
                confidence.generator == "boltz2"
                and {c.chain_id: c.sequence() for c in s.chains} != confidence.token_sequences
            ):
                raise ValueError("Boltz-2 coordinates differ from producer chain sequences")
        row = {
            **entry,
            "descriptor_status": "available" if ternary else "incomplete_construct",
            "structure_sha256": hashlib.sha256(Path(entry["structure"]).read_bytes()).hexdigest(),
            "confidence_sha256": confidence.sha256 if confidence else None,
            "confidence_token_mapping_status": confidence.token_mapping_status if confidence else "unavailable",
            "detail_confidence_sha256": confidence.detail_sha256 if confidence else None,
            "producer_structure_sha256": confidence.producer_structure_sha256
            if confidence
            else None,
            **(confidence.global_values if confidence else {}),
        }
        cm = ContactMap.from_structure(s, cutoff=cutoff, count_atoms=True)
        ct = confidence_on_contacts(cm.contacts, s, confidence)
        contacts.append(ct.with_columns(pl.lit(s.pdb_id).alias("complex.id")))
        pairs = confidence_chain_pairs(confidence or PredictionConfidence("", "unavailable"))
        counts = ct.group_by("chain.id.from", "chain.id.to").len(name="directed_contact_rows")
        pairs = pairs.join(
            counts, on=["chain.id.from", "chain.id.to"], how="left", validate="1:1"
        ).with_columns(
            pl.col("directed_contact_rows").fill_null(0), pl.lit(s.pdb_id).alias("complex.id")
        )
        interfaces.append(pairs)
        residues.append(residue_annotation(s).with_columns(pl.lit(s.pdb_id).alias("complex.id")))
        for interface in ("tcr_peptide", "tcr_mhc", "peptide_mhc"):
            row[f"contacts.{interface}"] = cm.interface(interface).height
        models.append(row)
        angles.append(
            {
                "complex.id": s.pdb_id,
                "angle_status": "available" if ternary else "incomplete_construct",
                **asdict(docking_angles(s)),
                **asdict(tcr_placement(s)),
            }
            if ternary
            else {"complex.id": s.pdb_id, "angle_status": "incomplete_construct"}
        )
        if ternary:
            complete.append((s.pdb_id, s))
    features = (
        pl.DataFrame(
            recognition_table(
                complete,
                annotate=False,
                threads=threads,
                organism=organism,
                include=("placement", "interface", "topology", "energetics"),
            ),
            infer_schema_length=None,
        )
        if complete
        else pl.DataFrame(schema={"complex.id": pl.String})
    )
    model_table = pl.DataFrame(models, infer_schema_length=None)
    joined = model_table.join(features, on="complex.id", how="left", validate="1:1")
    if "error" in joined.columns and joined["error"].is_not_null().any():
        raise RuntimeError(
            f"descriptor extraction failed: {joined.filter(pl.col('error').is_not_null()).select('complex.id', 'error').to_dicts()}"
        )
    assessment = pl.DataFrame(schema={"complex.id": pl.String, "assessment_status": pl.String})
    if historical_assessment:
        from .score import score_table

        ternary_table = joined.filter(pl.col("descriptor_status") == "available")
        assessment = (
            score_table(ternary_table)
            if ternary_table.height
            else pl.DataFrame(schema={"complex.id": pl.String})
        )
        if ternary_table.height:
            assessment = assessment.with_columns(
                pl.when(pl.col("assessment_available"))
                .then(pl.lit("available"))
                .otherwise(pl.lit("missing_descriptor_values"))
                .alias("assessment_status")
            )
        if ternary_table.height and "iptm" in ternary_table.columns:
            confident = ternary_table.filter(pl.col("iptm").is_not_null() & pl.col("iptm").is_finite())
            if confident.height:
                evidence = score_table(confident, iptm=confident["iptm"].to_numpy()).select(
                    "complex.id", "confidence_residual", "binder_iptm",
                    "confidence_residual_available", "binder_iptm_available",
                )
                assessment = assessment.join(evidence, on="complex.id", how="left", validate="1:1")
    joined = joined.join(
        assessment.select("complex.id", "assessment_status")
        if "assessment_status" in assessment.columns
        else pl.DataFrame(schema={"complex.id": pl.String, "assessment_status": pl.String}),
        on="complex.id",
        how="left",
        validate="1:1",
    ).with_columns(
        pl.col("assessment_status").fill_null(
            pl.when(pl.col("descriptor_status") == "incomplete_construct")
            .then(pl.lit("incomplete_construct"))
            .otherwise(pl.lit("historical_not_requested"))
        )
    )
    return {
        "models": joined,
        "residues": pl.concat(residues, how="diagonal_relaxed"),
        "contacts": pl.concat(contacts, how="diagonal_relaxed"),
        "assessment": assessment,
        "angles": pl.DataFrame(angles, infer_schema_length=None),
        "chains": pl.DataFrame(chain_rows, infer_schema_length=None),
        "interfaces": pl.concat(interfaces, how="diagonal_relaxed"),
    }


def inspect_models_command(
    manifest: Path,
    out: Path,
    cutoff: float = 5.0,
    organism: str = "human",
    threads: int = 1,
    historical_assessment: bool = False,
) -> None:
    """Inspect an AF3/TCRmodel2/native TSV manifest; write models, residues and contacts TSVs."""
    from .provenance import stamp

    table = pl.read_csv(manifest, separator="\t", infer_schema_length=None)
    for column in ("structure", "confidence", "detail_confidence", "producer_structure"):
        if column in table.columns:
            table = table.with_columns(
                pl.col(column).map_elements(
                    lambda p: str((manifest.parent / p).resolve()) if p else None,
                    return_dtype=pl.String,
                )
            )
    results = inspect_predictions(
        table, cutoff=cutoff, organism=organism, threads=threads,
        historical_assessment=historical_assessment,
    )
    out.mkdir(parents=True, exist_ok=True)
    for name, result in results.items():
        path = out / f"{name}.tsv"
        result.write_csv(path, separator="\t")
        stamp(
            path,
            command="tcren inspect-models",
            columns=result.columns,
            extra={
                "manifest": str(manifest.resolve()),
                "cutoff_A": cutoff,
                "threads": threads,
                "historical_assessment": historical_assessment,
                "rows": result.height,
            },
        )
