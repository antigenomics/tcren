"""Literal AF3 input auditing without inferring missing biological domains."""

import hashlib
import csv
import io
import json
from pathlib import Path


def af3_to_boltz_input(source, output_dir):
    """Export literal proteins and retained unpaired MSAs for a Boltz control.

    No templates or pairing are transferred. CSV key -1 explicitly prevents
    accidental pairing through reused accession headers. The query is retained
    first; all alignment rows and lowercase insertions are preserved.
    """
    data = json.loads(Path(source).read_text())
    expected = {e["protein"]["id"]: e["protein"]["sequence"] for e in data["sequences"]}
    audit = audit_af3_input(source, expected)
    output = Path(output_dir).resolve()
    output.mkdir(exist_ok=False, parents=True)
    entities, provenance = [], []
    for entry in data["sequences"]:
        protein = entry["protein"]
        cid, seq = protein["id"], protein["sequence"]
        rows = []
        for line in (protein.get("unpairedMsa") or "").splitlines():
            if line.startswith(">"):
                rows.append("")
            elif line.strip():
                if not rows:
                    raise ValueError(f"MSA has no header: {cid}")
                rows[-1] += line.strip()
        if not rows:
            rows = [seq]
        if any(len("".join(c for c in row if not c.islower())) != len(seq) for row in rows):
            raise ValueError(f"MSA alignment width differs from query: {cid}")
        msa = output / f"{cid}.csv"
        with msa.open("w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(["key", "sequence"])
            writer.writerows((-1, row) for row in rows)
        entities.append({"protein": {"id": cid, "sequence": seq, "msa": str(msa)}})
        provenance.append(
            {
                "chain": cid,
                "rows": len(rows),
                "msa_sha256": hashlib.sha256(msa.read_bytes()).hexdigest(),
            }
        )
    target = output / "input.yaml"
    target.write_text(json.dumps({"version": 1, "sequences": entities}, indent=2))
    report = {
        "source_audit": audit,
        "output": str(target),
        "msas": provenance,
        "pairing": "all CSV keys -1",
        "templates": "none transferred",
        "interpretation": "retained AF3 unpaired search outputs; no fresh MSA search",
    }
    (output / "provenance.json").write_text(json.dumps(report, indent=2))
    return report


def reference_complete_pmhc(source, reference, target, *, beta2m_sequence):
    """Explicit resolved-ectodomain/B2M control; preserve observed TCR sequences.

    The reference MHC must begin with the exact original MHC sequence. B2M is
    an explicitly supplied canonical mature sequence; reference variants are
    recorded rather than silently copied. No TCR constants are inferred.
    """
    from .structures import read_chains

    data = json.loads(Path(source).read_text())
    proteins = {e["protein"]["id"]: e["protein"] for e in data["sequences"]}
    mhc = proteins["D"]["sequence"]
    chains = read_chains(reference)
    candidates = [
        (cid, seq) for cid, seq, _ in chains if seq.startswith(mhc) and len(seq) > len(mhc)
    ]
    if (
        len(candidates) != 1
        or not beta2m_sequence
        or any(c not in "ACDEFGHIKLMNPQRSTVWY" for c in beta2m_sequence)
    ):
        raise ValueError(
            "Require exact original MHC prefix and explicit standard-amino-acid B2M sequence"
        )
    b2m_variants = [
        {
            "chain": cid,
            "substitutions": [
                {"position": i + 1, "canonical": a, "reference": b}
                for i, (a, b) in enumerate(zip(beta2m_sequence, seq))
                if a != b
            ],
        }
        for cid, seq, _ in chains
        if len(seq) == len(beta2m_sequence)
    ]
    if Path(target).exists():
        raise FileExistsError(target)
    proteins["D"].clear()
    proteins["D"].update({"id": "D", "sequence": candidates[0][1]})
    data["sequences"].append({"protein": {"id": "E", "sequence": beta2m_sequence}})
    data["name"] = "tcrpmhc_reference_completed_pmhc"
    data["modelSeeds"] = [42]
    Path(target).write_text(json.dumps(data))
    return {
        "reference": str(reference),
        "reference_sha256": hashlib.sha256(Path(reference).read_bytes()).hexdigest(),
        "source_sha256": hashlib.sha256(Path(source).read_bytes()).hexdigest(),
        "mhc_reference_chain": candidates[0][0],
        "original_mhc_residues": len(mhc),
        "resolved_reference_mhc_residues": len(candidates[0][1]),
        "beta2m_residues": len(beta2m_sequence),
        "equal_length_reference_chain_variants": b2m_variants,
        "interpretation": "explicit reference-completed pMHC input control, not the original four-chain construct",
    }


def audit_af3_input(path, expected_sequences):
    """Validate chain identity, MSA query sequences and template index bounds."""
    raw = Path(path).read_bytes()
    data = json.loads(raw)
    if data.get("dialect") != "alphafold3":
        raise ValueError("Expected native alphafold3 input dialect")
    observed = {}
    rows = []
    for entity in data["sequences"]:
        protein = entity["protein"]
        cid = protein["id"]
        if (not isinstance(cid, str) or not cid or cid in observed
                or cid in {".", ".."} or Path(cid).name != cid
                or any(character in cid for character in ("\\", "\x00", ":"))):
            raise ValueError("Expected unique safe single-component protein identifiers")
        seq = protein["sequence"]
        observed[cid] = seq
        row = {"chain": cid, "residues": len(seq)}
        for kind in ("unpairedMsa", "pairedMsa"):
            msa = protein.get(kind)
            row[kind + "_provided"] = msa is not None
            row[kind + "_rows"] = 0 if not msa else msa.count(">")
            if msa:
                lines = msa.splitlines()
                query = []
                for line in lines[1:]:
                    if line.startswith(">"):
                        break
                    query.append(line.strip())
                query = "".join(query)
                query = "".join(c for c in query if c.isupper())
                if query != seq:
                    raise ValueError(f"MSA query differs from chain {cid}/{kind}")
        templates = protein.get("templates", [])
        row["templates"] = len(templates)
        row["template_mapped_residues"] = []
        row["template_details"] = []
        for template in templates:
            query = template["queryIndices"]
            target = template["templateIndices"]
            if len(query) != len(target) or len(set(query)) != len(query):
                raise ValueError(f"Invalid template mapping for {cid}")
            if any(not isinstance(i, int) or isinstance(i, bool) for i in query + target):
                raise ValueError(f"Template indices must be integers for {cid}")
            if any(i < 0 or i >= len(seq) for i in query) or any(i < 0 for i in target):
                raise ValueError(f"Out-of-range template mapping for {cid}")
            from Bio.PDB.MMCIF2Dict import MMCIF2Dict

            if not template.get("mmcif"):
                raise ValueError("Processed-input audit requires inline template mmCIF")
            cif = MMCIF2Dict(io.StringIO(template["mmcif"]))
            numbers = list(map(int, cif.get("_entity_poly_seq.num", [])))
            entities = set(cif.get("_entity_poly_seq.entity_id", []))
            if len(entities) != 1 or numbers != list(range(1, len(numbers) + 1)):
                raise ValueError(f"Expected one sequential template polymer for {cid}")
            if any(i >= len(numbers) for i in target):
                raise ValueError(f"Template target index exceeds polymer length for {cid}")
            resolved = {int(i) - 1 for i in cif.get("_atom_site.label_seq_id", []) if i.isdigit()}
            row["template_mapped_residues"].append(len(query))
            row["template_details"].append(
                {
                    "entry_id": cif.get("_entry.id", [None])[0],
                    "description": cif.get("_entity.pdbx_description", [None])[0],
                    "polymer_residues": len(numbers),
                    "mapped_residues": len(query),
                    "mapped_resolved_residues": len(set(target) & resolved),
                    "query_coverage": len(query) / len(seq),
                    "template_index_max": max(target) if target else None,
                }
            )
        rows.append(row)
    if observed != expected_sequences:
        raise ValueError("Input chains/sequences differ from the explicit expected construct")
    return {
        "input": str(path),
        "sha256": hashlib.sha256(raw).hexdigest(),
        "model_seeds": data["modelSeeds"],
        "chains": rows,
    }


def disable_af3_pairing(source, target, expected_sequences, *, name, seeds):
    """Write an explicit controlled arm, preserving unpaired MSAs and templates."""
    audit = audit_af3_input(source, expected_sequences)
    data = json.loads(Path(source).read_text())
    if Path(target).exists():
        raise FileExistsError(target)
    for entry in data["sequences"]:
        entry["protein"]["pairedMsa"] = ""
    data["name"] = name
    data["modelSeeds"] = list(seeds)
    Path(target).write_text(json.dumps(data))
    return {
        "source_audit": audit,
        "output_audit": audit_af3_input(target, expected_sequences),
        "change": "pairedMsa set to empty; only name/seeds otherwise changed",
        "interpretation": "controlled hypothesis test; reused explicit search outputs, not fresh searches",
    }


def remove_af3_templates(source, target, expected_sequences, *, chains, name):
    """Controlled template ablation; preserve every MSA, sequence and seed."""
    audit = audit_af3_input(source, expected_sequences)
    data = json.loads(Path(source).read_text())
    if Path(target).exists():
        raise FileExistsError(target)
    if not set(chains) <= expected_sequences.keys():
        raise ValueError("Unknown ablation chain")
    for entry in data["sequences"]:
        if entry["protein"]["id"] in chains:
            entry["protein"]["templates"] = []
    data["name"] = name
    Path(target).write_text(json.dumps(data))
    return {
        "source_audit": audit,
        "output_audit": audit_af3_input(target, expected_sequences),
        "change": "templates removed for explicit chains",
        "chains": list(chains),
    }


def prepare_af3_search_input(source, target, sequences, *, name, seeds, reuse_roles=('D',),
                             single_sequence_roles=('C',)):
    """Reuse exact-sequence prepared proteins; search changed chains explicitly.

    This prepares a partial-search input, not inference-ready features. New
    protein chains omit both MSAs/templates so upstream searches them together.
    Short peptide roles explicitly use a single sequence. Pairing ablation is
    applied to the completed prepared input using disable_af3_pairing.
    """
    import copy

    data = json.loads(Path(source).read_text())
    proteins = {row['protein']['id']: row['protein'] for row in data['sequences']}
    audit = audit_af3_input(source, {k:v['sequence'] for k,v in proteins.items()})
    if (data.get('userCCD') or data.get('bondedAtomPairs') or
            any(row.get('modifications') for row in proteins.values())):
        raise ValueError('Partial replay search requires unmodified unbonded proteins')
    if not sequences or any(not isinstance(v,str) or not v for v in sequences.values()):
        raise ValueError('Explicit nonempty protein sequences required')
    if not set(reuse_roles) | set(single_sequence_roles) <= set(sequences):
        raise ValueError('Unknown reuse/single-sequence role')
    if set(reuse_roles) & set(single_sequence_roles):
        raise ValueError('Reuse and single-sequence roles must be disjoint')
    if Path(target).exists():
        raise FileExistsError(target)
    entities, decisions = [], []
    for role, sequence in sequences.items():
        if role in single_sequence_roles:
            protein = dict(id=role,sequence=sequence,unpairedMsa='',pairedMsa='',templates=[])
            decision = 'explicit_single_sequence'
        elif role in reuse_roles and role in proteins and proteins[role]['sequence'] == sequence:
            protein = copy.deepcopy(proteins[role])
            if any(protein.get(key) is None for key in ('unpairedMsa','pairedMsa','templates')):
                raise ValueError('Reused protein must have fully prepared MSAs and templates')
            decision = 'reused_exact_sequence_prepared_chain'
        else:
            protein = dict(id=role,sequence=sequence)
            decision = 'fresh_msa_template_search_required'
        entities.append(dict(protein=protein))
        decisions.append(dict(role=role,decision=decision))
    output = dict(dialect='alphafold3',version=3,name=name,modelSeeds=list(seeds),sequences=entities)
    Path(target).write_text(json.dumps(output))
    return dict(source_audit=audit,output_audit=audit_af3_input(target,sequences),decisions=decisions,
                scope='Partial fresh search with declared exact-sequence prepared-chain reuse; no model coordinates reused')
