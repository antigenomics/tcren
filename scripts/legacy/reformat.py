"""Reformat a native TCR:pMHC PDB into the isalgo/vdjdb_structure_models layout.

2026-07-19

Structure work (parse, canonical orientation, chain A/B/C/D/E rename, annotation) is done with
built-in ``tcren``; the coordinate / skeleton-plot / contact *format* is ported from the legacy
``tcr-structures-visualization`` pipeline (``produce_plots_pipline/``). One native structure →
the files the HF dataset ships plus a metadata row:

* ``aligned_<pdbid>.pdb``                     (pdb_files_native.tgz)  — tcren canonical PDB
* ``<pdbid>_aa_coordinates.tsv``              (coordinates_aa.tgz)    — Cα table, chains mapped
* ``<hash>_<pdbid>_aa_contacts.tsv``          (contacts_aa.tgz)       — CDR3α/β–peptide ≤5 Å
* ``<hash>.svg`` / ``<hash>_simplified.svg``  (maps, deduped by hash) — skeleton plots
* one metadata row (45 cols, ``is_native=True``)

The ``tcr_pmhc_hash`` is ``sha256`` over the 9 VDJdb identity fields (see ``HASH_KEYS``); it is
the dataset's identity key, so metadata fields come from the joined VDJdb record when available
and fall back to tcren annotation otherwise. Run ``python reformat.py`` for the self-check.
"""

from __future__ import annotations

import hashlib
import os
import sys
import warnings
from pathlib import Path

warnings.simplefilter("ignore")

# ---------------------------------------------------------------------------------------------
# Paths (env-overridable). All legacy assets are local — no GitLab / network at run time.
# ---------------------------------------------------------------------------------------------
LEGACY_REPO = Path(os.environ.get(
    "TCREN_LEGACY_REPO", "/Users/mikesh/vcs/code/tcr-structures-visualization"))
LEGACY_PIPELINE = LEGACY_REPO / "produce_plots_pipline"
PCA_PATH = Path(os.environ.get("TCREN_LEGACY_PCA", LEGACY_PIPELINE / "pca_all_structures.sav"))
NATIVE2026_DIR = Path(os.environ.get(
    "TCREN_NATIVE2026", os.path.expanduser("~/hf/tcren_structures/Native2026")))

# The 9 identity fields, concatenated (no separator, UTF-8) → sha256 = tcr_pmhc_hash.
HASH_KEYS = ["cdr3.alpha", "v.alpha", "j.alpha", "cdr3.beta", "v.beta", "j.beta",
             "mhc.a", "mhc.b", "antigen.epitope"]

# Legacy chain-letter → dataset chain name (tcren's canonical rename matches these letters).
CHAIN_FICT = {"A": "TCR_alpha", "B": "TCR_beta", "C": "peptide", "D": "MHC"}
AA3TO1 = {
    "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F", "GLY": "G", "HIS": "H",
    "ILE": "I", "LYS": "K", "LEU": "L", "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q",
    "ARG": "R", "SER": "S", "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y",
}


def _clean_allele(allele: str) -> str:
    """Drop a trailing ``:UniProt`` accession arda appends (``H2-Aa:P01910`` → ``H2-Aa``).

    HLA alleles keep their ``:NN`` field (``HLA-A*02:01``): a UniProt tail starts with a letter,
    an HLA field is digits.
    """
    if not allele or ":" not in allele:
        return allele
    head, _, tail = allele.rpartition(":")
    return head if (tail and tail[0].isalpha() and tail.isalnum()) else allele


def tcr_hash(fields: dict) -> str:
    """SHA-256 of the 9 identity fields concatenated in ``HASH_KEYS`` order (verified vs dataset)."""
    joined = "".join(str(fields[k]) for k in HASH_KEYS)
    return hashlib.sha256(joined.encode("utf-8")).hexdigest()


# ---------------------------------------------------------------------------------------------
# Lazy legacy-module + PCA loading (matplotlib Agg + monospace font, as the legacy pipeline).
# ---------------------------------------------------------------------------------------------
_LEGACY = {}


def _legacy():
    if _LEGACY:
        return _LEGACY
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams["font.family"] = "monospace"
    import pickle
    sys.path.insert(0, str(LEGACY_PIPELINE))
    import coordinates  # noqa: E402  (legacy module)
    import plotting  # noqa: E402  (legacy module)
    _LEGACY.update(coordinates=coordinates, plotting=plotting, plt=plt,
                   pca=pickle.load(open(PCA_PATH, "rb")))
    return _LEGACY


# ---------------------------------------------------------------------------------------------
# tcren: canonical orientation + identity-field annotation
# ---------------------------------------------------------------------------------------------
def _native_path(pdbid: str) -> Path:
    for ext in (".pdb", ".pdb.gz", ".cif", ".cif.gz"):
        p = NATIVE2026_DIR / f"{pdbid}{ext}"
        if p.exists():
            return p
    raise FileNotFoundError(f"{pdbid} not found under {NATIVE2026_DIR}")


def annotate_and_orient(pdbid: str):
    """Parse → annotate (chains + MHC) → canonicalize. Returns ``(oriented, CanonResult, fields)``.

    Identity fields are read from the annotated structure *before* ``canonicalize_structure``,
    because tcren's ``apply_transform`` clears region markup.
    """
    from tcren.structure import parse_structure
    from tcren.annotation import classify_chains
    from tcren.mhc import annotate_mhc
    from tcren.orient import canonicalize_structure
    structure = parse_structure(str(_native_path(pdbid)), pdb_id=pdbid)
    classify_chains(structure)
    annotate_mhc(structure)
    fields = annotate_fields(structure)
    oriented, res = canonicalize_structure(structure)
    return oriented, res, fields


def _chain_by_type(structure, ctype: str):
    return next((c for c in structure.chains if c.chain_type == ctype), None)


def _region(chain, name: str) -> str:
    if not chain or not getattr(chain, "regions", None):
        return ""
    for rm in chain.regions:
        rt = getattr(rm.region_type, "name", None) or str(rm.region_type)
        if rt.upper() == name:
            return rm.sequence or ""
    return ""


def _cdr3_junction(chain) -> str:
    """VDJdb-style CDR3 junction: conserved Cys104 + IMGT CDR3 + conserved Phe/Trp118."""
    cdr3 = _region(chain, "CDR3")
    if not cdr3:
        return ""
    fr3, fr4 = _region(chain, "FR3"), _region(chain, "FR4")
    return (fr3[-1] if fr3 else "") + cdr3 + (fr4[0] if fr4 else "")


def annotate_fields(structure) -> dict:
    """Extract the 9 identity fields from a tcren-annotated structure (region markup must be live).

    The metadata source for structures absent from VDJdb, and the cross-check of the VDJdb-derived
    hash on joinable structures (see 10_validate.py). CDR3s are junction-form to match VDJdb.
    """
    def _seq(chain):
        return "".join(r.aa for r in chain.residues if getattr(r, "aa", None)) if chain else ""

    def _vj(chain):
        info = getattr(chain, "allele_info", None) if chain else None
        if info and ":" in info:
            v, j = info.split(":", 1)
            return v, j
        return "", ""

    tra, trb = _chain_by_type(structure, "TRA"), _chain_by_type(structure, "TRB")
    pep = _chain_by_type(structure, "PEPTIDE")
    mhca, mhcb = _chain_by_type(structure, "MHCa"), _chain_by_type(structure, "MHCb")
    va, ja = _vj(tra)
    vb, jb = _vj(trb)
    mhc_a = _clean_allele((getattr(mhca, "allele_info", "") or "") if mhca else "")
    # class I: mhc.b is B2M; class II: the MHCb allele.
    mhc_b = "B2M" if (mhca and getattr(mhca, "chain_supertype", "") == "MHCI") else (
        _clean_allele((getattr(mhcb, "allele_info", "") or "") if mhcb else ""))
    return {
        "cdr3.alpha": _cdr3_junction(tra), "v.alpha": va, "j.alpha": ja,
        "cdr3.beta": _cdr3_junction(trb), "v.beta": vb, "j.beta": jb,
        "mhc.a": mhc_a, "mhc.b": mhc_b, "antigen.epitope": _seq(pep),
    }


# ---------------------------------------------------------------------------------------------
# File generation (ports the legacy coordinate / plot / contact writers)
# ---------------------------------------------------------------------------------------------
def _write_aligned_pdb(oriented, pdbid: str, out_dir: Path) -> Path:
    from tcren.structure.io import write_structure
    out = out_dir / f"aligned_{pdbid}.pdb"
    write_structure(oriented, out)
    return out


def _write_coords_tsv(aligned_pdb: Path, pdbid: str, out_dir: Path) -> Path:
    """Cα table ``<pdbid>_aa_coordinates.tsv`` (chains mapped A→TCR_alpha…), no index column."""
    import pandas as pd
    rows = []
    for line in aligned_pdb.read_text().splitlines():
        if not line.startswith("ATOM  ") or line[12:16].strip() != "CA":
            continue
        rows.append([line[21].strip() or "-", line[17:20].strip(), line[22:27].strip(),
                     round(float(line[30:38]), 3), round(float(line[38:46]), 3),
                     round(float(line[46:54]), 3)])
    df = pd.DataFrame(rows, columns=["Chain", "Residue", "ResNum", "X", "Y", "Z"])
    df["Chain"] = df["Chain"].map(CHAIN_FICT).fillna(df["Chain"])
    df["Residue"] = df["Residue"].map(AA3TO1).fillna(df["Residue"])
    out = out_dir / f"{pdbid}_aa_coordinates.tsv"
    df.to_csv(out, sep="\t", index=False)
    return out


def _write_plots_and_contacts(aligned_pdb: Path, pdbid: str, hash_: str,
                              cdr3a: str, cdr3b: str, epitope: str, out_dir: Path,
                              make_map: bool) -> dict:
    """Skeleton SVGs (deduped by hash upstream) + ``<hash>_<pdbid>_aa_contacts.tsv`` (no index).

    Returns ``{"contacts": path, "svg": path|None, "simplified": path|None, "num_contacts": int}``.
    """
    import pandas as pd
    lg = _legacy()
    coords_ca, coords_all, atom_to_ca = lg["coordinates"].extract_coords_from_pdb_by_seq(
        str(aligned_pdb), cdr3a, cdr3b, epitope)
    if not coords_ca:
        raise ValueError(f"{pdbid}: no CDR3α/β/peptide residues matched by sequence")
    coords_ca = lg["coordinates"].apply_pca(coords_ca, lg["pca"])

    # The legacy plotter keys every output off its ``pdb_filename`` base; use the hash so the
    # SVGs land as <hash>.svg / <hash>_simplified.svg, then adjust the contacts name + index.
    base = hash_
    lg["plotting"].plot_combined_residue_graph_pca(
        coords_ca, coords_all, atom_to_ca, max_distance=5.0,
        pdb_filename=base, save_dir=str(out_dir), simplified=False)
    lg["plotting"].plot_combined_residue_graph_pca(
        coords_ca, coords_all, atom_to_ca, max_distance=5.0,
        pdb_filename=base, save_dir=str(out_dir), simplified=True)

    # Contacts: legacy writes <hash>_aa_contacts.tsv WITH a pandas index; re-emit without index
    # under the native name <hash>_<pdbid>_aa_contacts.tsv, matching the shipped format.
    raw = out_dir / f"{base}_aa_contacts.tsv"
    contacts = pd.read_csv(raw, sep="\t", index_col=0)
    contacts = contacts[["chain_from", "aa_from", "res_num_from", "aa_to", "res_num_to", "chain_to"]]
    out_contacts = out_dir / f"{hash_}_{pdbid}_aa_contacts.tsv"
    contacts.to_csv(out_contacts, sep="\t", index=False)
    raw.unlink(missing_ok=True)
    (out_dir / f"{base}_contacts.txt").unlink(missing_ok=True)

    num_contacts = int(contacts[(contacts.chain_to == "peptide") |
                                (contacts.chain_from == "peptide")].shape[0])

    svg = out_dir / f"{hash_}.svg"
    simplified = out_dir / f"{hash_}_simplified.svg"
    if not make_map:  # hash already has a (predicted) map — drop the duplicate.
        svg.unlink(missing_ok=True)
        simplified.unlink(missing_ok=True)
        svg = simplified = None
    return {"contacts": out_contacts, "svg": svg, "simplified": simplified,
            "num_contacts": num_contacts}


def angles(aligned_pdb: Path) -> dict:
    """STCRpy scanning/pitch angles (rounded 4). Returns ``{}`` if stcrpy is unavailable."""
    try:
        import stcrpy
    except Exception:
        return {}
    try:
        tcr = stcrpy.load_TCR(str(aligned_pdb))
        return {"scanning_angle": round(float(tcr.get_scanning_angle()), 4),
                "pitch_angle": round(float(tcr.get_pitch_angle()), 4)}
    except Exception as e:  # noqa: BLE001 — angles are optional; never fail the structure over them
        print(f"  [angles] {aligned_pdb.name}: {type(e).__name__}: {e}", file=sys.stderr)
        return {}


# ---------------------------------------------------------------------------------------------
# Metadata: VDJdb join (joinable) + tcren-annotation fallback (unjoinable) → the 45-column row
# ---------------------------------------------------------------------------------------------
ANNOTATED_TSV = Path(os.environ.get(
    "TCREN_VDJDB_ANNOTATED", LEGACY_REPO / "utils_scripts_and_notebooks/vdjdb_structures_annotated.tsv"))

# The 45 published columns, in order.
PUBLISHED_COLUMNS = [
    "idx", "cdr3.alpha", "v.alpha", "j.alpha", "cdr3.beta", "v.beta", "d.beta", "j.beta",
    "species", "mhc.a", "mhc.b", "mhc.class", "antigen.epitope", "antigen.gene", "antigen.species",
    "reference.id", "method.identification", "method.frequency", "method.singlecell",
    "method.sequencing", "method.verification", "meta.study.id", "meta.cell.subset",
    "meta.subject.cohort", "meta.subject.id", "meta.replica.id", "meta.clone.id",
    "meta.epitope.id", "meta.tissue", "meta.donor.MHC", "meta.donor.MHC.method",
    "meta.structure.id", "cdr3fix.alpha", "cdr3fix.beta", "vdjdb.score", "tcr_pmhc_hash",
    "num_contacts", "ranking_confidence", "plddt", "ptm", "iptm", "tcr_pmhc_iptm",
    "scanning_angle", "pitch_angle", "is_native",
]
# VDJdb "full table" columns copied verbatim from the annotated TSV (published cols 2-35).
VDJDB_COLUMNS = PUBLISHED_COLUMNS[1:35]  # cdr3.alpha … vdjdb.score
# Model-only quality columns — always blank for natives.
MODEL_SCORE_COLUMNS = ["ranking_confidence", "plddt", "ptm", "iptm", "tcr_pmhc_iptm"]
_SPECIES = {"human": "HomoSapiens", "mouse": "MusMusculus",
            "Human": "HomoSapiens", "Mouse": "MusMusculus"}


def load_vdjdb_index(annotated_tsv: Path = ANNOTATED_TSV) -> dict:
    """Map ``meta.structure.id`` (lowercased) → the VDJdb-annotated row dict (first wins if dup)."""
    import csv
    index = {}
    with open(annotated_tsv, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            sid = (row.get("meta.structure.id") or "").strip().lower()
            if sid and sid not in index:
                index[sid] = row
    return index


def metadata_row(pdbid: str, hash_: str, *, vdjdb_rec: dict | None, fields: dict,
                 oriented, num_contacts, angles: dict) -> dict:
    """Build the 45-column native metadata row (``is_native=True``, model scores blank).

    ``vdjdb_rec`` (joinable) supplies the VDJdb columns verbatim; otherwise they are filled from
    tcren ``fields`` (identity) with VDJdb-provenance columns left blank.
    """
    row = {c: "" for c in PUBLISHED_COLUMNS}
    if vdjdb_rec is not None:
        for c in VDJDB_COLUMNS:
            row[c] = vdjdb_rec.get(c, "") or ""
        row["idx"] = vdjdb_rec.get("idx", "") or ""
    else:  # tcren-annotation fallback
        row.update(fields)  # cdr3/v/j/mhc.a/mhc.b/antigen.epitope
        mhca = _chain_by_type(oriented, "MHCa")
        row["mhc.class"] = getattr(mhca, "chain_supertype", "") if mhca else ""
        row["species"] = _SPECIES.get(getattr(oriented, "complex_species", ""),
                                      getattr(oriented, "complex_species", "") or "")
        row["meta.structure.id"] = pdbid
    row["tcr_pmhc_hash"] = hash_
    row["num_contacts"] = "" if num_contacts is None else float(num_contacts)
    row["scanning_angle"] = angles.get("scanning_angle", "")
    row["pitch_angle"] = angles.get("pitch_angle", "")
    row["is_native"] = True
    return row


def native_num_contacts(oriented) -> int:
    """TCR–peptide residue-pair count (tcren ``ContactMap``, cutoff 5.0) — the ``num_contacts`` def."""
    from tcren.contactmap import ContactMap
    return int(ContactMap.from_structure(oriented, cutoff=5.0)
               .interface("tcr_peptide", tcr_regions="all").height)


def process(pdbid: str, out_dir: Path, *, vdjdb_rec: dict | None, make_map: bool = True,
            do_angles: bool = True) -> dict:
    """Full per-structure pipeline → writes the files under ``out_dir``; returns the metadata row.

    ``vdjdb_rec`` joins to a VDJdb record (hash + metadata from it); ``None`` uses tcren annotation.
    """
    oriented, res, fields = annotate_and_orient(pdbid)
    if not (fields["cdr3.alpha"] and fields["cdr3.beta"] and fields["antigen.epitope"]):
        missing = [k for k in ("cdr3.alpha", "cdr3.beta", "antigen.epitope") if not fields[k]]
        raise ValueError(f"{pdbid}: incomplete αβ TCR:pMHC (missing {', '.join(missing)})")
    if vdjdb_rec is not None:
        hash_ = (vdjdb_rec.get("TCR_hash") or vdjdb_rec.get("tcr_pmhc_hash")
                 or tcr_hash({k: vdjdb_rec.get(k, "") for k in HASH_KEYS}))
        # File matching still uses tcren's guaranteed-substring CDR3s / peptide.
        cdr3a, cdr3b, epitope = fields["cdr3.alpha"], fields["cdr3.beta"], fields["antigen.epitope"]
    else:
        hash_ = tcr_hash(fields)
        cdr3a, cdr3b, epitope = fields["cdr3.alpha"], fields["cdr3.beta"], fields["antigen.epitope"]

    apdb = _write_aligned_pdb(oriented, pdbid, out_dir)
    _write_coords_tsv(apdb, pdbid, out_dir)
    plots = _write_plots_and_contacts(apdb, pdbid, hash_, cdr3a, cdr3b, epitope, out_dir, make_map)
    ncontacts = native_num_contacts(oriented)
    ang = angles(apdb) if do_angles else {}
    return metadata_row(pdbid, hash_, vdjdb_rec=vdjdb_rec, fields=fields, oriented=oriented,
                        num_contacts=ncontacts, angles=ang)


def process_prealigned(pdbid: str, aligned_pdb, out_dir: Path, *, make_map: bool = True,
                       do_angles: bool = True) -> dict:
    """Complete a structure that already ships an aligned PDB but has no metadata/hash.

    Annotates the shipped (already-oriented) PDB — no re-orientation, and pdb/coords are NOT
    written (they already exist in the dataset). Emits only contacts + maps + the metadata row.
    Raises if the structure is not a complete αβ TCR:pMHC.
    """
    from tcren.structure import parse_structure
    from tcren.annotation import classify_chains
    from tcren.mhc import annotate_mhc
    aligned_pdb = Path(aligned_pdb)
    s = parse_structure(str(aligned_pdb), pdb_id=pdbid)
    classify_chains(s)
    annotate_mhc(s)
    fields = annotate_fields(s)
    if not (fields["cdr3.alpha"] and fields["cdr3.beta"] and fields["antigen.epitope"]):
        missing = [k for k in ("cdr3.alpha", "cdr3.beta", "antigen.epitope") if not fields[k]]
        raise ValueError(f"{pdbid}: not a complete αβ TCR:pMHC (missing {', '.join(missing)})")
    hash_ = tcr_hash(fields)
    _write_plots_and_contacts(aligned_pdb, pdbid, hash_, fields["cdr3.alpha"], fields["cdr3.beta"],
                              fields["antigen.epitope"], out_dir, make_map)
    ncontacts = native_num_contacts(s)
    ang = angles(aligned_pdb) if do_angles else {}
    return metadata_row(pdbid, hash_, vdjdb_rec=None, fields=fields, oriented=s,
                        num_contacts=ncontacts, angles=ang)


def _self_check():
    # Hash contract: reproduce the published e65469… from idx-1159 fields.
    rec = {"cdr3.alpha": "CAMREGGSNYQLIW", "v.alpha": "TRAV14/DV4*01", "j.alpha": "TRAJ33*01",
           "cdr3.beta": "CASSMIPDMNTEAFF", "v.beta": "TRBV19*01", "j.beta": "TRBJ1-1*01",
           "mhc.a": "HLA-B*08:01", "mhc.b": "B2M", "antigen.epitope": "RPIIRPATL"}
    want = "e65469dacd705ce72895eba2a2429496b764a2b1757f9f31f1d4c69ba7b14319"
    got = tcr_hash(rec)
    assert got == want, f"hash mismatch: {got} != {want}"
    print("OK  tcr_hash reproduces published e65469…")
    assert PCA_PATH.exists(), f"missing PCA model: {PCA_PATH}"
    assert LEGACY_PIPELINE.exists(), f"missing legacy pipeline: {LEGACY_PIPELINE}"
    print(f"OK  legacy assets present ({PCA_PATH.name}, {LEGACY_PIPELINE.name}/)")


if __name__ == "__main__":
    _self_check()
