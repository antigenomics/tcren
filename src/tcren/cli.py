"""Command-line interface for tcren.

Commands are grouped in ``tcren --help``:

Scoring & prediction
    * ``tcren score`` — end-to-end candidate-epitope scoring (drop-in for ``run_TCRen.R``).
    * ``tcren rank`` — percentile-rank a peptide's energy against a random pMHC background.
    * ``tcren ddg`` — ΔΔG of peptide mutations (fast virtual-matrix path; alanine scan / neoantigens).
    * ``tcren binder`` — TCR binder vs non-binder from AF-orthogonal interface geometry.
    * ``tcren energy`` — DOPE atom-level interface interaction energy (the ΔΔG ``e_native`` scorer).
    * ``tcren mechanics`` — interface mechanics (stiffness / rupture / coupling) — the koff proxies.
    * ``tcren pipeline`` — full pipeline → per-interface energies (TCRen + MJ) + total.

Annotation & contacts
    * ``tcren annotate`` — chain typing + region markup (TCR CDR/FR, MHC groove, peptide; ``--pseudo``).
    * ``tcren contacts`` — annotated residue-pair contact table for an interface.

Orientation & refinement
    * ``tcren superimpose`` — orient structure(s) onto the canonical database by MHC.
    * ``tcren refine`` — potential-guided peptide-pose refinement (DOPE MC; optional ``--substitute``).

Reference data & potentials
    * ``tcren orient`` — build a canonical database from native complexes.
    * ``tcren derive-potential`` — derive a TCRen potential from a contact-map table.
    * ``tcren fetch-data`` / ``fetch-recent`` — fetch reference sets / recent RCSB TCR-pMHC entries.
    * ``tcren build-mhc-ref`` — build the IMGT/HLA + mouse MHC allele reference.

Info
    * ``tcren info`` — version + dependency check.
    * ``tcren paper …`` — Nat Comput Sci 2022 reproduction helpers.
"""

from __future__ import annotations

from pathlib import Path

import polars as pl
import typer

from . import __version__
from .annotation import classify_chains
from .contactmap import TCR_REGIONS, ContactMap
from .potential import Potential, derive_tcren, derive_tcren_loo, keskin, mj, tcren
from .scoring import score_peptides
from .structure import iter_structures, parse_structure

app = typer.Typer(
    add_completion=True,  # `tcren --install-completion` for bash/zsh/fish; --show-completion to print
    help="Structure-based TCR–epitope recognition: score epitopes, rank binders, ΔΔG, "
         "orient/refine poses, and derive potentials from TCR:pMHC structures.",
)

# Help panels grouping the commands in `tcren --help` (rich_help_panel below).
_P_SCORE = "Scoring & prediction"
_P_ANNOT = "Annotation & contacts"
_P_ORIENT = "Orientation & refinement"
_P_DATA = "Reference data & potentials"
_P_INFO = "Info"
paper_app = typer.Typer(add_completion=False, help="Nat Comput Sci 2022 reproduction.")
app.add_typer(paper_app, name="paper")


@paper_app.command("bootstrap")
def paper_bootstrap(
    structures: bool = typer.Option(True, "--structures/--no-structures"),
    canonical: bool = typer.Option(False, "--canonical", help="also fetch the Canonical2026 set"),
) -> None:
    """Fetch HF structure sets into notebooks/data/<Set>/ (gitignored; non-structure inputs
    are already committed under natcompsci2022/data_legacy/)."""
    from .paper import bootstrap as run

    summary = run(structures=structures, canonical=canonical)
    for k, v in summary.items():
        typer.echo(f"{k}: {v}")

_BUNDLED_POTENTIALS = {"tcren": tcren, "mj": mj, "keskin": keskin}


def _load_potential(spec: str | None) -> Potential:
    """Resolve a potential from ``None`` (bundled TCRen), a bundled name, or a CSV path."""
    if spec is None:
        return tcren()
    if spec in _BUNDLED_POTENTIALS:
        return _BUNDLED_POTENTIALS[spec]()
    p = Path(spec)
    if p.exists():
        return Potential.from_csv(p)
    raise typer.BadParameter(
        f"potential not recognised: {spec!r} (use a bundled name "
        f"{sorted(_BUNDLED_POTENTIALS)} or an existing CSV path)"
    )


def _read_candidates(path: Path) -> list[str]:
    lines = [line.strip() for line in path.read_text().splitlines() if line.strip()]
    return [line for line in lines if line.lower() != "peptide"]


@app.command(rich_help_panel=_P_INFO)
def info() -> None:
    """Show version and dependency availability."""
    typer.echo(f"tcren {__version__}")
    try:
        import arda  # noqa: PLC0415

        arda_status = f"available ({Path(arda.__file__).parent})"
    except ImportError:
        arda_status = "NOT available — run: pip install arda-mapper (or bash setup.sh)"
    typer.echo(f"arda: {arda_status}")
    typer.echo(f"bundled TCRen potential: {tcren().matrix.height} pairs")


_REGION_CHAINS = {
    "tcr": {"TRA", "TRB", "TRD", "TRG"},
    "mhc": {"MHCa", "MHCb", "B2M"},
    "peptide": {"PEPTIDE"},
}


@app.command(rich_help_panel=_P_ANNOT)
def annotate(
    structures: Path = typer.Option(..., "-s", "--structures", help="structure file, directory, or .tar.gz (.pdb/.cif/.pdb.gz/.cif.gz)"),
    out: Path = typer.Option("markup.csv", "-o", "--out", help="output residue-markup CSV"),
    regions: str = typer.Option("all", "--regions", help="which chains to annotate: all|tcr|mhc|peptide"),
    pseudo: bool = typer.Option(False, "--pseudo", help="also mark MHC pseudosequence (MPS) residues"),
    organism: str = typer.Option("human", "--organism"),
) -> None:
    """Annotate chains and emit a per-residue region-markup table.

    Covers TCR (CDR/FR), MHC groove (helices/floor) and peptide in one pass — ``--regions``
    restricts the output to one chain class. ``--pseudo`` additionally marks the NetMHCpan MHC
    pseudosequence residues (region ``MPS``). MHC groove + ``MPS`` require MHC annotation, which
    runs automatically when needed.
    """
    from .contacts.table import residue_annotation

    if regions not in ("all", "tcr", "mhc", "peptide"):
        raise typer.BadParameter("--regions must be one of all|tcr|mhc|peptide")
    want_mhc = pseudo or regions in ("all", "mhc")
    keep = None if regions == "all" else _REGION_CHAINS[regions]

    frames = []
    for pid, s in iter_structures(structures, importer=parse_structure):
        classify_chains(s, organism=organism)
        if want_mhc:
            from .mhc import annotate_mhc
            annotate_mhc(s)
        if pseudo:
            from .mhc import annotate_pseudo
            annotate_pseudo(s)
        df = residue_annotation(s).with_columns(pl.lit(pid).alias("pdb.id"))
        if keep is not None:
            df = df.filter(pl.col("chain.type").is_in(list(keep)))
        frames.append(df)
    pl.concat(frames).write_csv(str(out))
    typer.echo(f"wrote {out}")


@app.command(rich_help_panel=_P_ANNOT)
def contacts(
    structures: Path = typer.Option(..., "-s", "--structures", help="structure file, directory, or .tar.gz (.pdb/.cif/.pdb.gz/.cif.gz)"),
    out: Path = typer.Option("contacts.csv", "-o", "--out"),
    cutoff: float = typer.Option(5.0, "--cutoff"),
    interface: str = typer.Option("tcr_peptide", "--interface", help="tcr_peptide|tcr_mhc|peptide_mhc|all"),
    regions: str = typer.Option("all", "--regions", help="TCR regions on the TCR side: all|cdr|cdr+fr (default: all)"),
    organism: str = typer.Option("human", "--organism"),
) -> None:
    """Compute and emit an annotated contact table."""
    if regions not in TCR_REGIONS:
        raise typer.BadParameter("--regions must be one of all|cdr|cdr+fr")
    frames = []
    for _pid, s in iter_structures(structures, importer=parse_structure):
        classify_chains(s, organism=organism)
        cm = ContactMap.from_structure(s, cutoff=cutoff)
        frames.append(
            cm.contacts if interface == "all" else cm.interface(interface, tcr_regions=regions)
        )
    pl.concat(frames).write_csv(str(out))
    typer.echo(f"wrote {out}")


@app.command(rich_help_panel=_P_DATA)
def orient(
    structures: Path = typer.Option(..., "-s", "--structures", help="PDB/CIF file or directory of native complexes"),
    out: Path = typer.Option("oriented", "-o", "--out", help="output dir for oriented structures"),
    metadata: Path = typer.Option("orient_metadata.csv", "--metadata"),
    organism: str = typer.Option("human", "--organism"),
    reference_id: str = typer.Option(None, "--reference", help="force a reference complex id"),
    force_pca: bool = typer.Option(False, "--force-pca", help="skip native superposition"),
    threads: int = typer.Option(None, "--threads", "-t", help="threads for alignment/IO (default: all cores)"),
    mmcif: bool = typer.Option(False, "--mmCIF", help="write mmCIF (.cif) instead of PDB"),
    compress: bool = typer.Option(False, "--compress", help="gzip the output (.gz)"),
) -> None:
    """Build a canonical database: orient native TCR-pMHC complexes into the common MHC frame.

    Derives the per-class canonical frame and writes every complex into it (A–E chains). This is
    how the bundled ``Canonical2026`` set is produced; use ``superimpose`` to bring a *new*
    structure into an existing canonical database.
    """
    from .orient import run_folder

    run_folder(structures, out, metadata=metadata, organism=organism,
               reference_id=reference_id, force_pca=force_pca, threads=threads,
               mmcif=mmcif, compress=compress)


@app.command(rich_help_panel=_P_ORIENT)
def superimpose(
    structures: str = typer.Option(..., "-s", "--structures", help="structure file, directory, .tar.gz, or a glob ('data/*.pdb')"),
    out: Path = typer.Option("superimposed", "-o", "--out", help="output directory, or a single structure file (one input)"),
    db: Path = typer.Option(None, "--db", help="canonical database dir (default: data/Canonical2026, fetched at install)"),
    organism: str = typer.Option("human", "--organism"),
    threads: int = typer.Option(None, "--threads", "-t", help="threads for the alignment/write (default: all cores)"),
    mmcif: bool = typer.Option(False, "--mmCIF", help="write mmCIF (.cif) instead of PDB"),
    compress: bool = typer.Option(False, "--compress", help="gzip the output (.gz)"),
) -> None:
    """Superimpose structure(s) onto a canonical database by MHC.

    Detects each input's MHC chains, class, and species, then superposes its conserved groove Cα
    onto *every* database structure of the same class and species and averages the transforms into
    one consensus placement. The database defaults to ``data/Canonical2026`` (populated at install).

    ``-s`` accepts a file, directory, ``.tar.gz``, or a shell glob. ``-o`` is an output directory,
    or — for a single input — a structure file whose extension must match ``--mmCIF``/``--compress``.
    Annotation is one batched mmseqs call; ``-t`` threads the alignment + write.
    """
    from .orient import run_superimpose

    try:
        run_superimpose(structures, out, db_dir=db, organism=organism,
                        threads=threads, mmcif=mmcif, compress=compress)
    except ValueError as exc:
        raise typer.BadParameter(str(exc)) from exc


@app.command("derive-potential", rich_help_panel=_P_DATA)
def derive_potential(
    contact_maps: Path | None = typer.Option(None, "-i", "--contact-maps", help="contact-map CSV"),
    out: Path = typer.Option("TCRen_potential.csv", "-o", "--out"),
    summary: Path | None = typer.Option(None, "--summary", help="summary CSV with a nonred flag"),
    nonred: bool = typer.Option(False, "--nonred", help="restrict to non-redundant structures"),
    structure_dir: Path | None = typer.Option(
        None, "--structure-dir",
        help="folder of PDBs to assemble contacts from (PDBs→contacts) when no -i CSV is given",
    ),
    redundancy_t: float | None = typer.Option(
        None, "--redundancy-t",
        help="non-redundancy clustering cutoff over αβ structures (off by default; "
             "requires markup, available with --structure-dir)",
    ),
    variant: str = typer.Option("classic", "--variant", help="classic|am"),
    pseudocount: int = typer.Option(1, "--pseudocount"),
    loo: bool = typer.Option(False, "--loo", help="emit leave-one-out potentials instead"),
) -> None:
    """Derive a TCRen potential from observed contacts.

    Provide contacts either as a precomputed ``-i`` CSV or as a ``--structure-dir`` of
    PDBs (assembled via ``annotate_structure_set``); pass exactly one. With a structure
    directory, ``--redundancy-t`` additionally restricts derivation to one representative
    per non-redundant cluster of αβ complexes (PDBs→contacts→cluster→derive in one call).
    """
    if (contact_maps is None) == (structure_dir is None):
        raise typer.BadParameter("pass exactly one of -i/--contact-maps or --structure-dir")

    markup = None
    if structure_dir is not None:
        from .paper import annotate_structure_set
        contacts, markup = annotate_structure_set(structure_dir)
    else:
        contacts = pl.read_csv(contact_maps)

    include = None
    if nonred:
        if summary is None:
            raise typer.BadParameter("--nonred requires --summary")
        include = pl.read_csv(summary).filter(pl.col("nonred"))["pdb.id"].to_list()
    if redundancy_t is not None:
        from .potential import alphabeta_ids, nonredundant_ids
        if markup is None:
            raise typer.BadParameter("--redundancy-t requires markup (use --structure-dir)")
        ab = alphabeta_ids(contacts)
        include = nonredundant_ids(markup.filter(pl.col("pdb.id").is_in(ab)), t=redundancy_t)

    if loo:
        ids = include or contacts["pdb.id"].unique().to_list()
        derive_tcren_loo(contacts, ids, variant=variant, pseudocount=pseudocount).write_csv(str(out))
    else:
        pot = derive_tcren(contacts, include=include, variant=variant, pseudocount=pseudocount)
        pot.to_csv(out)
    typer.echo(f"wrote {out}")


@app.command("fetch-data", rich_help_panel=_P_DATA)
def fetch_data(
    canonical: bool = typer.Option(True, "--canonical/--no-canonical", help="also fetch Canonical2026"),
) -> None:
    """Populate ``data/`` with the reference structure sets from the HF dataset.

    Run once at install (``setup.sh`` does this). Fetches ``Native2026`` (orientation
    references) and, by default, ``Canonical2026`` (the default ``superimpose`` database) into
    ``$TCREN_DATA_DIR`` / repo ``data/``. Skips folders already present.
    """
    from .paper.bootstrap import fetch_hf_structures
    from .paths import data_dir

    folders = ("Native2026",) + (("Canonical2026",) if canonical else ())
    summary = fetch_hf_structures(data_dir(), folders=folders)
    for k, v in summary.items():
        typer.echo(f"{k}: {v} structures")


@app.command("build-mhc-ref", rich_help_panel=_P_DATA)
def build_mhc_ref(
    species: str = typer.Option("human,mouse", "--species", help="comma-separated"),
    force_download: bool = typer.Option(False, "--force-download"),
) -> None:
    """Download and curate the MHC allele reference (IMGT/HLA + UniProt mouse)."""
    from .mhc import reference

    fasta = reference.build(
        species=tuple(s.strip() for s in species.split(",")), force_download=force_download
    )
    typer.echo(f"MHC reference written to {fasta}")


@app.command("fetch-recent", rich_help_panel=_P_DATA)
def fetch_recent(
    dest: Path = typer.Option(None, "--dest", help="output dir (default: data/pdb_recent)"),
    discover: bool = typer.Option(False, "--discover", help="also RCSB-search new TCR-pMHC entries"),
    after: str = typer.Option(None, "--after", help="discovery: release date >= YYYY-MM-DD"),
    organism: str = typer.Option("human", "--organism"),
) -> None:
    """Download recent TCR-pMHC structures from RCSB into data/pdb_recent.

    Seeds with the Native2026 ids; with --discover also full-text-searches RCSB for new
    entries. Each is pulled as mmCIF (.cif.gz; handles extended PDB ids), annotated, and kept
    only if it has all 5 required chains (MHCa + b2m/MHCb + peptide + TCR pair).
    """
    from .recent import discover_similar, fetch_ids, native2026_ids

    ids = native2026_ids()
    if discover:
        have = set(ids)
        ids = ids + [i for i in discover_similar(after_date=after) if i not in have]
    summary = fetch_ids(ids, dest=dest, organism=organism)
    for k, v in summary.items():
        typer.echo(f"{k}: {v}")


@app.command(rich_help_panel=_P_SCORE)
def score(
    structures: Path = typer.Option(..., "-s", "--structures", help="structure file, directory, or .tar.gz (.pdb/.cif/.pdb.gz/.cif.gz)"),
    candidates: Path = typer.Option(..., "-c", "--candidates", help="candidate epitopes file"),
    potential: str | None = typer.Option(None, "-p", "--potential", help="potential CSV (default: bundled TCRen)"),
    out: Path = typer.Option("candidate_epitopes_TCRen.csv", "-o", "--out"),
    interface: str = typer.Option("tcr_peptide", "--interface"),
    regions: str = typer.Option("all", "--regions", help="TCR regions on the TCR side: all|cdr|cdr+fr (default: all)"),
    organism: str = typer.Option("human", "--organism"),
    cutoff: float = typer.Option(5.0, "--cutoff"),
) -> None:
    """Score candidate epitopes against input structures (end-to-end pipeline)."""
    if regions not in TCR_REGIONS:
        raise typer.BadParameter("--regions must be one of all|cdr|cdr+fr")
    pot = _load_potential(potential)
    cands = _read_candidates(candidates)
    frames = []
    for _pid, s in iter_structures(structures, importer=parse_structure):
        classify_chains(s, organism=organism)
        cm = ContactMap.from_structure(s, cutoff=cutoff)
        frames.append(score_peptides(cm, cands, pot, interface=interface, tcr_regions=regions))
    result = pl.concat(frames) if frames else pl.DataFrame()
    result.write_csv(str(out))
    typer.echo(f"The ranked list of candidate epitopes can be found in {out}")


@app.command("ddg", rich_help_panel=_P_SCORE)
def ddg_cmd(
    structures: Path = typer.Option(..., "-s", "--structures", help="structure file, directory, or .tar.gz (.pdb/.cif/.pdb.gz/.cif.gz)"),
    native: str = typer.Option(..., "--native", help="native peptide sequence"),
    alanine_scan: bool = typer.Option(False, "--alanine-scan", help="ΔΔG of every position mutated to alanine"),
    mutant: list[str] = typer.Option(None, "--mutant", help="mutant peptide(s); repeat for several (neoantigen mode)"),
    potential: str | None = typer.Option(None, "-p", "--potential", help="potential CSV (default: bundled TCRen)"),
    out: Path = typer.Option("ddg.csv", "-o", "--out"),
    interface: str = typer.Option("tcr_peptide", "--interface", help="tcr_peptide|tcr_mhc|peptide_mhc"),
    regions: str = typer.Option("all", "--regions", help="TCR regions on the TCR side: all|cdr|cdr+fr (default: all)"),
    organism: str = typer.Option("human", "--organism"),
    cutoff: float = typer.Option(5.0, "--cutoff"),
) -> None:
    """ΔΔG of peptide mutations (fast virtual-matrix path; no atoms move).

    Re-scores the mutant sequence on the native contact map; ``ddG = E(native) - E(mutant)``
    (positive => destabilising). Use ``--alanine-scan`` for a per-position scan, or one or more
    ``--mutant`` for specific neoantigen substitutions.
    """
    if regions not in TCR_REGIONS:
        raise typer.BadParameter("--regions must be one of all|cdr|cdr+fr")
    if alanine_scan == bool(mutant):
        raise typer.BadParameter("pass exactly one of --alanine-scan or --mutant")
    from .ddg import alanine_scan as run_scan, neoantigen_ddg

    pot = _load_potential(potential)
    frames = []
    for _pid, s in iter_structures(structures, importer=parse_structure):
        classify_chains(s, organism=organism)
        cm = ContactMap.from_structure(s, cutoff=cutoff)
        if alanine_scan:
            df = run_scan(cm, native, pot, interface=interface, tcr_regions=regions)
        else:
            df = neoantigen_ddg(cm, native, mutant, pot, interface=interface, tcr_regions=regions)
        frames.append(df.with_columns(pl.lit(cm.pdb_id).alias("complex.id")))
    pl.concat(frames).write_csv(str(out))
    typer.echo(f"wrote {out}")


@app.command(rich_help_panel=_P_SCORE)
def rank(
    structures: Path = typer.Option(..., "-s", "--structures", help="structure file, directory, or .tar.gz (.pdb/.cif/.pdb.gz/.cif.gz)"),
    candidates: Path = typer.Option(None, "-c", "--candidates", help="peptides to rank; default: each structure's native peptide"),
    potential: str | None = typer.Option(None, "-p", "--potential", help="potential CSV (default: bundled TCRen)"),
    out: Path = typer.Option("rank.csv", "-o", "--out"),
    interface: str = typer.Option("tcr_peptide", "--interface", help="tcr_peptide|tcr_mhc|peptide_mhc"),
    regions: str = typer.Option("all", "--regions", help="TCR regions on the TCR side: all|cdr|cdr+fr (default: all)"),
    background: int = typer.Option(1000, "--background", help="number of random background peptides"),
    background_source: Path = typer.Option(None, "--background-source", help="FASTA/text of epitopes to sample the background from (default: uniform-random)"),
    seed: int = typer.Option(0, "--seed"),
    organism: str = typer.Option("human", "--organism"),
    cutoff: float = typer.Option(5.0, "--cutoff"),
) -> None:
    """Percentile-rank peptides' TCRen energy against a random pMHC background.

    For each structure, scores the supplied candidate peptides (or the structure's own
    peptide when ``-c`` is omitted) together with ``--background`` random peptides of the
    same length and reports ``rank_pct`` — the fraction of background scoring at least as
    well (lower energy = better binder, so a small ``rank_pct`` means a strong binder).
    """
    if regions not in TCR_REGIONS:
        raise typer.BadParameter("--regions must be one of all|cdr|cdr+fr")
    from .scoring_rank import percentile_rank

    from .structure.model import PEPTIDE_TYPE

    pot = _load_potential(potential)
    cands = _read_candidates(candidates) if candidates is not None else None
    src = str(background_source) if background_source is not None else None
    rows = []
    for _pid, s in iter_structures(structures, importer=parse_structure):
        classify_chains(s, organism=organism)
        cm = ContactMap.from_structure(s, cutoff=cutoff)
        if cands is not None:
            peptides = cands
        else:
            native = next((c.sequence() for c in s.chains if c.chain_type == PEPTIDE_TYPE), None)
            if native is None:
                raise typer.BadParameter(f"no peptide chain in {cm.pdb_id}; pass -c/--candidates")
            peptides = [native]
        for pep in peptides:
            bg = None
            if src is not None:
                from .scoring_rank import background_peptides
                bg = background_peptides(len(pep), n=background, seed=seed, source=src)
            res = percentile_rank(cm, pep, pot, interface=interface, n_background=background,
                                  seed=seed, tcr_regions=regions, background=bg)
            rows.append({"complex.id": cm.pdb_id, **res})
    pl.DataFrame(rows).write_csv(str(out))
    typer.echo(f"wrote {out}")


@app.command(rich_help_panel=_P_SCORE)
def binder(
    structures: Path = typer.Option(..., "-s", "--structures", help="TCR-pMHC model file, directory, or .tar.gz"),
    out: Path = typer.Option("binder.csv", "-o", "--out", help="per-structure descriptor + P(bind) table"),
    organism: str = typer.Option("human", "--organism"),
    cutoff: float = typer.Option(5.0, "--cutoff", help="contact cutoff (Å)"),
    features_only: bool = typer.Option(False, "--features-only", help="emit the 5 descriptors, skip P(bind)"),
) -> None:
    """Predict TCR binder/non-binder from AF-orthogonal interface geometry (native _geom + frozen model).

    Scores each complex from interface size, dual-chain balance, H-bonds, buried ΔSASA and the
    CDR1/2-vs-CDR3α TCRen potential — signal that beats AlphaFold/TCRmodel2 confidence for ranking
    candidate TCRs against a fixed pMHC. All descriptors are computed natively (no PyRosetta/Biopython
    SASA/sklearn). Low ``p_bind`` = unlikely binder.
    """
    from .binder import binder_features, binder_score

    rows = []
    for _pid, s in iter_structures(structures, importer=parse_structure):
        feats = binder_features(s, organism=organism, cutoff=cutoff)
        row = {"complex.id": s.pdb_id, **feats}
        if not features_only:
            row["p_bind"] = binder_score(feats)
        rows.append(row)
    pl.DataFrame(rows).write_csv(str(out))
    typer.echo(f"wrote {out}")


@app.command(rich_help_panel=_P_SCORE)
def pipeline(
    structures: Path = typer.Option(..., "-s", "--structures", help="structure file, directory, or .tar.gz"),
    out: Path = typer.Option("pipeline_scores.csv", "-o", "--out", help="per-structure interface-score table"),
    no_superimpose: bool = typer.Option(False, "--no-superimpose", help="skip canonical orientation"),
    db: Path = typer.Option(None, "--db", help="canonical database dir (default: data/Canonical2026)"),
    organism: str = typer.Option("human", "--organism"),
    cutoff: float = typer.Option(5.0, "--cutoff"),
    tcr_peptide_potential: str = typer.Option(None, "--tcr-peptide-potential", help="potential for the TCR↔peptide interface: bundled name (tcren|mj|keskin) or CSV path (default: tcren)"),
    tcr_mhc_potential: str = typer.Option(None, "--tcr-mhc-potential", help="potential for the TCR↔MHC interface: bundled name or CSV path (default: mj)"),
    peptide_mhc_potential: str = typer.Option(None, "--peptide-mhc-potential", help="potential for the peptide↔MHC interface: bundled name or CSV path (default: mj)"),
    regions: str = typer.Option("all", "--regions", help="TCR regions on the TCR side: all|cdr|cdr+fr (default: all)"),
) -> None:
    """Run the full pipeline and write per-interface energies for each structure.

    structure → annotate (alleles + chains) → superimpose → resmarkup / canonical Cα / contacts
    → score (TCRen for TCR↔peptide, MJ for TCR↔MHC and peptide↔MHC) + total.

    Each interface's potential can be overridden with a bundled name (``tcren``/``mj``/
    ``keskin``) or a CSV path; an unset option keeps the default family for that interface.
    """
    from .pipeline import run as run_pipeline, score_row

    if regions not in TCR_REGIONS:
        raise typer.BadParameter("--regions must be one of all|cdr|cdr+fr")
    potentials = {
        "tcr_peptide": tcr_peptide_potential,
        "tcr_mhc": tcr_mhc_potential,
        "peptide_mhc": peptide_mhc_potential,
    }
    rows = []
    for _pid, s in iter_structures(structures, importer=parse_structure):
        try:
            res = run_pipeline(s, organism=organism, superimpose=not no_superimpose,
                               db_dir=db, cutoff=cutoff,
                               potentials=potentials, tcr_regions=regions)
            rows.append(score_row(res))
        except Exception as exc:  # noqa: BLE001 - keep the batch resilient
            rows.append({"pdb.id": s.pdb_id, "total": None,
                         "error": f"{type(exc).__name__}: {str(exc)[:80]}"})
    pl.DataFrame(rows).write_csv(str(out))
    typer.echo(f"wrote {out}")


@app.command(rich_help_panel=_P_SCORE)
def energy(
    structures: str = typer.Option(..., "-s", "--structures", help="structure file, directory, .tar.gz, or glob"),
    out: Path = typer.Option("energy.csv", "-o", "--out", help="per-structure DOPE interface-energy table"),
    relax: bool = typer.Option(False, "--relax", help="also report the energy after DOPE refinement + the gap"),
    shell: float = typer.Option(12.0, "--shell", help="partner atoms within this many Å of the peptide (DOPE range)"),
    organism: str = typer.Option("human", "--organism"),
    seed: int = typer.Option(0, "--seed", help="refinement seed (with --relax)"),
) -> None:
    """DOPE atom-level interaction energy across the peptide↔partner interface (the ``_relax`` kernel).

    Sums the DOPE potential over peptide↔partner heavy-atom pairs — the interface ΔG contribution of the
    peptide (lower = more favourable). With ``--relax`` it also reports the energy after a rigid-body DOPE
    refinement (:func:`tcren.refine_peptide`) and the relaxation ``gap`` = e_native − e_relax. This is the
    single-structure scorer behind the ΔΔG benchmark (``e_native``/``e_relax``).
    """
    from .refine import interface_energy, refine_peptide
    from .structure.io import import_structure

    rows = []
    for pid, s in iter_structures(structures, importer=import_structure):
        try:
            classify_chains(s, organism=organism)
            row = {"pdb.id": pid, "e_native": interface_energy(s, shell=shell)}
            if relax:
                row["e_relax"] = interface_energy(refine_peptide(s, seed=seed)[0], shell=shell)
                row["gap"] = row["e_native"] - row["e_relax"]
            rows.append(row)
        except Exception as exc:  # noqa: BLE001 - keep the batch resilient
            rows.append({"pdb.id": pid, "e_native": None,
                         "error": f"{type(exc).__name__}: {str(exc)[:80]}"})
    pl.DataFrame(rows).write_csv(str(out))
    typer.echo(f"wrote {out}")


@app.command(rich_help_panel=_P_SCORE)
def mechanics(
    structures: str = typer.Option(..., "-s", "--structures", help="structure file, directory, .tar.gz, or glob"),
    out: Path = typer.Option("mechanics.csv", "-o", "--out", help="per-structure interface-mechanics table"),
    cutoff: float = typer.Option(8.0, "--cutoff", help="heavy-atom contact cutoff (Å) defining a spring"),
    weight: str = typer.Option("invdist2", "--weight", help="spring stiffness model: unit|count|invdist2"),
    direction: str = typer.Option("tensile", "--direction", help="rupture pull: tensile|shear|auto"),
    break_strain: float = typer.Option(0.5, "--break-strain", help="fractional extension at which a spring breaks"),
    organism: str = typer.Option("human", "--organism"),
) -> None:
    """Interface mechanics — the koff proxies: stiffness tensor + steered rupture + coupling residues.

    Treats the TCR↔pMHC contact map as a network of breakable springs and reports, per structure:
    ``n_spring``, ``S_tot``/``K_tens``/``K_shear``/``aniso`` (stiffness tensor), ``rupture_force``/
    ``rupture_work`` (steered unbinding), and ``couple_pep``/``couple_total`` (coupling residues).
    Validated on ATLAS: the tensile stiffness / rupture resistance track the dissociation off-rate
    (koff) far better than the equilibrium ΔG/Kd (Bell–Evans; the TCR is a mechanosensor).
    """
    from .mechanics import WEIGHTS, coupling_residues, rupture, stiffness_tensor

    if weight not in WEIGHTS:
        raise typer.BadParameter(f"--weight must be one of {'|'.join(WEIGHTS)}")
    if direction not in ("tensile", "shear", "auto"):
        raise typer.BadParameter("--direction must be one of tensile|shear|auto")
    rows = []
    for pid, s in iter_structures(structures, importer=parse_structure):
        try:
            classify_chains(s, organism=organism)
            row = {"pdb.id": pid, **stiffness_tensor(s, cutoff=cutoff, weight=weight)}
            row.update(rupture(s, direction=direction, cutoff=cutoff, weight=weight,
                               break_strain=break_strain))
            row.update(coupling_residues(s))
            rows.append(row)
        except Exception as exc:  # noqa: BLE001 - keep the batch resilient
            rows.append({"pdb.id": pid, "K_tens": None,
                         "error": f"{type(exc).__name__}: {str(exc)[:80]}"})
    pl.DataFrame(rows).write_csv(str(out))
    typer.echo(f"wrote {out}")


@app.command(rich_help_panel=_P_ORIENT)
def refine(
    structures: str = typer.Option(..., "-s", "--structures", help="structure file, directory, .tar.gz, or glob"),
    out: Path = typer.Option("refined", "-o", "--out", help="output directory for refined structures"),
    substitute: str = typer.Option(None, "--substitute", help="thread this peptide onto the backbone first"),
    organism: str = typer.Option("human", "--organism"),
    n_steps: int = typer.Option(2000, "--steps", help="Monte-Carlo steps"),
    restraint_w: float = typer.Option(0.5, "--restraint", help="harmonic restraint to the input pose"),
    seed: int = typer.Option(0, "--seed"),
    mmcif: bool = typer.Option(False, "--mmCIF", help="write mmCIF (.cif) instead of PDB"),
    compress: bool = typer.Option(False, "--compress", help="gzip the output (.gz)"),
) -> None:
    """Potential-guided rigid-body refinement of the peptide pose (knowledge-based, not physics).

    Optionally ``--substitute`` a new equal-length peptide first, then run a Monte-Carlo refinement
    scored by the DOPE atom-level statistical potential (restrained to the input pose; independent of
    the TCRen/MJ scoring potentials). Writes one structure per input and prints the final DOPE
    energy. (For physics-grade relaxation use Rosetta FlexPepDock externally.)
    """
    from .refine import refine_peptide, substitute_peptide
    from .structure.io import import_structure, structure_output_path, write_structure

    out.mkdir(parents=True, exist_ok=True)
    rows = []
    for pid, s in iter_structures(structures, importer=import_structure):
        try:
            classify_chains(s, organism=organism)
            if substitute:
                s = substitute_peptide(s, substitute)
            oriented, energy = refine_peptide(s, restraint_w=restraint_w, n_steps=n_steps, seed=seed)
            write_structure(oriented, structure_output_path(out, pid, mmcif, compress))
            rows.append({"pdb.id": pid, "energy": energy})
        except Exception as exc:  # noqa: BLE001 - keep the batch resilient
            rows.append({"pdb.id": pid, "energy": None,
                         "error": f"{type(exc).__name__}: {str(exc)[:80]}"})
    pl.DataFrame(rows).write_csv(str(out / "refine_energies.csv"))
    typer.echo(f"refined {sum(r.get('energy') is not None for r in rows)}/{len(rows)} -> {out}")


if __name__ == "__main__":
    app()
