"""Command-line interface for tcren.

Commands are grouped in ``tcren --help``:

Scoring & prediction
    * ``tcren score`` — end-to-end candidate-epitope scoring (drop-in for ``run_TCRen.R``).
    * ``tcren rank`` — percentile-rank a peptide's energy against a random pMHC background.
    * ``tcren ddg`` — ΔΔG of peptide mutations (fast virtual-matrix path; alanine scan / neoantigens).
    * ``tcren binder`` — TCR binder vs non-binder from AF-orthogonal interface geometry.
    * ``tcren energy`` — DOPE atom-level interface interaction energy (the ΔΔG ``e_native`` scorer).
    * ``tcren mechanics`` — interface mechanics (stiffness / rupture / coupling) — the koff proxies.
    * ``tcren scoring`` — per-interface contact energies Φ (``--delta`` for ΔΦ, ``--geometry`` for Q).
    * ``tcren surface`` — pMHC surface topology: height/hydropathy/charge maps + epitope comparison.
    * ``tcren cpl`` — combinatorial-peptide-library response matrix from one template structure.
    * ``tcren recognize`` — every interface descriptor + joint P(real), one row per structure.

Annotation & contacts
    * ``tcren annotate`` — chain typing + region markup (TCR CDR/FR, MHC groove, peptide; ``--pseudo``).
    * ``tcren contacts`` — annotated residue-pair contact table for an interface.

Orientation & refinement
    * ``tcren superimpose`` — orient structure(s) onto the canonical database by MHC.
    * ``tcren refine`` — potential-guided peptide-pose refinement (DOPE MC; optional ``--substitute``).
    * ``tcren substitute-tcr`` — graft a donor TCR onto a host pMHC (a chimeric complex).

Reference data & potentials
    * ``tcren orient`` — build a canonical database from native complexes.
    * ``tcren shuffle`` — wrong-TCR-on-real-pMHC decoys, the negatives for a recognition model.
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
from .potential import Potential, derive_tcren, derive_tcren_loo, tcren, tcren2
from .scoring import score_peptides
from .structure import iter_structures, parse_structure

app = typer.Typer(
    add_completion=True,  # `tcren --install-completion` for bash/zsh; --show-completion to print
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

from .pipeline import _BUNDLED_POTENTIALS  # the one map; see pipeline.py


def _load_potential(spec: str | None) -> Potential:
    """Resolve a potential from ``None`` (TCRen2, the default), a bundled name, or a CSV path."""
    if spec is None:
        return tcren2()
    if spec in _BUNDLED_POTENTIALS:
        return _BUNDLED_POTENTIALS[spec]()
    p = Path(spec)
    if p.exists():
        return Potential.from_csv(p)
    raise typer.BadParameter(
        f"potential not recognised: {spec!r} (use a bundled name "
        f"{sorted(_BUNDLED_POTENTIALS)} or an existing CSV path)"
    )


def _iter_typed(structures: Path, organism: str = "human"):
    """Yield chain-typed structures, annotating a whole set in one batch.

    ``classify_chains`` per structure spawns one ``mmseqs easy-search`` per structure, each
    building a temporary database from a handful of query sequences; the process startup dominates
    and the cost is roughly an order of magnitude. :func:`~tcren.paper.iter_annotated_set` sends
    every chain of every structure in a single call per organism, which is why anything resolving
    to more than one structure -- a directory, a glob, a manifest -- goes through it. A single
    file has nothing to batch and takes the direct path.
    """
    from .structure import structure_paths
    try:
        many = len(structure_paths(structures)) > 1
    except Exception:  # .tar.gz and anything structure_paths cannot enumerate
        many = False
    if many:
        from .paper.helpers import iter_annotated_set
        yield from iter_annotated_set(structures)
        return
    for _pid, s in iter_structures(structures, importer=parse_structure):
        classify_chains(s, organism=organism)
        yield s


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
    typer.echo(f"default potential TCRen2: {tcren2().matrix.height} pairs "
               f"(karnaukhov2022: {tcren().matrix.height})")


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
    for s in _iter_typed(structures, organism):
        pid = s.pdb_id
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
    for s in _iter_typed(structures, organism):
        cm = ContactMap.from_structure(s, cutoff=cutoff)
        frames.append(
            (cm.contacts if interface == "all" else cm.interface(interface, tcr_regions=regions))
            .with_columns(pl.lit(s.pdb_id).alias("pdb.id"))
        )
    pl.concat(frames).write_csv(str(out))
    typer.echo(f"wrote {out}")


@app.command(rich_help_panel=_P_DATA)
def orient(
    structures: Path = typer.Option(..., "-s", "--structures", help="PDB/CIF file or directory of native complexes"),
    out: Path = typer.Option("oriented", "-o", "--out", help="output dir for oriented structures"),
    metadata: Path = typer.Option(None, "--metadata", help="metadata table (default: <out>/orient_metadata.json, the format `superimpose` reads; .csv writes CSV)"),
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


@app.command(rich_help_panel=_P_DATA)
def shuffle(
    structures: Path = typer.Option(..., "-s", "--structures", help="dir of ORIENTED TCR-pMHC (co-framed; run `orient`/`superimpose` first)"),
    out: Path = typer.Option("shuffled", "-o", "--out", help="output dir for decoy complexes"),
    n: int = typer.Option(10, "--n", help="decoys generated per pMHC"),
    seed: int = typer.Option(0, "--seed"),
    within_class: bool = typer.Option(True, "--within-class/--any-class", help="graft only same-MHC-class TCRs"),
    organism: str = typer.Option("human", "--organism"),
    compress: bool = typer.Option(False, "--compress", help="gzip the output (.gz)"),
) -> None:
    """Generate wrong-TCR-on-real-pMHC decoys (a Shuffled set) for recognition models.

    Keeps each oriented complex's pMHC intact and grafts on ``n`` different complexes' TCRs (a within-MHC-class
    derangement, so no decoy reproduces a real pairing). Real (label 1) vs these decoys (label 0) trains a
    label-free TCR-recognition classifier. Inputs must be co-framed — run ``tcren orient`` first.
    """
    from .shuffle import run_shuffle

    written = run_shuffle(structures, out, n=n, seed=seed, within_class=within_class,
                          organism=organism, compress=compress)
    typer.echo(f"wrote {written} decoys to {out}")


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
    nonred: bool = typer.Option(False, "--nonred", help="NOT USED FOR TCRen2. restrict to non-redundant structures"),
    structure_dir: Path | None = typer.Option(
        None, "--structure-dir",
        help="folder of PDBs to assemble contacts from (PDBs→contacts) when no -i CSV is given",
    ),
    redundancy_t: float | None = typer.Option(
        None, "--redundancy-t",
        help="NOT USED FOR TCRen2. non-redundancy clustering cutoff over αβ structures (off by default; "
             "requires markup, available with --structure-dir)",
    ),
    variant: str = typer.Option("classic", "--variant", help="classic|am"),
    pseudocount: int = typer.Option(1, "--pseudocount"),
    balance: str | None = typer.Option(
        None, "--balance",
        help="down-weight PDB redundancy: epitope|tcr|both. Each structure is weighted by "
             "the mean of 1/(structures sharing its value) over the chosen axes, so a "
             "re-solved complex counts once while a novel receptor on a common epitope "
             "still counts (requires --structure-dir)",
    ),
    loo: bool = typer.Option(False, "--loo", help="NOT USED FOR TCRen2. emit leave-one-out potentials instead"),
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

    if markup is not None:
        # HARD RULE: tcren is for alpha-beta TCR : peptide-MHC and nothing else. A structure
        # missing either CDR3 or the peptide is not in scope, and it is not merely useless --
        # `--balance` skips a structure with a null on any axis, and `derive_tcren` then defaults
        # it to weight 1.0, the maximum. On Native2026 that admitted 12 files (3 pMHC-only, 8
        # single-chain, one gamma-delta), three near-duplicate pairs among them at full weight
        # each. There is no flag to turn this off.
        before = markup.height
        markup = markup.filter(
            pl.col("cdr3a").is_not_null() & pl.col("cdr3b").is_not_null()
            & pl.col("peptide").is_not_null()
        )
        if markup.height < before:
            typer.echo(f"alpha-beta TCR:pMHC only: {markup.height} of {before} structures kept")
        keep = set(markup["pdb.id"])
        include = sorted(keep) if include is None else [i for i in include if i in keep]

    weights = None
    if balance is not None:
        from .potential import balanced_weights
        axes = {"epitope": (("peptide",),),
                "tcr": (("cdr3a", "cdr3b"),),
                "both": (("peptide",), ("cdr3a", "cdr3b"))}.get(balance)
        if axes is None:
            raise typer.BadParameter("--balance must be one of: epitope, tcr, both")
        if markup is None:
            raise typer.BadParameter("--balance requires markup (use --structure-dir)")
        weights = balanced_weights(markup, axes=axes)

    if loo:
        ids = include or contacts["pdb.id"].unique().to_list()
        derive_tcren_loo(contacts, ids, variant=variant, pseudocount=pseudocount).write_csv(str(out))
    else:
        pot = derive_tcren(contacts, include=include, variant=variant,
                           pseudocount=pseudocount, weights=weights)
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


def _soft_score(structure, candidates, potential, interface, cutoff, contact_probabilities):
    """Score candidates over rotamer-averaged contact probabilities instead of a hard contact map.

    Built as its own map because the probabilities carry pairs the input pose never contacted, which
    a ContactMap by construction cannot hold.
    """
    import numpy as np

    p = contact_probabilities(structure, interface, cutoff=cutoff)
    matrix, index = potential.as_matrix()
    pep = next(c for c in structure.chains if c.chain_type == "PEPTIDE")
    pos = {r.seq_index: i for i, r in enumerate(pep.residues)}
    to_pos = np.array([pos.get(i, -1) for i in p["residue.index.to"].to_list()], dtype=np.int64)
    fixed = np.array([index.get(a, -1) for a in p["residue.aa.from"].to_list()], dtype=np.int64)
    w = p["p"].to_numpy()

    rows = []
    for peptide in candidates:
        if len(peptide) != len(pep.residues):
            continue
        sub = np.array([index.get(peptide[q], -1) if 0 <= q < len(peptide) else -1 for q in to_pos],
                       dtype=np.int64)
        ok = (fixed >= 0) & (sub >= 0)
        rows.append({"complex.id": structure.pdb_id, "peptide": peptide,
                     "score": float(np.nansum(matrix[fixed[ok], sub[ok]] * w[ok]))})
    return (pl.DataFrame(rows, schema={"complex.id": pl.Utf8, "peptide": pl.Utf8,
                                       "score": pl.Float64})
            .with_columns(pl.lit(potential.name).alias("potential"))
            .select("complex.id", "peptide", "potential", "score").sort("complex.id", "score"))


def _score_weights(structure, cm, interface, regions, drop_untyped, position_scheme):
    """Combine the optional contact-type and peptide-position weights, or None when both are off."""
    import numpy as np

    if not drop_untyped and position_scheme == "uniform":
        return None
    w = np.ones(cm.interface(interface, tcr_regions=regions).height, dtype=np.float64)
    if drop_untyped:
        from .contact_types import classify_contacts, stacked_pairs, type_weights
        typed = classify_contacts(cm.interface(interface, tcr_regions=regions), "v2",
                                  stacked_pairs(structure))
        w *= type_weights(typed)
    if position_scheme != "uniform":
        from .scoring import peptide_positions, position_weights
        w *= position_weights(peptide_positions(cm, structure, interface, regions), position_scheme)
    return w


@app.command(rich_help_panel=_P_SCORE)
def score(
    structures: Path = typer.Option(..., "-s", "--structures", help="structure file, directory, or .tar.gz (.pdb/.cif/.pdb.gz/.cif.gz)"),
    candidates: Path = typer.Option(..., "-c", "--candidates", help="candidate epitopes file"),
    potential: str | None = typer.Option(None, "-p", "--potential", help="potential: bundled name (tcren2|karnaukhov2022|mj|keskin) or CSV path (default: tcren2)"),
    out: Path = typer.Option("candidate_epitopes_TCRen.csv", "-o", "--out"),
    interface: str = typer.Option("tcr_peptide", "--interface"),
    regions: str = typer.Option("all", "--regions", help="TCR regions on the TCR side: all|cdr|cdr+fr (default: all)"),
    organism: str = typer.Option("human", "--organism"),
    cutoff: float = typer.Option(5.0, "--cutoff"),
    intra_weight: float = typer.Option(0.0, "--intra-weight", help="weight of the intra-peptide term: score = interface energy + w x the candidate's contact energy with itself (0 = off, the default)"),
    drop_untyped: bool = typer.Option(False, "--drop-untyped", help="ignore contacts that are only proximity (no salt bridge / h-bond / stacking / hydrophobic / polar chemistry)"),
    position_scheme: str = typer.Option("uniform", "--position-weights", help="weight contacts by where they sit on the peptide: uniform|central|tcr_facing"),
    soft: bool = typer.Option(False, "--soft", help="rotamer-averaged contacts: replace the hard 5 A cutoff with a Boltzmann-weighted contact probability over chi rotamers"),
) -> None:
    """Score candidate epitopes against input structures (end-to-end pipeline).

    ``--intra-weight`` adds the term the interface sum omits: a candidate threaded onto the
    template's peptide conformation also pays for the contacts that conformation makes it have
    with **itself** (5 Å, sequence separation >= 3, MJ). It is off by default; a class-I 9-mer is
    extended and makes zero to two such contacts, so it separates candidates only where the peptide
    is genuinely bulged or self-packed.

    ``--drop-untyped``, ``--position-weights`` and ``--soft`` all reweight the same sum, and all
    default to off so the score is unchanged unless asked. The first uses the chemical typing to
    ignore pairs that are within 5 Å but make no interaction; the second says a contact under the
    CDR3 loops in the middle of the peptide is not worth the same as one at an anchor the TCR never
    touches; the third replaces the hard cutoff with a contact *probability* averaged over side-chain
    rotamers, which is what stops a single wrong χ1 from moving the energy by more than the energy
    itself (measured ``|ΔΦ|`` 0.524 → 0.054 under a deliberately wrong rotamer).
    """
    if regions not in TCR_REGIONS:
        raise typer.BadParameter("--regions must be one of all|cdr|cdr+fr")
    from .scoring import POSITION_SCHEMES
    if position_scheme not in POSITION_SCHEMES:
        raise typer.BadParameter(f"--position-weights must be one of {'|'.join(POSITION_SCHEMES)}")
    pot = _load_potential(potential)
    cands = _read_candidates(candidates)
    frames = []
    for s in _iter_typed(structures, organism):
        cm = ContactMap.from_structure(s, cutoff=cutoff, peptide_internal=bool(intra_weight))
        if soft:
            from .mhc import annotate_mhc
            from .rotamers import contact_probabilities
            annotate_mhc(s)
            frames.append(_soft_score(s, cands, pot, interface, cutoff, contact_probabilities))
            continue
        w = _score_weights(s, cm, interface, regions, drop_untyped, position_scheme)
        frames.append(score_peptides(cm, cands, pot, interface=interface, tcr_regions=regions,
                                     intra_weight=intra_weight, weights=w))
    result = pl.concat(frames) if frames else pl.DataFrame()
    result.write_csv(str(out))
    typer.echo(f"The ranked list of candidate epitopes can be found in {out}")


@app.command("ddg", rich_help_panel=_P_SCORE)
def ddg_cmd(
    structures: Path = typer.Option(..., "-s", "--structures", help="structure file, directory, or .tar.gz (.pdb/.cif/.pdb.gz/.cif.gz)"),
    native: str = typer.Option(..., "--native", help="native peptide sequence"),
    alanine_scan: bool = typer.Option(False, "--alanine-scan", help="ΔΔG of every position mutated to alanine"),
    mutant: list[str] = typer.Option(None, "--mutant", help="mutant peptide(s); repeat for several (neoantigen mode)"),
    potential: str | None = typer.Option(None, "-p", "--potential", help="potential: bundled name (tcren2|karnaukhov2022|mj|keskin) or CSV path (default: tcren2)"),
    out: Path = typer.Option("ddg.csv", "-o", "--out"),
    interface: str = typer.Option("tcr_peptide", "--interface", help="tcr_peptide|tcr_mhc|peptide_mhc"),
    regions: str = typer.Option("all", "--regions", help="TCR regions on the TCR side: all|cdr|cdr+fr (default: all)"),
    organism: str = typer.Option("human", "--organism"),
    cutoff: float = typer.Option(5.0, "--cutoff"),
) -> None:
    """ΔΔG of peptide mutations (fast virtual-matrix path; no atoms move).

    Re-scores the mutant sequence on the native contact map; ``ddG = E(native) - E(mutant)``
    (positive => STABILISING: the mutant scores lower, which is the better binder). Use
    ``--alanine-scan`` for a per-position scan, or one or more
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


@app.command("cpl", rich_help_panel=_P_SCORE)
def cpl_cmd(
    structures: Path = typer.Option(..., "-s", "--structures", help="template TCR:pMHC structure, directory, or .tar.gz (.pdb/.cif/.pdb.gz/.cif.gz)"),
    peptide: str = typer.Option(None, "--peptide", help="peptide to thread substitutions off (default: the template's own)"),
    position: int = typer.Option(None, "--position", help="1-based peptide position; restrict the output to this position"),
    mutation: str = typer.Option(None, "--mutation", help="one-letter residue; with --position, report just that cell"),
    to_mixture: bool = typer.Option(False, "--to-mixture", help="report the cost of giving a position up to the equimolar 1/20 mixture"),
    reference: str = typer.Option("both", "--reference", help="equimolar|wild_type|both (default: both)"),
    potential: str | None = typer.Option(None, "-p", "--potential", help="TCR:peptide potential: bundled name or CSV path (default: tcren2)"),
    mhc_potential: str | None = typer.Option(None, "--mhc-potential", help="peptide:MHC potential (default: Miyazawa-Jernigan)"),
    out: Path = typer.Option("cpl_matrix.csv", "-o", "--out"),
    regions: str = typer.Option("all", "--regions", help="TCR regions on the TCR side: all|cdr|cdr+fr (default: all)"),
    organism: str = typer.Option("human", "--organism"),
    cutoff: float = typer.Option(5.0, "--cutoff"),
) -> None:
    """Predict a combinatorial-peptide-library response matrix from a template TCR:pMHC structure.

    One row per (peptide position, amino acid) cell. Every cell carries BOTH peptide-bearing
    interfaces summed -- TCRen over TCR:peptide plus Miyazawa-Jernigan over peptide:MHC -- because
    the assay reads activation, which needs the peptide presented as well as the receptor engaged.

    Two reference states are reported, and a cell is only meaningful against one of them:
    ``effect_equimolar`` scores a residue against the 1/20 mixture, which is the CPL background and
    the right axis to compare with a measured matrix; ``effect_wild_type`` scores it against the
    residue the template carries, which is the mutation-scan / neoantigen question. Positive is
    favourable on both.

    Three narrower questions, all from the same matrix:

    \b
      --position 5                 every substitution at position 5, best first
      --position 5 --mutation W    just that one cell
      --position 5 --to-mixture    the cost of giving position 5 up to the 1/20 mixture
    """
    if reference not in ("equimolar", "wild_type", "both"):
        raise typer.BadParameter("--reference must be equimolar|wild_type|both")
    if regions not in TCR_REGIONS:
        raise typer.BadParameter("--regions must be one of all|cdr|cdr+fr")
    if (mutation or to_mixture) and position is None:
        raise typer.BadParameter("--mutation and --to-mixture both need --position")
    if mutation and to_mixture:
        raise typer.BadParameter("pass at most one of --mutation or --to-mixture")
    from .cpl import equimolar_effect, mutation_effect, position_scan, response_matrix
    from .mhc import annotate_mhc
    from .potential import mj

    pot = _load_potential(potential)
    mhc_pot = _load_potential(mhc_potential) if mhc_potential else mj()
    frames = []
    for _pid, s in iter_structures(structures, importer=parse_structure):
        classify_chains(s, organism=organism)
        # the peptide:MHC interface comes out EMPTY without this second pass, which would silently
        # reduce every anchor cell to zero rather than failing
        annotate_mhc(s)
        cm = ContactMap.from_structure(s, cutoff=cutoff)
        rm = response_matrix(cm, peptide, tcr_potential=pot, mhc_potential=mhc_pot,
                             tcr_regions=regions)
        if to_mixture:
            df = pl.DataFrame({"pos": [position], "aa": [rm.wild_type_at(position)],
                               "effect_equimolar": [equimolar_effect(rm, position)]})
        elif mutation:
            df = pl.DataFrame({"pos": [position], "aa": [mutation.upper()], **{
                f"effect_{r}": [mutation_effect(rm, position, mutation, reference=r)]
                for r in (("equimolar", "wild_type") if reference == "both" else (reference,))}})
        elif position is not None:
            df = position_scan(rm, position,
                               reference="equimolar" if reference == "both" else reference)
        else:
            df = rm.to_frame(None if reference == "both" else reference)
        frames.append(df.with_columns(pl.lit(rm.structure_id).alias("complex.id")))
    pl.concat(frames, how="diagonal_relaxed").write_csv(str(out))
    typer.echo(f"wrote {out}")


@app.command(rich_help_panel=_P_SCORE)
def rank(
    structures: Path = typer.Option(..., "-s", "--structures", help="structure file, directory, or .tar.gz (.pdb/.cif/.pdb.gz/.cif.gz)"),
    candidates: Path = typer.Option(None, "-c", "--candidates", help="peptides to rank; default: each structure's native peptide"),
    potential: str | None = typer.Option(None, "-p", "--potential", help="potential: bundled name (tcren2|karnaukhov2022|mj|keskin) or CSV path (default: tcren2)"),
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
def surface(
    structures: Path = typer.Option(..., "-s", "--structures", help="pMHC / TCR-pMHC file, directory, or .tar.gz"),
    out: Path = typer.Option("surface.csv", "-o", "--out", help="per-structure surface-topology descriptors"),
    organism: str = typer.Option("human", "--organism"),
    grid: str = typer.Option("64x32", "--grid", help="map cells as <n_y>x<n_x> (along groove x across)"),
    scale: str = typer.Option("kd", "--scale", help="hydropathy scale: kd (Kyte-Doolittle) or mj"),
    channel: str = typer.Option("h", "--channel", help="channel for --svg and --compare: h, phobic, charge"),
    region: str = typer.Option(None, "--region", help="restrict --compare to one source, e.g. peptide"),
    compare: Path = typer.Option(None, "--compare", help="also write the pairwise map-distance matrix here"),
    cells: Path = typer.Option(None, "--cells", help="also write the long per-cell table here"),
    svg: Path = typer.Option(None, "--svg", help="directory to write one SVG map per structure"),
) -> None:
    """Map the pMHC surface a TCR sees — height + hydropathy + charge over the groove.

    Emits, per structure, the scalars that say how *featured* the presented surface is: ``relief``
    (height spread over the peptide's footprint), ``peak_to_valley``, ``frac_above_ridge`` (how much
    peptide surface clears the MHC helix crests) and the mean/central hydropathy. A flat,
    MHC-dominated landscape — a "featureless" epitope — scores low on all of them.

    The groove frame is refit from each structure, so maps are comparable without prealigning the
    inputs: ``--compare`` writes the pairwise Manhattan map distance, which clusters structures of
    the same epitope together.
    """
    from .surface import surface_distance, surface_map, surface_table

    try:
        n_y, n_x = (int(v) for v in grid.lower().split("x"))
    except ValueError as exc:
        raise typer.BadParameter(f"--grid must look like 64x32, got {grid!r}") from exc

    maps, rows = [], []
    for _pid, s in iter_structures(structures, importer=parse_structure):
        try:
            if all(c.chain_type is None for c in s.chains):
                classify_chains(s, organism=organism, autodetect_species=True)
            from .mhc import annotate_mhc
            annotate_mhc(s)
            maps.append(surface_map(s, grid=(n_y, n_x), scale=scale))
        except Exception as exc:  # noqa: BLE001 - keep the batch resilient, report per structure
            rows.append({"structure.id": s.pdb_id, "error": f"{type(exc).__name__}: {exc}"})

    table = surface_table(maps)
    if rows:
        table = pl.concat([table, pl.DataFrame(rows)], how="diagonal")
    table.write_csv(str(out))
    typer.echo(f"wrote {out} ({len(maps)} mapped, {len(rows)} failed)")

    if cells is not None and maps:
        pl.concat([m.to_frame() for m in maps], how="vertical").write_csv(str(cells))
        typer.echo(f"wrote {cells}")
    if compare is not None and maps:
        ids, dist = surface_distance(maps, channel=channel, region=region)
        pl.DataFrame({"structure.id": ids, **{i: dist[:, k] for k, i in enumerate(ids)}}
                     ).write_csv(str(compare))
        typer.echo(f"wrote {compare}")
    if svg is not None and maps:
        from .viz.surface2d import render_surface_map
        svg.mkdir(parents=True, exist_ok=True)
        for m in maps:
            (svg / f"{m.structure_id}_{channel}.svg").write_text(render_surface_map(m, channel))
        typer.echo(f"wrote {len(maps)} SVG maps to {svg}")


@app.command(rich_help_panel=_P_SCORE)
def recognize(
    structures: str = typer.Option(None, "-s", "--structures", help="TCR-pMHC structure file, directory, .tar.gz, or glob"),
    features_table: Path = typer.Option(None, "--features", help="score a table already written by `tcren features` instead of re-reading the structures; emits Q, P_native and the channel posteriors G/T/E"),
    out: Path = typer.Option("recognize.tsv", "-o", "--out", help="per-structure descriptors + P(real) table (TSV)"),
    organism: str = typer.Option("human", "--organism"),
    features_only: bool = typer.Option(False, "--features-only", help="emit the descriptors, skip P(real)"),
    full: bool = typer.Option(False, "--full", help="add the 18 CDR3-local frame descriptors (the FramePose strain layer) and the intra-peptide term F_pep_int/n_pep_int"),
    scores: bool = typer.Option(False, "--scores", help="LEGACY (v1 reproduction): the fit-free q_bind + s_strain and the fitted p_bind + p_forced; implies --full. Prefer --features for Q/P_native"),
    mechanics: bool = typer.Option(False, "--mechanics", help="also append the koff proxies from `tcren mechanics` (stiffness tensor, steered rupture, coupling residues) — same table, no second annotation pass"),
    threads: int = typer.Option(1, "-t", "--threads", help="concurrent annotation batches for a multi-structure run (0 = all cores); cohort scores stay computed over the whole set"),
    autodetect_species: bool = typer.Option(True, "--autodetect-species/--no-autodetect-species", help="also search mouse to catch a mis-declared organism; --no- halves the annotation cost"),
) -> None:
    """Full interface descriptor table + joint P(real) for each TCR-pMHC complex (one TSV row per PDB).

    One row per structure with the complete recognition feature set
    (``tcren.recognition.RECOGNITION_FEATURES``): docking geometry (pitch, crossing, the 6
    TCRdock rigid-body params), per-interface energies ``F_{tcr_pep,tcr_mhc,pep_mhc}`` and
    poly-alanine ``dF``, CDR-loop energies ``F_{cdr12,cdr3a,cdr3b}``, contact-type tallies, ΔSASA ``burial`` and the MHC-class
    indicator — plus ``p_real`` (the distribution-aware Bayesian logistic) and ``p_real_bn`` (the
    Gaussian BN): the joint probability the complex is a genuine recognition interface rather than a
    wrong-TCR shuffle. ``--full`` also emits the 18 CDR3-local frame descriptors (the FramePose strain
    layer) and the intra-peptide term ``F_pep_int``/``n_pep_int`` — the peptide's contact energy with
    **itself**, which the three interface energies omit. ``--scores`` is kept for v1 reproduction:
    it adds the fit-free ``q_bind`` (binder-ID; the directional-decorrelated interface-quality
    score, calibrated on the native crystal reference so it is defined per structure and
    transfers) and ``s_strain`` (forced-pose), alongside the fitted ``p_bind`` / ``p_forced``. The
    recommended scores are ``Q`` and ``P_native``, which come from ``--features``.
    ``--features-only`` skips the models. Output is TSV.

    ``--mechanics`` appends the koff proxies ``tcren mechanics`` reports — stiffness tensor, steered
    rupture, coupling residues — to these same rows. Prefer it to running the two commands: they
    need the identical annotated structure, so the second command repeats the parse and both mmseqs
    searches to produce a second table in a different format (CSV) under a different key
    (``pdb.id``) that then has to be joined. One flag costs about a sixth of the descriptor pass and
    returns one table.

    Complementary scorer on the same inputs: ``tcren ddg`` (per-mutation alanine/neoantigen ΔΔF).

    Examples::

        tcren features  -s models/ -o feats.tsv                      # the descriptor pass, once
        tcren recognize --features feats.tsv -o scores.tsv           # Q, P_native and its channels
        tcren recognize -s models/ -o out.tsv                        # descriptors + p_real, no feature file
        tcren recognize -s models/ --mechanics -t 0 -o out.tsv       # + the spring-network terms
    """
    from .recognition import recognition_table
    from .structure.io import import_structure

    if (structures is None) == (features_table is None):
        raise typer.BadParameter("pass exactly one of -s/--structures or --features")
    if features_table is not None:
        _score_feature_table(features_table, out)
        return

    full = full or scores                                             # p_forced needs the CDR3-frame feats
    items = list(iter_structures(structures, importer=import_structure))
    import os as _os
    rows = recognition_table(items, organism=organism, full=full, scores=scores,
                             with_p_real=not features_only, autodetect_species=autodetect_species,
                             mechanics=mechanics,
                             threads=threads if threads > 0 else (_os.cpu_count() or 1))
    table = pl.DataFrame(rows)
    table.write_csv(str(out), separator="\t")
    typer.echo(f"wrote {out} ({len(rows)} rows)")


@app.command(rich_help_panel=_P_SCORE)
def scoring(
    structures: list[str] = typer.Option(..., "-s", "--structures", help="structure file(s), directory, .tar.gz, glob, or a .txt manifest of paths; repeatable and comma-separable"),
    out: Path = typer.Option("scores.csv", "-o", "--out", help="per-structure interface-score table"),
    no_superimpose: bool = typer.Option(False, "--no-superimpose", help="skip canonical orientation"),
    db: Path = typer.Option(None, "--db", help="canonical database dir (default: data/Canonical2026)"),
    organism: str = typer.Option("human", "--organism"),
    cutoff: float = typer.Option(5.0, "--cutoff", help="heavy-atom contact distance threshold (Å)"),
    tcr_peptide_potential: str = typer.Option(None, "--tcr-peptide-potential", help="potential for the TCR↔peptide interface: bundled name (tcren2|karnaukhov2022|mj|keskin) or CSV path (default: tcren2)"),
    tcr_mhc_potential: str = typer.Option(None, "--tcr-mhc-potential", help="potential for the TCR↔MHC interface: bundled name or CSV path (default: mj)"),
    peptide_mhc_potential: str = typer.Option(None, "--peptide-mhc-potential", help="potential for the peptide↔MHC interface: bundled name or CSV path (default: mj)"),
    regions: str = typer.Option("all", "--regions", help="TCR regions on the TCR side: all|cdr|cdr+fr (default: all)"),
    contact_weight: str = typer.Option("residue", "--contact-weight", help="residue (default, one per contacting pair) or atomic (weight by heavy-atom-pair count)"),
    intra_weight: float = typer.Option(0.0, "--intra-weight", help="weight of the intra-peptide term: report F_pep_int and add w x it to F_total (0 = off, the default)"),
    delta: bool = typer.Option(False, "--delta", help="also report the poly-alanine-referenced ΔΦ per interface and ΔΦ total"),
    reference_aa: str = typer.Option("A", "--reference-aa", help="reference residue for --delta (default: alanine)"),
    geometry: bool = typer.Option(False, "--geometry", help="also report the interface-geometry descriptors and the decorrelated quality score Q"),
    skip_errors: bool = typer.Option(False, "--skip-errors", help="drop structures that fail instead of writing an error row"),
    threads: int = typer.Option(1, "-t", "--threads", help="worker processes for a multi-structure run (0 = all cores); each also gets mmseqs threads"),
) -> None:
    """Score structures: per-interface contact energies Φ (and ΔΦ, and interface geometry).

    This is **scoring only** — it reads structures and writes numbers. The preparation steps
    (canonicalisation, region mapping, Cα / contact / atom-distance matrices) are separate
    commands: ``tcren annotate``, ``tcren superimpose``, ``tcren contacts``.

    Columns ``F_tcr_pep``, ``F_tcr_mhc``, ``F_pep_mhc`` are the three interface terms
    Φ_TP, Φ_TM, Φ_PM; ``F_total`` is their sum Φ. With ``--delta`` each also gets its
    poly-alanine-referenced counterpart ``dF_*`` (ΔΦ_TP, ΔΦ_TM≡0, ΔΦ_PM) and ``dF_total`` = ΔΦ.
    The names match ``tcren recognize``, so the two tables join on ``pdb.id``.
    ΔΦ is the score to use across candidates that each carry their **own** generated pose,
    where raw Φ partly reads the pose geometry rather than the peptide sequence.

    ``--intra-weight w`` adds the term the three interface sums omit — ``F_pep_int``, the peptide's
    contact energy with **itself** (5 Å, sequence separation >= 3, MJ) — and folds ``w x F_pep_int``
    into ``F_total``. The energy is reported raw, so the term and the weight stay separable.

    ``--geometry`` appends the interface descriptors (buried surface ``burial``, peptide
    coverage ``n_pep_contacted``, ``chain_balance``, ``n_hbond``, docking ``pitch``/``crossing``)
    and ``Q`` — the directional, decorrelated interface-quality score, standardised against the
    native-crystal reference so it is defined for a single structure (:func:`tcren.q_score`).
    For the complete descriptor catalogue plus P(real), use ``tcren recognize``.

    Each interface's potential can be overridden with a bundled name (``tcren2``,
    ``karnaukhov2022``, ``mj``, ``keskin``) or a CSV path; an unset option keeps the default
    family for that interface.

    Examples::

        tcren scoring -s complex.pdb.gz -o scores.csv
        tcren scoring -s a.pdb.gz -s b.pdb.gz --delta          # repeat -s, or comma-separate
        tcren scoring -s 'models/*.pdb.gz' --delta --geometry  # quote the glob
        tcren scoring -s models/ --delta -t 8                  # a directory, 8 workers
        tcren scoring -s models.txt --delta                    # one path per line
        tcren scoring -s models.tar.gz --regions cdr           # CDR contacts only

    Scoring a cohort is embarrassingly parallel and dominated by the per-structure mmseqs
    annotation, so ``-t`` is worth setting for anything above a handful of structures
    (``-t 0`` uses every core).
    """
    from .pipeline import run as run_pipeline, score_row
    from .structure.io import resolve_sources

    if regions not in TCR_REGIONS:
        raise typer.BadParameter("--regions must be one of all|cdr|cdr+fr")
    if contact_weight not in ("residue", "atomic"):
        raise typer.BadParameter("--contact-weight must be residue or atomic")
    potentials = {
        "tcr_peptide": tcr_peptide_potential,
        "tcr_mhc": tcr_mhc_potential,
        "peptide_mhc": peptide_mhc_potential,
    }
    kw = dict(organism=organism, superimpose=not no_superimpose, db_dir=db, cutoff=cutoff,
              potentials=potentials, tcr_regions=regions, contact_weight=contact_weight,
              reference_aa=reference_aa if delta else None, intra_weight=intra_weight)
    def score_struct(s, typed=False):
        return score_row(run_pipeline(s, typed=typed, **kw))

    # Chain typing is one mmseqs search per structure and dominated a cohort run. It is done
    # here for the whole cohort in a handful of searches (one per organism) and the scorer is
    # told not to repeat it -- threads only ever hid that cost behind more concurrent mmseqs
    # processes, each of which defaults to every core.
    structs = [s for src in resolve_sources(structures) for s in _iter_typed(src, organism)]

    rows, failed, first_error = [], 0, None
    def one(s):
        try:
            return score_struct(s, typed=True)
        except Exception as exc:  # noqa: BLE001 - keep the batch resilient
            return {"pdb.id": s.pdb_id, "F_total": None,
                    "error": f"{type(exc).__name__}: {str(exc)[:80]}"}

    if threads == 1:
        results = map(one, structs)
    else:
        import os as _os
        from concurrent.futures import ThreadPoolExecutor
        n = threads if threads > 0 else (_os.cpu_count() or 1)
        # What is left after annotation is contact-map construction and the energy sums, which
        # spend their time in numpy and release the GIL.
        ex = ThreadPoolExecutor(max_workers=n)
        results = ex.map(one, structs)
    for r in results:
        if r.get("F_total") is None and "error" in r:
            failed += 1
            first_error = first_error or r["error"]
            if skip_errors:
                continue
        rows.append(r)
    if not rows:
        raise typer.BadParameter(f"no structures scored from {list(structures)}")
    table = pl.DataFrame(rows, strict=False)

    if geometry:
        # Reuse the recognition descriptors verbatim rather than recomputing geometry here:
        # `tcren recognize` stays the one definition of every descriptor.
        from .cohort import Q_FEATURES_GEOM, q_score
        from .recognition import recognition_table
        from .structure.io import import_structure

        items = [it for src in resolve_sources(structures)
                 for it in iter_structures(src, importer=import_structure, on_error="skip")]
        import os as _os2
        geo = pl.DataFrame(recognition_table(
            items, organism=organism, with_p_real=False,
            threads=threads if threads > 0 else (_os2.cpu_count() or 1)))
        geo = geo.with_columns(pl.Series("Q", q_score(geo, features=Q_FEATURES_GEOM)))
        keep = ["complex.id", *Q_FEATURES_GEOM, "pitch", "crossing", "Q"]
        geo = geo.select([c for c in keep if c in geo.columns]).rename({"complex.id": "pdb.id"})
        table = table.join(geo, on="pdb.id", how="left")

    table.write_csv(str(out))
    typer.echo(f"wrote {out} ({table.height} rows"
               + (f", {failed} failed" if failed else "") + ")")
    if failed:  # an error buried in a column is an error nobody reads
        # --skip-errors drops the failed rows, so the message cannot come from the table.
        typer.secho(f"first failure: {first_error}", fg="red", err=True)
    if failed == len(structs):
        # Every structure failing is an environment fault, not a data one -- a missing reference,
        # a broken install. Exiting 0 with a table that has no score columns makes the caller die
        # later on a missing column, several stages downstream of the real cause.
        raise typer.Exit(1)


@app.command(rich_help_panel=_P_SCORE, hidden=True)
def pipeline(ctx: typer.Context) -> None:
    """Deprecated alias for ``tcren scoring`` (this command never ran the full pipeline)."""
    raise typer.BadParameter(
        "`tcren pipeline` is now `tcren scoring` — it scores structures, it does not run the "
        "preparation pipeline (see `tcren annotate`, `tcren superimpose`, `tcren contacts`)."
    )


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
    threads: int = typer.Option(1, "-t", "--threads", help="worker processes for the scoring stage (0 = all cores); the mmseqs searches always use all cores"),
    autodetect_species: bool = typer.Option(True, "--autodetect-species/--no-autodetect-species", help="also search mouse to catch a mis-declared organism; --no- halves the annotation cost"),
) -> None:
    """Interface mechanics — the koff proxies: stiffness tensor + steered rupture + coupling residues.

    Treats the TCR↔pMHC contact map as a network of breakable springs and reports, per structure:
    ``n_spring``, ``S_tot``/``K_tens``/``K_shear``/``aniso`` (stiffness tensor), ``rupture_force``/
    ``rupture_work`` (steered unbinding), and ``couple_pep``/``couple_total`` (coupling residues).
    Validated on ATLAS: the tensile stiffness / rupture resistance track the dissociation off-rate
    (koff) far better than the equilibrium ΔG/Kd (Bell–Evans; the TCR is a mechanosensor).
    """
    from .mechanics import WEIGHTS

    if weight not in WEIGHTS:
        raise typer.BadParameter(f"--weight must be one of {'|'.join(WEIGHTS)}")
    if direction not in ("tensile", "shear", "auto"):
        raise typer.BadParameter("--direction must be one of tensile|shear|auto")
    # Stage 1: annotate the whole set in one mmseqs search per organism. One call per structure
    # rebuilds the arda index every time -- that is the difference between ~2/s and ~8/s on a
    # cohort, and it is why this command used to be the slow one.
    items = list(iter_structures(structures, importer=parse_structure))
    _annotate_set([s for _, s in items], organism=organism,
                  autodetect_species=autodetect_species)
    # Stage 2: the spring network, which is where the time goes and is pure Python/numpy.
    # `threads` > 1 runs it in that many worker processes; the stages never overlap.
    work = [(pid, s, cutoff, weight, direction, break_strain) for pid, s in items]
    if threads != 1 and len(work) > 1:
        import os as _os
        from concurrent.futures import ProcessPoolExecutor
        n = min(threads if threads > 0 else (_os.cpu_count() or 1), len(work))
        with ProcessPoolExecutor(max_workers=n) as ex:
            rows = list(ex.map(_mechanics_one, work, chunksize=max(1, len(work) // (n * 4))))
    else:
        rows = [_mechanics_one(w) for w in work]
    pl.DataFrame(rows).write_csv(str(out))
    typer.echo(f"wrote {out}")


def _mechanics_one(args) -> dict:
    """One annotated structure -> one mechanics row. Module-level so it pickles to a worker."""
    from .mechanics import interface_mechanics

    pid, s, cutoff, weight, direction, break_strain = args
    try:
        return {"pdb.id": pid, **interface_mechanics(
            s, cutoff=cutoff, weight=weight, direction=direction, break_strain=break_strain)}
    except Exception as exc:  # noqa: BLE001 - keep the batch resilient
        return {"pdb.id": pid, "K_tens": None, "error": f"{type(exc).__name__}: {str(exc)[:80]}"}


def _annotate_set(structs, *, organism: str, autodetect_species: bool) -> None:
    """Chain-type and MHC-annotate a whole set in place: one mmseqs search per organism, all cores.

    mmseqs parallelises internally, so it is handed the whole set once rather than N concurrent
    batches. The previous shape ran a thread pool over 64-structure batches, and since arda leaves
    mmseqs at its all-cores default, N batches asked for N x cpu_count and the run thrashed. The
    MHC search was worse off still: it was called with no thread count at all, so it ran
    single-threaded over every batch.

    The MHC pass must stay AFTER chain typing, or the MHC interfaces come out silently empty.
    """
    import os as _os

    from .annotation import classify_chains
    from .annotation.arda_adapter import _import_arda
    from .mhc import annotate_mhc_batch
    from .paper.helpers import annotate_batch

    structs = [s for s in structs if s is not None]
    if not structs:
        return
    cores = _os.cpu_count() or 1
    orgs = (organism, "mouse") if autodetect_species else (organism,)
    recs = annotate_batch(structs, _import_arda(), organisms=orgs, threads=cores)
    for i, s in enumerate(structs):
        try:
            classify_chains(s, organism=organism, autodetect_species=autodetect_species,
                            precomputed_records=recs[i])
        except Exception:  # noqa: BLE001 - unannotatable chains stay unset
            pass
    annotate_mhc_batch(structs, threads=cores)


@app.command(rich_help_panel=_P_ORIENT)
def refine(
    structures: str = typer.Option(..., "-s", "--structures", help="structure file, directory, .tar.gz, or glob"),
    out: Path = typer.Option("refined", "-o", "--out", help="output directory for refined structures"),
    substitute: str = typer.Option(None, "--substitute", help="thread this peptide onto the backbone first"),
    organism: str = typer.Option("human", "--organism"),
    n_steps: int = typer.Option(2000, "--steps", help="Monte-Carlo steps"),
    restraint_w: float = typer.Option(0.5, "--restraint", help="harmonic restraint to the input pose"),
    seed: int = typer.Option(0, "--seed"),
    repack: bool = typer.Option(False, "--repack", help="also place each peptide side chain in the chi rotamer DOPE prefers (native _relax packer)"),
    max_chi: int = typer.Option(2, "--max-chi", help="chi angles sampled per residue when repacking"),
    mmcif: bool = typer.Option(False, "--mmCIF", help="write mmCIF (.cif) instead of PDB"),
    compress: bool = typer.Option(False, "--compress", help="gzip the output (.gz)"),
) -> None:
    """Potential-guided rigid-body refinement of the peptide pose (knowledge-based, not physics).

    Optionally ``--substitute`` a new equal-length peptide first, then run a Monte-Carlo refinement
    scored by the DOPE atom-level statistical potential (restrained to the input pose; independent of
    the TCRen/MJ scoring potentials). Writes one structure per input and prints the final DOPE
    energy.

    ``--repack`` adds the side-chain half: the MC moves the peptide rigidly and leaves every χ where
    it found it, so a full-atom model whose side chains a predictor placed keeps them. The packer
    re-samples χ discretely, which is what a *local* minimiser structurally cannot do — measured on
    five crystals with χ1 deliberately rotated 120°, it recovers side-chain RMSD from 4.13 Å to
    2.36 Å in 6 ms, where OpenMM's restrained minimisation returns 4.13 Å (unchanged) in 3.1 s,
    because gradient descent cannot cross a torsional barrier. It rotates the side chains a model
    *has*; it cannot rebuild ones ``--substitute`` stripped.
    """
    from .refine import refine_peptide, substitute_peptide
    from .rotamers import repack as repack_sidechains
    from .structure.io import import_structure, structure_output_path, write_structure

    out.mkdir(parents=True, exist_ok=True)
    rows = []
    for pid, s in iter_structures(structures, importer=import_structure):
        try:
            classify_chains(s, organism=organism)
            if substitute:
                s = substitute_peptide(s, substitute)
            oriented, energy = refine_peptide(s, restraint_w=restraint_w, n_steps=n_steps, seed=seed)
            row = {"pdb.id": pid, "energy": energy}
            if repack:
                oriented, report = repack_sidechains(oriented, max_chi=max_chi)
                row["repack.energy"] = float(report["energy"].sum())
                row["repack.p_best.min"] = float(report["p_best"].min())
            write_structure(oriented, structure_output_path(out, pid, mmcif, compress))
            rows.append(row)
        except Exception as exc:  # noqa: BLE001 - keep the batch resilient
            rows.append({"pdb.id": pid, "energy": None,
                         "error": f"{type(exc).__name__}: {str(exc)[:80]}"})
    pl.DataFrame(rows).write_csv(str(out / "refine_energies.csv"))
    typer.echo(f"refined {sum(r.get('energy') is not None for r in rows)}/{len(rows)} -> {out}")


@app.command("substitute-tcr", rich_help_panel=_P_ORIENT)
def substitute_tcr_cmd(
    host: Path = typer.Option(..., "--host", help="host complex — keeps its peptide + MHC"),
    donor: Path = typer.Option(..., "--donor", help="donor complex — its TCR is grafted on"),
    out: Path = typer.Option("chimera.pdb", "-o", "--out",
                             help="output structure (format from the suffix: .pdb/.cif/.pdb.gz/.cif.gz)"),
    by: str = typer.Option("mhc", "--by",
                           help="superposition anchor: mhc (donor keeps native docking) | tcr (inherits host pose)"),
    organism: str = typer.Option("human", "--organism"),
) -> None:
    """Graft the donor TCR onto the host pMHC → a chimeric TCR:pMHC complex.

    Keeps the host peptide + MHC and the donor TCR. ``--by mhc`` superposes the donor MHC groove onto
    the host groove (the donor TCR keeps its native docking geometry); ``--by tcr`` superposes the
    donor TCR onto the host TCR (the donor TCR inherits the host's docking pose). Both inputs are
    chain-typed automatically (and, for ``--by mhc``, MHC-annotated).
    """
    if by not in ("mhc", "tcr"):
        raise typer.BadParameter("--by must be 'mhc' or 'tcr'")
    from .orient import substitute_tcr
    from .structure.io import import_structure, write_structure

    h = import_structure(host)
    classify_chains(h, organism=organism)
    d = import_structure(donor)
    classify_chains(d, organism=organism)
    if by == "mhc":
        from .mhc import annotate_mhc
        annotate_mhc(h)
        annotate_mhc(d)
    try:
        chimera = substitute_tcr(h, d, by=by)
    except ValueError as exc:
        raise typer.BadParameter(str(exc)) from exc
    write_structure(chimera, out)
    typer.echo(f"grafted {d.pdb_id} TCR onto {h.pdb_id} pMHC (by {by}) -> {out}")


def _score_feature_table(path: Path, out: Path) -> None:
    """`tcren recognize --features`: turn a `tcren features` table into cohort scores.

    The whole point of splitting the two commands is that the expensive pass runs once. Everything
    here is arithmetic over an existing table -- no structure is parsed and nothing is annotated.
    """
    from .cohort import P_NATIVE_CHANNELS, Q_FEATURES_GEOM, p_native, q_score

    sep = "," if path.suffix.lower() == ".csv" else "\t"
    t = pl.read_csv(path, separator=sep, infer_schema_length=None)
    if "complex.id" not in t.columns:
        raise typer.BadParameter(f"--features needs a 'complex.id' column; got {list(t.columns)[:6]}")
    scores = t.select("complex.id")

    try:
        scores = scores.with_columns(pl.Series("Q", q_score(t, features=Q_FEATURES_GEOM)))
    except KeyError as exc:
        typer.echo(f"  Q skipped: {exc.args[0].split(';')[0]}")
    # Score only the channels this table actually carries. `tcren features -i topology` is a
    # legitimate call and must still yield `T`; asking for all three would raise on the absent
    # geometry columns and drop the whole score.
    from .cohort import _channel_columns
    have = tuple(c for c in P_NATIVE_CHANNELS
                 if sum(n in t.columns for n in _channel_columns(c)) >= 2)
    if not have:
        typer.echo("  P_native skipped: no channel has two usable columns in this table")
    try:
        if not have:
            raise ValueError("no usable channel")
        v, models = p_native(t, channels=have, return_model=True)
        if len(have) == 1:
            models = {have[0]: models}
        scores = scores.with_columns(pl.Series("P_native", v))
        # each channel on its own as well: `T` is the shape score the paper reports, and a caller
        # who wants to know WHY a structure scored as it did needs the parts, not just the sum
        for ch in models:
            scores = scores.with_columns(
                pl.Series({"geometry": "G", "topology": "T", "energetics": "E"}.get(ch, ch),
                          p_native(t, channels=(ch,), rule="flat")))
        typer.echo(f"  P_native: {'sum of per-channel log-odds over ' if len(have) > 1 else ''}"
                   f"{', '.join(models)} ({len(have)} of {len(P_NATIVE_CHANNELS)} channels present)")
        for ch, model in models.items():
            typer.echo(f"    {ch}: EM over {len(model.feature_names)} features, "
                       f"{len(model.loglik_)} rounds, "
                       + ("converged" if model.converged_ else "hit the round cap"))
            typer.echo("      class weight per feature (sign and size EM chose, no labels used):")
            for j, nm in enumerate(model.feature_names):
                typer.echo(f"        {nm:18s} {model.nodes_[j]['beta'][-2]:+.3f}")
    except ValueError as exc:
        typer.echo(f"  P_native skipped: {exc}")

    scores.write_csv(str(out), separator="\t")
    typer.echo(f"wrote {out} ({scores.height} rows, {len(scores.columns) - 1} scores)")


@app.command(rich_help_panel=_P_SCORE)
def features(
    structures: str = typer.Option(..., "-s", "--structures", help="TCR-pMHC structure file, directory, glob, or .tar.gz"),
    out: Path = typer.Option("features.tsv", "-o", "--out", help="per-structure descriptor table (TSV)"),
    include: str = typer.Option("placement,interface,topology,energetics", "-i", "--include", help="comma-separated feature families: placement, interface, topology, energetics, kinetics"),
    all_families: bool = typer.Option(False, "--all", help="every family, kinetics included"),
    organism: str = typer.Option("human", "--organism"),
    radii: str = typer.Option("7,8", "--radii", help="Calpha thresholds for the footprint flag complex (topology family)"),
    threads: int = typer.Option(1, "-t", "--threads", help="worker processes for featurisation (0 = all cores); annotation is always one batched call"),
    autodetect_species: bool = typer.Option(True, "--autodetect-species/--no-autodetect-species", help="also search mouse to catch a mis-declared organism; --no- halves the annotation cost"),
) -> None:
    """Raw per-structure descriptors, one row per structure, in the four (+1) feature families.

    This command emits **features only** — no model, no probability, no cohort score. Its companion
    is ``tcren recognize``, which turns a feature table into scores and can read this file back with
    ``--features`` instead of re-reading the structures.

    The families are split by what each quantity is invariant under, which is also the axis along
    which they carry independent evidence:

    \b

    * ``placement`` -- where the receptor sits in the groove frame: docking angles, the TCRdock
      rigid-body parameters, ride height / shift / offset, and the CDR3 loop frames.
      Frame-DEPENDENT.
    * ``interface`` -- how much contact there is and of what chemical kind: buried area, contact
      counts and types, hydrogen bonds, clashes, chain and loop balance.
    * ``topology`` -- the SHAPE of the contact set, free of its size: coverage entropy and Hill
      numbers, the footprint's Betti numbers and persistence entropy, the canonical germline/CDR3
      preference.
    * ``energetics`` -- statistical-potential interface energies F and their poly-alanine
      references dF.
    * ``kinetics`` -- the interface as a spring network: stiffness, rupture, coupling residues.
      Off by default (it is the most expensive family); add it with --all or -i.

    Only what you ask for is computed: ``-i topology`` never builds the energies and ``-i placement``
    never runs the spring network. Whatever the selection, the whole set is annotated in **one**
    arda call per organism and **one** mmseqs MHC search, never one per structure.

    Every emitted column is catalogued in ``tcren.recognition.DESCRIPTORS``, so the families are a
    partition of the table rather than a label on it.

    Examples::

        tcren features -s models/ -o feats.tsv                       # the four default families
        tcren features -s models/ -i topology -o shape.tsv           # footprint shape alone
        tcren features -s models/ --all -t 0 -o feats.tsv            # everything, all cores
        tcren recognize --features feats.tsv -o scores.tsv           # score without re-reading structures
    """
    import os as _os

    from .recognition import FAMILIES, recognition_table
    from .structure.io import import_structure

    fams = list(FAMILIES) if all_families else [f.strip() for f in include.split(",") if f.strip()]
    unknown = [f for f in fams if f not in FAMILIES]
    if unknown:
        raise typer.BadParameter(f"unknown feature families {unknown}; expected {list(FAMILIES)}")
    try:
        rr = tuple(float(v) for v in radii.replace(" ", "").split(",") if v)
    except ValueError as exc:
        raise typer.BadParameter(f"--radii must be comma-separated numbers, got {radii!r}") from exc

    items = list(iter_structures(structures, importer=import_structure))
    rows = recognition_table(items, organism=organism, include=fams, radii=rr,
                             autodetect_species=autodetect_species,
                             threads=threads if threads > 0 else (_os.cpu_count() or 1))
    table = pl.DataFrame(rows)
    table.write_csv(str(out), separator="\t")
    n_err = int(table["error"].is_not_null().sum()) if "error" in table.columns else 0
    typer.echo(f"features [{','.join(fams)}]: {table.height} structures, {len(table.columns) - 1} "
               f"descriptors -> {out}" + (f"  ({n_err} failed)" if n_err else ""))


@app.command(rich_help_panel=_P_SCORE, hidden=True)
def footprint(
    structures: str = typer.Option(..., "-s", "--structures", help="TCR-pMHC structure file, directory, glob, or .tar.gz"),
    out: Path = typer.Option("footprint.tsv", "-o", "--out", help="per-structure coverage + topology table (TSV)"),
    organism: str = typer.Option("human", "--organism"),
    cutoff: float = typer.Option(5.0, "--cutoff", help="heavy-atom contact threshold (A)"),
    radii: str = typer.Option("7,8", "--radii", help="Calpha thresholds for the footprint flag complex; b0 is most informative at 7 A and b1 at 8 A"),
    group: str = typer.Option(None, "--group", help="column to fit T within (e.g. epitope); needs --meta"),
    meta: Path = typer.Option(None, "--meta", help="TSV/CSV with a 'pdb.id' column plus --group, joined before scoring"),
    score: bool = typer.Option(False, "--score", help="append T, the cohort-fitted shape-channel posterior (p_native over the topology family)"),
) -> None:
    """Footprint shape: how a receptor's contacts are DISTRIBUTED, not what they score.

    Superseded by ``tcren features -i topology``, which emits the same columns from the shared
    feature pass; kept working, and hidden from the command list.

    One TSV row per structure with the coverage measures -- normalised Shannon entropy
    ``H_cell`` and the Hill numbers ``D1``/``D2`` over the 6 CDR loops x {peptide, MHC}, plus
    ``D2_pep24`` on the finer partition that splits the peptide into thirds -- the canonical
    docking preference (``L_canon``, ``p_germ_mhc``, ``p_cdr3_pep``), the alpha/beta contact
    imbalance, and the footprint's topology (``fp_b0_*`` patches, ``fp_b1_*`` holes, the Euler
    characteristic, and the H0 persistence entropy).

    None of these is an energy and none needs a potential, a reference structure or a fitted
    parameter. They are invariant under rigid motion, so the inputs do **not** have to be
    oriented -- only chain-typed with CDR region markup, which this command does for you in one
    batched annotation pass over the whole set.

    ``--score`` adds ``T``, the shape channel's posterior: the same latent-class Bayes network
    ``tcren recognize`` fits, restricted to these descriptors. It is fitted on the cohort being
    scored, so it needs a set and is undefined for a single input. Use ``--group`` with ``--meta``
    to fit within epitope rather than across the whole table.

    Complementary scorer on the same inputs: ``tcren recognize`` (energies + P(real)).
    """
    from .footprint import footprint_batch

    try:
        rr = tuple(float(v) for v in radii.replace(" ", "").split(",") if v)
    except ValueError as exc:
        raise typer.BadParameter(f"--radii must be comma-separated numbers, got {radii!r}") from exc
    if not rr:
        raise typer.BadParameter("--radii needs at least one value")

    table = footprint_batch(structures, cutoff=cutoff, radii=rr, organism=organism)
    if table.height == 0:
        raise typer.BadParameter(f"no structures scored from {structures!r}")
    if meta is not None:
        sep = "," if meta.suffix.lower() == ".csv" else "\t"
        m = pl.read_csv(meta, separator=sep, infer_schema_length=None)
        if "pdb.id" not in m.columns:
            raise typer.BadParameter(f"--meta needs a 'pdb.id' column; got {list(m.columns)}")
        table = table.join(m, on="pdb.id", how="left")
    if score or group:
        # 2.12.0 removed the fp_score z-sum in favour of the fitted channel posterior; this is the
        # replacement its changelog names. Fitted per group where one is given, because the EM is a
        # cohort operation and pooling epitopes fits a mixture across them.
        from .cohort import p_native  # noqa: PLC0415

        if group and group not in table.columns:
            raise typer.BadParameter(f"--group {group!r} is not a column; pass it via --meta")
        if group:
            table = pl.concat([
                part.with_columns(pl.Series("T", p_native(part, channels=("topology",))))
                for (_key,), part in table.group_by([group], maintain_order=True)
            ])
        else:
            table = table.with_columns(pl.Series("T", p_native(table, channels=("topology",))))
    table.write_csv(out, separator="\t")
    typer.echo(f"footprint: {table.height} structures, {len(table.columns)} columns -> {out}")


def main() -> None:
    """Console-script entry point.

    Every command takes ``-s`` as a free-form spec (file, directory, glob, manifest, archive), so
    Typer cannot check it exists and a typo'd path surfaced as an 80-line Rich traceback out of
    Biopython. One line is what the user needs; ``typer.BadParameter`` covers the rest.
    """
    try:
        app()
    except OSError as exc:
        typer.secho(f"Error: {exc}", fg=typer.colors.RED, err=True)
        raise SystemExit(1) from None


if __name__ == "__main__":
    main()
