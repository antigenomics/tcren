"""Command-line interface for tcren.

Subcommands:

* ``tcren info`` — environment / dependency check.
* ``tcren annotate`` — chain typing + region markup (TCR/MHC/peptide; ``--regions`` to filter,
  ``--pseudo`` for MHC pseudosequence residues) for input structures.
* ``tcren contacts`` — annotated contact table for input structures.
* ``tcren derive-potential`` — derive a TCRen potential from a contact-map table.
* ``tcren score`` — end-to-end candidate scoring (drop-in for ``run_TCRen.R``).
"""

from __future__ import annotations

from pathlib import Path

import polars as pl
import typer

from . import __version__
from .annotation import classify_chains
from .contactmap import ContactMap
from .potential import Potential, derive_tcren, derive_tcren_loo, tcren
from .scoring import score_peptides
from .structure import iter_structures, parse_structure

app = typer.Typer(add_completion=False, help="Structure-based TCR–epitope recognition scoring.")
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

def _load_potential(spec: str | None) -> Potential:
    if spec is None:
        return tcren()
    p = Path(spec)
    if p.exists():
        return Potential.from_csv(p)
    raise typer.BadParameter(f"potential file not found: {spec}")


def _read_candidates(path: Path) -> list[str]:
    lines = [line.strip() for line in path.read_text().splitlines() if line.strip()]
    return [line for line in lines if line.lower() != "peptide"]


@app.command()
def info() -> None:
    """Show version and dependency availability."""
    typer.echo(f"tcren {__version__}")
    try:
        import arda  # noqa: PLC0415

        arda_status = f"available ({Path(arda.__file__).parent})"
    except ImportError:
        arda_status = "NOT available — run: bash setup.sh (installs arda@2.0.1)"
    typer.echo(f"arda: {arda_status}")
    typer.echo(f"bundled TCRen potential: {tcren().matrix.height} pairs")


_REGION_CHAINS = {
    "tcr": {"TRA", "TRB", "TRD", "TRG"},
    "mhc": {"MHCa", "MHCb", "B2M"},
    "peptide": {"PEPTIDE"},
}


@app.command()
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


@app.command()
def contacts(
    structures: Path = typer.Option(..., "-s", "--structures", help="structure file, directory, or .tar.gz (.pdb/.cif/.pdb.gz/.cif.gz)"),
    out: Path = typer.Option("contacts.csv", "-o", "--out"),
    cutoff: float = typer.Option(5.0, "--cutoff"),
    interface: str = typer.Option("tcr_peptide", "--interface", help="tcr_peptide|tcr_mhc|peptide_mhc|all"),
    organism: str = typer.Option("human", "--organism"),
) -> None:
    """Compute and emit an annotated contact table."""
    frames = []
    for _pid, s in iter_structures(structures, importer=parse_structure):
        classify_chains(s, organism=organism)
        cm = ContactMap.from_structure(s, cutoff=cutoff)
        frames.append(cm.contacts if interface == "all" else cm.interface(interface))
    pl.concat(frames).write_csv(str(out))
    typer.echo(f"wrote {out}")


@app.command()
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


@app.command()
def superimpose(
    structures: Path = typer.Option(..., "-s", "--structures", help="structure file, directory, or .tar.gz to orient"),
    out: Path = typer.Option("superimposed", "-o", "--out", help="output dir for oriented structures"),
    db: Path = typer.Option(None, "--db", help="canonical database dir (default: data/Canonical2026, fetched at install)"),
    organism: str = typer.Option("human", "--organism"),
    mmcif: bool = typer.Option(False, "--mmCIF", help="write mmCIF (.cif) instead of PDB"),
    compress: bool = typer.Option(False, "--compress", help="gzip the output (.gz)"),
) -> None:
    """Superimpose structure(s) onto a canonical database by MHC.

    Detects each input's MHC chains, class, and species, then superposes its conserved groove Cα
    onto *every* database structure of the same class and species and averages the transforms into
    one consensus placement. The database defaults to ``data/Canonical2026`` (populated at install).
    """
    from .orient import run_superimpose

    run_superimpose(structures, out, db_dir=db, organism=organism, mmcif=mmcif, compress=compress)


@app.command("derive-potential")
def derive_potential(
    contact_maps: Path = typer.Option(..., "-i", "--contact-maps", help="contact-map CSV"),
    out: Path = typer.Option("TCRen_potential.csv", "-o", "--out"),
    summary: Path | None = typer.Option(None, "--summary", help="summary CSV with a nonred flag"),
    nonred: bool = typer.Option(False, "--nonred", help="restrict to non-redundant structures"),
    variant: str = typer.Option("classic", "--variant", help="classic|am"),
    pseudocount: int = typer.Option(1, "--pseudocount"),
    loo: bool = typer.Option(False, "--loo", help="emit leave-one-out potentials instead"),
) -> None:
    """Derive a TCRen potential from observed contacts."""
    contacts = pl.read_csv(contact_maps)
    include = None
    if nonred:
        if summary is None:
            raise typer.BadParameter("--nonred requires --summary")
        include = pl.read_csv(summary).filter(pl.col("nonred"))["pdb.id"].to_list()
    if loo:
        ids = include or contacts["pdb.id"].unique().to_list()
        derive_tcren_loo(contacts, ids, variant=variant, pseudocount=pseudocount).write_csv(str(out))
    else:
        pot = derive_tcren(contacts, include=include, variant=variant, pseudocount=pseudocount)
        pot.to_csv(out)
    typer.echo(f"wrote {out}")


@app.command("fetch-data")
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


@app.command("build-mhc-ref")
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


@app.command("fetch-recent")
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


@app.command()
def score(
    structures: Path = typer.Option(..., "-s", "--structures", help="structure file, directory, or .tar.gz (.pdb/.cif/.pdb.gz/.cif.gz)"),
    candidates: Path = typer.Option(..., "-c", "--candidates", help="candidate epitopes file"),
    potential: str | None = typer.Option(None, "-p", "--potential", help="potential CSV (default: bundled TCRen)"),
    out: Path = typer.Option("candidate_epitopes_TCRen.csv", "-o", "--out"),
    interface: str = typer.Option("tcr_peptide", "--interface"),
    organism: str = typer.Option("human", "--organism"),
    cutoff: float = typer.Option(5.0, "--cutoff"),
) -> None:
    """Score candidate epitopes against input structures (end-to-end pipeline)."""
    pot = _load_potential(potential)
    cands = _read_candidates(candidates)
    frames = []
    for _pid, s in iter_structures(structures, importer=parse_structure):
        classify_chains(s, organism=organism)
        cm = ContactMap.from_structure(s, cutoff=cutoff)
        frames.append(score_peptides(cm, cands, pot, interface=interface))
    result = pl.concat(frames) if frames else pl.DataFrame()
    result.write_csv(str(out))
    typer.echo(f"The ranked list of candidate epitopes can be found in {out}")


if __name__ == "__main__":
    app()
