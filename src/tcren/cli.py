"""Command-line interface for tcren.

Subcommands:

* ``tcren info`` — environment / dependency check.
* ``tcren annotate`` — chain typing + region markup for input structures.
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
        arda_status = "NOT available — install with: pip install -e ../arda"
    typer.echo(f"arda: {arda_status}")
    typer.echo(f"bundled TCRen potential: {tcren().matrix.height} pairs")


@app.command()
def annotate(
    structures: Path = typer.Option(..., "-s", "--structures", help="structure file, directory, or .tar.gz (.pdb/.cif/.pdb.gz/.cif.gz)"),
    out: Path = typer.Option("markup.csv", "-o", "--out", help="output residue-markup CSV"),
    organism: str = typer.Option("human", "--organism"),
) -> None:
    """Annotate chains and emit a per-residue region-markup table."""
    from .contacts.table import residue_annotation

    frames = []
    for pid, s in iter_structures(structures, importer=parse_structure):
        classify_chains(s, organism=organism)
        frames.append(residue_annotation(s).with_columns(pl.lit(pid).alias("pdb.id")))
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
    structures: Path = typer.Option(..., "-s", "--structures", help="PDB/CIF file or directory"),
    out: Path = typer.Option("oriented", "-o", "--out", help="output dir for oriented PDBs"),
    metadata: Path = typer.Option("orient_metadata.csv", "--metadata"),
    organism: str = typer.Option("human", "--organism"),
    reference_id: str = typer.Option(None, "--reference", help="force a reference complex id"),
    force_pca: bool = typer.Option(False, "--force-pca", help="skip native superposition"),
    workers: int = typer.Option(1, "--workers", "-j", help="parallel worker processes"),
    push_to_hub: str = typer.Option(None, "--push-to-hub", help="HF dataset repo id to upload to"),
    hub_folder: str = typer.Option("Native2026", "--hub-folder"),
) -> None:
    """Canonicalize TCR-pMHC structures into the common MHC frame (A-E chains)."""
    from .orient import run_folder

    run_folder(structures, out, metadata=metadata, organism=organism,
               reference_id=reference_id, force_pca=force_pca, workers=workers)
    if push_to_hub:
        from .orient.hub import push_oriented

        push_oriented(out, push_to_hub, folder=hub_folder)
        typer.echo(f"pushed {out} -> {push_to_hub}/{hub_folder}")


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


@app.command()
def mhc(
    structures: Path = typer.Option(..., "-s", "--structures", help="structure file, directory, or .tar.gz (.pdb/.cif/.pdb.gz/.cif.gz)"),
    out: Path = typer.Option("mhc_calls.csv", "-o", "--out"),
    organism: str = typer.Option("human", "--organism"),
) -> None:
    """Map MHC chains to allele / class / role for input structures."""
    from .mhc import map_mhc

    rows = []
    for pid, s in iter_structures(structures, importer=parse_structure):
        classify_chains(s, organism=organism)
        for call in map_mhc(s):
            rows.append(
                {
                    "pdb.id": pid,
                    "chain.id": call.chain_id,
                    "chain.role": call.chain_role,
                    "mhc.class": call.mhc_class,
                    "allele": call.allele,
                    "identity": call.identity,
                }
            )
    pl.DataFrame(rows).write_csv(str(out))
    typer.echo(f"wrote {out}")


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
