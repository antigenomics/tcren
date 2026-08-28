"""The contact map, read as a map — predicted contact frequency and peptide residue importance.

A marimo notebook — reactive, so changing the structure or the grouping redraws only what depends
on it, and plain Python, so it diffs and imports like any other module.

    marimo edit notebooks/potts_contact_map.py     # explore and change the code
    marimo run  notebooks/potts_contact_map.py     # read-only app, controls still live

This is the released path, not a re-implementation: the same call `tcren potts map` makes. It needs
only a structure file, and defaults to the crystal shipped with the tests.

The point is that a fitted contact model predicts *where* a receptor touches a peptide, not only
how favourably. Two readings come out of the same marginals: the loop x position grid, which is the
observable a molecular-dynamics trajectory reports as a contact frequency, and its collapse onto
the peptide alone, which is a residue-importance profile.
"""

import marimo

__generated_with = "0.23.15"
app = marimo.App(width="medium", app_title="tcren · potts contact map")


@app.cell
def _():
    import marimo as mo
    return (mo,)


@app.cell
def _(mo):
    mo.md(
        """
        # The contact map, read as a map

        `contact_probabilities` gives $p_j$, the coupled model's marginal for one
        receptor-residue : peptide-residue pair. An experiment measures something coarser — whether
        a whole CDR loop touches a given peptide position at all. The residues of a loop are
        distinct pairs with different marginals, so the count of simultaneous contacts is
        Poisson-binomially distributed and has no closed form. The event "at least one" does:

        $$P(N \\ge 1) = 1 - \\prod_j \\left(1 - p_j\\right)$$

        That is what `contact_map` computes, accumulated in $\\log(1-p)$ so a twelve-residue loop
        does not underflow. These are **frequencies, not energies** — dimensionless, in $[0, 1]$,
        carrying no $k_\\mathrm{B}T$.
        """
    )
    return


@app.cell
def _(mo):
    from pathlib import Path

    import numpy as np
    import polars as pl

    from tcren.annotation import classify_chains
    from tcren.potts import PottsModel, available_pairs, contact_map
    from tcren.structure import parse_structure

    default = Path(__file__).resolve().parents[1] / "tests" / "assets" / "pdb" / "1ao7.pdb"
    path = mo.ui.text(value=str(default), label="structure (PDB)", full_width=True)
    path
    return (Path, PottsModel, available_pairs, classify_chains, contact_map, np,
            parse_structure, path, pl)


@app.cell
def _(PottsModel, available_pairs, classify_chains, mo, parse_structure, path):
    _p = path.value
    structure = parse_structure(_p, pdb_id=str(_p).rsplit("/", 1)[-1].split(".")[0])
    classify_chains(structure, organism="human", autodetect_species=True)
    pairs = available_pairs(structure)
    model = PottsModel.bundled()
    mo.md(
        f"**{pairs['pdb.id'][0]}** — {pairs.height} available residue pairs, "
        f"{int(pairs['sigma'].sum())} of them in contact at {model.cutoff:.0f} Å, over "
        f"{pairs['pos.par'].n_unique()} peptide positions and "
        f"{pairs['region.rec'].n_unique()} receptor regions."
    )
    return model, pairs


@app.cell
def _(contact_map, model, pairs):
    # One sampler run, reused by both panels: `by="pair"` is the ungrouped table the two grids
    # aggregate, so drawing them from one call keeps them exactly consistent.
    loops = contact_map(pairs, model, by="loop")
    positions = contact_map(pairs, model, by="position")
    return loops, positions


@app.cell
def _(mo):
    mo.md(
        """
        ## The map, and what the structure did

        Left: the model's predicted frequency for each (CDR loop, peptide position) cell. Right:
        the same cells, coloured by whether this structure actually made a contact there. A cell
        the model calls likely and the structure leaves empty is the interesting kind of
        disagreement — the pose could have made that contact and did not.
        """
    )
    return


@app.cell
def _(loops, np, pl):
    import matplotlib.pyplot as plt

    _reg = sorted(loops["region.rec"].unique().to_list())
    _pos = sorted(loops["pos.par"].unique().to_list())
    _pred = np.full((len(_reg), len(_pos)), np.nan)
    _obs = np.full((len(_reg), len(_pos)), np.nan)
    for _r in loops.iter_rows(named=True):
        _i, _j = _reg.index(_r["region.rec"]), _pos.index(_r["pos.par"])
        _pred[_i, _j] = _r["p_any"]
        _obs[_i, _j] = _r["observed"]

    _fig, _ax = plt.subplots(1, 2, figsize=(10, 3.2), constrained_layout=True)
    for _a, _m, _t in ((_ax[0], _pred, "predicted  $P(N \\geq 1)$"),
                       (_ax[1], _obs, "observed  (0/1)")):
        _im = _a.imshow(_m, aspect="auto", vmin=0, vmax=1, cmap="viridis")
        _a.set_xticks(range(len(_pos)), [str(p + 1) for p in _pos])
        _a.set_yticks(range(len(_reg)), _reg, fontsize=8)
        _a.set_xlabel("peptide position")
        _a.set_title(_t, fontsize=10)
        _fig.colorbar(_im, ax=_a, shrink=0.85)
    _fig
    return (plt,)


@app.cell
def _(mo):
    mo.md(
        """
        ## Peptide residue importance

        Collapsing the loops leaves one number per peptide position: how engaged the model expects
        that residue to be. `p_any` is the probability of *any* receptor contact; `p_expected` is
        the expected *number* of them, so it keeps rising after `p_any` has saturated and is the
        better read-out at the centre of the peptide. Neither has seen a residue identity — this is
        the footprint before any energy is scored.
        """
    )
    return


@app.cell
def _(np, plt, positions):
    _d = positions.sort("pos.par")
    _x = np.arange(_d.height)
    _fig, _ax = plt.subplots(figsize=(10, 3.0), constrained_layout=True)
    _ax.bar(_x - 0.2, _d["p_any"].to_numpy(), width=0.4, label="$P(N \\geq 1)$")
    _ax.bar(_x + 0.2, _d["p_expected"].to_numpy() / max(_d["p_expected"].max(), 1e-9),
            width=0.4, label="expected count (scaled)")
    _ax.plot(_x, _d["n_observed"].to_numpy() / max(_d["n_observed"].max(), 1),
             "k.--", lw=0.8, ms=6, label="observed count (scaled)")
    _ax.set_xticks(_x, [f"{a}{p + 1}" for a, p in zip(_d["aa.par"], _d["pos.par"])])
    _ax.set_xlabel("peptide residue")
    _ax.set_ylabel("engagement")
    _ax.legend(fontsize=8, frameon=False)
    _fig
    return


@app.cell
def _(mo, positions):
    _top = positions.sort("p_expected", descending=True).head(3)
    mo.md(
        "The three positions the model expects to be most engaged are "
        + ", ".join(f"**{r['aa.par']}{r['pos.par'] + 1}** ({r['p_expected']:.2f} expected "
                    f"contacts, {r['n_observed']} observed)" for r in _top.iter_rows(named=True))
        + ". Reproduce this table from the command line with "
        "`tcren potts map -s <structure> --by position`."
    )
    return


if __name__ == "__main__":
    app.run()
