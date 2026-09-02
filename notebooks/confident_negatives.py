"""Confident negatives — reading a generator's confidence together with the coordinates.

A marimo notebook — reactive, so moving the confidence slider redraws only what depends on it, and
plain Python, so it diffs and imports like any other module.

    marimo edit notebooks/confident_negatives.py     # explore and change the code
    marimo run  notebooks/confident_negatives.py     # read-only app, controls still live

This is the released path, not a re-implementation: `tcren.reliability.af_band` and
`tcren.reliability.s_score`, the same calls `tcren assess` makes. It runs on the shipped
native-crystal reference, so it needs no data of your own.

The point is the last section. A co-folding model returns a confident complex for any pair, binding
or not — so at a FIXED confidence the generator has already said everything it can say, and any
remaining separation has to come from the coordinates. Move the slider to the top of the confidence
range and watch the structural reading stay spread out while the confidence-only one collapses to a
single number.
"""

import marimo

__generated_with = "0.23.15"
app = marimo.App(width="medium", app_title="tcren · confident negatives")


@app.cell
def _():
    import marimo as mo
    return (mo,)


@app.cell
def _(mo):
    mo.md(
        """
        # Confident negatives

        `af_band` says how often a confidence band is wrong: it looks a reported confidence up in a
        frozen table of benchmark deciles and returns the observed non-binder fraction with a Wilson
        interval. It does not say *which* of the models in that band to distrust.

        `S` does, and it reads no label anywhere:

        $$S = \\frac{Q}{\\sigma_Q} + \\frac{T}{\\sigma_T}
          + \\frac{\\Pi - \\mu_\\Pi}{\\sigma_\\Pi}$$

        Three fit-free directional blocks — interface geometry $Q$, footprint shape $T$ and the
        contact energy read against the partition function $\\Pi$ — each divided by its own spread
        over the 374 Native2026 crystals, so all three carry equal weight in native-sd units.
        Nothing is fitted at score time, so $S$ is defined for a single structure.
        """
    )
    return


@app.cell
def _():
    import numpy as np

    from tcren.reliability import (af_band, available_bands, inversion_flag,
                                   reliability_reference, s_score, screening_yield)

    ref = reliability_reference()
    energy = np.asarray(ref["neg_energy"], float)
    S = s_score(ref, energy=energy)
    inversion = inversion_flag(ref, energy=energy)
    return (S, af_band, available_bands, energy, inversion, np, ref, s_score,
            screening_yield)


@app.cell
def _(available_bands, mo):
    which = mo.ui.dropdown(options=available_bands(), value="tcrvdb|ipTM",
                           label="frozen band table")
    conf = mo.ui.slider(0.50, 0.99, 0.01, value=0.88, label="the generator's confidence",
                        show_value=True)
    mo.hstack([which, conf])
    return conf, which


@app.cell
def _(af_band, conf, mo, which):
    b = af_band([conf.value], reference=which.value)[0]
    mo.md(
        f"""
        **{which.value}**, band {b['band']} — reported confidence
        {b['lo']:.3f} to {b['hi']:.3f}, {b['n']} benchmark models in it.

        | quantity | value |
        |---|---|
        | non-binders in this band | {b['p_nonbinder']:.1%} [{b['ci_lo']:.1%}, {b['ci_hi']:.1%}] |
        | what `S` still separates *inside* it | {b['s_roc_in_band']:.3f} ROC-AUC |

        The band is a decile of the benchmark's own confidence distribution, never a threshold
        scanned for an effect. A confidence outside the table's range clamps to the end band.
        """
    )
    return (b,)


@app.cell
def _(S, mo, np):
    ok = np.isfinite(S)
    mo.md(
        f"""
        ## What the structure adds at one confidence

        Over the {int(ok.sum())} native crystals the shipped reference carries, `S` spans
        **{np.nanmin(S):+.2f} to {np.nanmax(S):+.2f} native-sd units**, with an interquartile range
        of **{np.subtract(*np.nanpercentile(S, [75, 25])):.2f}**.

        The confidence-only reading is the same number for every one of them — it cannot be
        otherwise, since it never saw a structure. That flat line is the whole point.
        """
    )
    return (ok,)


@app.cell
def _(S, b, inversion, np, ok):
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(6.5, 3.6))
    sc = ax.scatter(S[ok], inversion[ok], c=inversion[ok], cmap="RdYlBu_r", s=22,
                    alpha=0.8, linewidths=0)
    ax.axvline(float(np.nanmedian(S)), ls="--", lw=1.2, color="0.35")
    ax.set_xlabel("$S$ (native-sd units, higher = more like a real complex)")
    ax.set_ylabel("inversion flag (energy minus the two shape blocks)")
    ax.set_title(f"the confidence alone says {b['p_nonbinder']:.0%} non-binder "
                 f"for every one of these")
    fig.colorbar(sc, ax=ax, label="inversion flag")
    fig.tight_layout()
    fig
    return ax, fig, plt, sc


@app.cell
def _(mo):
    mo.md(
        """
        A large positive **inversion flag** is the energy vouching for a structure the footprint does
        not, which is the pattern to distrust: a generator fakes favourable contacts far more easily
        than a well-formed footprint. It ranks and triages; it is not calibrated, and nothing in this
        package is — every out-of-fold-fitted read-out was removed at 2.28.0.
        """
    )
    return


@app.cell
def _(S, mo, screening_yield):
    y = screening_yield(S, budget=0.10, prevalence=0.48)
    mo.md(
        f"""
        ## What testing the top decile implies

        `screening_yield` returns the cut and nothing else: **{y['n_tested']} structures** above
        `S` = {y['threshold']:.2f}, the {y['rank_cut']:.1%} rank, and **{y['expected_hits']:.1f}
        expected hits** at a stated prevalence of 0.48. Enrichment over random is deliberately not
        returned — it needs labels the function does not have, and a NaN there would read like a
        measurement.
        """
    )
    return (y,)


@app.cell
def _(mo):
    mo.md(
        """
        ## Where this is validated, and where it is not

        On the balanced 22-cohort VDJdb panel the **top ipTM decile is 26.2 %** [18.7, 35.5]
        **non-binders**, and it is also the band where `S` reads highest. The panel is reported
        cohort by cohort and stratified by template coverage, because where no receptor has been
        co-crystallized with the peptide, nothing works — the generator's own confidence included,
        and it does not fall to warn you.

        For the two-class read-outs on the same input — `binder_score`, `confidence_residual` and
        the five channel scores — see `tcren assess` and the `score_vdjdb_panel` notebook.
        """
    )
    return


if __name__ == "__main__":
    app.run()
