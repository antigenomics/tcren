"""Confident negatives — reading a generator's confidence together with the coordinates.

A marimo notebook — reactive, so moving the confidence slider redraws only what depends on it, and
plain Python, so it diffs and imports like any other module.

    marimo edit notebooks/confident_negatives.py     # explore and change the code
    marimo run  notebooks/confident_negatives.py     # read-only app, controls still live

This is the released path, not a re-implementation: the same call `tcren diagnose` makes. Needs a
`tcren features` table carrying the geometry, topology and Potts columns; the notebook falls back
to the shipped native-crystal reference so it runs with no data of your own.

The point is the last section. A co-folding model returns a confident complex for any pair, binding
or not — so at a FIXED confidence the generator has already said everything it can say, and any
remaining separation has to come from the coordinates. Move the slider to the top of the confidence
range and watch the corrected probabilities stay spread out while the confidence-only reading
collapses to a single number.
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

        `af_band` says how often a confidence band is wrong. It does not say what to believe
        instead. `correct_confidence` does:

        $$\\mathrm{logit}\\,P(\\mathrm{binder}) = b_0 + b_c\\,z(c) + b_S\\,S_{\\mathrm{free}}
          + b_N\\,N$$

        Four frozen coefficients, fitted out of fold on the benchmarks. **This is the one shipped
        read-out that is not fit-free** — the structural terms take no label, but their weighting
        was learned from one.
        """
    )
    return


@app.cell
def _():
    import numpy as np

    from tcren.reliability import (available_corrections, correct_confidence, moments,
                                   reliability_reference)

    ref = reliability_reference()
    energy = np.asarray(ref["neg_energy"], float)
    contacts = np.asarray(ref["n_contacts"], float)
    return (available_corrections, contacts, correct_confidence, energy, moments, np,
            ref, reliability_reference)


@app.cell
def _(available_corrections, mo):
    which = mo.ui.dropdown(options=available_corrections(), value="tcrvdb|ipTM",
                           label="frozen correction")
    conf = mo.ui.slider(0.50, 0.99, 0.01, value=0.88, label="the generator's confidence",
                        show_value=True)
    mo.hstack([which, conf])
    return conf, which


@app.cell
def _(mo, moments, which):
    c = moments()["corrections"][which.value]
    mo.md(
        f"""
        **{which.value}** — fitted over {c['n_folds']} folds, {c['fold_scheme']}.

        | term | coefficient |
        |---|---|
        | intercept | {c['b0']:+.1f} |
        | the generator's confidence | {c['b_conf']:+.1f} |
        | $S_{{\\mathrm{{free}}}}$ | {c['b_s_free']:+.1f} |
        | observed contact count | {c['b_n_contacts']:+.1f} |

        Out of fold this reads {c['roc_corrected_out_of_fold']:.4f} macro ROC-AUC against the raw
        confidence's {c['roc_raw']:.4f}.
        """
    )
    return


@app.cell
def _(conf, contacts, correct_confidence, energy, np, ref, which):
    out = correct_confidence(ref, np.full(len(energy), conf.value), reference=which.value,
                             energy=energy, contacts=contacts)
    ok = np.isfinite(out["p_corrected"]) & np.isfinite(out["s_free"])
    return ok, out


@app.cell
def _(mo, np, ok, out):
    d = out["delta_logit"][ok]
    mo.md(
        f"""
        ## What the structure adds

        Over {int(ok.sum())} native crystals at this one confidence, the coordinates move the
        log-odds by **{np.mean(d):+.2f} nats on average**, range
        **[{np.min(d):+.2f}, {np.max(d):+.2f}]**. The structure argues *against*
        **{int((d < 0).sum())}** of them.

        The confidence-only reading is the same number for every row — it cannot be otherwise,
        since it never saw a structure. That flat line is the whole problem the correction exists
        to fix.
        """
    )
    return


@app.cell
def _(np, ok, out):
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(6.5, 3.6))
    s, p, d = out["s_free"][ok], out["p_corrected"][ok], out["delta_logit"][ok]
    sc = ax.scatter(s, p, c=d, cmap="RdYlBu", s=22, alpha=0.75, linewidths=0)
    ax.axhline(float(np.nanmean(out["p_confidence"][ok])), ls="--", lw=1.2, color="0.35")
    ax.set_xlabel("$S_{free}$ (native-sd units)")
    ax.set_ylabel("probability")
    ax.set_ylim(0, 1)
    ax.set_title("dashed: the confidence alone, identical for every structure")
    fig.colorbar(sc, ax=ax, label="$\\Delta$ log-odds from the structure")
    fig.tight_layout()
    fig
    return ax, d, fig, p, plt, s, sc


@app.cell
def _(mo):
    mo.md(
        """
        ## Where this is validated, and where it is not

        Leave-one-epitope-out on the balanced VDJdb panel, the correction adds **+0.051** macro
        ROC-AUC to ipTM and **+0.068** to pLDDT over the 6 cohorts whose epitope has a solved
        complex to template on (*n* = 284) — and **subtracts about 0.04** over the 16 that do not
        (*n* = 743).

        That is the same template covariate everything else in this framework divides under. Where
        no receptor has been co-crystallized with the peptide, nothing works, the generator's own
        confidence included. Read the correction as an improvement for epitopes with structural
        precedent, and report the template covariate beside it rather than inferring it.
        """
    )
    return


if __name__ == "__main__":
    app.run()
