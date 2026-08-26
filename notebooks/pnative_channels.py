"""P(native) — the three descriptor channels, and why they are read separately.

A marimo notebook — reactive, so changing the structure set or the channel selection redraws only
what depends on it, and plain Python, so it diffs and imports like any other module.

    marimo edit notebooks/pnative_channels.py     # explore and change the code
    marimo run  notebooks/pnative_channels.py     # read-only app, controls still live

This is the released scoring path, not a re-implementation of it: the same two calls a shell user
makes as ``tcren features`` then ``tcren recognize --features``. Needs a directory of TCR:pMHC
structures, the MHC allele reference (``tcren build-mhc-ref``, once) and the canonical database
(``tcren fetch-data``, or ``$TCREN_DATA_DIR`` pointed at a local copy). Install with
``pip install "tcren[marimo]"``.

The point of the notebook is the last section. The energetics channel is the one that changes sign
between cohorts — favourable contacts read backwards on an epitope the generator could copy from a
template — and the reason P(native) does not care is that each channel's class coefficient is
fitted rather than asserted. That is visible here as the sign of one number.
"""

import marimo

__generated_with = "0.23.15"
app = marimo.App(width="medium", app_title="tcren · P(native) channels")


@app.cell
def _():
    import marimo as mo
    return (mo,)


@app.cell
def _(mo):
    mo.md(
        r"""
        # tcren · P(native), channel by channel

        A generative structure model builds a confident complex for any receptor–peptide pair,
        binding or not. What separates a real interface from a manufactured one is not one number
        but several, and they are not interchangeable: where the receptor sits, what shape its
        contact set has, and what the contact chemistry is worth are invariant under different
        things.

        `tcren` splits its descriptor catalogue on exactly that criterion and fits **one
        latent-class Bayesian network per channel**. The class node is unobserved, so the fit takes
        no binder label from the cohort it scores, and the channel posteriors combine by adding
        log-odds:

        \[
        \operatorname{logit} P_{\mathrm{native}}(x)
          = \sum_{c \in \{G,\,T,\,E\}} \operatorname{logit} P_c(x) - (C-1)\operatorname{logit}\pi .
        \]

        Adding log-odds is the exact posterior only across channels that are conditionally
        independent given the class, which is why *placement* and *interface* — the most dependent
        pair measured — are pooled into one geometry network rather than counted twice.
        """
    )
    return


@app.cell
def _(mo):
    from tcren import paths

    default_dir = paths.data_dir() / "Canonical2026"
    struct_dir = mo.ui.text(
        value=str(default_dir),
        label="structure directory (a cohort — one epitope, or a mixed set)",
        full_width=True,
    )
    limit = mo.ui.slider(20, 400, value=120, step=20, label="structures to score")
    mo.vstack([struct_dir, limit])
    return limit, paths, struct_dir


@app.cell
def _(limit, mo, struct_dir):
    from pathlib import Path

    from tcren.structure.io import STRUCTURE_SUFFIXES

    d = Path(struct_dir.value).expanduser()
    files = sorted(
        p for p in d.glob("*") if any(str(p).endswith(s) or str(p).endswith(s + ".gz")
                                      for s in STRUCTURE_SUFFIXES)
    )[: limit.value]
    mo.stop(
        not files,
        mo.md(
            f"**No structures under `{d}`.** Run `tcren fetch-data` to populate the canonical "
            "database, or point the box above at your own directory."
        ),
    )
    mo.md(f"`{len(files)}` structures from `{d.name}`.")
    return (files,)


@app.cell
def _(mo):
    mo.md(
        r"""
        ## One pass over the coordinates

        `recognition_table` is what `tcren features` calls. `include` selects descriptor
        **families**; the four below are the ones the three channels draw on. Everything is read
        from the coordinates: nothing here consumes a confidence emitted by the structure
        generator, and nothing consumes a binding label.
        """
    )
    return


@app.cell
def _(files, mo):
    import polars as pl

    from tcren.recognition import recognition_table

    with mo.status.spinner(title="annotating and featurising…"):
        feats = pl.DataFrame(
            recognition_table(
                files,
                include=("placement", "interface", "topology", "energetics"),
                with_p_real=False,
                threads=0,
            )
        )
    mo.md(f"`{feats.height}` rows × `{feats.width}` descriptors.")
    return feats, pl


@app.cell
def _(mo):
    mo.md(
        r"""
        ## The channels, fitted one at a time

        Each channel is fitted on the cohort being scored. `P_NATIVE_POOL` records which descriptor
        families feed which channel, and `P_NATIVE_ORIENT` the feature whose mean orients the
        mixture — a finite mixture is identified only up to permutation of its components, so
        without an orientation rule two runs of the same data can disagree about which side is
        native.
        """
    )
    return


@app.cell
def _(feats, mo, pl):
    from tcren.cohort import P_NATIVE_CHANNELS, P_NATIVE_ORIENT, P_NATIVE_POOL, p_native

    with mo.status.spinner(title="fitting each channel by EM…"):
        per_channel = {
            ch: p_native(feats, channels=(ch,)) for ch in P_NATIVE_CHANNELS
        }
        combined = p_native(feats)

    table = pl.DataFrame(
        {
            "channel": list(P_NATIVE_CHANNELS) + ["P_native"],
            "families": [", ".join(P_NATIVE_POOL[c]) for c in P_NATIVE_CHANNELS] + ["all three"],
            "orients on": [P_NATIVE_ORIENT[c] for c in P_NATIVE_CHANNELS] + ["—"],
            "mean posterior": [float(per_channel[c].mean()) for c in P_NATIVE_CHANNELS]
            + [float(combined.mean())],
        }
    )
    table
    return P_NATIVE_CHANNELS, combined, per_channel


@app.cell
def _(mo):
    mo.md(
        r"""
        ## The one that changes sign

        In a genuine complex the geometry and energetics channels are tied, because a larger,
        better-packed interface holds more favourable contacts. A pose the generator manufactured
        breaks that tie, and the correlation below is where it shows: positive on a cohort whose
        interfaces the physics produced, negative on one the generator copied from a near-native
        template and then filled with the wrong sequence.

        This is a **diagnostic, not a gate**. Nothing downstream is multiplied by it. The reason
        `P_native` survives a cohort where the energy runs backwards is that expectation-maximization
        gives that channel a negative class coefficient and leaves the other two alone.
        """
    )
    return


@app.cell
def _(P_NATIVE_CHANNELS, combined, mo, per_channel, pl):
    import itertools

    import numpy as np

    rows = []
    for a, b in itertools.combinations(P_NATIVE_CHANNELS, 2):
        r = float(np.corrcoef(per_channel[a], per_channel[b])[0, 1])
        rows.append({"channel A": a, "channel B": b, "Pearson R": r})
    mo.vstack(
        [
            pl.DataFrame(rows),
            mo.md(
                f"`P_native` over `{len(combined)}` structures: "
                f"mean `{combined.mean():.3f}`, "
                f"interquartile range `{np.percentile(combined, 25):.3f}`–"
                f"`{np.percentile(combined, 75):.3f}`.\n\n"
                "A cohort of one has no cohort to fit, so `P_native` is defined on a **set**. "
                "For a score defined on a single structure, use `q_score`, which is standardised "
                "against the native-crystal reference and carries no fitted coefficient."
            ),
        ]
    )
    return


if __name__ == "__main__":
    app.run()
