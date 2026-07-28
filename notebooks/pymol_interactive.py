"""Interactive PyMOL render explorer for canonically-oriented TCR–pMHC complexes.

A marimo notebook — reactive, so changing a control re-renders only what depends on it, and
plain Python, so it diffs and imports like any other module.

    marimo edit notebooks/pymol_interactive.py     # explore and change the code
    marimo run  notebooks/pymol_interactive.py     # read-only app, controls still live

Every ray-traced panel comes from `tcren.viz.pymol`, the same library code the static gallery
notebook and any script would use. Renders take a second or two each, so `render_cached` is
content-addressed: a combination you have already looked at comes back instantly.

Needs a `pymol` binary on PATH for the ray-traced sections and `py3Dmol` for the live 3D one;
each section says so and degrades to a message rather than a traceback.
"""

import marimo

__generated_with = "0.23.15"
app = marimo.App(width="medium", app_title="tcren · PyMOL render explorer")


@app.cell
def _():
    import marimo as mo
    return (mo,)


@app.cell
def _(mo):
    mo.md(
        r"""
        # tcren · PyMOL render explorer

        Ray-traced views of **canonically-oriented** TCR–pMHC complexes, driven by
        `tcren.viz.pymol`. Each panel carries a thin, arrow-headed **axis gizmo** naming the
        canonical frame, because an oriented structure is not interpretable unless the reader
        can tell which way it points:

        | axis | label | direction |
        |---|---|---|
        | x | `width` | groove width, across the cleft (α1↔α2) |
        | y | `N→C` | groove axis, toward the peptide C-terminus |
        | z | `TCR` | docking normal, MHC floor → TCR |

        An axis pointing at the viewer foreshortens to a dot, and its label drops to the lower
        left of it — the convention for an axis normal to the page.
        """
    )
    return


@app.cell
def _(mo):
    import hashlib
    import json
    import shutil
    from pathlib import Path

    # marimo is reactive, so a name may be defined by exactly one cell. These are imported once
    # here and passed into the cells that need them rather than re-imported per cell.
    from tcren.annotation import classify_chains
    from tcren.mhc import annotate_mhc
    from tcren.structure import import_structure
    from tcren.viz.pymol import (
        CANONICAL_AXES,
        CHAIN_COLOURS,
        CORNERS,
        PALETTES,
        groove_scene,
        importance_scene,
        interface_scene,
        overlay_scene,
        render,
        residue_importance,
    )

    HERE = Path(__file__).resolve().parent
    CANON = HERE.parent / "data" / "Canonical2026"
    NATIVE = HERE.parent / "data" / "Native2026"
    META = HERE.parent / "data" / "orient_metadata.json"
    CACHE = Path("/tmp/tcren_marimo_renders")
    CACHE.mkdir(exist_ok=True)

    HAVE_PYMOL = shutil.which("pymol") is not None
    HAVE_DATA = CANON.is_dir() and any(CANON.glob("*.pdb.gz"))

    def render_cached(scene: str, **kw):
        """Ray-trace `scene`, reusing an earlier identical render.

        Keyed on the scene text plus every render option, so the cache can never serve a panel
        that was made with different settings — the usual failure of a "cache by structure id"
        shortcut.
        """
        key = hashlib.sha1(
            (scene + json.dumps(kw, sort_keys=True, default=str)).encode()
        ).hexdigest()[:16]
        png = CACHE / f"{key}.png"
        if not png.exists():
            render(scene, png, **kw)
        return png

    mo.md(
        f"""
        **Environment** — pymol on PATH: `{HAVE_PYMOL}` · Canonical2026 present: `{HAVE_DATA}`
        {"" if HAVE_DATA else "· run `tcren fetch-data` to populate `data/`"}
        """
    )
    return (
        CACHE,
        CANON,
        CANONICAL_AXES,
        CHAIN_COLOURS,
        CORNERS,
        HAVE_DATA,
        HAVE_PYMOL,
        META,
        NATIVE,
        PALETTES,
        Path,
        annotate_mhc,
        classify_chains,
        groove_scene,
        importance_scene,
        import_structure,
        interface_scene,
        json,
        overlay_scene,
        render_cached,
        residue_importance,
    )


@app.cell
def _(CANON, HAVE_DATA, META, json):
    import collections

    # Structures grouped by MHC class x species, so the pickers below offer real choices.
    by_group = collections.defaultdict(list)
    all_ids = []
    if HAVE_DATA and META.exists():
        for _r in json.load(open(META)):
            if _r.get("status") == "ok" and (CANON / f"{_r['pdb.id']}.pdb.gz").exists():
                by_group[(_r["mhc.class"], _r["species"])].append(_r["pdb.id"])
                all_ids.append(_r["pdb.id"])
    all_ids = sorted(all_ids)
    GROUPS = [("MHCI", "Human"), ("MHCII", "Human"), ("MHCI", "Mouse"), ("MHCII", "Mouse")]
    return GROUPS, all_ids, by_group


@app.cell
def _(CANONICAL_AXES, CHAIN_COLOURS, mo):
    mo.accordion(
        {
            "What the axis labels mean": mo.md(
                "\n".join(
                    [f"- **`{a.short}`** ({a.letter}, *{a.label}*) — {a.definition}"
                     for a in CANONICAL_AXES]
                )
                + "\n\nSame three directions the docking-geometry literature uses (SwiftTCR, "
                "TCR3d); the principal-component *ranking* differs because "
                "`tcren.orient.frame` fits the whole complex where those fit the MHC groove "
                "alone."
            ),
            "Chain roles after `tcren orient`": mo.md(
                "\n".join(f"- **{cid}** — {role} (`{colour}`)"
                          for cid, (role, colour) in CHAIN_COLOURS.items())
            ),
        }
    )
    return


@app.cell
def _(mo):
    mo.md(r"""## 1 · Scene explorer — the three view families""")
    return


@app.cell
def _(all_ids, mo):
    pick_id = mo.ui.dropdown(
        options=all_ids or ["<no structures>"],
        value=("1ao7" if "1ao7" in all_ids else (all_ids[0] if all_ids else "<no structures>")),
        label="structure",
    )
    pick_scene = mo.ui.radio(
        options=["groove (top-down)", "interface (side-on)", "overlay (ensemble, side-on)"],
        value="groove (top-down)",
        label="scene",
    )
    pick_surface = mo.ui.checkbox(label="molecular surface (histo.fyi look; groove scene only)")
    pick_size = mo.ui.slider(400, 1400, step=100, value=800, label="pixels", show_value=True)
    mo.hstack([mo.vstack([pick_id, pick_scene]),
               mo.vstack([pick_surface, pick_size])], justify="start", gap=2)
    return pick_id, pick_scene, pick_size, pick_surface


@app.cell
def _(
    CANON,
    HAVE_DATA,
    HAVE_PYMOL,
    NATIVE,
    by_group,
    classify_chains,
    groove_scene,
    import_structure,
    interface_scene,
    mo,
    overlay_scene,
    pick_id,
    pick_scene,
    pick_size,
    pick_surface,
    render_cached,
):
    def cdr_resi(pid):
        """CDR residue numbers, read off the pre-orientation structure.

        PDB numbering survives orientation, so the same numbers select chain A/B in the
        oriented file. Returns empty lists if the structure cannot be annotated, which
        `interface_scene` renders as "no loops" rather than failing.
        """
        out = {"TRA": [], "TRB": []}
        try:
            s = import_structure(NATIVE / f"{pid}.pdb.gz", pdb_id=pid)
            classify_chains(s, organism="human")
            for c in s.chains:
                if c.chain_type in out:
                    for reg in c.regions:
                        if reg.region_type.startswith("CDR"):
                            out[c.chain_type] += [r.pdb_index for r in reg.residues]
        except Exception:
            pass
        return out

    if not (HAVE_PYMOL and HAVE_DATA):
        scene_out = mo.callout(
            mo.md("Needs a `pymol` binary on PATH and `data/Canonical2026` "
                  "(`tcren fetch-data`)."),
            kind="warn",
        )
    else:
        _pid = pick_id.value
        if pick_scene.value.startswith("groove"):
            _scene = groove_scene(_pid, CANON, surface=pick_surface.value)
        elif pick_scene.value.startswith("interface"):
            _scene = interface_scene(_pid, CANON, cdr_resi(_pid))
        else:
            # The ensemble the picked structure belongs to, so the choice above still steers it.
            _ids = next((v for v in by_group.values() if _pid in v), [_pid])
            _scene = overlay_scene(_ids, CANON)
        _png = render_cached(_scene, size=(pick_size.value, pick_size.value))
        scene_out = mo.vstack([
            mo.image(_png, width=pick_size.value),
            mo.md(f"`{pick_scene.value}` · **{_pid}** · {pick_size.value}px"),
        ])
    scene_out
    return (cdr_resi,)


@app.cell
def _(mo):
    mo.md(
        r"""
        ## 2 · Camera — watch the gizmo follow

        The triad reports the **world** axes, so turning the camera turns it. Swing the view and
        the labels track: an axis rotating toward the viewer shrinks to a dot and its label drops
        to the lower left of it.
        """
    )
    return


@app.cell
def _(mo):
    azim = mo.ui.slider(-180, 180, step=15, value=0, label="azimuth °",
                        debounce=True, show_value=True)
    elev = mo.ui.slider(-90, 90, step=15, value=0, label="elevation °",
                        debounce=True, show_value=True)
    roll = mo.ui.slider(-180, 180, step=15, value=0, label="roll °",
                        debounce=True, show_value=True)
    mo.vstack([azim, elev, roll])
    return azim, elev, roll


@app.cell
def _(
    CANON,
    HAVE_DATA,
    HAVE_PYMOL,
    azim,
    elev,
    groove_scene,
    mo,
    pick_id,
    render_cached,
    roll,
):
    if not (HAVE_PYMOL and HAVE_DATA):
        cam_out = mo.md("*(needs pymol + data)*")
    else:
        # `cmd.turn` after the preset's own set_view/zoom. The gizmo is built from the rotation
        # PyMOL actually ends up with, read back by `probe_rotation`, so it stays correct for
        # any combination of turns rather than only the presets.
        _scene = groove_scene(pick_id.value, CANON) + (
            f'\ncmd.turn("y", {azim.value})'
            f'\ncmd.turn("x", {elev.value})'
            f'\ncmd.turn("z", {roll.value})'
        )
        cam_out = mo.vstack([
            mo.image(render_cached(_scene, size=(700, 700)), width=700),
            mo.md(f"azimuth **{azim.value}°** · elevation **{elev.value}°** · roll **{roll.value}°**"),
        ])
    cam_out
    return


@app.cell
def _(mo):
    mo.md(r"""## 3 · Gizmo styling — corner, size, palette, labels""")
    return


@app.cell
def _(CORNERS, PALETTES, mo):
    g_corner = mo.ui.dropdown(options=sorted(CORNERS), value="bottom-left", label="corner")
    g_scale = mo.ui.slider(0.08, 0.40, step=0.01, value=0.19, label="tile / image width",
                           debounce=True, show_value=True)
    g_palette = mo.ui.radio(options=sorted(PALETTES), value="mono", label="palette")
    g_short = mo.ui.switch(value=True, label="short labels")
    g_on = mo.ui.switch(value=True, label="gizmo on")
    mo.hstack([mo.vstack([g_corner, g_palette]), mo.vstack([g_scale, g_short, g_on])],
              justify="start", gap=2)
    return g_corner, g_on, g_palette, g_scale, g_short


@app.cell
def _(
    CANON,
    HAVE_DATA,
    HAVE_PYMOL,
    g_corner,
    g_on,
    g_palette,
    g_scale,
    g_short,
    groove_scene,
    mo,
    pick_id,
    render_cached,
):
    if not (HAVE_PYMOL and HAVE_DATA):
        style_out = mo.md("*(needs pymol + data)*")
    else:
        style_out = mo.image(
            render_cached(
                groove_scene(pick_id.value, CANON),
                size=(700, 700),
                gizmo=g_on.value,
                corner=g_corner.value,
                gizmo_scale=g_scale.value,
                palette=g_palette.value,
                short_labels=g_short.value,
            ),
            width=700,
        )
    style_out
    return


@app.cell
def _(mo):
    mo.md(
        r"""
        ## 4 · Which residues carry the score

        Φ is a sum over residue–residue contacts, so it decomposes exactly: a residue's share is
        the sum of `φ(a_i, a_j)` over the contacts it makes across the interface. The total says
        how large the score is; this says what it is **made of**.

        `phi` is the energy share, **blue favourable / red unfavourable** on a ramp centred on
        zero so those words keep their meaning. `n_contacts` is the geometric share — how much of
        the interface the residue physically occupies. A residue can be large on one and small on
        the other, and that difference is usually the interesting part.
        """
    )
    return


@app.cell
def _(mo):
    imp_by = mo.ui.radio(options=["phi", "n_contacts"], value="phi", label="colour by")
    imp_regions = mo.ui.multiselect(
        options=["CDR1", "CDR2", "CDR3", "PEPTIDE"], value=["CDR3", "PEPTIDE"], label="regions"
    )
    mo.hstack([imp_by, imp_regions], justify="start", gap=2)
    return imp_by, imp_regions


@app.cell
def _(
    CANON,
    HAVE_DATA,
    HAVE_PYMOL,
    NATIVE,
    annotate_mhc,
    classify_chains,
    imp_by,
    imp_regions,
    importance_scene,
    import_structure,
    mo,
    pick_id,
    render_cached,
    residue_importance,
):
    if not (HAVE_PYMOL and HAVE_DATA):
        imp_out = mo.md("*(needs pymol + data)*")
    else:
        _s = import_structure(NATIVE / f"{pick_id.value}.pdb.gz", pdb_id=pick_id.value)
        classify_chains(_s, organism="human")
        annotate_mhc(_s)
        _imp = residue_importance(_s)
        _regions = tuple(imp_regions.value) or ("CDR3", "PEPTIDE")
        _spectrum = "blue_white_red" if imp_by.value == "phi" else "white_red"
        _png = render_cached(
            importance_scene(pick_id.value, CANON, _imp, by=imp_by.value,
                             regions=_regions, spectrum=_spectrum),
            size=(760, 760),
        )
        # The table beside the render, so a colour can be traced back to a number.
        _shown = _imp.filter(_imp["region.type"].is_in(list(_regions)))
        imp_out = mo.hstack(
            [mo.image(_png, width=560),
             mo.vstack([mo.md(f"**{pick_id.value}** · {len(_shown)} residues shown"),
                        mo.ui.table(_shown.to_dicts(), page_size=12, selection=None)])],
            justify="start", gap=1, align="start",
        )
    imp_out
    return


@app.cell
def _(mo):
    mo.md(
        r"""
        ## 5 · Live 3D — drag to rotate

        The sections above ray-trace, which is what a figure needs but not what exploring needs.
        This one is **3Dmol.js**: drag to rotate, scroll to zoom, no re-render. Same canonical
        orientation, the groove as a translucent surface with the peptide as sticks and the CDR
        loops traced — `tcren.viz.pocket3d.view_pocket_cdr`.
        """
    )
    return


@app.cell
def _(NATIVE, annotate_mhc, classify_chains, import_structure, mo, pick_id):
    try:
        import py3Dmol  # noqa: F401

        from tcren.viz import view_pocket_cdr

        _s = import_structure(NATIVE / f"{pick_id.value}.pdb.gz", pdb_id=pick_id.value)
        classify_chains(_s, organism="human")
        annotate_mhc(_s)
        _view = view_pocket_cdr(_s, surface=True, width=760, height=520)
        # `_make_html` is the standalone document 3Dmol.js needs; an iframe keeps its CDN script
        # and canvas from touching the rest of the page.
        live_out = mo.iframe(_view._make_html(), height=540)
    except ImportError:
        live_out = mo.callout(
            mo.md("Needs `py3Dmol` — `pip install 'tcren[viz]'`."), kind="warn"
        )
    except Exception as exc:  # noqa: BLE001 - a structure that will not annotate is not fatal here
        live_out = mo.callout(mo.md(f"Could not build the 3D view: `{exc}`"), kind="warn")
    live_out
    return


@app.cell
def _(mo):
    mo.md(
        r"""
        ## 6 · Contact sheet — one representative per MHC class × species

        Four ray-traces, so this one is behind a button rather than reacting to every change
        above.
        """
    )
    return


@app.cell
def _(mo):
    sheet_go = mo.ui.run_button(label="render the contact sheet")
    sheet_scene = mo.ui.radio(options=["groove", "interface", "overlay"], value="groove",
                              label="scene")
    mo.hstack([sheet_scene, sheet_go], justify="start", gap=2)
    return sheet_go, sheet_scene


@app.cell
def _(
    CANON,
    GROUPS,
    HAVE_DATA,
    HAVE_PYMOL,
    by_group,
    cdr_resi,
    groove_scene,
    interface_scene,
    mo,
    overlay_scene,
    render_cached,
    sheet_go,
    sheet_scene,
):
    if not sheet_go.value:
        sheet_out = mo.md("*(press the button)*")
    elif not (HAVE_PYMOL and HAVE_DATA):
        sheet_out = mo.md("*(needs pymol + data)*")
    else:
        _cells = []
        for _g in GROUPS:
            _ids = by_group.get(_g, [])
            if not _ids:
                continue
            if sheet_scene.value == "groove":
                _sc = groove_scene(_ids[0], CANON)
            elif sheet_scene.value == "interface":
                _sc = interface_scene(_ids[0], CANON, cdr_resi(_ids[0]))
            else:
                _sc = overlay_scene(_ids, CANON)
            _cells.append(mo.vstack([
                mo.image(render_cached(_sc, size=(520, 520)), width=380),
                mo.md(f"**{_g[0]} · {_g[1]}** — {_ids[0]} (n={len(_ids)})"),
            ]))
        sheet_out = mo.vstack([mo.hstack(_cells[:2], gap=1), mo.hstack(_cells[2:], gap=1)])
    sheet_out
    return


@app.cell
def _(CACHE, mo):
    mo.md(
        f"""
        ---
        Renders are cached under `{CACHE}`, keyed on the scene text and every render option, so
        revisiting a combination is instant and a changed option can never serve a stale panel.
        Delete that directory to force a re-render.
        """
    )
    return


if __name__ == "__main__":
    app.run()
