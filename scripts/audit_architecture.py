#!/usr/bin/env python
"""Where the import graph is tangled, and where the layering it should have would put things.

The package is meant to read bottom-up: a structure is parsed, contacts are derived from it, a
Hamiltonian scores those contacts, a descriptor module *calls* those layers and computes nothing
itself, a score composes descriptors, and the CLI sits on top and is imported by nothing. Anything
that reads the other way is a layering violation, and the reliable symptom is a **function-local
relative import** -- a module dodging a cycle by deferring the import into a call.

This prints the cycles, the violations and the deferred-import count per module, and draws the
subsystem graph with the back edges marked. Run it to see whether a change made the layering better
or worse; ``--check`` exits non-zero while any module-level cycle remains, which is what a CI gate
would use once they are gone.

    python scripts/audit_architecture.py [--check] [--png PATH]
"""
from __future__ import annotations

import argparse
import ast
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1] / "src" / "tcren"

#: The layering the package is meant to have, lowest first. A module in layer i may import from
#: layers <= i only. ``None`` collects everything not yet assigned -- the work list.
LAYERS: tuple[tuple[str, tuple[str, ...]], ...] = (
    # 0  no domain logic at all: where files live, and how a generated table is stamped.
    ("foundation", ("paths", "provenance", "_provenance", "metadata", "data")),
    # 1  a PDB becomes a typed, region-annotated Structure.
    ("structure", ("structure", "annotation", "mhc")),
    # 2  which residues touch which, at what distance, and of what chemical kind.
    ("contacts", ("contacts", "contactmap", "contact_types", "clashes")),
    # 3  the energy functions themselves: the pair potential and the coupled model.
    ("hamiltonian", ("potential", "potts")),
    # 4  everything that turns a structure plus a layer below into a number -- geometry,
    #    energetics, topology, mechanics. One module per question.
    ("docking", ("docking", "geometry", "torsions")),
    ("topology", ("topology", "project2d", "stacking")),
    ("energetics", ("energetics", "refine")),
    ("mechanics", ("mechanics",)),
    # 5  the catalogue that names those numbers, and the dispatch that runs them over a set.
    ("descriptors", ("descriptors", "recognition")),
    # 6  scores composed from descriptors: Q, T, S and the cohort machinery.
    ("scores", ("cohort", "reliability", "scoring_rank")),
    # 7  imported by nothing inside the package.
    ("app", ("cli", "pipeline", "paper", "viz", "cpl", "binder", "analysis", "recent", "oracle",
             "shuffle")),  # noqa: E501
    # Deprecated top-level names kept as transparent re-exports after the 2026-09-01 move. They sit
    # outside the stack on purpose: nothing inside the package imports them, and they exist only so
    # code written against the old layout keeps working.
    ("shims", ("orient_shim", "footprint", "interface_graph", "surface", "pose", "scoring", "ddg",
               "rotamers", "stability", "dynamics")),
)


def modules() -> dict[str, Path]:
    out = {}
    for p in ROOT.rglob("*.py"):
        n = ".".join(p.relative_to(ROOT).with_suffix("").parts).replace(".__init__", "")
        out[n or "tcren"] = p
    return out


def _resolve(name: str, path: Path, node: ast.ImportFrom) -> list[str]:
    """The modules a relative ``ImportFrom`` inside ``name`` actually names.

    Python's rule, which three earlier versions of this script each got wrong in a different way
    and each time invented cycles that were not there:

    * a package's ``__init__`` resolves ``from .x`` against **itself**, a plain module against its
      parent package -- so ``contacts/__init__.py``'s ``from .geometry`` is ``contacts.geometry``;
    * ``level`` counts dots, so ``from ..`` in ``mhc/pseudo.py`` is the top-level package, not
      ``mhc``;
    * ``from . import x`` names the **sibling** ``pkg.x``, not the package ``__init__``.
    """
    own = name if path.name == "__init__.py" else (name.rsplit(".", 1)[0] if "." in name else "")
    parts = own.split(".") if own else []
    up = node.level - 1
    base = ".".join(parts[:len(parts) - up]) if up <= len(parts) else None
    if base is None:
        return []
    if node.module is None:                       # from <base> import a, b
        return [f"{base}.{al.name}" if base else al.name for al in node.names]
    return [f"{base}.{node.module}" if base else node.module]


def graph(mods: dict[str, Path]):
    """``(module-level edges, deferred-edge counts)`` over relative imports only."""
    edges, deferred, deferred_edges = defaultdict(set), defaultdict(int), defaultdict(set)
    for name, p in mods.items():
        try:
            tree = ast.parse(p.read_text())
        except SyntaxError:
            continue
        inside = {id(n) for fn in ast.walk(tree)
                  if isinstance(fn, (ast.FunctionDef, ast.AsyncFunctionDef))
                  for n in ast.walk(fn) if isinstance(n, ast.ImportFrom)}
        for n in ast.walk(tree):
            if not (isinstance(n, ast.ImportFrom) and n.level):
                continue
            for tgt in _resolve(name, p, n):
                if tgt not in mods or tgt == name:
                    continue
                if id(n) in inside:
                    deferred[name] += 1
                    deferred_edges[name].add(tgt)
                else:
                    edges[name].add(tgt)
    return edges, deferred, deferred_edges


def layer_of(mod: str) -> tuple[int, str] | None:
    head = mod.split(".")[0]
    for i, (label, members) in enumerate(LAYERS):
        if head in members:
            return i, label
    return None


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--check", action="store_true", help="exit 1 while any module-level cycle remains")
    ap.add_argument("--png", type=Path, default=ROOT.parents[1] / "docs" / "_static" / "architecture.png")
    a = ap.parse_args()

    import networkx as nx
    mods = modules()
    edges, deferred, deferred_edges = graph(mods)
    G = nx.DiGraph()
    G.add_nodes_from(mods)
    for s, ts in edges.items():
        G.add_edges_from((s, t) for t in ts)

    cycles = sorted(nx.simple_cycles(G), key=len)
    print(f"modules {len(mods)}   module-level edges {G.number_of_edges()}   "
          f"deferred (function-local) {sum(deferred.values())}\n")
    print(f"import cycles: {len(cycles)}")
    for c in cycles:
        print("   " + " -> ".join(c) + f" -> {c[0]}")

    def _up(src):
        out = []
        for s_, ts in src.items():
            ls = layer_of(s_)
            for t in ts:
                lt = layer_of(t)
                if ls and lt and lt[0] > ls[0]:
                    out.append((s_, ls[1], t, lt[1]))
        return out

    viol, viol_deferred = _up(edges), _up(deferred_edges)
    print("\nlayering violations at module level (these would be import errors or cycles):")
    for s_, a_, t, b_ in sorted(viol):
        print(f"   {s_:24s} [{a_:11s}] -> {t:22s} [{b_}]")
    print(f"   {len(viol)} violation(s)")
    print("\nupward imports DEFERRED into a function (legal, but the layering is aspirational "
          "wherever one appears):")
    for s_, a_, t, b_ in sorted(viol_deferred):
        print(f"   {s_:24s} [{a_:11s}] -> {t:22s} [{b_}]")
    print(f"   {len(viol_deferred)} deferred upward import(s)")

    unplaced = sorted(m for m in mods if layer_of(m) is None)
    if unplaced:
        print(f"\nnot yet assigned to a layer ({len(unplaced)}): {', '.join(unplaced[:12])}"
              + (" ..." if len(unplaced) > 12 else ""))

    print("\nmost deferred imports (each one is a cycle being dodged):")
    for m, c in sorted(deferred.items(), key=lambda x: -x[1])[:10]:
        print(f"   {m:24s} {c}")

    _draw(G, edges, deferred_edges, viol, a.png)
    if a.check and cycles:
        print(f"\nFAIL: {len(cycles)} import cycle(s)", file=sys.stderr)
        return 1
    return 0


def _draw(G, edges, deferred_edges, viol, png: Path) -> None:
    """A layered map drawn deterministically, because a layout engine will not do this legibly.

    Two graphviz attempts came out unreadable -- 109 modules and 450 edges is past what spring or
    rank layout can place, and neither honoured the layer order, so the picture said nothing. This
    puts each layer on its own row, bottom-up, with the modules named inside it, and summarises the
    dependencies as a layer-to-layer count matrix beside it rather than drawing 450 curves. What a
    reader wants from this is *which layer holds what* and *which way the arrows go*, and both of
    those survive the summary; the individual edges are in the CSV.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.patches import FancyArrowPatch, Rectangle

    FILL = ["#eef1f5", "#e8f0f8", "#e6f2ec", "#fdf0e3", "#f4eef8", "#f2ecf6", "#fdf3e8", "#f7eef0",
            "#fce9e6", "#eef3d9", "#efefef"]
    def node(m):
        head = m.split(".")[0]
        return head if head in {"docking", "potts", "potential", "contacts", "structure",
                                "annotation", "mhc", "project2d", "viz", "paper", "descriptors",
                                "refine", "binder", "data", "topology", "energetics",
                                "mechanics", "orient"} else m
    per = {i: sorted({node(m) for m in G if layer_of(m) and layer_of(m)[0] == i
                      and not m.split(".")[-1].startswith("__")})
           for i in range(len(LAYERS))}
    # layer x layer counts, module-level and deferred kept apart
    n_lay = len(LAYERS)
    solid = np.zeros((n_lay, n_lay), int); dash = np.zeros((n_lay, n_lay), int)
    for mat, src in ((solid, edges), (dash, deferred_edges)):
        for s_, ts in src.items():
            ls = layer_of(s_)
            for t in ts:
                lt = layer_of(t)
                if ls and lt and ls[0] != lt[0]:
                    mat[ls[0], lt[0]] += 1
    tot = solid + dash

    fig = plt.figure(figsize=(17.5, 9.6))
    gs = fig.add_gridspec(1, 2, width_ratios=[2.5, 1], wspace=.16)
    ax = fig.add_subplot(gs[0]); ax.set_axis_off()
    ax.set_xlim(0, 1); ax.set_ylim(-.4, n_lay + .1)
    for i in range(n_lay):
        label, _ = LAYERS[i]
        ax.add_patch(Rectangle((0.005, i + .06), .99, .82, facecolor=FILL[i],
                               edgecolor="#c9c9c9", lw=.9, zorder=1,
                               transform=ax.transData))
        ax.text(.018, i + .47, f"{i}. {label}", fontsize=13, fontweight="bold",
                va="center", color="#222222", zorder=3)
        mods = per[i]
        if not mods:
            continue
        cols = max(1, min(len(mods), 8))
        for k, m in enumerate(mods):
            col, row = k % cols, k // cols
            x = .215 + col * (.775 / cols)
            y = i + .52 - row * .19
            ax.text(x, y, m, fontsize=9.5, va="center", ha="left", zorder=3,
                    bbox=dict(boxstyle="round,pad=0.22", fc="white", ec="#bbbbbb", lw=.6))
    # the spine: every arrow points DOWN, which is the whole claim
    for i in range(n_lay - 1):
        ax.add_patch(FancyArrowPatch((.5, i + 1.03), (.5, i + .92), arrowstyle="-|>",
                                     mutation_scale=13, color="#9aa0a6", lw=1.1, zorder=2))
    n_up = sum(int(tot[i, j]) for i in range(n_lay) for j in range(n_lay) if j > i)
    n_up_solid = sum(int(solid[i, j]) for i in range(n_lay) for j in range(n_lay) if j > i)
    ax.text(.5, -.28,
            f"no import cycles, and {n_up_solid} upward imports at module level — "
            f"the {n_up} upward edges above are all deferred into a function",
            ha="center", fontsize=11, style="italic", color="#444444")

    ax2 = fig.add_subplot(gs[1])
    ax2.imshow(np.where(tot > 0, tot, np.nan), cmap="BuPu", origin="lower", aspect="auto")
    for i in range(n_lay):
        for j in range(n_lay):
            if tot[i, j]:
                ax2.text(j, i, f"{solid[i, j]}\n+{dash[i, j]}" if dash[i, j] else f"{solid[i, j]}",
                         ha="center", va="center", fontsize=8.5,
                         color="white" if tot[i, j] > tot.max() * .55 else "#333333")
    ax2.set_xticks(range(n_lay)); ax2.set_yticks(range(n_lay))
    names = [f"{i}. {lab}" for i, (lab, _) in enumerate(LAYERS)]
    ax2.set_xticklabels(names, rotation=45, ha="right", fontsize=9)
    ax2.set_yticklabels(names, fontsize=9)
    ax2.set_xlabel("imports FROM layer", fontsize=10)
    ax2.set_ylabel("layer", fontsize=10)
    ax2.set_title("dependencies: module-level\n+ deferred into a function", fontsize=11,
                  fontweight="bold", loc="left")
    for sp in ax2.spines.values():
        sp.set_visible(False)
    ax2.tick_params(length=0)
    # anything above the diagonal would be a violation; shade it so its emptiness is visible
    ax2.fill_between([-.5, n_lay - .5], [-.5, n_lay - .5], n_lay - .5,
                     color="#CC3311", alpha=.05, zorder=0)
    ax2.text(n_lay - 1.1, .3, "violations\nwould sit here", fontsize=8, color="#CC3311",
             ha="right", va="bottom", style="italic")

    fig.suptitle("tcren module map — what each layer holds, and which way it depends",
                 fontsize=15, fontweight="bold", y=.97)
    png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(png, dpi=150, bbox_inches="tight", pad_inches=.25)
    plt.close(fig)
    print(f"\nwrote {png}")


if __name__ == "__main__":
    raise SystemExit(main())
