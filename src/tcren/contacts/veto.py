"""Contact graph diagnostics with a supplied development cell-effect reference."""

from __future__ import annotations

from collections.abc import Mapping
from math import isfinite

import polars as pl

from ..topology.graph import _contact_graph_counts

_GRAPHS = ("full", "kept", "vetoed")
_INTEGER_FIELDS = ("b0", "b1", "n_edges", "n_nodes")
_FRACTIONS = ("kept_edge_frac", "kept_b1_frac", "veto_b1_frac", "core_size")


def contact_veto_features(
    opportunities: pl.DataFrame, effects: Mapping[str, float],
) -> pl.DataFrame:
    """Summarize full, positive-effect and vetoed residue contact graphs.

    Supply directed rows from :func:`~tcren.contacts.opportunities.native_opportunities`
    and development estimates with keys such as ``tcr_peptide/A|G``. Each graph
    has only its incident vertices. Missing cell effects count as vetoes and
    are counted separately. A structure without contacts has no graph row.
    """
    if not effects or not all(isfinite(value) for value in effects.values()):
        raise ValueError("supply a nonempty finite development cell-effect reference")
    cells = pl.concat_str(
        "interface", pl.concat_str("aa.from", "aa.to", separator="|"), separator="/")
    contacts = opportunities.filter("contacted").with_columns(
        cells.replace_strict(dict(effects), default=None, return_dtype=pl.Float64).alias("effect"),
    )
    rows = []
    for (id_,), group in contacts.group_by("pdb.id"):
        kept, vetoed, unknown = set(), set(), set()
        for a, i, b, j, effect in group.select(
            "chain.id.from", "residue.index.from", "chain.id.to", "residue.index.to", "effect",
        ).iter_rows():
            edge = tuple(sorted(((a, i), (b, j))))
            (kept if effect is not None and effect > 0 else vetoed).add(edge)
            if effect is None:
                unknown.add(edge)
        row = {"id": id_, "n_unknown_edges": len(unknown)}
        for name, edges in zip(_GRAPHS, (kept | vetoed, kept, vetoed), strict=True):
            n_edges, n_nodes, b0, b1 = _contact_graph_counts(edges)
            row.update({f"{name}_b0": b0, f"{name}_b1": b1,
                        f"{name}_n_edges": n_edges, f"{name}_n_nodes": n_nodes,
                        f"{name}_cyclo_frac": b1 / max(n_edges, 1)})
        row.update(
            kept_edge_frac=row["kept_n_edges"] / row["full_n_edges"],
            kept_b1_frac=row["kept_b1"] / max(row["full_b1"], 1),
            veto_b1_frac=row["vetoed_b1"] / max(row["full_b1"], 1),
            core_size=row["kept_n_nodes"] / row["full_n_nodes"],
        )
        rows.append(row)
    schema = {
        "id": pl.String, "n_unknown_edges": pl.Int64,
        **{f"{g}_{c}": pl.Int64 for g in _GRAPHS for c in _INTEGER_FIELDS},
        **{f"{g}_cyclo_frac": pl.Float64 for g in _GRAPHS},
        **dict.fromkeys(_FRACTIONS, pl.Float64),
    }
    return pl.DataFrame(rows, schema=schema).sort("id")
