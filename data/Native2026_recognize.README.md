# Native2026_recognize.tsv

`tcren recognize --full` over the 374 `Native2026` reference crystals: 374 rows, 63 columns, keyed
on `complex.id`. Regenerate with

```bash
tcren recognize --full --scores -s "$TCREN_DATA_DIR/Native2026" -o Native2026_recognize.tsv
```

This retained annotation table records the reference set, including `mhc_class_bin`:
280 class-I structures and 94 class-II structures.

Everything in it is **computed**, not measured. Preserve it as a historical
output; it does not establish results for a changed descriptor implementation or
potential. Run the command above only when a fresh calculation is requested,
and keep the new output separate from this retained table.
