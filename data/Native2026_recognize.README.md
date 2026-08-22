# Native2026_recognize.tsv

`tcren recognize --full` over the 374 `Native2026` reference crystals: 374 rows, 63 columns, keyed
on `complex.id`. Regenerate with

```bash
tcren recognize --full --scores -s "$TCREN_DATA_DIR/Native2026" -o Native2026_recognize.tsv
```

It is kept in the repository because it is the only per-structure annotation table for the reference
set that is small enough to commit, and several analyses need one column of it without wanting to
re-featurise 374 structures — `mhc_class_bin` in particular, which is how the reference set is split
into class I (280) and class II (94).

Everything in it is **computed**, not measured. It is a convenience cache of this library's own
output, so it must be regenerated whenever the features or the potential change; nothing should
treat it as an input.
