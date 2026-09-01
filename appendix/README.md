# The descriptor whitepaper

`descriptors.pdf` documents every descriptor `tcren` computes: the reduction chain the catalogue
sits on, the nine operators that generate it, the formulas and their derivations, the MHC class
I/II conditioning, and the full catalogue with units and `STATUS` flags.

```
make          # regenerate generated/, then build descriptors.pdf
make check    # exit 1 if generated/ is stale against the catalogue
```

`descriptors.tex` is written by hand. Everything under `generated/` is emitted by
`scripts/gen_appendix.py` from `tcren.recognition`, so a descriptor added to the catalogue reaches
the whitepaper by re-running the generator rather than by being remembered. Do not edit
`generated/*`.
