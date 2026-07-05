# Third-party notices

`tcren` is distributed under GPL-3.0-or-later. It includes code derived from the following third-party
projects, whose licenses are reproduced below.

## TCRdock (Bradley lab)

`src/tcren/orient/tcrdock_geometry.py` is a native reimplementation of the TCR:pMHC docking-geometry
computation from **TCRdock** (https://github.com/phbradley/TCRdock), commit
`c5a7af42eeb0c2a4492a4d4fe803f1f9aafb6193` (2024-03-04) — specifically the rigid-body "docking geometry"
of `tcrdock/docking_geometry.py`, the MHC/TCR symmetry stubs of `tcrdock/mhc_util.py` and
`tcrdock/tcr_util.py`, and the geometry helpers of `tcrdock/superimpose.py` and `tcrdock/geom_util.py`.
No TCRdock source is copied verbatim; the algorithm is ported to tcren's own structure model and annotation.
The class-I template sequence and β-sheet core positions, and the conserved IMGT TCR core positions, are
constants taken from TCRdock.

Reference: Bradley P. *Structure-based prediction of T cell receptor:peptide-MHC interactions.* eLife 2023;12:e82813.

TCRdock is licensed MIT:

```
MIT License

Copyright (c) 2022 Philip Bradley

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```
