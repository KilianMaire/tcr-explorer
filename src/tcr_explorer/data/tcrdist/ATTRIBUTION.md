# tcrdist germline CDR reference table

`alphabeta_gammadelta_db.tsv` is the V and J gene reference table used to look up
the germline CDR1, CDR2, and CDR2.5 (pMHC) loops for each TCR V gene. TCR Explorer
uses it to reproduce the authoritative tcrdist metric offline, together with the
`pwseqdist` distance engine.

## Source and license

The table is vendored verbatim from tcrdist3.

- **tcrdist3**, Copyright (c) Koshlan Mayer-Blackwell, released under the MIT License.
  <https://github.com/kmayerb/tcrdist3> (`tcrdist/db/alphabeta_gammadelta_db.tsv`).
  Mayer-Blackwell K. et al. TCR meta-clonotypes for biomarker discovery with tcrdist3.
  eLife, 2021.

The MIT License permits redistribution provided the copyright notice and permission
notice are retained. The full MIT terms are reproduced below.

The CDR sequences in the table derive from the IMGT reference directory. IMGT is
distributed under CC BY 4.0 (Lefranc M-P. et al., IMGT, the international ImMunoGeneTics
information system, <https://www.imgt.org>). The same attribution applies to the
bundled germline in `../germline/ATTRIBUTION.md`.

## MIT License (tcrdist3)

```
MIT License

Copyright (c) Koshlan Mayer-Blackwell

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
