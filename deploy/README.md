---
title: TCR Explorer
emoji: 🧬
colorFrom: indigo
colorTo: blue
sdk: docker
app_port: 7860
pinned: false
license: mit
---

# TCR Explorer (public demo)

A germline-only public demo of [TCR Explorer](https://github.com/KilianMaire/tcr-explorer),
a federated T cell receptor analysis tool.

This instance serves the germline features that need no external record data:

- allele assignment (IMGT V and J, to the allele, with tie sets),
- full chain reconstruction from a V gene, a J gene, and a CDR3,
- germline gene alignment,
- CDR1 and CDR2 loops.

The record databases (VDJdb, IEDB, McPAS, TCR3d) are not included here. Their
licenses do not all permit public redistribution, so this demo ships none of
them. Record retrieval and neighbour search therefore return nothing on this
instance and say so.

For the full tool, with records and the tcrdist similarity engine, install it
locally:

```bash
pip install "tcr-explorer[tcrdist]"
tcr-explorer-refresh          # downloads the record databases into a local folder
```

Source (MIT): https://github.com/KilianMaire/tcr-explorer .
Package: https://pypi.org/project/tcr-explorer/ .
DOI: https://doi.org/10.5281/zenodo.21204936 .

The bundled IMGT germline reference is used under CC BY 4.0 (see the package
attribution).
