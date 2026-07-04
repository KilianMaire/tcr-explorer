# Bundled IMGT germline: attribution and license

This directory holds the IMGT germline (TCR V, D, J, C region FASTA and region
motif tables) for human and mouse, bundled with TCR Explorer so germline
features work offline.

## Data: IMGT, CC BY 4.0

The germline sequences and annotations are from IMGT, the international
ImMunoGeneTics information system (IMGT/GENE-DB).

- Release: 20268-7 (produced 2026-02-25, per each species' `data-production-date.tsv`).
- License: Creative Commons Attribution 4.0 International (CC BY 4.0), per
  <https://www.imgt.org/about/termsofuse.php>. Redistribution, including
  commercial, is permitted with attribution.
- Please cite: Lefranc M-P. et al. IMGT, the international ImMunoGeneTics
  information system. <https://www.imgt.org>

## Reformatting: stitchr / IMGTgeneDL (MIT)

The FASTA layout was produced by IMGTgeneDL and stitchr (Jamie Heather), both
under the MIT license. The underlying sequences are unchanged IMGT data; only
the file organisation is stitchr's. This is the "indicate if modified" note that
CC BY 4.0 asks for.

- stitchr: <https://github.com/JamieHeather/stitchr>

## Refreshing

The bundled copy updates with TCR Explorer releases. To pull a newer IMGT
germline yourself, run `tcr-explorer-refresh --germline` (needs IMGT/GENE-DB
reachable); it writes into the user data dir, which the tool prefers over this
bundled copy.
