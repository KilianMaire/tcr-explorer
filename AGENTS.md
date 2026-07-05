# Working with TCR Explorer (instructions for AI assistants)

This file tells an AI coding or research assistant (Codex, Cursor, Cline, and similar MCP-capable tools) how to drive TCR Explorer to answer a user's T cell receptor questions against real immunology data. Read it before you touch the tool.

TCR Explorer retrieves known TCR records (VDJdb, IEDB, McPAS, TCR3d), assigns germline V and J genes to the allele level, reconstructs full membrane bound chains, builds per receptor dossiers, and finds similar receptors. The same pure functions back a REST API, a web query box, and an MCP server, so you can drive everything through one interface.

## First, make the data ready

The package bundles the IMGT germline but ships no record datasets. Before any tool can return results, the user runs one command that downloads the four record sources (VDJdb, IEDB, McPAS, TCR3d) into a local folder and builds a single harmonized index. The germline is already present, so this step does not touch IMGT.

```bash
pip install "tcr-explorer[tcrdist]"   # recommended; drop [tcrdist] for a leaner install
tcr-explorer-refresh                  # first run: download plus build (a few minutes, roughly 60 MB)
```

The `tcrdist` extra gives authoritative similarity scoring. Either way the tool does both automatically: it uses tcrdist when the extra is present and falls back to the bundled BLOSUM CDR3 distance otherwise, reporting which one it used in each result's `engine` field. Paired similarity needs the extra (it is tcrdist only).

If a tool reports that the data is not downloaded yet (a `*_index_unavailable` warning, or a `DataNotReadyError`), do not try to work around it. Tell the user to run `tcr-explorer-refresh` once, then retry. A query never downloads on its own. This is deliberate.

If a tool returns a `records_index_stale` warning, the local index is older than the staleness window (30 days by default). Results are still valid. Suggest the user re-run `tcr-explorer-refresh` to pull fresh data, but do not force it.

## Connect over MCP (recommended)

For a desktop or CLI MCP client, add an MCP server that runs the console entry point:

```json
{"mcpServers":{"tcr-explorer":{"command":"uvx","args":["--from","tcr-explorer[tcrdist]","tcr-explorer-mcp"]}}}
```

If `uvx` is unavailable, run `python -m tcr_explorer.mcp_server`. You can also run it straight from the public repo: `uvx --from "git+https://github.com/KilianMaire/tcr-explorer[tcrdist]" tcr-explorer-mcp`.

Some web only assistants do not run local stdio MCP servers the same way. For a paste and go flow, use a desktop or CLI MCP client, or call the REST API directly (see below).

## The tools you have

All are read only and run in one process against the local index.

| Tool | Use it for |
|------|-----------|
| `retrieve_tcr_records` | Look up known receptors by CDR3, gene, epitope, or record id in the harmonized index |
| `assign_tcr_alleles` | Assign a nucleotide or amino acid sequence (CDR3, region, or full chain) to V and J alleles, with per region identity and co optimal ties |
| `get_tcr_dossier` | Build a full per receptor dossier (assignment, records, neighbours, context) in one call |
| `find_similar_tcrs` | Find receptors near a single chain query (alpha or beta) by tcrdist when the tcrdist extra is installed, otherwise BLOSUM CDR3 distance, filtered by V family and species |
| `find_similar_paired_tcrs` | Find alpha/beta paired receptors near a paired query by paired tcrdist (alpha plus beta). Needs the tcrdist extra and both query V genes in the reference table |
| `align_tcr_genes` | Align a sequence against germline genes |
| `ask_tcr` | Route a single free text or sequence input to the right tool when you are unsure which to call |

When you do not know which tool fits, pass the raw user input to `ask_tcr`. It routes a CDR3, a full chain, a gene name, a record id, or a phrase to the correct tool deterministically, no LLM required.

## How to answer a TCR question

1. Confirm the data is ready. If a tool returns an unavailable warning, stop and ask the user to run `tcr-explorer-refresh`.
2. Identify what the user gave you: a bare CDR3, a full chain, a V or J gene, an epitope, or a record id. If in doubt, hand it to `ask_tcr`.
3. Call the narrowest tool that answers the question. Use `get_tcr_dossier` when the user wants a full picture of one receptor.
4. Read the `warnings` field on every result and pass its meaning to the user. These warnings are honest limits, not noise.

## Honesty rules (do not break these)

- The tool refuses to call a V allele from a bare CDR3, because the CDR3 alone does not determine the V gene. If the user wants a V call, ask for the full chain or region. Do not invent one.
- Similarity scores are within query relative: only `similarity == 1.0` (identical) is absolute. For cross query comparisons, threshold on the raw `distance` field, not on similarity.
- The similarity result reports which engine scored it in the `engine` field. When `engine` is `tcrdist`, the distances are the authoritative tcrdist metric (installed via the optional `tcrdist` extra) and are absolute, so cross query comparisons are valid on the `distance` field. When `engine` is `blosum_cdr3`, the tool used the bundled BLOSUM CDR3 distance and says so in a `tcrdist_unavailable` warning (suggest `pip install tcr-explorer[tcrdist]` for authoritative scoring). Never present a BLOSUM distance as tcrdist. If a `tcrdist_candidates_skipped` warning appears, some reference candidates had a V gene absent from the tcrdist table and were skipped (their germline loops are unknown), so they are not in the results.
- The records index does not contain MHC allele sequences. Those come from the optional hla and mhc proxies, reached only when the user starts them (`docker-compose up`). Do not claim MHC sequence data from the records index.
- Report what the data says, including empty results. If there are no records for a CDR3, say so. Do not fill the gap with a guess.

## REST API (alternative to MCP)

Start the web app and API in one process:

```bash
PYTHONPATH=src uvicorn tcr_explorer.api:app --port 8000
```

Then either open the query box at `http://localhost:8000/ui`, or call the endpoints. The unified box routes any single input:

```bash
curl -s -X POST http://localhost:8000/v1/tcr/query \
  -H "Content-Type: application/json" \
  -d '{"query":"CASSLGGAGGTDTQYF","species":"human"}'
```

Named endpoints: `/v1/tcr/assign`, `/v1/tcr/records`, `/reconstruct`, and `/predict/cdr`. See `README.md` for request shapes.

## Data sources to cite

When you report results, cite the sources the user relied on. Full references are in `README.md`. In short: VDJdb (Goncharov 2022), IEDB (Vita 2025, CC BY 4.0), McPAS-TCR (Tickotsky 2017), TCR3d (Lin 2025), IMGT germline (Lefranc, CC BY 4.0). The four record datasets are downloaded on the user's machine and not redistributed; the IMGT germline is bundled with the package under CC BY 4.0 (release 20268-7, attribution in `src/tcr_explorer/data/germline/ATTRIBUTION.md`).
