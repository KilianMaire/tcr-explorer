# TCR Explorer

A federated tool for T cell receptor sequence analysis. It retrieves known TCR records (VDJdb, IEDB, McPAS, TCR3d), assigns germline V and J genes down to the allele level, reconstructs full membrane bound chains, builds per receptor dossiers, and finds similar receptors. The same pure functions back a web UI, a REST API, and an MCP server.

---

## What you get

On first use, `tcr-explorer-refresh` downloads every dataset (VDJdb, IEDB, McPAS, TCR3d, plus the IMGT germline) from each source's official endpoint into a local folder and builds a harmonized records index. After that, everything below runs **in a single process from that local index**. No API keys, and no internet after the initial download until you refresh.

- **Records retrieval** across VDJdb, IEDB, McPAS, and TCR3d.
- **Germline assignment** of any TCR sequence (nucleotide or amino acid, CDR3, region, or full chain) to V and J alleles, with per region identity, co optimal ties reported, CDR3 extraction, and an honest refusal to call a V allele from a bare CDR3.
- **Chain reconstruction** from V, J, and CDR3 (or from a CDR3 alone by inferring V and J from the database).
- **Dossiers and similarity** for a receptor.
- **A unified query box** at `/ui`: one central box routes any input to the right tool.

Optional live enrichment for MHC allele sequences (EBI IMGT/HLA and IPD-MHC), which the snapshot does not contain, runs as a separate service group and is never required to use the tool. See [Optional: live enrichment stack](#optional-live-enrichment-stack).

---

## Requirements

- Python 3.11+

Internet access is needed once, for the initial `tcr-explorer-refresh` download. After that the core runs offline against the local index until you refresh.

---

## Install

```bash
python -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

---

## Run it

Pick one of two front doors. Both run as a **single process**.

### As an MCP server (recommended)

Point your own assistant at TCR Explorer over MCP and ask questions in plain language. See [Connect your assistant](#connect-your-assistant). Nothing else to start.

### As a web app and REST API

```bash
PYTHONPATH=src uvicorn tcr_explorer.api:app --port 8000
```

Then open the query box at <http://localhost:8000/ui>, or call the REST API directly.

```bash
curl http://localhost:8000/health
# {"status":"ok"}
```

That is the whole tool. The five separate servers described in older revisions are not required.

---

## Connect your assistant

TCR Explorer ships an MCP server (`tcr_explorer.mcp_server:main`, console entry point `tcr-explorer-mcp`) so an assistant can query the databases and run the analyses directly. On first use, run `tcr-explorer-refresh` once. It downloads every dataset (VDJdb, IEDB, McPAS, TCR3d, plus the IMGT germline) from each source's official endpoint into a local folder and builds the index (a few minutes). The package itself ships no data; each user fetches it directly from the source. Re-run `tcr-explorer-refresh` any time to update.

This targets Claude Desktop and Claude Code, which run a local stdio MCP server from the config below. ChatGPT does not run local stdio MCP servers the same way, so use Claude for the paste and go flow.

Add this to your assistant's MCP configuration:

```json
{"mcpServers":{"tcr-explorer":{"command":"uvx","args":["--from","tcr-explorer","tcr-explorer-mcp"]}}}
```

Or paste this prompt into your assistant to have it set the connection up for you:

```
Set up the TCR Explorer MCP server so you can answer T cell receptor questions against real immunology databases. First install it (`pip install tcr-explorer`, or use `uvx --from tcr-explorer`), then run `tcr-explorer-refresh` once in a terminal to download the datasets into a local folder (a few minutes). Then add an MCP server named tcr-explorer that runs `uvx --from tcr-explorer tcr-explorer-mcp` (if uvx is unavailable, run `python -m tcr_explorer.mcp_server`). It exposes these read only tools: retrieve_tcr_records, assign_tcr_alleles, get_tcr_dossier, find_similar_tcrs, align_tcr_genes, and ask_tcr. If a tool reports the data is not downloaded yet, tell me to run `tcr-explorer-refresh`. After adding it, confirm the connection and suggest three example questions I can ask.
```

Until the package is published to PyPI, the git form works instead: `uvx --from git+<your-repo-url> tcr-explorer-mcp`.

The exposed read only tools are `retrieve_tcr_records`, `assign_tcr_alleles`, `get_tcr_dossier`, `find_similar_tcrs`, `align_tcr_genes`, and `ask_tcr`.

---

## REST API

All of these run in process, no external services.

### Unified query box

**POST** `/v1/tcr/query` routes a single input (a CDR3, a full chain, a gene name, a record id, or a phrase) to the right tool.

```bash
curl -s -X POST http://localhost:8000/v1/tcr/query \
  -H "Content-Type: application/json" \
  -d '{"query":"CASSLGGAGGTDTQYF","species":"human"}'
```

### Germline assignment

**POST** `/v1/tcr/assign` assigns a TCR sequence to V and J alleles.

```bash
curl -s -X POST http://localhost:8000/v1/tcr/assign \
  -H "Content-Type: application/json" \
  -d '{"sequence":"CASSLGGAGGTDTQYF","species":"human"}'
```

### Records retrieval

**POST** `/v1/tcr/records` searches the harmonized records index.

```bash
curl -s -X POST http://localhost:8000/v1/tcr/records \
  -H "Content-Type: application/json" \
  -d '{"cdr3":"CASSLGGAGGTDTQYF","species":"human","limit":20}'
```

### Chain reconstruction

**POST** `/reconstruct` builds a full membrane bound chain. Provide V, J, and CDR3, or a CDR3 alone (V and J are inferred from the database).

```bash
curl -s -X POST http://localhost:8000/reconstruct \
  -H "Content-Type: application/json" \
  -d '{"cdr3_aa":"CASSLGGAGGTDTQYF","species":"human"}'
```

### CDR1 / CDR2 prediction

**GET** `/predict/cdr` returns germline CDR1 and CDR2 for a TCR V gene from IMGT germline data.

```bash
curl "http://localhost:8000/predict/cdr?v_gene=TRBV12-3&species=human"
```

---

## Optional: live enrichment stack

The core reads from the locally downloaded, harmonized index and needs nothing else. The one thing the snapshot does not contain is **MHC allele sequences**. If you want live lookup of those from EBI IMGT/HLA and IPD-MHC, start the optional enrichment stack with **one command**, never five terminals:

```bash
docker-compose up
```

This launches the hla tool server and the main API together on a shared network. When the hla/mhc servers are absent, the `/search` endpoint still responds; it just cannot serve allele sequences.

TCR **record** search is served in process by `POST /v1/tcr/records` from the local index. The `/search` endpoint is now scoped to MHC allele sequences only: `{"source": "hla"}` or `{"source": "mhc"}`. Any other source returns HTTP 400 pointing at `/v1/tcr/records`.

---

## Environment variables

All are optional. Defaults give you the single process core.

| Variable | Default | Description |
|----------|---------|-------------|
| `TCR_EXPLORER_DATA` | platformdirs user data dir | Local folder where `tcr-explorer-refresh` downloads datasets and builds the index |
| `TCR_EXPLORER_MAX_AGE_DAYS` | `30` | Age after which the local index is flagged stale (a refresh is suggested, never forced) |
| `RECORDS_INDEX_PATH` | `<data dir>/records_index.parquet` | Override the records index path directly |
| `HLA_SERVER_URL` | `http://127.0.0.1:8101` | HLA allele sequence proxy (optional live stack) |
| `MHC_SERVER_URL` | `http://127.0.0.1:8105` | IPD-MHC allele sequence proxy (optional live stack) |
| `LLM_BASE_URL` | *(empty)* | OpenAI compatible endpoint for the free text `ask` path (falls back to a heuristic parser when unset) |
| `LLM_MODEL` | `local-model` | Model id for the `ask` path |

---

## Run tests

```bash
PYTHONPATH=src pytest tests/ -v
```

---

## Architecture

The core is a set of single source pure functions (records retrieval, germline assignment, reconstruction, dossiers, similarity) that read from the locally downloaded index. The REST API, the MCP server, and the web query box all call these same functions in one process.

```
        MCP server        REST API + /ui query box
             \                    /
              \                  /
        single source pure functions
        (records, assign, reconstruct,
         dossier, similar)
                   |
        local index (records index parquet,
        IMGT germline)
```

The optional hla and mhc allele sequence proxies (EBI IMGT/HLA, IPD-MHC) sit beside this core and are reached over HTTP only when running. They are launched with `docker-compose up` and are never required for the core.

IMGT (IMGT/HLA, IMGT/GENE-DB, IMGT germline, IMGT numbering) is a data source cited throughout. TCR Explorer is an independent tool and is not affiliated with IMGT.
