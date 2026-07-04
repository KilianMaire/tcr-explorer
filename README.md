# TCR Explorer

A federated tool for T cell receptor sequence analysis. It retrieves known TCR records (VDJdb, IEDB, McPAS, TCR3d), assigns germline V and J genes down to the allele level, reconstructs full membrane bound chains, builds per receptor dossiers, and finds similar receptors. The same pure functions back a web UI, a REST API, and an MCP server.

---

## What you get

Everything below runs **in a single process from vendored data** (a harmonized records index plus IMGT germline). No external services, no API keys, no internet required for the core.

- **Records retrieval** across VDJdb, IEDB, McPAS, and TCR3d.
- **Germline assignment** of any TCR sequence (nucleotide or amino acid, CDR3, region, or full chain) to V and J alleles, with per region identity, co optimal ties reported, CDR3 extraction, and an honest refusal to call a V allele from a bare CDR3.
- **Chain reconstruction** from V, J, and CDR3 (or from a CDR3 alone by inferring V and J from the database).
- **Dossiers and similarity** for a receptor.
- **A unified query box** at `/ui`: one central box routes any input to the right tool.

Optional live enrichment (fresh EBI IMGT/HLA, NCBI Entrez, live VDJdb and IEDB proxying, plus ML scoring backends) runs as a separate service group and is never required to use the tool. See [Optional: live enrichment stack](#optional-live-enrichment-stack).

---

## Requirements

- Python 3.11+

No internet access is needed for the core features. The vendored data ships with the package.

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

TCR Explorer ships an MCP server (`tcr_explorer.mcp_server:main`, console entry point `tcr-explorer-mcp`) so an assistant can query the databases and run the analyses directly.

Add this to your assistant's MCP configuration:

```json
{"mcpServers":{"tcr-explorer":{"command":"uvx","args":["--from","tcr-explorer","tcr-explorer-mcp"]}}}
```

Or paste this prompt into your assistant to have it set the connection up for you:

```
Set up the TCR Explorer MCP server so you can answer T cell receptor questions against real immunology databases. Add an MCP server named tcr-explorer that runs `uvx --from tcr-explorer tcr-explorer-mcp` (if uvx is unavailable, `pip install tcr-explorer` then run `python -m tcr_explorer.mcp_server`). It exposes these read only tools: retrieve_tcr_records, assign_tcr_alleles, get_tcr_dossier, find_similar_tcrs, align_tcr_genes, and ask_tcr. After adding it, confirm the connection and suggest three example questions I can ask.
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

The core reads from a vendored, harmonized snapshot. If you also want **live** proxying to EBI IMGT/HLA and NCBI Entrez, live VDJdb and IEDB queries, or the ML scoring backends, start the enrichment services. Do this with **one command**, never five terminals:

```bash
docker-compose up
```

This launches the tool servers and the main API together on a shared network. When these services are absent, the endpoints that would use them degrade gracefully to the vendored data rather than failing.

The legacy federated `/search` endpoint (`{"source":"hla|tcr|vdjdb|iedb", ...}`) is served by this stack. Whether the live federation still earns its place, given the vendored snapshot, is an open design question and is not needed for day to day use.

---

## Environment variables

All are optional. Defaults give you the single process core.

| Variable | Default | Description |
|----------|---------|-------------|
| `RECORDS_INDEX_PATH` | `data/records_index.parquet` | Vendored harmonized records index |
| `HLA_SERVER_URL` | `http://127.0.0.1:8101` | HLA tool server (optional live stack) |
| `TCR_SERVER_URL` | `http://127.0.0.1:8102` | TCR tool server (optional live stack) |
| `VDJDB_SERVER_URL` | `http://127.0.0.1:8103` | VDJdb tool server (optional live stack) |
| `IEDB_SERVER_URL` | `http://127.0.0.1:8104` | IEDB tool server (optional live stack) |
| `NCBI_API_KEY` | *(empty)* | NCBI API key, raises the rate limit to 10 req/s |
| `LLM_BASE_URL` | *(empty)* | OpenAI compatible endpoint for the free text `ask` path (falls back to a heuristic parser when unset) |
| `LLM_MODEL` | `local-model` | Model id for the `ask` path |

---

## Run tests

```bash
PYTHONPATH=src pytest tests/ -v
```

---

## Architecture

The core is a set of single source pure functions (records retrieval, germline assignment, reconstruction, dossiers, similarity) that read from vendored data. The REST API, the MCP server, and the web query box all call these same functions in one process.

```
        MCP server        REST API + /ui query box
             \                    /
              \                  /
        single source pure functions
        (records, assign, reconstruct,
         dossier, similar)
                   |
        vendored data (records index parquet,
        IMGT germline)
```

The optional live enrichment services (EBI IMGT/HLA, NCBI Entrez, live VDJdb and IEDB, ML scoring) sit beside this core and are reached over HTTP only when running. They are launched together with `docker-compose up` and are never required for the core.

IMGT (IMGT/HLA, IMGT/GENE-DB, IMGT germline, IMGT numbering) is a data source cited throughout. TCR Explorer is an independent tool and is not affiliated with IMGT.
