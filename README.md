# TCR Explorer

A federated tool for T cell receptor analysis. It retrieves known TCR records (VDJdb, IEDB, McPAS, TCR3d), assigns germline V and J genes down to the allele level, reconstructs full membrane bound chains, builds per receptor dossiers, and finds similar receptors. The same pure functions back a web UI, a REST API, and an MCP server, so an assistant can drive the whole tool.

## How the data works

The package ships **no datasets**. On first use you run `tcr-explorer-refresh` once. It downloads every dataset (VDJdb, IEDB, McPAS, TCR3d, plus the IMGT germline) from each source's own official endpoint into a local folder, then harmonizes them into a single records index. After that, everything runs in one process against that local index, offline, until you refresh again to pull fresh data.

This design keeps the data current and means the tool never redistributes third party datasets. Each user fetches them directly from the source under that source's own terms. See [Data sources](#data-sources).

## Requirements

- Python 3.11 or newer.
- Internet access for the initial `tcr-explorer-refresh` (a few minutes, roughly 60 MB). Offline afterward until you refresh.

## Install

```bash
pip install tcr-explorer
```

Or from a checkout:

```bash
python -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
pip install -e .
```

## First run

```bash
tcr-explorer-refresh
```

This downloads the datasets, builds the index into a local data folder (a platform specific user data directory, or wherever `TCR_EXPLORER_DATA` points), and fetches the IMGT germline. Re-run it any time to update. If a tool is used before the first refresh, it returns a clear message asking you to run this command.

## Use it

Two front doors, both a single process.

### As an MCP server (recommended)

Point your own assistant at TCR Explorer over MCP and ask questions in plain language. This targets Claude Desktop and Claude Code, which run a local stdio MCP server. ChatGPT does not run local stdio MCP servers the same way, so use Claude for the paste and go flow. See [Connect your assistant](#connect-your-assistant).

### As a web app and REST API

```bash
PYTHONPATH=src uvicorn tcr_explorer.api:app --port 8000
```

Open the query box at <http://localhost:8000/ui>, or call the REST API directly. Health check:

```bash
curl http://localhost:8000/health   # {"status":"ok"}
```

## Connect your assistant

TCR Explorer ships an MCP server (console entry point `tcr-explorer-mcp`). Add this to your assistant's MCP configuration:

```json
{"mcpServers":{"tcr-explorer":{"command":"uvx","args":["--from","tcr-explorer","tcr-explorer-mcp"]}}}
```

Or paste this prompt into Claude to have it set the connection up for you:

```
Set up the TCR Explorer MCP server so you can answer T cell receptor questions against real immunology databases. First install it (pip install tcr-explorer, or use uvx --from tcr-explorer), then run tcr-explorer-refresh once in a terminal to download the datasets into a local folder (a few minutes). Then add an MCP server named tcr-explorer that runs `uvx --from tcr-explorer tcr-explorer-mcp` (if uvx is unavailable, run python -m tcr_explorer.mcp_server). It exposes these read only tools: retrieve_tcr_records, assign_tcr_alleles, get_tcr_dossier, find_similar_tcrs, align_tcr_genes, and ask_tcr. If a tool reports the data is not downloaded yet, tell me to run tcr-explorer-refresh. After adding it, confirm the connection and suggest three example questions I can ask.
```

Until the package is on PyPI, the git form works: `uvx --from git+<your-repo-url> tcr-explorer-mcp`.

The read only MCP tools are `retrieve_tcr_records`, `assign_tcr_alleles`, `get_tcr_dossier`, `find_similar_tcrs`, `align_tcr_genes`, and `ask_tcr`.

## REST API

All of these run in process against the local index.

### Unified query box

**POST** `/v1/tcr/query` routes a single input (a CDR3, a full chain, a gene name, a record id, or a phrase) to the right tool.

```bash
curl -s -X POST http://localhost:8000/v1/tcr/query \
  -H "Content-Type: application/json" \
  -d '{"query":"CASSLGGAGGTDTQYF","species":"human"}'
```

### Germline assignment

**POST** `/v1/tcr/assign` assigns a TCR sequence (nucleotide or amino acid, CDR3, region, or full chain) to V and J alleles, with per region identity, co optimal ties, CDR3 extraction, and an honest refusal to call a V allele from a bare CDR3.

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

**POST** `/reconstruct` builds a full membrane bound chain. Provide V, J, and CDR3, or a CDR3 alone (V and J are inferred from the records that carry the same CDR3).

```bash
curl -s -X POST http://localhost:8000/reconstruct \
  -H "Content-Type: application/json" \
  -d '{"cdr3_aa":"CASSLGGAGGTDTQYF","species":"human"}'
```

### CDR1 and CDR2 prediction

**GET** `/predict/cdr` returns germline CDR1 and CDR2 for a TCR V gene from IMGT germline data.

```bash
curl "http://localhost:8000/predict/cdr?v_gene=TRBV12-3&species=human"
```

## Optional: MHC allele sequences

The records index does not contain MHC allele sequences. If you want live lookup of those from EBI IMGT/HLA and IPD-MHC, start the optional hla and mhc proxies with one command:

```bash
docker-compose up
```

The `/search` endpoint is scoped to these two sources: `{"source": "hla"}` or `{"source": "mhc"}`. Any other source returns HTTP 400 pointing at `/v1/tcr/records`, which is where TCR record search lives.

## Environment variables

All optional.

| Variable | Default | Description |
|----------|---------|-------------|
| `TCR_EXPLORER_DATA` | platform user data dir | Local folder where `tcr-explorer-refresh` downloads datasets and builds the index |
| `TCR_EXPLORER_MAX_AGE_DAYS` | `30` | Age after which the local index is flagged stale (a refresh is suggested in query warnings, never forced) |
| `RECORDS_INDEX_PATH` | `<data dir>/records_index.parquet` | Override the records index path directly |
| `HLA_SERVER_URL` | `http://127.0.0.1:8101` | HLA allele sequence proxy (optional) |
| `MHC_SERVER_URL` | `http://127.0.0.1:8105` | IPD-MHC allele sequence proxy (optional) |
| `LLM_BASE_URL` | *(empty)* | OpenAI compatible endpoint for the free text `ask` path (falls back to a heuristic parser when unset) |
| `LLM_MODEL` | `local-model` | Model id for the `ask` path |

## Data sources

TCR Explorer downloads and cites the following. It does not redistribute them; each is fetched on your machine from its official endpoint. Please cite the ones you use.

- **VDJdb**. Goncharov M. et al. VDJdb in the pandemic era: a compendium of T cell receptors specific for SARS-CoV-2. Nature Methods, 2022. <https://github.com/antigenomics/vdjdb-db>
- **IEDB** (CC BY 4.0). Vita R. et al. The Immune Epitope Database (IEDB): 2024 update. Nucleic Acids Research, 2025. <https://www.iedb.org>
- **McPAS-TCR**. Tickotsky N. et al. McPAS-TCR: a manually curated catalogue of pathology associated T cell receptor sequences. Bioinformatics, 2017. <https://friedmanlab.weizmann.ac.il/McPAS-TCR/>
- **TCR3d**. Lin V. et al. TCR3d 2.0: expanding the T cell receptor structure database. Nucleic Acids Research, 2025. <https://tcr3d.ibbr.umd.edu>
- **IMGT germline** (CC BY 4.0). Lefranc M-P. et al. IMGT, the international ImMunoGeneTics information system. Fetched via stitchr. <https://www.imgt.org>

## Run tests

```bash
PYTHONPATH=src pytest tests/ -v
```

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
        local index (built by tcr-explorer-refresh
        from the downloaded datasets + IMGT germline)
```

IMGT (IMGT/HLA, IMGT/GENE-DB, IMGT germline, IMGT numbering) is a data source cited throughout. TCR Explorer is an independent tool and is not affiliated with IMGT.
