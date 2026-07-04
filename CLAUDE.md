# TCR Explorer

Federated TCR analysis tool: records retrieval (VDJdb, IEDB, McPAS, TCR3d), germline V/J allele assignment, chain reconstruction, dossiers, similarity. Single source pure functions back a web UI, a REST API, and an MCP server.

## Quick Start (single process, no servers)

The core runs in one process from vendored data (`data/records_index.parquet` plus IMGT germline). No external services, no keys, no internet needed.

```bash
pip install -r requirements.txt
pip install -r requirements-dev.txt
# Web UI + REST API:
PYTHONPATH=src uvicorn tcr_explorer.api:app --port 8000 --reload   # then open /ui
# Or the MCP server (recommended front door):
python -m tcr_explorer.mcp_server
```

## Architecture

The core feature modules run entirely in process from vendored data. They do NOT open any socket:

- `src/tcr_explorer/api.py` — FastAPI app (web `/ui` + REST), port 8000
- `src/tcr_explorer/mcp_server.py` — MCP server (same pure functions)
- `records.py`, `tcr_align.py`, `germline_db.py`, `reconstructor.py`, `query_router.py`, `dossier.py` — the pure core, reading from `data/` parquet + germline

Optional live enrichment (reached over HTTP only when running, degrades gracefully when absent):

- `servers/hla_server.py` — EBI IMGT/HLA allele sequence proxy (port 8101)
- `servers/mhc_server.py` — IPD-MHC allele sequence proxy (port 8105)
- `servers/batman_server.py`, `tempo_server.py`, `structural_server.py` — ML scoring / structure backends (audit pending)

Never hand start these in separate terminals. Launch the whole optional stack with one command: `docker-compose up`.

## Running Tests

```bash
make test
# or
PYTHONPATH=src:. python -m pytest tests/ -v --tb=short
```

## Optional live stack (one command, not five terminals)

```bash
docker-compose up
```

## Conventions

- Core = single source pure functions backing REST + MCP + UI, reading vendored `data/` (no live dependency)
- The only record path is the vendored parquet (`retrieve_records`, `/v1/tcr/records`, MCP, `/ui`); there is no SQLite index or ingest anymore
- `/search` is scoped to hla/mhc allele sequences (400 for any other source); the hla/mhc servers under `servers/` are OPTIONAL live enrichment reached via httpx
- `query_router.route_query` deterministically routes any `/ui` input to records/assign/dossier/similar/ask (no LLM required)
- Tests use `conftest.py` fixtures for TestClient
