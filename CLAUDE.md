# TCR Explorer

Distributed API for searching IMGT (HLA, TCR), VDJdb, IEDB, and IPD-MHC data.

## Quick Start

```bash
pip install -r requirements.txt
pip install -r requirements-dev.txt
PYTHONPATH=src uvicorn tcr_explorer.api:app --port 8000 --reload
```

## Architecture

- `src/tcr_explorer/api.py` — Main FastAPI gateway (port 8000)
- `servers/hla_server.py` — EBI IMGT/HLA proxy (port 8101)
- `servers/tcr_server.py` — NCBI Entrez TCR proxy (port 8102)
- `servers/vdjdb_server.py` — VDJdb integration (port 8103)
- `servers/iedb_server.py` — IEDB integration (port 8104)
- `servers/mhc_server.py` — IPD-MHC proxy (port 8105)

## Running Tests

```bash
make test
# or
PYTHONPATH=src:. python -m pytest tests/ -v --tb=short
```

## Docker

```bash
docker-compose up -d
```

## Conventions

- All tool servers are independent FastAPI apps (one per data source)
- Main API routes requests to tool servers via httpx
- NL query heuristic in `nl_query.py` detects source from free text
- Tests use `conftest.py` fixtures for TestClient
