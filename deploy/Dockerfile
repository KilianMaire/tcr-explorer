# TCR Explorer, public IEDB + germline instance (Hugging Face Spaces, Docker SDK).
#
# This image serves three surfaces from the same pure functions:
#   /ui           the human web UI
#   REST + /openapi.json  for chat assistants that call a REST API (GPT Actions)
#   /mcp          the MCP server over streamable HTTP, for remote MCP clients
#
# Data: only IEDB is served (CC BY 4.0, redistributable with attribution). The
# IEDB records index is built once AT IMAGE BUILD TIME and baked into the image,
# so container cold starts are fast and IEDB is never committed to the repo or
# the PyPI package. VDJdb, McPAS and TCR3d are NOT served (their licenses do not
# permit redistribution). The IMGT germline ships inside the package (CC BY 4.0).
#
# The package is installed from a pinned public GitHub commit tarball (no git
# binary needed). To update, bump the commit SHA below and rebuild.

FROM python:3.11-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    TCR_EXPLORER_DATA=/index \
    TCR_EXPLORER_MCP_HTTP=1 \
    TCR_EXPLORER_DEMO=1 \
    TCR_EXPLORER_RATE_LIMIT_PER_MIN=120

# Install the package (pins the tcrdist extra for the authoritative metric).
RUN pip install "tcr-explorer[tcrdist] @ https://github.com/KilianMaire/tcr-explorer/archive/db45c273de5fa53bcff42b86970b0e917c5ee663.tar.gz"

# Build the IEDB-only records index into $TCR_EXPLORER_DATA at build time.
# (A build path avoids IEDB being fetched on every container start; /index is
# used instead of /data so a Space persistent-storage mount cannot shadow it.)
RUN python -c "import os; from pathlib import Path; \
from tcr_explorer.data_sources import download_iedb; \
from tcr_explorer.records_build import build_index; \
raw=Path('/tmp/raw'); raw.mkdir(parents=True, exist_ok=True); \
r=download_iedb(raw); assert r.ok, ('IEDB download failed: '+str(r.error)); \
os.makedirs('/index', exist_ok=True); \
m=build_index(str(raw), '/index/records_index.parquet', '/index/records_index.meta.json'); \
print('IEDB index rows:', m['rows_total'])"

# Hugging Face Spaces routes to the port declared as app_port in README.md.
EXPOSE 7860

CMD ["uvicorn", "tcr_explorer.api:app", "--host", "0.0.0.0", "--port", "7860"]
