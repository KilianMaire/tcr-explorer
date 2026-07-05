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

# Bake the pre-built IEDB-only records index. It is built once with
# deploy/build_iedb_index.py and copied into the Space repo next to this
# Dockerfile (see deploy/DEPLOY.md), so the build needs no network to iedb.org
# (which times out from datacenter build runners) and cold starts are fast.
# /index is used instead of /data so a Space persistent-storage mount cannot
# shadow it.
RUN mkdir -p /index
COPY records_index.parquet /index/records_index.parquet
COPY records_index.meta.json /index/records_index.meta.json

# Hugging Face Spaces routes to the port declared as app_port in README.md.
EXPOSE 7860

CMD ["uvicorn", "tcr_explorer.api:app", "--host", "0.0.0.0", "--port", "7860"]
