# TCR Explorer, public germline-only demo (Hugging Face Spaces, Docker SDK).
#
# This image serves the web UI and REST API with germline features only:
# allele assignment, chain reconstruction, alignment, and CDR loops. It never
# runs `tcr-explorer-refresh`, so no record database (VDJdb, IEDB, McPAS, TCR3d)
# is downloaded or served, which keeps the public instance within the data
# licenses. The bundled IMGT germline (CC BY 4.0) ships inside the package.
#
# The package is installed straight from a pinned public GitHub commit tarball
# (no git binary needed), so the build is reproducible. To update the demo, bump
# the commit SHA below to a newer tcr-explorer commit.

FROM python:3.11-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    TCR_EXPLORER_DEMO=1 \
    TCR_EXPLORER_DATA=/tmp/tcr-explorer-data

RUN pip install "https://github.com/KilianMaire/tcr-explorer/archive/c1b8b881bd43c4f71d215c5ecf8bee49af8ae6fe.tar.gz"

# Hugging Face Spaces routes to the port declared as app_port in README.md.
EXPOSE 7860

CMD ["uvicorn", "tcr_explorer.api:app", "--host", "0.0.0.0", "--port", "7860"]
