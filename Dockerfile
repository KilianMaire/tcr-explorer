# =============================================================================
# TCRpredictor — Multi-stage Dockerfile
# =============================================================================
# Stage 1 (base):  System deps + production Python packages
# Stage 2 (dev):   Add dev/test dependencies on top of base
# Stage 3 (prod):  Slim runtime with only production files
# =============================================================================

# ---------------------------------------------------------------------------
# Stage 1: base — system dependencies + production Python packages
# ---------------------------------------------------------------------------
FROM python:3.11-slim AS base

# Prevent Python from buffering stdout/stderr and writing .pyc files
ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

# Install system dependencies required by scientific/bioinformatics packages
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        curl \
        git \
        libhdf5-dev \
        libffi-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Install production Python dependencies (layer cached unless requirements change)
COPY requirements.txt ./
RUN pip install --no-cache-dir -r requirements.txt

# ---------------------------------------------------------------------------
# Stage 2: dev — development and test dependencies
# ---------------------------------------------------------------------------
FROM base AS dev

COPY requirements-dev.txt requirements-ml.txt ./
RUN pip install --no-cache-dir -r requirements-dev.txt \
    && pip install --no-cache-dir -r requirements-ml.txt

# Copy entire project for development/testing
COPY . .

ENV PYTHONPATH=src:.

# Default command for dev: run tests
CMD ["pytest", "tests/", "-v", "--tb=short"]

# ---------------------------------------------------------------------------
# Stage 3: prod — slim runtime image
# ---------------------------------------------------------------------------
FROM python:3.11-slim AS prod

# Prevent Python from buffering stdout/stderr and writing .pyc files
ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

# Install minimal runtime system dependencies (curl for healthcheck)
RUN apt-get update && apt-get install -y --no-install-recommends \
        curl \
        libhdf5-103-1 \
    && rm -rf /var/lib/apt/lists/*

# Create non-root user for security
RUN groupadd --gid 1000 appuser \
    && useradd --uid 1000 --gid appuser --shell /bin/bash --create-home appuser

WORKDIR /app

# Copy installed Python packages from base stage
COPY --from=base /usr/local/lib/python3.11/site-packages /usr/local/lib/python3.11/site-packages
COPY --from=base /usr/local/bin /usr/local/bin

# Copy application source code (only what is needed at runtime)
COPY src/ ./src/
COPY servers/ ./servers/
COPY pipeline/ ./pipeline/
COPY models/ ./models/
COPY structural/ ./structural/
COPY mlops/ ./mlops/
COPY examples/ ./examples/
COPY requirements.txt ./

# Set PYTHONPATH so imports work correctly
ENV PYTHONPATH=src:.

# Create directories that may be needed at runtime
RUN mkdir -p /app/batman_cache /app/data/models/tier1 \
    && chown -R appuser:appuser /app

# Switch to non-root user
USER appuser

# Expose common service ports
EXPOSE 8000 8101 8102 8103 8104 8105 8106 8110 8120 8501

# Health check — default assumes main API gateway on port 8000
HEALTHCHECK --interval=30s --timeout=5s --start-period=30s --retries=3 \
    CMD curl -f http://localhost:8000/health || exit 1

# Default command: run the main API gateway
CMD ["uvicorn", "tcr_explorer.api:app", "--host", "0.0.0.0", "--port", "8000"]
