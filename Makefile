.PHONY: help install install-dev test test-cov lint format serve serve-all clean

# Default target
.DEFAULT_GOAL := help

# Show help
help: ## Show this help message
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | \
		awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-20s\033[0m %s\n", $$1, $$2}'

# ──────────────────────────── Setup ────────────────────────────

install: ## Install production dependencies
	pip install -r requirements.txt

install-dev: ## Install dev + production dependencies
	pip install -r requirements-dev.txt

# ──────────────────────────── Testing ──────────────────────────

test: ## Run test suite
	PYTHONPATH=src:. python -m pytest tests/ \
		--ignore=tests/test_cdr_enricher.py \
		--ignore=tests/test_nl_query.py \
		-v --tb=short

test-cov: ## Run tests with coverage report
	PYTHONPATH=src:. python -m pytest tests/ \
		--ignore=tests/test_cdr_enricher.py \
		--ignore=tests/test_nl_query.py \
		--cov=src/ --cov=servers/ --cov-report=html --cov-report=term-missing

# ──────────────────────────── Code quality ─────────────────────

lint: ## Run linters (ruff)
	ruff check src/ tests/ servers/

format: ## Auto-format code (ruff)
	ruff format src/ tests/ servers/

# ──────────────────────────── Serving ──────────────────────────

serve: ## Start main API server (port 8000)
	PYTHONPATH=src:. uvicorn tcr_explorer.api:app --host 0.0.0.0 --port 8000 --reload

serve-all: ## Start all services with docker-compose
	docker-compose up -d

# ──────────────────────────── Cleanup ──────────────────────────

clean: ## Remove generated files and caches
	rm -rf __pycache__ .pytest_cache htmlcov .coverage
	find . -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null || true
	find . -type f -name "*.pyc" -delete 2>/dev/null || true
