"""Tests for .env.example completeness, Dockerfile validity, and .dockerignore presence.

Validates:
- .env.example has no duplicate keys
- All environment variables referenced in Python source exist in .env.example
- Dockerfile exists, is non-empty, and has valid multi-stage build syntax
- .dockerignore exists
"""
from __future__ import annotations

import re
from pathlib import Path

import pytest

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
PROJECT_ROOT = Path(__file__).resolve().parent.parent
ENV_EXAMPLE = PROJECT_ROOT / ".env.example"
DOCKERFILE = PROJECT_ROOT / "Dockerfile"
DOCKERIGNORE = PROJECT_ROOT / ".dockerignore"

# Source directories to scan for env var references
SOURCE_DIRS = [
    PROJECT_ROOT / "src",
    PROJECT_ROOT / "servers",
]

# Pattern to extract env var names from os.getenv("X", ...) or os.environ.get("X", ...)
# Also matches os.environ["X"] for completeness
_ENV_REF_PATTERN = re.compile(
    r"""os\.(?:getenv|environ\.get|environ\[)\s*\(\s*["']([A-Z_][A-Z0-9_]*)["']"""
)

# Env vars that are set internally (e.g. PYTHONPATH in docker-compose) or are
# standard Python/system vars — not expected in .env.example
_SKIP_VARS = {
    "PYTHONPATH",
    "PATH",
    "HOME",
    "USER",
    "LANG",
    "LC_ALL",
    "TZ",
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _parse_env_example() -> list[str]:
    """Return list of keys defined in .env.example (preserving order)."""
    keys: list[str] = []
    for line in ENV_EXAMPLE.read_text().splitlines():
        stripped = line.strip()
        # Skip blank lines and comments
        if not stripped or stripped.startswith("#"):
            continue
        # Keys are of the form: KEY=value or KEY= (with optional inline comment)
        match = re.match(r"^([A-Z_][A-Z0-9_]*)=", stripped)
        if match:
            keys.append(match.group(1))
    return keys


def _collect_env_vars_from_code() -> set[str]:
    """Scan Python source files and return all referenced env var names."""
    env_vars: set[str] = set()
    for src_dir in SOURCE_DIRS:
        if not src_dir.exists():
            continue
        for py_file in src_dir.rglob("*.py"):
            content = py_file.read_text(errors="replace")
            for match in _ENV_REF_PATTERN.finditer(content):
                var_name = match.group(1)
                if var_name not in _SKIP_VARS:
                    env_vars.add(var_name)
    return env_vars


# ---------------------------------------------------------------------------
# Tests: .env.example
# ---------------------------------------------------------------------------


class TestEnvExample:
    """Validate .env.example structure and completeness."""

    def test_env_example_exists(self):
        """The .env.example file must exist at the project root."""
        assert ENV_EXAMPLE.exists(), f"Missing {ENV_EXAMPLE}"

    def test_env_example_not_empty(self):
        """The .env.example file must not be empty."""
        assert ENV_EXAMPLE.stat().st_size > 0, ".env.example is empty"

    def test_no_duplicate_keys(self):
        """Each environment variable key must appear exactly once."""
        keys = _parse_env_example()
        seen: dict[str, int] = {}
        duplicates: list[str] = []
        for key in keys:
            if key in seen:
                duplicates.append(key)
            seen[key] = seen.get(key, 0) + 1

        assert not duplicates, (
            f"Duplicate keys in .env.example: {', '.join(duplicates)}"
        )

    def test_all_code_env_vars_present(self):
        """Every os.getenv / os.environ.get referenced in source code must
        have a corresponding entry in .env.example."""
        env_example_keys = set(_parse_env_example())
        code_vars = _collect_env_vars_from_code()

        missing = code_vars - env_example_keys
        assert not missing, (
            f"Environment variables referenced in code but missing from "
            f".env.example: {', '.join(sorted(missing))}"
        )

    def test_keys_are_uppercase_underscore(self):
        """All keys must follow the UPPER_CASE_WITH_UNDERSCORES convention."""
        keys = _parse_env_example()
        for key in keys:
            assert re.match(r"^[A-Z][A-Z0-9_]*$", key), (
                f"Key '{key}' does not follow UPPER_CASE convention"
            )

    def test_minimum_key_count(self):
        """Sanity check: .env.example should have a reasonable number of keys."""
        keys = _parse_env_example()
        assert len(keys) >= 20, (
            f"Expected at least 20 keys in .env.example, found {len(keys)}"
        )


# ---------------------------------------------------------------------------
# Tests: Dockerfile
# ---------------------------------------------------------------------------


class TestDockerfile:
    """Validate Dockerfile existence and basic syntax."""

    def test_dockerfile_exists(self):
        """The Dockerfile must exist at the project root."""
        assert DOCKERFILE.exists(), f"Missing {DOCKERFILE}"

    def test_dockerfile_not_empty(self):
        """The Dockerfile must not be empty."""
        assert DOCKERFILE.stat().st_size > 0, "Dockerfile is empty"

    def test_dockerfile_has_from(self):
        """The Dockerfile must contain at least one FROM instruction."""
        content = DOCKERFILE.read_text()
        from_lines = [
            line for line in content.splitlines()
            if re.match(r"^\s*FROM\s+", line, re.IGNORECASE)
        ]
        assert len(from_lines) >= 1, "Dockerfile has no FROM instruction"

    def test_dockerfile_is_multistage(self):
        """The Dockerfile must use multi-stage build (multiple FROM + AS)."""
        content = DOCKERFILE.read_text()
        stage_lines = re.findall(
            r"^\s*FROM\s+.+\s+AS\s+\w+", content, re.MULTILINE | re.IGNORECASE
        )
        assert len(stage_lines) >= 2, (
            f"Expected at least 2 build stages (FROM ... AS ...), "
            f"found {len(stage_lines)}"
        )

    def test_dockerfile_has_healthcheck(self):
        """The Dockerfile should contain a HEALTHCHECK instruction."""
        content = DOCKERFILE.read_text()
        assert "HEALTHCHECK" in content, "Dockerfile missing HEALTHCHECK"

    def test_dockerfile_has_nonroot_user(self):
        """The Dockerfile should create and switch to a non-root user."""
        content = DOCKERFILE.read_text()
        assert "USER" in content, "Dockerfile missing USER instruction (non-root)"

    def test_dockerfile_has_workdir(self):
        """The Dockerfile should set a WORKDIR."""
        content = DOCKERFILE.read_text()
        assert "WORKDIR" in content, "Dockerfile missing WORKDIR instruction"

    def test_dockerfile_sets_pythonpath(self):
        """The Dockerfile should set PYTHONPATH in at least one stage."""
        content = DOCKERFILE.read_text()
        assert "PYTHONPATH" in content, "Dockerfile missing PYTHONPATH"

    def test_dockerfile_has_cmd_or_entrypoint(self):
        """The Dockerfile should define a default CMD or ENTRYPOINT."""
        content = DOCKERFILE.read_text()
        has_cmd = bool(re.search(r"^\s*CMD\s+", content, re.MULTILINE))
        has_entrypoint = bool(re.search(r"^\s*ENTRYPOINT\s+", content, re.MULTILINE))
        assert has_cmd or has_entrypoint, (
            "Dockerfile missing CMD or ENTRYPOINT"
        )

    def test_dockerfile_copies_requirements(self):
        """The Dockerfile should COPY requirements files for caching."""
        content = DOCKERFILE.read_text()
        assert "requirements" in content.lower(), (
            "Dockerfile does not reference any requirements file"
        )

    def test_dockerfile_python311(self):
        """The Dockerfile should use Python 3.11."""
        content = DOCKERFILE.read_text()
        assert "python:3.11" in content, (
            "Dockerfile does not use python:3.11 base image"
        )


# ---------------------------------------------------------------------------
# Tests: .dockerignore
# ---------------------------------------------------------------------------


class TestDockerignore:
    """Validate .dockerignore existence and content."""

    def test_dockerignore_exists(self):
        """The .dockerignore file must exist at the project root."""
        assert DOCKERIGNORE.exists(), f"Missing {DOCKERIGNORE}"

    def test_dockerignore_not_empty(self):
        """The .dockerignore file must not be empty."""
        assert DOCKERIGNORE.stat().st_size > 0, ".dockerignore is empty"

    def test_dockerignore_excludes_git(self):
        """The .dockerignore should exclude the .git directory."""
        content = DOCKERIGNORE.read_text()
        assert ".git" in content, ".dockerignore should exclude .git"

    def test_dockerignore_excludes_env(self):
        """The .dockerignore should exclude .env files (secrets)."""
        content = DOCKERIGNORE.read_text()
        assert ".env" in content, ".dockerignore should exclude .env"

    def test_dockerignore_excludes_pycache(self):
        """The .dockerignore should exclude __pycache__."""
        content = DOCKERIGNORE.read_text()
        assert "__pycache__" in content, ".dockerignore should exclude __pycache__"
