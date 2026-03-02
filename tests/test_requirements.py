"""Validate requirements files: format, no duplicates, no missing packages."""
from __future__ import annotations

import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# Matches a PEP 508-ish requirement line (package name with optional extras).
# Examples: "fastapi==0.116.1", "uvicorn[standard]==0.35.0", "numpy>=1.24.0"
_REQ_RE = re.compile(
    r"^([A-Za-z0-9]([A-Za-z0-9._-]*[A-Za-z0-9])?)"  # package name
    r"(\[[A-Za-z0-9,._-]+\])?"                        # optional extras [standard]
    r"(==|>=|<=|!=|~=|>|<)"                            # version operator
    r"[A-Za-z0-9.*+!,><=~\s]+$"                       # version specifier(s)
)

_INCLUDE_RE = re.compile(r"^-r\s+\S+")   # e.g. "-r requirements.txt"
_COMMENT_RE = re.compile(r"^\s*#")        # comment line


def _strip_inline_comment(line: str) -> str:
    """Remove trailing inline comments (e.g. ``pkg==1.0  # note``)."""
    # Split on `` # `` (space-hash-space) which is the pip convention.
    # Also handle tab-separated comments.
    idx = line.find("  #")
    if idx != -1:
        return line[:idx].rstrip()
    return line


def _parse_requirements(path: Path) -> list[tuple[str, str]]:
    """Return list of (normalised_name, raw_line) from a requirements file.

    Skips blank lines, comments, and ``-r`` include directives.
    """
    entries: list[tuple[str, str]] = []
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line or _COMMENT_RE.match(line) or _INCLUDE_RE.match(line):
            continue
        spec = _strip_inline_comment(line)
        m = _REQ_RE.match(spec)
        assert m, f"Malformed requirement line in {path.name}: {line!r}"
        # PEP 503 normalisation: lower-case, underscores/dots to hyphens
        name = re.sub(r"[-_.]+", "-", m.group(1)).lower()
        entries.append((name, line))
    return entries


# ---------------------------------------------------------------------------
# Tests for requirements.txt
# ---------------------------------------------------------------------------

class TestRequirementsTxt:
    """Tests for the main requirements.txt file."""

    REQ_PATH = ROOT / "requirements.txt"

    def test_file_exists(self) -> None:
        assert self.REQ_PATH.is_file(), "requirements.txt not found at project root"

    def test_no_duplicates(self) -> None:
        entries = _parse_requirements(self.REQ_PATH)
        seen: dict[str, str] = {}
        duplicates: list[str] = []
        for name, line in entries:
            if name in seen:
                duplicates.append(f"{name} (first: {seen[name]!r}, duplicate: {line!r})")
            seen[name] = line
        assert not duplicates, f"Duplicate packages: {duplicates}"

    def test_all_lines_valid_format(self) -> None:
        """Every non-blank, non-comment line must match requirement format."""
        for raw in self.REQ_PATH.read_text().splitlines():
            line = raw.strip()
            if not line or _COMMENT_RE.match(line) or _INCLUDE_RE.match(line):
                continue
            spec = _strip_inline_comment(line)
            assert _REQ_RE.match(spec), f"Malformed line: {line!r}"

    def test_no_trailing_whitespace(self) -> None:
        for i, raw in enumerate(self.REQ_PATH.read_text().splitlines(), 1):
            assert raw == raw.rstrip(), f"Trailing whitespace on line {i}: {raw!r}"

    def test_numpy_present(self) -> None:
        """numpy is imported across many modules and must be listed."""
        entries = _parse_requirements(self.REQ_PATH)
        names = {n for n, _ in entries}
        assert "numpy" in names, "numpy missing from requirements.txt"

    def test_no_unused_fair_esm(self) -> None:
        """fair-esm is not imported anywhere; ESM-2 is loaded via transformers."""
        entries = _parse_requirements(self.REQ_PATH)
        names = {n for n, _ in entries}
        assert "fair-esm" not in names, "fair-esm should not be in requirements.txt"


# ---------------------------------------------------------------------------
# Tests for requirements-dev.txt
# ---------------------------------------------------------------------------

class TestRequirementsDevTxt:
    """Tests for the dev requirements file."""

    DEV_PATH = ROOT / "requirements-dev.txt"

    def test_file_exists(self) -> None:
        assert self.DEV_PATH.is_file(), "requirements-dev.txt not found"

    def test_includes_base(self) -> None:
        """requirements-dev.txt must include -r requirements.txt."""
        text = self.DEV_PATH.read_text()
        assert "-r requirements.txt" in text, "Missing '-r requirements.txt' include"

    def test_pytest_listed(self) -> None:
        entries = _parse_requirements(self.DEV_PATH)
        names = {n for n, _ in entries}
        assert "pytest" in names, "pytest missing from requirements-dev.txt"

    def test_pytest_asyncio_listed(self) -> None:
        entries = _parse_requirements(self.DEV_PATH)
        names = {n for n, _ in entries}
        assert "pytest-asyncio" in names, "pytest-asyncio missing from requirements-dev.txt"

    def test_no_duplicates(self) -> None:
        entries = _parse_requirements(self.DEV_PATH)
        seen: dict[str, str] = {}
        duplicates: list[str] = []
        for name, line in entries:
            if name in seen:
                duplicates.append(f"{name} (first: {seen[name]!r}, dup: {line!r})")
            seen[name] = line
        assert not duplicates, f"Duplicate packages: {duplicates}"


# requirements-ml.txt removed (migrated to imgt-ml repo)
