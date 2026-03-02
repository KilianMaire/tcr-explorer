"""Tests for the project Makefile."""

from __future__ import annotations

import pathlib
import subprocess

import pytest

ROOT = pathlib.Path(__file__).resolve().parent.parent
MAKEFILE = ROOT / "Makefile"


# ── 1. Makefile exists and is valid ──────────────────────────────────────────

def test_makefile_exists():
    """Makefile must exist at the project root."""
    assert MAKEFILE.is_file(), f"Makefile not found at {MAKEFILE}"


def test_makefile_is_valid():
    """Makefile must be parseable by make (dry-run on 'help')."""
    result = subprocess.run(
        ["make", "-n", "help"],
        cwd=str(ROOT),
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, (
        f"make -n help failed:\nstdout: {result.stdout}\nstderr: {result.stderr}"
    )


# ── 2. All key targets are defined ──────────────────────────────────────────

EXPECTED_TARGETS = [
    "help",
    "install",
    "install-dev",
    "test",
    "test-cov",
    "lint",
    "format",
    "serve",
    "serve-all",
    "clean",
]


@pytest.mark.parametrize("target", EXPECTED_TARGETS)
def test_target_defined(target: str):
    """Each expected target must be recognised by make."""
    result = subprocess.run(
        ["make", "-n", target],
        cwd=str(ROOT),
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, (
        f"Target '{target}' not found or invalid:\n{result.stderr}"
    )


# ── 3. Help target works and contains descriptions ──────────────────────────

def test_help_output_contains_descriptions():
    """Running `make help` should list targets with '##' descriptions."""
    result = subprocess.run(
        ["make", "help"],
        cwd=str(ROOT),
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, (
        f"make help failed:\nstdout: {result.stdout}\nstderr: {result.stderr}"
    )

    output = result.stdout

    # At minimum, these key targets should appear in help output
    must_appear = [
        "install",
        "test",
        "serve",
        "clean",
    ]
    for target in must_appear:
        assert target in output, (
            f"Expected target '{target}' not found in help output:\n{output}"
        )

    # Each line that shows a target must also show a description
    lines_with_targets = [
        line for line in output.strip().splitlines() if line.strip()
    ]
    assert len(lines_with_targets) >= 10, (
        f"Help output has too few lines ({len(lines_with_targets)}); "
        f"expected at least 10 documented targets."
    )
