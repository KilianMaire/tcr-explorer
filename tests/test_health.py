"""Tests for /health endpoint on the hla and mhc tool servers."""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

# Ensure src/ is importable (the servers import cdr_enricher).
_src = Path(__file__).resolve().parent.parent / "src"
if str(_src) not in sys.path:
    sys.path.insert(0, str(_src))

from fastapi.testclient import TestClient  # noqa: E402


@pytest.fixture()
def hla_client():
    from servers.hla_server import app
    return TestClient(app)


@pytest.fixture()
def mhc_client():
    from servers.mhc_server import app
    return TestClient(app)


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("fixture_name,expected_server", [
    ("hla_client", "hla"),
    ("mhc_client", "mhc"),
])
def test_health_returns_ok_and_server_name(fixture_name, expected_server, request):
    client = request.getfixturevalue(fixture_name)
    resp = client.get("/health")
    assert resp.status_code == 200
    data = resp.json()
    assert data["status"] == "ok"
    assert data["server"] == expected_server


def test_health_response_keys(hla_client):
    """Health response must contain exactly 'status' and 'server'."""
    data = hla_client.get("/health").json()
    assert set(data.keys()) == {"status", "server"}
