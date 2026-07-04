"""The reduced /search endpoint serves only hla/mhc; TCR record search moved to /v1/tcr/records."""
from fastapi.testclient import TestClient

from tcr_explorer.api import app

client = TestClient(app)


def test_search_rejects_vdjdb_source():
    r = client.post("/search", json={"source": "vdjdb", "limit": 5})
    assert r.status_code == 400
    assert "/v1/tcr/records" in r.json()["detail"]


def test_search_rejects_tcr_source():
    r = client.post("/search", json={"source": "tcr", "limit": 5})
    assert r.status_code == 400


def test_search_rejects_iedb_source():
    r = client.post("/search", json={"source": "iedb", "limit": 5})
    assert r.status_code == 400
