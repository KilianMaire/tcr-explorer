"""Integration tests for MHC source routing in the main /search handler."""
from __future__ import annotations

import sys
from pathlib import Path
from unittest.mock import AsyncMock, patch

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from fastapi.testclient import TestClient


# ---------------------------------------------------------------------------
# Fake responses
# ---------------------------------------------------------------------------

_MHC_RECORD = {
    "source": "mhc",
    "gene_name": "Mamu-A1",
    "species": "other",
    "allele_name": "Mamu-A1*026:01:01:01",
    "sequence": "ATGGCGGTC",
    "region": "coding",
    "metadata": {
        "accession": "NHP01224",
        "mhc_class": "I",
        "organism_common": "rhesus monkey",
        "organism_group": "NHP",
        "backend": "ebi-ipd-mhc",
    },
}

_MHC_RESP = {"total": 1, "records": [_MHC_RECORD], "limit": 50, "offset": 0}
_EMPTY_RESP = {"total": 0, "records": [], "limit": 50, "offset": 0}


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture()
def api_client():
    from imgt_app.api import app
    return TestClient(app)


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

def test_mhc_source_routes_to_mhc_client(api_client):
    """source='mhc' should route to mhc_client.search."""
    with patch("imgt_app.api.mhc_client") as mock_mhc:
        mock_mhc.search = AsyncMock(return_value=_MHC_RESP)
        resp = api_client.post("/search", json={
            "source": "mhc",
            "gene_name": "Mamu-A1",
        })

    assert resp.status_code == 200
    data = resp.json()
    assert data["total"] == 1
    rec = data["records"][0]
    assert rec["source"] == "mhc"
    assert rec["gene_name"] == "Mamu-A1"
    assert rec["metadata"]["organism_group"] == "NHP"
    mock_mhc.search.assert_called_once()


def test_mhc_source_not_called_for_hla(api_client):
    """source='hla' should NOT call mhc_client."""
    with (
        patch("imgt_app.api.hla_client") as mock_hla,
        patch("imgt_app.api.mhc_client") as mock_mhc,
    ):
        mock_hla.search = AsyncMock(return_value=_EMPTY_RESP)
        mock_mhc.search = AsyncMock(return_value=_EMPTY_RESP)
        resp = api_client.post("/search", json={"source": "hla", "gene_name": "HLA-A"})

    assert resp.status_code == 200
    mock_mhc.search.assert_not_called()


def test_mhc_source_not_called_for_vdjdb(api_client):
    """source='vdjdb' should NOT call mhc_client."""
    with (
        patch("imgt_app.api.vdjdb_client") as mock_vdjdb,
        patch("imgt_app.api.iedb_client") as mock_iedb,
        patch("imgt_app.api.mhc_client") as mock_mhc,
    ):
        mock_vdjdb.search = AsyncMock(return_value=_EMPTY_RESP)
        mock_iedb.search = AsyncMock(return_value=_EMPTY_RESP)
        mock_mhc.search = AsyncMock(return_value=_EMPTY_RESP)
        resp = api_client.post("/search", json={"source": "vdjdb"})

    assert resp.status_code == 200
    mock_mhc.search.assert_not_called()


def test_mhc_server_failure_raises(api_client):
    """If mhc_client raises, the exception propagates (uncaught for non-vdjdb sources)."""
    with patch("imgt_app.api.mhc_client") as mock_mhc:
        mock_mhc.search = AsyncMock(side_effect=Exception("MHC server down"))
        with pytest.raises(Exception, match="MHC server down"):
            api_client.post("/search", json={
                "source": "mhc",
                "gene_name": "Mamu-A1",
            })


def test_mhc_source_accepted_in_request():
    """Pydantic should accept source='mhc' in SearchRequest."""
    from imgt_app.models import SearchRequest
    req = SearchRequest(source="mhc", gene_name="Mamu-A1")
    assert req.source == "mhc"
