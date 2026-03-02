"""Integration tests for IEDB enrichment wired into the /search handler."""
from __future__ import annotations

import sys
from pathlib import Path
from unittest.mock import AsyncMock, patch

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from fastapi.testclient import TestClient


# ---------------------------------------------------------------------------
# Fake responses returned by mocked tool-server clients
# ---------------------------------------------------------------------------

_VDJDB_RECORD = {
    "source": "vdjdb",
    "gene_name": "TRBV19",
    "species": "human",
    "sequence": "CASSIRSSYEQYF",
    "region": "cdr3",
    "antigen_epitope": "GILGFVFTL",
    "metadata": {},
}

_IEDB_RECORD = {
    "source": "iedb",
    "gene_name": "HLA-A*02:01",
    "species": "human",
    "sequence": "GILGFVFTL",
    "region": "epitope",
    "metadata": {
        "mhc_class": "I",
        "antigen_name": "Matrix protein 1",
        "source_organism": "Influenza A virus",
        "assay_type": "cytokine release",
        "effector_cell_type": "CD8+",
        "qualitative_measure": "Positive",
    },
}

_VDJDB_RESP = {"total": 1, "records": [_VDJDB_RECORD], "limit": 50, "offset": 0}
_IEDB_RESP  = {"total": 1, "records": [_IEDB_RECORD],  "limit": 50, "offset": 0}
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

def test_vdjdb_result_enriched_with_iedb_hits(api_client):
    from imgt_app.models import SearchResponse
    _local_empty = SearchResponse(total=0, records=[], limit=50, offset=0)
    with (
        patch("imgt_app.api.vdjdb_client") as mock_vdjdb,
        patch("imgt_app.api.iedb_client") as mock_iedb,
        patch("imgt_app.api.index") as mock_index,
    ):
        mock_vdjdb.search = AsyncMock(return_value=_VDJDB_RESP)
        mock_iedb.search  = AsyncMock(return_value=_IEDB_RESP)
        mock_index.search.return_value = _local_empty

        resp = api_client.post("/search", json={"source": "vdjdb", "species": "human"})

    assert resp.status_code == 200
    data = resp.json()
    assert data["total"] == 1
    records = data["records"]
    enriched = [r for r in records if r["sequence"] == "CASSIRSSYEQYF"]
    assert len(enriched) == 1
    hits = enriched[0]["iedb_hits"]
    assert hits is not None
    assert len(hits) == 1
    assert hits[0]["epitope_sequence"] == "GILGFVFTL"
    assert hits[0]["mhc_allele"] == "HLA-A*02:01"


def test_orphan_iedb_hit_creates_phantom(api_client):
    from imgt_app.models import SearchResponse
    _local_empty = SearchResponse(total=0, records=[], limit=50, offset=0)
    with (
        patch("imgt_app.api.vdjdb_client") as mock_vdjdb,
        patch("imgt_app.api.iedb_client") as mock_iedb,
        patch("imgt_app.api.index") as mock_index,
    ):
        mock_vdjdb.search = AsyncMock(return_value=_EMPTY_RESP)
        mock_iedb.search  = AsyncMock(return_value=_IEDB_RESP)
        mock_index.search.return_value = _local_empty

        resp = api_client.post("/search", json={"source": "vdjdb"})

    assert resp.status_code == 200
    data = resp.json()
    assert data["total"] == 0
    records = data["records"]
    phantoms = [r for r in records if r["sequence"] == ""]
    assert len(phantoms) == 1
    assert phantoms[0]["gene_name"] == "GILGFVFTL"
    assert phantoms[0]["iedb_hits"] is not None


def test_iedb_failure_returns_vdjdb_records_without_hits(api_client):
    """If IEDB is down, VDJdb records are returned as-is (no error raised)."""
    from imgt_app.models import SearchResponse
    _local_empty = SearchResponse(total=0, records=[], limit=50, offset=0)
    with (
        patch("imgt_app.api.vdjdb_client") as mock_vdjdb,
        patch("imgt_app.api.iedb_client") as mock_iedb,
        patch("imgt_app.api.index") as mock_index,
    ):
        mock_vdjdb.search = AsyncMock(return_value=_VDJDB_RESP)
        mock_iedb.search  = AsyncMock(side_effect=Exception("IEDB unreachable"))
        mock_index.search.return_value = _local_empty

        resp = api_client.post("/search", json={"source": "vdjdb"})

    assert resp.status_code == 200
    records = resp.json()["records"]
    assert len(records) == 1
    assert records[0]["sequence"] == "CASSIRSSYEQYF"
    assert records[0]["iedb_hits"] is None


def test_non_vdjdb_source_not_enriched(api_client):
    """Searching HLA does not trigger IEDB enrichment."""
    with (
        patch("imgt_app.api.hla_client") as mock_hla,
        patch("imgt_app.api.iedb_client") as mock_iedb,
    ):
        mock_hla.search = AsyncMock(return_value={
            "total": 1,
            "records": [{"source": "hla", "gene_name": "HLA-A*02:01",
                         "sequence": "ATCG", "species": "human", "metadata": {}}],
            "limit": 50, "offset": 0,
        })
        mock_iedb.search = AsyncMock(return_value=_IEDB_RESP)

        resp = api_client.post("/search", json={"source": "hla", "gene_name": "HLA-A"})

    assert resp.status_code == 200
    # IEDB must NOT have been called
    mock_iedb.search.assert_not_called()
