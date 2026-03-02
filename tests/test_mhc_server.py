"""Tests for the IPD-MHC tool server (servers/mhc_server.py)."""
from __future__ import annotations

import asyncio
import sys
from pathlib import Path
from unittest.mock import AsyncMock, patch, MagicMock

import httpx
import pytest

_src = Path(__file__).resolve().parent.parent / "src"
if str(_src) not in sys.path:
    sys.path.insert(0, str(_src))

from fastapi.testclient import TestClient

from servers.mhc_server import (
    _detect_project,
    _build_query_prefix,
    _allele_to_record,
    app,
)


def _run(coro):
    """Run an async coroutine synchronously (no pytest-asyncio needed)."""
    return asyncio.get_event_loop().run_until_complete(coro)


# ---------------------------------------------------------------------------
# Unit tests: helper functions
# ---------------------------------------------------------------------------

class TestDetectProject:
    def test_mamu_prefix(self):
        assert _detect_project("Mamu-A1") == "MHC"

    def test_bola_prefix(self):
        assert _detect_project("BoLA-2") == "MHC"

    def test_patr_prefix(self):
        assert _detect_project("Patr-A") == "MHC"

    def test_dla_prefix(self):
        assert _detect_project("DLA-88") == "MHC"

    def test_sla_prefix(self):
        assert _detect_project("SLA-1") == "MHC"

    def test_rt1_prefix(self):
        assert _detect_project("RT1-A") == "MHC"

    def test_hla_prefix(self):
        assert _detect_project("HLA-A") == "HLA"

    def test_bare_locus_with_star(self):
        assert _detect_project("A*02") == "HLA"

    def test_unknown_defaults_to_mhc(self):
        assert _detect_project("UNKNOWN-X") == "MHC"

    def test_case_insensitive(self):
        assert _detect_project("mamu-b") == "MHC"
        assert _detect_project("BOLA-1") == "MHC"


class TestBuildQueryPrefix:
    def test_mamu_kept_as_is(self):
        assert _build_query_prefix("Mamu-A1") == "Mamu-A1"

    def test_bola_kept_as_is(self):
        assert _build_query_prefix("BoLA-2") == "BoLA-2"

    def test_hla_prefix_stripped(self):
        assert _build_query_prefix("HLA-A") == "A"

    def test_hla_drb1_stripped(self):
        assert _build_query_prefix("HLA-DRB1") == "DRB1"

    def test_bare_locus_unchanged(self):
        assert _build_query_prefix("DRB1") == "DRB1"

    def test_whitespace_stripped(self):
        assert _build_query_prefix("  Mamu-A1  ") == "Mamu-A1"


class TestAlleleToRecord:
    def test_nhp_allele(self):
        item = {
            "accession": "NHP01224",
            "name": "Mamu-A1*026:01:01:01",
            "class": "I",
            "locus": "A1",
            "organism": {
                "commonName": "rhesus monkey",
                "scientificName": "Macaca mulatta",
                "group": "NHP",
                "taxon": 9544,
            },
            "sequence": {
                "coding": "ATGGCGGTC",
                "genomic": "ATGGCGGTCAAAA",
                "protein": "MAVMAPRTL",
            },
            "date_assigned": "2010-01-01",
            "date_modified": "2021-06-15",
            "status": "public",
        }
        rec = _allele_to_record(item, "MHC")
        assert rec["source"] == "mhc"
        assert rec["species"] == "other"
        assert rec["gene_name"] == "Mamu-A1"
        assert rec["allele_name"] == "Mamu-A1*026:01:01:01"
        assert rec["region"] == "coding"
        assert rec["sequence"] == "ATGGCGGTC"
        assert rec["metadata"]["accession"] == "NHP01224"
        assert rec["metadata"]["mhc_class"] == "I"
        assert rec["metadata"]["organism_common"] == "rhesus monkey"
        assert rec["metadata"]["organism_group"] == "NHP"
        assert rec["metadata"]["taxon"] == 9544
        assert rec["metadata"]["backend"] == "ebi-ipd-mhc"

    def test_hla_allele_gets_hla_prefix(self):
        item = {
            "accession": "HLA00001",
            "name": "A*01:01:01:01",
            "class": "I",
            "locus": "A",
            "organism": {
                "commonName": "human",
                "scientificName": "Homo sapiens",
                "group": "HLA",
                "taxon": 9606,
            },
            "sequence": {"coding": "ATCGATCG", "genomic": "", "protein": "MK"},
        }
        rec = _allele_to_record(item, "HLA")
        assert rec["gene_name"] == "HLA-A"
        assert rec["species"] == "human"

    def test_missing_sequence(self):
        item = {
            "accession": "NHP99999",
            "name": "Mamu-B*099:01",
            "organism": {},
        }
        rec = _allele_to_record(item, "MHC")
        assert rec["sequence"] == ""
        assert rec["metadata"]["protein_seq"] == ""

    def test_missing_organism(self):
        item = {
            "accession": "NHP00001",
            "name": "Mamu-A*01:01",
        }
        rec = _allele_to_record(item, "MHC")
        assert rec["species"] == "other"
        assert rec["metadata"]["organism_common"] is None


# ---------------------------------------------------------------------------
# Integration tests: /health and /search endpoints
# ---------------------------------------------------------------------------

@pytest.fixture()
def client():
    return TestClient(app)


class TestHealthEndpoint:
    def test_health_ok(self, client):
        resp = client.get("/health")
        assert resp.status_code == 200
        data = resp.json()
        assert data["status"] == "ok"
        assert data["server"] == "mhc"


class TestSearchEndpoint:
    def test_empty_gene_returns_empty(self, client):
        resp = client.post("/search", json={"gene_name": "", "limit": 10})
        assert resp.status_code == 200
        assert resp.json()["total"] == 0
        assert resp.json()["records"] == []

    def test_missing_gene_returns_empty(self, client):
        resp = client.post("/search", json={"limit": 10})
        assert resp.status_code == 200
        assert resp.json()["total"] == 0

    @patch("servers.mhc_server._fetch_ebi")
    def test_search_returns_records(self, mock_fetch, client):
        mock_fetch.return_value = [
            {
                "source": "mhc",
                "species": "other",
                "gene_name": "Mamu-A1",
                "allele_name": "Mamu-A1*026:01:01:01",
                "region": "coding",
                "sequence": "ATGGCG",
                "metadata": {"accession": "NHP01224", "backend": "ebi-ipd-mhc"},
            }
        ]
        resp = client.post("/search", json={"gene_name": "Mamu-A1", "limit": 10})
        assert resp.status_code == 200
        data = resp.json()
        assert data["total"] == 1
        assert data["records"][0]["gene_name"] == "Mamu-A1"
        assert data["records"][0]["source"] == "mhc"

    @patch("servers.mhc_server._fetch_ebi")
    def test_search_respects_limit(self, mock_fetch, client):
        mock_fetch.return_value = [
            {"source": "mhc", "species": "other", "gene_name": "Mamu-A1",
             "allele_name": f"Mamu-A1*{i:03d}:01", "region": "coding",
             "sequence": "ATG", "metadata": {}}
            for i in range(5)
        ]
        resp = client.post("/search", json={"gene_name": "Mamu-A1", "limit": 2})
        assert resp.status_code == 200
        assert resp.json()["total"] == 2

    @patch("servers.mhc_server._fetch_ebi")
    def test_search_region_filter(self, mock_fetch, client):
        mock_fetch.return_value = [
            {"source": "mhc", "species": "other", "gene_name": "Mamu-A1",
             "allele_name": "Mamu-A1*001:01", "region": "coding",
             "sequence": "ATG", "metadata": {}},
            {"source": "mhc", "species": "other", "gene_name": "Mamu-A1",
             "allele_name": "Mamu-A1*002:01", "region": "genomic",
             "sequence": "ATG", "metadata": {}},
        ]
        resp = client.post("/search", json={
            "gene_name": "Mamu-A1", "region": "coding", "limit": 10
        })
        assert resp.status_code == 200
        data = resp.json()
        assert data["total"] == 1
        assert data["records"][0]["region"] == "coding"


# ---------------------------------------------------------------------------
# Tests for _fetch_ebi with mocked httpx
# ---------------------------------------------------------------------------

def _mock_client_for(list_resp, *single_resps):
    """Build a mocked httpx.AsyncClient context manager."""
    mc = AsyncMock()
    mc.get = AsyncMock(side_effect=[list_resp, *single_resps])
    mc.__aenter__ = AsyncMock(return_value=mc)
    mc.__aexit__ = AsyncMock(return_value=False)
    return mc


def _ok_response(json_data):
    r = MagicMock()
    r.status_code = 200
    r.raise_for_status = MagicMock()
    r.json.return_value = json_data
    return r


class TestFetchEbi:
    def test_fetch_ebi_success(self):
        """Mock the full EBI API flow: list -> parallel individual fetches."""
        from servers.mhc_server import _fetch_ebi, _CACHE
        _CACHE.clear()

        list_resp = _ok_response({
            "data": [{"accession": "NHP01224", "name": "Mamu-A1*026:01:01:01"}],
            "meta": {"total": 1},
        })
        single_resp = _ok_response({
            "accession": "NHP01224",
            "name": "Mamu-A1*026:01:01:01",
            "class": "I",
            "locus": "A1",
            "organism": {
                "commonName": "rhesus monkey",
                "scientificName": "Macaca mulatta",
                "group": "NHP",
                "taxon": 9544,
            },
            "sequence": {"coding": "ATGGCGGTC", "genomic": "", "protein": "MAV"},
        })

        mc = _mock_client_for(list_resp, single_resp)
        with patch("servers.mhc_server.httpx.AsyncClient", return_value=mc):
            results = _run(_fetch_ebi("Mamu-A1", "", 10))

        assert len(results) == 1
        assert results[0]["gene_name"] == "Mamu-A1"
        assert results[0]["sequence"] == "ATGGCGGTC"
        assert results[0]["metadata"]["organism_group"] == "NHP"

    def test_fetch_ebi_empty_list(self):
        """When EBI returns no alleles, result is empty."""
        from servers.mhc_server import _fetch_ebi, _CACHE
        _CACHE.clear()

        list_resp = _ok_response({"data": [], "meta": {"total": 0}})
        mc = _mock_client_for(list_resp)
        with patch("servers.mhc_server.httpx.AsyncClient", return_value=mc):
            results = _run(_fetch_ebi("UNKNOWN-X", "", 10))

        assert results == []

    def test_fetch_ebi_network_error(self):
        """Network errors return empty list, don't raise."""
        from servers.mhc_server import _fetch_ebi, _CACHE
        _CACHE.clear()

        mc = AsyncMock()
        mc.get = AsyncMock(side_effect=httpx.ConnectError("timeout"))
        mc.__aenter__ = AsyncMock(return_value=mc)
        mc.__aexit__ = AsyncMock(return_value=False)

        with patch("servers.mhc_server.httpx.AsyncClient", return_value=mc):
            results = _run(_fetch_ebi("Mamu-A1", "", 10))

        assert results == []

    def test_fetch_ebi_seq_contains_filter(self):
        """sequence_contains filters results post-fetch."""
        from servers.mhc_server import _fetch_ebi, _CACHE
        _CACHE.clear()

        list_resp = _ok_response({
            "data": [
                {"accession": "NHP001", "name": "Mamu-A1*001:01"},
                {"accession": "NHP002", "name": "Mamu-A1*002:01"},
            ],
            "meta": {"total": 2},
        })

        single1 = _ok_response({
            "accession": "NHP001", "name": "Mamu-A1*001:01",
            "organism": {"commonName": "rhesus monkey", "group": "NHP"},
            "sequence": {"coding": "ATGGCGGTCAAA", "genomic": "", "protein": ""},
        })
        single2 = _ok_response({
            "accession": "NHP002", "name": "Mamu-A1*002:01",
            "organism": {"commonName": "rhesus monkey", "group": "NHP"},
            "sequence": {"coding": "ATGCCCTTTAGG", "genomic": "", "protein": ""},
        })

        mc = _mock_client_for(list_resp, single1, single2)
        with patch("servers.mhc_server.httpx.AsyncClient", return_value=mc):
            results = _run(_fetch_ebi("Mamu-A1", "GCGGTC", 10))

        assert len(results) == 1
        assert "GCGGTC" in results[0]["sequence"]


# ---------------------------------------------------------------------------
# Cache tests
# ---------------------------------------------------------------------------

class TestCache:
    def test_cache_hit_avoids_second_request(self):
        """Second call with same params should use cache, not hit API."""
        from servers.mhc_server import _fetch_ebi, _CACHE
        _CACHE.clear()

        list_resp = _ok_response({
            "data": [{"accession": "NHP01224", "name": "Mamu-A1*026:01"}],
            "meta": {"total": 1},
        })
        single_resp = _ok_response({
            "accession": "NHP01224", "name": "Mamu-A1*026:01",
            "organism": {"commonName": "rhesus monkey", "group": "NHP"},
            "sequence": {"coding": "ATGGCG", "genomic": "", "protein": ""},
        })

        mc = _mock_client_for(list_resp, single_resp)
        with patch("servers.mhc_server.httpx.AsyncClient", return_value=mc):
            r1 = _run(_fetch_ebi("Mamu-A1", "", 10))
            # Second call should use cache
            r2 = _run(_fetch_ebi("Mamu-A1", "", 10))

        assert len(r1) == 1
        assert len(r2) == 1
        # Only 2 get() calls total (list + single from first call)
        assert mc.get.call_count == 2
