"""Tests for pagination (limit / offset) on SearchRequest, SearchResponse,
SearchIndex, and the _merge_results helper in api.py."""
from __future__ import annotations

import json
import sys
import tempfile
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from tcr_explorer.models import GeneRecord, SearchRequest, SearchResponse


# ---------------------------------------------------------------------------
# Helper — build a GeneRecord with minimal required fields
# ---------------------------------------------------------------------------

def _rec(gene: str, seq: str = "ATGC") -> GeneRecord:
    return GeneRecord(source="hla", species="human", gene_name=gene, sequence=seq)


# ---------------------------------------------------------------------------
# Model validation
# ---------------------------------------------------------------------------

def test_search_request_default_offset():
    req = SearchRequest()
    assert req.offset == 0


def test_search_request_offset_accepted():
    req = SearchRequest(offset=10, limit=5)
    assert req.offset == 10
    assert req.limit == 5


def test_search_request_offset_negative_rejected():
    with pytest.raises(Exception):
        SearchRequest(offset=-1)


def test_search_request_offset_too_large_rejected():
    with pytest.raises(Exception):
        SearchRequest(offset=10001)


def test_search_response_default_pagination_fields():
    resp = SearchResponse(total=0, records=[])
    assert resp.limit == 50
    assert resp.offset == 0


def test_search_response_pagination_fields_set():
    resp = SearchResponse(total=100, records=[], limit=10, offset=20)
    assert resp.limit == 10
    assert resp.offset == 20


# ---------------------------------------------------------------------------
# _merge_results pagination logic
# ---------------------------------------------------------------------------

def test_merge_results_no_offset():
    """With offset=0, result should be the first `limit` records."""
    from tcr_explorer.api import _merge_results

    local = SearchResponse(total=3, records=[_rec("A"), _rec("B"), _rec("C")])
    remote = SearchResponse(total=2, records=[_rec("D"), _rec("E")])
    req = SearchRequest(limit=3, offset=0)

    result = _merge_results(local, remote, req)
    assert result.total == 5
    assert len(result.records) == 3
    assert [r.gene_name for r in result.records] == ["A", "B", "C"]
    assert result.limit == 3
    assert result.offset == 0


def test_merge_results_with_offset():
    """Offset skips the first N combined records."""
    from tcr_explorer.api import _merge_results

    local = SearchResponse(total=3, records=[_rec("A"), _rec("B"), _rec("C")])
    remote = SearchResponse(total=2, records=[_rec("D"), _rec("E")])
    req = SearchRequest(limit=2, offset=2)

    result = _merge_results(local, remote, req)
    assert result.total == 5
    assert len(result.records) == 2
    assert [r.gene_name for r in result.records] == ["C", "D"]
    assert result.limit == 2
    assert result.offset == 2


def test_merge_results_offset_beyond_available():
    """Offset larger than available records returns empty list."""
    from tcr_explorer.api import _merge_results

    local = SearchResponse(total=2, records=[_rec("A"), _rec("B")])
    remote = SearchResponse(total=0, records=[])
    req = SearchRequest(limit=10, offset=5)

    result = _merge_results(local, remote, req)
    assert result.total == 2
    assert result.records == []


# ---------------------------------------------------------------------------
# SearchIndex pagination
# ---------------------------------------------------------------------------

def test_search_index_fetch_limit_is_inflated():
    """SQLite search fetches limit+offset rows so merge has enough to slice."""
    from tcr_explorer.search_index import SearchIndex

    with tempfile.TemporaryDirectory() as tmpdir:
        idx = SearchIndex(f"{tmpdir}/test.db")
        # Insert 10 records
        records = [_rec(f"GENE-{i:02d}") for i in range(10)]
        idx.upsert_many(records)

        # Page 1: offset=0, limit=3
        req1 = SearchRequest(limit=3, offset=0)
        resp1 = idx.search(req1)
        assert resp1.total == 10
        assert len(resp1.records) == 3

        # Page 2: offset=3, limit=3 — index returns 6 rows (limit+offset), caller slices
        req2 = SearchRequest(limit=3, offset=3)
        resp2 = idx.search(req2)
        assert resp2.total == 10
        # The index returns 6 records (limit+offset=6) for the merge layer to slice
        assert len(resp2.records) == 6

        # Response carries pagination metadata
        assert resp2.limit == 3
        assert resp2.offset == 3
