"""Tests for pagination (limit / offset) fields on SearchRequest and SearchResponse."""
from __future__ import annotations

import sys
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


