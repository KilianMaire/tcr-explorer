"""Tests for batman enrichment helper in api.py."""
from __future__ import annotations
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).parent.parent / "src"))

import asyncio
import pytest
from unittest.mock import patch, AsyncMock, MagicMock
from imgt_app.models import GeneRecord, SearchRequest


def _run(coro):
    loop = asyncio.new_event_loop()
    try:
        return loop.run_until_complete(coro)
    finally:
        loop.close()


class TestCompositeScore:
    def test_all_scores_geometric_mean(self):
        from imgt_app.api import _composite_score
        result = _composite_score(0.8, 0.5, 0.4)
        assert result == pytest.approx((0.8 * 0.5 * 0.4) ** (1/3), rel=1e-3)

    def test_two_scores(self):
        from imgt_app.api import _composite_score
        result = _composite_score(0.8, 0.5, None)
        assert result == pytest.approx((0.8 * 0.5) ** 0.5, rel=1e-3)

    def test_one_score(self):
        from imgt_app.api import _composite_score
        assert _composite_score(0.9, None, None) == pytest.approx(0.9)

    def test_all_none_returns_none(self):
        from imgt_app.api import _composite_score
        assert _composite_score(None, None, None) is None


class TestEnrichWithBatman:
    def _make_vdjdb_record(self) -> GeneRecord:
        return GeneRecord(
            source="vdjdb",
            gene_name="TRBV19",
            sequence="CASSIRSSYEQYF",
            antigen_epitope="GILGFVFTL",
            metadata={
                "v_segm": "TRBV19", "j_segm": "TRBJ2-1",
                "mhc_a": "HLA-A*02:01", "mhc_class": "MHCI",
            },
        )

    @patch("imgt_app.api.batman_client")
    def test_enrichment_attaches_scores_to_records(self, mock_client):
        mock_client.post_score = AsyncMock(return_value={
            "batman_score": 0.75,
            "pmhc_score": 0.88,
            "tcrdist_score": 0.60,
            "cache_status": "hit",
        })
        from imgt_app.api import _enrich_with_batman
        records = [self._make_vdjdb_record()]
        req = SearchRequest(source="vdjdb", antigen_epitope="GILGFVFTL")
        enriched = _run(_enrich_with_batman(records, req))

        assert len(enriched) == 1
        r = enriched[0]
        assert r.batman_score == pytest.approx(0.75)
        assert r.pmhc_score == pytest.approx(0.88)
        assert r.composite_score is not None
        assert 0.0 <= r.composite_score <= 1.0

    @patch("imgt_app.api.batman_client")
    def test_enrichment_skips_records_without_epitope(self, mock_client):
        mock_client.post_score = AsyncMock(return_value={})
        from imgt_app.api import _enrich_with_batman
        record = GeneRecord(source="vdjdb", gene_name="TRBV19",
                            sequence="CASSIRSSYEQYF")
        enriched = _run(_enrich_with_batman([record], SearchRequest()))
        assert enriched[0].batman_score is None
        mock_client.post_score.assert_not_called()

    @patch("imgt_app.api.batman_client")
    def test_enrichment_graceful_on_batman_server_down(self, mock_client):
        import httpx
        mock_client.post_score = AsyncMock(side_effect=httpx.HTTPError("connection refused"))
        from imgt_app.api import _enrich_with_batman
        records = [self._make_vdjdb_record()]
        enriched = _run(_enrich_with_batman(records, SearchRequest()))
        assert len(enriched) == 1
        assert enriched[0].batman_score is None  # not raised

    @patch("imgt_app.api.batman_client")
    def test_composite_score_computed_from_available_scores(self, mock_client):
        mock_client.post_score = AsyncMock(return_value={
            "batman_score": 0.8,
            "pmhc_score": 0.5,
            "tcrdist_score": None,
            "cache_status": "hit",
        })
        from imgt_app.api import _enrich_with_batman, _composite_score
        assert _composite_score(0.8, 0.5, None) == pytest.approx((0.8 * 0.5) ** 0.5, rel=1e-3)
        assert _composite_score(None, None, None) is None
        assert _composite_score(0.9, None, None) == pytest.approx(0.9)
