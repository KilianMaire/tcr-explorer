"""Tests for batman enrichment helper in api.py."""
from __future__ import annotations
import sys
import pathlib
sys.path.insert(0, str(pathlib.Path(__file__).parent.parent / "src"))

import asyncio
import pytest


def _run(coro):
    loop = asyncio.new_event_loop()
    try:
        return loop.run_until_complete(coro)
    finally:
        loop.close()


class TestCompositeScore:
    def test_all_scores_geometric_mean(self):
        from tcr_explorer.api import _composite_score
        result = _composite_score(0.8, 0.5, 0.4)
        assert result == pytest.approx((0.8 * 0.5 * 0.4) ** (1/3), rel=1e-3)

    def test_two_scores(self):
        from tcr_explorer.api import _composite_score
        result = _composite_score(0.8, 0.5, None)
        assert result == pytest.approx((0.8 * 0.5) ** 0.5, rel=1e-3)

    def test_one_score(self):
        from tcr_explorer.api import _composite_score
        assert _composite_score(0.9, None, None) == pytest.approx(0.9)

    def test_all_none_returns_none(self):
        from tcr_explorer.api import _composite_score
        assert _composite_score(None, None, None) is None


