"""Tests for TEMPO binary wrapper."""
from __future__ import annotations
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).parent.parent / "servers"))

import pytest
from tempo.tempo_binary import (
    is_available,
    perc_rank_to_score,
    TempoPrediction,
)


class TestPercRankToScore:
    def test_rank_0_returns_1(self):
        assert perc_rank_to_score(0.0) == pytest.approx(1.0)

    def test_rank_100_returns_0(self):
        assert perc_rank_to_score(100.0) == pytest.approx(0.0)

    def test_rank_50_returns_half(self):
        assert perc_rank_to_score(50.0) == pytest.approx(0.5)

    def test_rank_over_100_clipped(self):
        assert perc_rank_to_score(150.0) == pytest.approx(0.0)

    def test_negative_rank_clipped(self):
        assert perc_rank_to_score(-10.0) == pytest.approx(1.0)


class TestIsAvailable:
    def test_returns_bool(self):
        result = is_available()
        assert isinstance(result, bool)


@pytest.mark.skipif(not is_available(), reason="TEMPO binary not installed")
class TestBinaryPredict:
    def test_predict_gilgfvftl(self):
        from tempo.tempo_binary import predict
        tcrs = [
            {
                "TRAV": "TRAV13-1", "TRAJ": "TRAJ28",
                "cdr3_TRA": "CAASITYSGAGSYQLTF",
                "TRBV": "TRBV11-2", "TRBJ": "TRBJ2-6",
                "cdr3_TRB": "CASSLRGEAFSGANVLTF",
            },
        ]
        results = predict(tcrs, "A0201_GILGFVFTL")
        assert len(results) == 1
        assert results[0].perc_rank is not None
        assert 0.0 <= results[0].perc_rank <= 100.0

    def test_predict_beta_only(self):
        from tempo.tempo_binary import predict
        tcrs = [
            {"TRBV": "TRBV19", "TRBJ": "TRBJ2-7", "cdr3_TRB": "CASSIRSSYEQYF"},
        ]
        results = predict(tcrs, "A0201_GILGFVFTL", chain="B")
        assert len(results) == 1
        assert results[0].perc_rank is not None

    def test_predict_invalid_epitope(self):
        from tempo.tempo_binary import predict
        tcrs = [{"TRBV": "TRBV19", "TRBJ": "TRBJ2-7", "cdr3_TRB": "CASSIRSSYEQYF"}]
        results = predict(tcrs, "NONEXISTENT_EPITOPE", chain="B")
        # Should return results with problem/error, not crash
        assert len(results) >= 1

    def test_predict_empty_list(self):
        from tempo.tempo_binary import predict
        results = predict([], "A0201_GILGFVFTL")
        assert results == []
