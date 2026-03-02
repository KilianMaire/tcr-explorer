"""Tests for BATMAN scorer: train() wrapper and distance→score normalization."""
from __future__ import annotations
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).parent.parent / "servers"))

import pytest
import numpy as np
import pandas as pd
from unittest.mock import patch, MagicMock
from batman.scorer import BATMANScorer, TrainingResult, normalize_distance


class TestNormalizeDistance:
    def test_zero_distance_returns_one(self):
        assert normalize_distance(0.0, max_dist=10.0) == pytest.approx(1.0)

    def test_max_distance_returns_zero(self):
        assert normalize_distance(10.0, max_dist=10.0) == pytest.approx(0.0)

    def test_half_distance_returns_half(self):
        assert normalize_distance(5.0, max_dist=10.0) == pytest.approx(0.5)

    def test_over_max_clipped_to_zero(self):
        assert normalize_distance(15.0, max_dist=10.0) == pytest.approx(0.0)

    def test_negative_clipped_to_one(self):
        assert normalize_distance(-1.0, max_dist=10.0) == pytest.approx(1.0)


class TestBATMANScorer:
    def _make_training_df(self) -> pd.DataFrame:
        return pd.DataFrame({
            "tcr": ["TCR1"] * 3,
            "index": ["GILGFVFTL"] * 3,
            "peptide": ["GILGFVFTL", "AILGFVFTL", "XILGFVFTL"],
            "activation": [2, 1, 0],
        })

    @patch("batman.scorer.train")
    def test_train_returns_training_result(self, mock_train):
        mock_weights = np.ones((1, 9))
        mock_matrix = np.eye(20)
        mock_train.return_value = (mock_weights, mock_matrix)

        scorer = BATMANScorer()
        df = self._make_training_df()
        result = scorer.train(df, steps=100, seed=42)

        assert isinstance(result, TrainingResult)
        assert result.tcr_id == "TCR1"
        assert result.index_peptide == "GILGFVFTL"
        np.testing.assert_array_equal(result.weights, mock_weights)
        mock_train.assert_called_once()

    @patch("batman.scorer.peptide2index")
    def test_score_returns_float_in_0_1(self, mock_p2i):
        mock_p2i.return_value = np.array([2.5])

        scorer = BATMANScorer()
        score = scorer.score(
            index_peptide="GILGFVFTL",
            candidate_peptide="NLVPMVATV",
            weights=np.ones((1, 9)),
            aa_matrix=np.eye(20),
        )
        assert 0.0 <= score <= 1.0

    @patch("batman.scorer.peptide2index")
    def test_score_same_peptide_is_high(self, mock_p2i):
        mock_p2i.return_value = np.array([0.0])

        scorer = BATMANScorer()
        score = scorer.score(
            index_peptide="GILGFVFTL",
            candidate_peptide="GILGFVFTL",
            weights=np.ones((1, 9)),
            aa_matrix=np.eye(20),
        )
        assert score == pytest.approx(1.0)

    @patch("batman.scorer.train")
    def test_train_raises_on_missing_columns(self, mock_train):
        scorer = BATMANScorer()
        bad_df = pd.DataFrame({"tcr": ["TCR1"], "peptide": ["GILG"]})
        with pytest.raises(ValueError, match="Missing columns"):
            scorer.train(bad_df)
        mock_train.assert_not_called()

    @patch("batman.scorer.train")
    def test_train_raises_on_non_consecutive_activation(self, mock_train):
        scorer = BATMANScorer()
        df = pd.DataFrame({
            "tcr": ["TCR1", "TCR1"],
            "index": ["GILGFVFTL", "GILGFVFTL"],
            "peptide": ["GILGFVFTL", "AILGFVFTL"],
            "activation": [0, 2],  # missing 1
        })
        with pytest.raises(ValueError, match="consecutive"):
            scorer.train(df)
        mock_train.assert_not_called()
