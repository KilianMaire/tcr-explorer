"""Tests for TCR similarity scoring module (tcrdist3)."""
from __future__ import annotations
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).parent.parent / "servers"))

import pytest
import numpy as np
import pandas as pd
from unittest.mock import patch, MagicMock
from batman.tcrdist import (
    normalize_tcrdist,
    build_reference_df,
    TCRdistScorer,
)


class TestNormalizeTcrdist:
    def test_zero_distance_returns_one(self):
        assert normalize_tcrdist(0.0) == pytest.approx(1.0)

    def test_max_distance_returns_zero(self):
        assert normalize_tcrdist(96.0) == pytest.approx(0.0)

    def test_half_returns_half(self):
        assert normalize_tcrdist(48.0) == pytest.approx(0.5)

    def test_over_max_clipped_to_zero(self):
        assert normalize_tcrdist(200.0) == pytest.approx(0.0)

    def test_negative_clipped_to_one(self):
        assert normalize_tcrdist(-5.0) == pytest.approx(1.0)

    def test_custom_max(self):
        assert normalize_tcrdist(50.0, max_dist=100.0) == pytest.approx(0.5)


class TestBuildReferenceDF:
    def test_valid_records_produce_dataframe(self):
        records = [
            {
                "sequence": "CASSIRSSYEQYF",
                "metadata": {"v_segm": "TRBV19", "j_segm": "TRBJ2-1",
                              "antigen_epitope": "GILGFVFTL"},
            },
            {
                "sequence": "CASSLDRGQETQYF",
                "metadata": {"v_segm": "TRBV12-3", "j_segm": "TRBJ2-5",
                              "antigen_epitope": "NLVPMVATV"},
            },
        ]
        df = build_reference_df(records)
        assert len(df) == 2
        assert "cdr3_b_aa" in df.columns
        assert "v_b_gene" in df.columns
        assert "j_b_gene" in df.columns
        assert "epitope" in df.columns

    def test_missing_v_gene_excluded(self):
        records = [
            {
                "sequence": "CASSIRSSYEQYF",
                "metadata": {"j_segm": "TRBJ2-1", "antigen_epitope": "GILGFVFTL"},
                # v_segm missing
            },
        ]
        df = build_reference_df(records)
        assert len(df) == 0

    def test_v_gene_normalized_to_tcrdist_format(self):
        """TRBV19 → TRBV19*01 (tcrdist3 requires *01 allele suffix)."""
        records = [
            {
                "sequence": "CASSIRSSYEQYF",
                "metadata": {"v_segm": "TRBV19", "j_segm": "TRBJ2-1",
                              "antigen_epitope": "GILGFVFTL"},
            },
        ]
        df = build_reference_df(records)
        assert df["v_b_gene"].iloc[0].endswith("*01")

    def test_empty_records_returns_empty_df(self):
        df = build_reference_df([])
        assert len(df) == 0


class TestTCRdistScorer:
    def _make_reference_df(self) -> pd.DataFrame:
        return pd.DataFrame({
            "cdr3_b_aa": ["CASSIRSSYEQYF", "CASSLDRGQETQYF"],
            "v_b_gene": ["TRBV19*01", "TRBV12-3*01"],
            "j_b_gene": ["TRBJ2-1*01", "TRBJ2-5*01"],
            "epitope": ["GILGFVFTL", "NLVPMVATV"],
            "count": [1, 1],
        })

    @patch("batman.tcrdist.TCRrep")
    def test_score_returns_float_in_0_1(self, mock_tcrep_cls):
        """score() returns a float in [0, 1]."""
        mock_tr = MagicMock()
        mock_tr.pw_beta = np.array([[0, 24], [24, 0]])
        mock_tr.clone_df = pd.DataFrame({
            "cdr3_b_aa": ["CASSIRSSYEQYF", "CASSLDRGQETQYF"],
        })
        mock_tcrep_cls.return_value = mock_tr

        scorer = TCRdistScorer()
        scorer._ref_df = self._make_reference_df()

        result = scorer.score(
            cdr3_b="CASSIRSSYEQYF",
            v_gene="TRBV19*01",
            j_gene="TRBJ2-1*01",
            epitope="GILGFVFTL",
        )
        assert result is not None
        assert 0.0 <= result <= 1.0

    @patch("batman.tcrdist.TCRrep")
    def test_score_identical_cdr3_returns_high(self, mock_tcrep_cls):
        """Query CDR3 identical to reference → distance 0 → score 1."""
        mock_tr = MagicMock()
        mock_tr.pw_beta = np.array([[0, 0], [0, 0]])
        mock_tr.clone_df = pd.DataFrame({"cdr3_b_aa": ["CASSIRSSYEQYF", "Q"]})
        mock_tcrep_cls.return_value = mock_tr

        scorer = TCRdistScorer()
        scorer._ref_df = self._make_reference_df()

        result = scorer.score("CASSIRSSYEQYF", "TRBV19*01", "TRBJ2-1*01", "GILGFVFTL")
        assert result == pytest.approx(1.0)

    def test_score_no_reference_for_epitope_returns_none(self):
        """No reference binders for the requested epitope → return None."""
        scorer = TCRdistScorer()
        scorer._ref_df = self._make_reference_df()
        result = scorer.score("CASSIRSSYEQYF", "TRBV19*01", "TRBJ2-1*01",
                              epitope="UNKNOWN_EPITOPE")
        assert result is None

    def test_score_empty_reference_returns_none(self):
        """Empty reference DB → return None."""
        scorer = TCRdistScorer()
        scorer._ref_df = pd.DataFrame(
            columns=["cdr3_b_aa", "v_b_gene", "j_b_gene", "epitope", "count"]
        )
        result = scorer.score("CASSIRSSYEQYF", "TRBV19*01", "TRBJ2-1*01", "GILGFVFTL")
        assert result is None

    def test_score_returns_none_when_tcrdist_not_installed(self):
        """score() returns None gracefully when tcrdist3 is not installed."""
        import batman.tcrdist as tcrdist_module
        original = tcrdist_module.TCRrep
        try:
            tcrdist_module.TCRrep = None
            scorer = TCRdistScorer()
            scorer._ref_df = self._make_reference_df()
            result = scorer.score("CASSIRSSYEQYF", "TRBV19*01", "TRBJ2-1*01", "GILGFVFTL")
            assert result is None
        finally:
            tcrdist_module.TCRrep = original

    @patch("batman.tcrdist.TCRrep")
    def test_score_returns_none_on_tcrdist_exception(self, mock_tcrep_cls):
        """score() returns None when TCRrep raises an exception."""
        mock_tcrep_cls.side_effect = ValueError("invalid CDR3")
        scorer = TCRdistScorer()
        scorer._ref_df = self._make_reference_df()
        result = scorer.score("BADSEQ!!!", "TRBV19*01", "TRBJ2-1*01", "GILGFVFTL")
        assert result is None


class TestBuildReferenceDFAlpha:
    def test_alpha_chain_records(self):
        records = [
            {
                "sequence": "CAVGDDKIIF",
                "metadata": {
                    "v_segm": "TRAV12-2",
                    "j_segm": "TRAJ30",
                    "antigen_epitope": "LLWNGPMAV",
                },
            },
        ]
        df = build_reference_df(records, chain="alpha")
        assert len(df) == 1
        assert "cdr3_a_aa" in df.columns
        assert "v_a_gene" in df.columns
        assert "j_a_gene" in df.columns

    def test_alpha_from_cdr3_a_metadata(self):
        records = [
            {
                "sequence": "CASSIRSSYEQYF",  # beta CDR3
                "metadata": {
                    "v_segm": "TRBV19",
                    "j_segm": "TRBJ2-7",
                    "antigen_epitope": "GILGFVFTL",
                    "cdr3_a": "CAGGGSQGNLIF",
                    "v_a_segm": "TRAV27",
                    "j_a_segm": "TRAJ42",
                },
            },
        ]
        df = build_reference_df(records, chain="alpha")
        assert len(df) == 1
        assert df["cdr3_a_aa"].iloc[0] == "CAGGGSQGNLIF"


class TestTCRdistScorerAlphaChain:
    @patch("batman.tcrdist.TCRrep")
    def test_score_alpha_chain(self, mock_tcrep_cls):
        mock_tr = MagicMock()
        mock_tr.pw_alpha = np.array([[0, 10], [10, 0]])
        mock_tcrep_cls.return_value = mock_tr

        scorer = TCRdistScorer()
        ref = [
            {
                "sequence": "CAVGDDKIIF",
                "metadata": {
                    "v_segm": "TRAV12-2",
                    "j_segm": "TRAJ30",
                    "antigen_epitope": "LLWNGPMAV",
                },
                "antigen_epitope": "LLWNGPMAV",
            },
        ]
        scorer.load_reference(ref, chain="alpha")
        score = scorer.score(
            cdr3="CAVGDDKIIF",
            v_gene="TRAV12-2",
            j_gene="TRAJ30",
            epitope="LLWNGPMAV",
            chain="alpha",
        )
        assert score is not None
        assert 0.0 <= score <= 1.0

    @patch("batman.tcrdist.TCRrep")
    def test_score_paired_alpha_beta(self, mock_tcrep_cls):
        # Mock alpha call
        mock_tr_alpha = MagicMock()
        mock_tr_alpha.pw_alpha = np.array([[0, 10], [10, 0]])
        # Mock beta call
        mock_tr_beta = MagicMock()
        mock_tr_beta.pw_beta = np.array([[0, 20], [20, 0]])

        mock_tcrep_cls.side_effect = [mock_tr_alpha, mock_tr_beta]

        scorer = TCRdistScorer()
        # Load beta reference
        beta_ref = [
            {
                "sequence": "CASSIRSSYEQYF",
                "metadata": {
                    "v_segm": "TRBV19",
                    "j_segm": "TRBJ2-7",
                    "antigen_epitope": "GILGFVFTL",
                },
                "antigen_epitope": "GILGFVFTL",
            },
        ]
        scorer.load_reference(beta_ref, chain="beta")
        # Load alpha reference
        alpha_ref = [
            {
                "sequence": "CAGGGSQGNLIF",
                "metadata": {
                    "v_segm": "TRAV27",
                    "j_segm": "TRAJ42",
                    "antigen_epitope": "GILGFVFTL",
                },
                "antigen_epitope": "GILGFVFTL",
            },
        ]
        scorer.load_reference(alpha_ref, chain="alpha")
        score = scorer.score_paired(
            cdr3_a="CAGGGSQGNLIF", v_a="TRAV27", j_a="TRAJ42",
            cdr3_b="CASSIRSSYEQYF", v_b="TRBV19", j_b="TRBJ2-7",
            epitope="GILGFVFTL",
        )
        assert score is not None
        assert 0.0 <= score <= 1.0


class TestTCRdistScorerMouse:
    @patch("batman.tcrdist.TCRrep")
    def test_score_mouse_beta(self, mock_tcrep_cls):
        mock_tr = MagicMock()
        mock_tr.pw_beta = np.array([[0, 15], [15, 0]])
        mock_tcrep_cls.return_value = mock_tr

        scorer = TCRdistScorer()
        ref = [
            {
                "sequence": "CASSQDWGAETLYF",
                "metadata": {
                    "v_segm": "TRBV13-1",
                    "j_segm": "TRBJ2-3",
                    "antigen_epitope": "SIINFEKL",
                },
                "antigen_epitope": "SIINFEKL",
            },
        ]
        scorer.load_reference(ref, chain="beta")
        score = scorer.score(
            cdr3="CASSQDWGAETLYF",
            v_gene="TRBV13-1",
            j_gene="TRBJ2-3",
            epitope="SIINFEKL",
            chain="beta",
            organism="mouse",
        )
        assert score is not None
        assert 0.0 <= score <= 1.0
        # Verify mouse organism was passed to TCRrep
        call_kwargs = mock_tcrep_cls.call_args
        assert call_kwargs[1]["organism"] == "mouse" or call_kwargs.kwargs["organism"] == "mouse"
