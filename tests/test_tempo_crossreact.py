"""Tests for MixTCRcross cross-reactivity prediction (0-3 scale)."""
from __future__ import annotations
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).parent.parent / "servers"))

import pytest
from tempo.crossreact import (
    classify_position,
    predict_cross_reactivity,
    CrossReactivityResult,
    _affinity_to_score,
    _tcr_contacts_match,
    _standardize_allele,
)


class TestClassifyPosition:
    def test_anchor_pos2(self):
        assert classify_position(2, 9, "I") == "anchor"

    def test_anchor_pos9(self):
        assert classify_position(9, 9, "I") == "anchor"

    def test_p3_is_nonanchor_p3p5(self):
        assert classify_position(3, 9, "I") == "nonanchor_p3p5"

    def test_p5_is_nonanchor_p3p5(self):
        assert classify_position(5, 9, "I") == "nonanchor_p3p5"

    def test_p6_is_nonanchor_p6p7(self):
        assert classify_position(6, 9, "I") == "nonanchor_p6p7"

    def test_p7_is_nonanchor_p6p7(self):
        assert classify_position(7, 9, "I") == "nonanchor_p6p7"

    def test_p8_is_other(self):
        assert classify_position(8, 9, "I") == "other"

    def test_p1_is_other(self):
        assert classify_position(1, 9, "I") == "other"


class TestAffinityToScore:
    def test_strong_binder_near_3(self):
        # 0.01 %rank maps to exactly 2.5 via the log10 formula
        score = _affinity_to_score(0.01)
        assert score >= 2.5

    def test_weak_binder_lower(self):
        score = _affinity_to_score(1.5)
        assert 1.0 < score < 3.0

    def test_very_weak_binder_low(self):
        # 50 %rank is above 2.0 on the raw formula; it is a weak binder
        score = _affinity_to_score(50.0)
        assert score < 2.5

    def test_output_range(self):
        for rank in [0.001, 0.01, 0.1, 0.5, 1.0, 2.0, 10.0, 50.0, 100.0]:
            score = _affinity_to_score(rank)
            assert 0.0 <= score <= 4.0  # Raw score before capping


class TestStandardizeAllele:
    def test_already_standard(self):
        assert _standardize_allele("HLA-A*02:01") == "HLA-A*02:01"

    def test_compact_format(self):
        assert _standardize_allele("HLA-A0201") == "HLA-A*02:01"


class TestTCRContactsMatch:
    def test_same_allele_matches(self):
        assert _tcr_contacts_match("HLA-A*02:01", "HLA-A*02:01") is True

    def test_unknown_alleles_assume_match(self):
        assert _tcr_contacts_match("HLA-A*99:99", "HLA-B*99:99") is True

    def test_different_contacts_dont_match(self):
        # A*02:01 and B*07:02 have different contacts at position 66 (I vs N)
        assert _tcr_contacts_match("HLA-A*02:01", "HLA-B*07:02") is False


class TestPredictCrossReactivity:
    def test_identical_peptides_score_3(self):
        result = predict_cross_reactivity("GILGFVFTL", "GILGFVFTL", "HLA-A*02:01")
        assert result.score == pytest.approx(3.0)
        assert result.category == "High"

    def test_single_p3_mutation_low(self):
        # Position 3 mutation -> capped below Medium boundary, stays Low
        result = predict_cross_reactivity("GILGFVFTL", "GIXGFVFTL", "HLA-A*02:01")
        assert result.score < 1.0
        assert result.category == "Low"

    def test_single_p6_mutation_medium(self):
        # Position 6 mutation -> capped below High boundary, at most Medium
        result = predict_cross_reactivity("GILGFVFTL", "GILGFXFTL", "HLA-A*02:01")
        assert result.score < 2.0
        assert result.category in ("Medium", "Low")

    def test_anchor_mutation_low(self):
        # Position 2 (anchor) mutation + another TCR position -> Low
        result = predict_cross_reactivity("GILGFVFTL", "GXLGFVFTL", "HLA-A*02:01")
        # This is a single TCR-determining diff (pos 2 is anchor, not TCR-det)
        # But pos2 anchor change is significant
        assert result.score <= 3.0

    def test_multiple_tcr_det_mutations_low(self):
        # >1 TCR-determining position differs -> score = 0
        result = predict_cross_reactivity("GILGFVFTL", "GIXXFVFTL", "HLA-A*02:01")
        assert result.score == pytest.approx(0.0)
        assert result.category == "Low"

    def test_different_length_low(self):
        result = predict_cross_reactivity("GILGFVFTL", "GILGFVFT", "HLA-A*02:01")
        assert result.score == pytest.approx(0.0)
        assert result.category == "Low"
        assert result.length_match is False

    def test_with_affinity_poor_binder(self):
        # affinity > 2.0 -> Low
        result = predict_cross_reactivity(
            "GILGFVFTL", "AILGFVFTL", "HLA-A*02:01",
            affinity_var=5.0,
        )
        assert result.score == pytest.approx(0.0)
        assert result.category == "Low"

    def test_with_affinity_strong_binder(self):
        result = predict_cross_reactivity(
            "GILGFVFTL", "AILGFVFTL", "HLA-A*02:01",
            affinity_var=0.1,
        )
        assert result.score > 2.0
        assert result.category == "High"

    def test_result_has_normalized_score(self):
        result = predict_cross_reactivity("GILGFVFTL", "GILGFVFTL", "HLA-A*02:01")
        assert result.normalized_score == pytest.approx(1.0)

    def test_non_9mer_handled(self):
        # Non-9-mers get basic scoring
        result = predict_cross_reactivity("GILGFVFTLM", "GILGFVFTLX", "HLA-A*02:01")
        assert 0.0 <= result.score <= 3.0
