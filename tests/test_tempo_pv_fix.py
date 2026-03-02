"""Tests for corrected TEMPO P(V)/P(J) computation."""
import math
import pytest
from servers.tempo.scorer import _safe_log_ratio

def test_safe_log_ratio_both_nonzero():
    result = _safe_log_ratio(0.15, 0.05)
    expected = math.log(0.15 / 0.05)
    assert abs(result - expected) < 0.01

def test_safe_log_ratio_p_greater_q():
    result = _safe_log_ratio(0.2, 0.1)
    assert result > 0

def test_safe_log_ratio_p_less_q():
    result = _safe_log_ratio(0.05, 0.2)
    assert result < 0

def test_safe_log_ratio_zero_q():
    result = _safe_log_ratio(0.1, 0.0)
    assert math.isfinite(result)

def test_safe_log_ratio_zero_p():
    result = _safe_log_ratio(0.0, 0.1)
    assert math.isfinite(result)

def test_safe_log_ratio_equal():
    result = _safe_log_ratio(0.1, 0.1)
    assert abs(result) < 0.01
