"""Tests for TSP (TCR Specificity Profile) computation."""
from __future__ import annotations
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).parent.parent / "servers"))

import pytest
import numpy as np
from tempo.baseline import BaselineModel
from tempo.tsp import compute_tsp, TSProfile


def _mini_baseline() -> BaselineModel:
    return BaselineModel(
        species="human", chain="beta",
        v_freq={"TRBV28": 0.10, "TRBV19": 0.08, "TRBV27": 0.04},
        j_freq={"TRBJ2-7": 0.06, "TRBJ1-1": 0.05},
        v_std={"TRBV28": 0.02, "TRBV19": 0.02, "TRBV27": 0.01},
        j_std={"TRBJ2-7": 0.01, "TRBJ1-1": 0.01},
        length_dist={("TRBV28", "TRBJ2-7"): {13: 0.5, 14: 0.5}},
        cdr3_freq={
            ("TRBV28", "TRBJ2-7", 13): np.full((20, 13), 1.0 / 20),
        },
    )


class TestComputeTSP:
    def test_returns_tsprofile(self):
        bl = _mini_baseline()
        tcrs = [
            {"v_gene": "TRBV28", "j_gene": "TRBJ2-7", "cdr3": "CASTPQTAYEQYF"},
            {"v_gene": "TRBV28", "j_gene": "TRBJ2-7", "cdr3": "CASSIRSSYEQYF"},
            {"v_gene": "TRBV19", "j_gene": "TRBJ2-7", "cdr3": "CASSQRSSYEQYF"},
        ]
        profile = compute_tsp(tcrs, bl)
        assert isinstance(profile, TSProfile)

    def test_v_enrichment_has_all_observed_genes(self):
        bl = _mini_baseline()
        tcrs = [
            {"v_gene": "TRBV28", "j_gene": "TRBJ2-7", "cdr3": "CASTPQTAYEQYF"},
            {"v_gene": "TRBV28", "j_gene": "TRBJ2-7", "cdr3": "CASSIRSSYEQYF"},
        ]
        profile = compute_tsp(tcrs, bl)
        assert "TRBV28" in profile.v_enrichment
        assert profile.v_enrichment["TRBV28"]["observed_freq"] == pytest.approx(1.0)

    def test_cdr3_length_distribution(self):
        bl = _mini_baseline()
        tcrs = [
            {"v_gene": "TRBV28", "j_gene": "TRBJ2-7", "cdr3": "CASTPQTAYEQYF"},  # len 13
            {"v_gene": "TRBV28", "j_gene": "TRBJ2-7", "cdr3": "CASTPQTAYEQYF"},
        ]
        profile = compute_tsp(tcrs, bl)
        assert 13 in profile.cdr3_length_dist
        assert profile.cdr3_length_dist[13] == pytest.approx(1.0)

    def test_empty_tcr_list_returns_empty_profile(self):
        bl = _mini_baseline()
        profile = compute_tsp([], bl)
        assert len(profile.v_enrichment) == 0
        assert len(profile.j_enrichment) == 0
