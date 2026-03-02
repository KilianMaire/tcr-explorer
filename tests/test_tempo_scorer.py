"""Tests for TEMPO log-likelihood scorer."""
from __future__ import annotations
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).parent.parent / "servers"))

import pytest
import numpy as np
from tempo.baseline import BaselineModel
from tempo.scorer import TempoScorer, TempoResult


def _mini_baseline() -> BaselineModel:
    v_freq = {"TRAV12-2": 0.15, "TRAV27": 0.05}
    j_freq = {"TRAJ30": 0.08, "TRAJ34": 0.04}
    v_std = {"TRAV12-2": 0.03, "TRAV27": 0.01}
    j_std = {"TRAJ30": 0.02, "TRAJ34": 0.01}
    length_dist = {
        ("TRAV12-2", "TRAJ30"): {10: 0.6, 11: 0.4},
    }
    cdr3_freq = {
        ("TRAV12-2", "TRAJ30", 10): np.full((20, 10), 1.0 / 20),
    }
    return BaselineModel(
        species="human", chain="alpha",
        v_freq=v_freq, j_freq=j_freq, v_std=v_std, j_std=j_std,
        length_dist=length_dist, cdr3_freq=cdr3_freq,
    )


class TestTempoScorer:
    def test_score_single_chain_returns_result(self):
        bl = _mini_baseline()
        scorer = TempoScorer(alpha_baseline=bl)
        result = scorer.score_single_chain(
            v_gene="TRAV12-2",
            j_gene="TRAJ30",
            cdr3="CAVGDDKIIF",  # length 10
            chain="alpha",
        )
        assert isinstance(result, TempoResult)
        assert isinstance(result.log_likelihood, float)

    def test_known_v_gene_higher_score_than_unknown(self):
        bl = _mini_baseline()
        scorer = TempoScorer(alpha_baseline=bl)
        known = scorer.score_single_chain("TRAV12-2", "TRAJ30", "CAVGDDKIIF", "alpha")
        unknown = scorer.score_single_chain("UNKNOWN_V", "TRAJ30", "CAVGDDKIIF", "alpha")
        assert known.log_likelihood > unknown.log_likelihood

    def test_score_paired_combines_chains(self):
        bl_a = _mini_baseline()
        bl_b = BaselineModel(
            species="human", chain="beta",
            v_freq={"TRBV28": 0.10},
            j_freq={"TRBJ2-7": 0.06},
            v_std={"TRBV28": 0.02},
            j_std={"TRBJ2-7": 0.01},
            length_dist={("TRBV28", "TRBJ2-7"): {13: 0.5}},
            cdr3_freq={("TRBV28", "TRBJ2-7", 13): np.full((20, 13), 1.0 / 20)},
        )
        scorer = TempoScorer(alpha_baseline=bl_a, beta_baseline=bl_b)
        result = scorer.score_paired(
            v_a="TRAV12-2", j_a="TRAJ30", cdr3_a="CAVGDDKIIF",
            v_b="TRBV28", j_b="TRBJ2-7", cdr3_b="CASTPQTAYEQYF",
        )
        assert isinstance(result.log_likelihood, float)

    def test_pseudocount_prevents_negative_infinity(self):
        bl = _mini_baseline()
        scorer = TempoScorer(alpha_baseline=bl)
        result = scorer.score_single_chain("TRAV12-2", "TRAJ30", "CWWWWWWWWF", "alpha")
        assert np.isfinite(result.log_likelihood)

    def test_cdr3_too_short_handled(self):
        bl = _mini_baseline()
        scorer = TempoScorer(alpha_baseline=bl)
        result = scorer.score_single_chain("TRAV12-2", "TRAJ30", "CDF", "alpha")
        assert np.isfinite(result.log_likelihood)
