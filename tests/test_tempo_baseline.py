"""Tests for TEMPO baseline TCR repertoire model."""
from __future__ import annotations
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).parent.parent / "servers"))

import pytest
import numpy as np
from tempo.baseline import BaselineModel, load_baseline


class TestBaselineModel:
    def _make_mini_baseline(self) -> BaselineModel:
        """Create a tiny baseline for testing with 2 V genes, 2 J genes, CDR3 len 10-11."""
        v_freq = {"TRAV12-2": 0.15, "TRAV27": 0.05}
        j_freq = {"TRAJ30": 0.08, "TRAJ34": 0.04}
        v_std = {"TRAV12-2": 0.03, "TRAV27": 0.01}
        j_std = {"TRAJ30": 0.02, "TRAJ34": 0.01}
        length_dist = {
            ("TRAV12-2", "TRAJ30"): {10: 0.6, 11: 0.4},
            ("TRAV12-2", "TRAJ34"): {10: 0.5, 11: 0.5},
            ("TRAV27", "TRAJ30"): {10: 0.4, 11: 0.6},
            ("TRAV27", "TRAJ34"): {10: 0.3, 11: 0.7},
        }
        cdr3_freq = {}
        for vj_key in length_dist:
            for length in length_dist[vj_key]:
                cdr3_freq[(*vj_key, length)] = np.full((20, length), 1.0 / 20)

        return BaselineModel(
            species="human",
            chain="alpha",
            v_freq=v_freq,
            j_freq=j_freq,
            v_std=v_std,
            j_std=j_std,
            length_dist=length_dist,
            cdr3_freq=cdr3_freq,
        )

    def test_baseline_model_has_required_fields(self):
        bl = self._make_mini_baseline()
        assert bl.species == "human"
        assert bl.chain == "alpha"
        assert "TRAV12-2" in bl.v_freq
        assert "TRAJ30" in bl.j_freq

    def test_q_v_returns_frequency(self):
        bl = self._make_mini_baseline()
        assert bl.q_v("TRAV12-2") == pytest.approx(0.15)
        assert bl.q_v("UNKNOWN_GENE") == pytest.approx(0.0)

    def test_q_j_returns_frequency(self):
        bl = self._make_mini_baseline()
        assert bl.q_j("TRAJ30") == pytest.approx(0.08)

    def test_q_length_given_vj(self):
        bl = self._make_mini_baseline()
        assert bl.q_length("TRAV12-2", "TRAJ30", 10) == pytest.approx(0.6)
        assert bl.q_length("TRAV12-2", "TRAJ30", 99) == pytest.approx(0.0)

    def test_q_cdr3_given_vjl_returns_matrix(self):
        bl = self._make_mini_baseline()
        mat = bl.q_cdr3("TRAV12-2", "TRAJ30", 10)
        assert mat.shape == (20, 10)
        assert mat[0, 0] == pytest.approx(1.0 / 20)

    def test_q_cdr3_unknown_vj_returns_uniform(self):
        bl = self._make_mini_baseline()
        mat = bl.q_cdr3("UNKNOWN", "UNKNOWN", 10)
        assert mat.shape == (20, 10)
        assert np.allclose(mat, 1.0 / 20)

    def test_v_zscore(self):
        bl = self._make_mini_baseline()
        z = bl.v_zscore("TRAV12-2", observed_freq=0.50)
        assert z > 10.0

    def test_j_zscore(self):
        bl = self._make_mini_baseline()
        z = bl.j_zscore("TRAJ30", observed_freq=0.08)
        assert z == pytest.approx(0.0)


class TestLoadBaseline:
    def test_load_baseline_raises_for_missing_species(self):
        with pytest.raises(FileNotFoundError):
            load_baseline("martian", "alpha")

    def test_load_baseline_raises_for_missing_chain(self):
        with pytest.raises(FileNotFoundError):
            load_baseline("human", "gamma")
