"""Integration tests for batman_server FastAPI endpoints."""
from __future__ import annotations
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).parent.parent / "servers"))

import pytest
import numpy as np
from unittest.mock import patch, MagicMock
from fastapi.testclient import TestClient


@pytest.fixture
def client():
    # Patch pybatman before importing batman_server
    with patch("batman.scorer.train") as mock_train, \
         patch("batman.scorer.peptide2index") as mock_p2i:
        mock_train.return_value = (np.ones((1, 9)), np.eye(20))
        mock_p2i.return_value = np.array([2.0])
        from batman_server import app
        yield TestClient(app)


class TestHealthEndpoint:
    def test_health_returns_ok(self, client):
        r = client.get("/health")
        assert r.status_code == 200
        assert r.json()["status"] == "ok"
        assert r.json()["server"] == "batman"


class TestScoreEndpoint:
    def test_score_with_cached_model(self, client):
        """POST /score with a pre-trained model in cache returns composite score."""
        payload = {
            "tcr_id": "TCR1",
            "index_peptide": "GILGFVFTL",
            "candidate_peptide": "NLVPMVATV",
            "hla_allele": "HLA-A*02:01",
        }
        r = client.post("/score", json=payload)
        # First call → cache miss → batman_score = None (training async)
        # OR returns cached result. Either way, status 200.
        assert r.status_code == 200
        data = r.json()
        assert "batman_score" in data
        assert "cache_status" in data

    def test_score_missing_required_field_returns_422(self, client):
        r = client.post("/score", json={"tcr_id": "TCR1"})
        assert r.status_code == 422

    def test_score_response_has_all_score_fields(self, client):
        payload = {
            "tcr_id": "TCR1",
            "index_peptide": "GILGFVFTL",
            "candidate_peptide": "GILGFVFTL",
            "hla_allele": "HLA-A*02:01",
        }
        r = client.post("/score", json=payload)
        assert r.status_code == 200
        data = r.json()
        # All score fields present (may be null if cache miss)
        assert "batman_score" in data
        assert "tcr_id" in data
        assert "candidate_peptide" in data

    def test_score_after_train_returns_hit(self, client):
        """Train a model then score with it → cache hit with non-null score."""
        from batman_server import _cache
        # Train
        train_payload = {
            "tcr_id": "TCR_INTEGRATION",
            "index_peptide": "GILGFVFTL",
            "activation_data": [
                {"peptide": "GILGFVFTL", "activation": 2},
                {"peptide": "AILGFVFTL", "activation": 1},
                {"peptide": "XILGFVFTL", "activation": 0},
            ],
        }
        r = client.post("/train", json=train_payload)
        assert r.status_code == 200

        try:
            # Score with the trained model
            score_payload = {
                "tcr_id": "TCR_INTEGRATION",
                "index_peptide": "GILGFVFTL",
                "candidate_peptide": "AILGFVFTL",
            }
            r = client.post("/score", json=score_payload)
            assert r.status_code == 200
            data = r.json()
            assert data["cache_status"] == "hit"
            assert data["batman_score"] is not None
            assert 0.0 <= data["batman_score"] <= 1.0
        finally:
            _cache.clear()


class TestTrainEndpoint:
    def test_train_with_valid_data(self, client):
        payload = {
            "tcr_id": "TCR1",
            "index_peptide": "GILGFVFTL",
            "activation_data": [
                {"peptide": "GILGFVFTL", "activation": 2},
                {"peptide": "AILGFVFTL", "activation": 1},
                {"peptide": "XILGFVFTL", "activation": 0},
            ],
        }
        r = client.post("/train", json=payload)
        assert r.status_code == 200
        data = r.json()
        assert data["status"] in ("trained", "queued")
        assert data["n_samples"] == 3

    def test_train_missing_index_peptide_returns_422(self, client):
        r = client.post("/train", json={"tcr_id": "TCR1", "activation_data": []})
        assert r.status_code == 422

    def test_train_empty_activation_data_returns_422(self, client):
        """Empty activation_data should return 422 (Pydantic validation)."""
        payload = {
            "tcr_id": "TCR1",
            "index_peptide": "GILGFVFTL",
            "activation_data": [],
        }
        r = client.post("/train", json=payload)
        assert r.status_code == 422


class TestTCRdistInScore:
    @patch("batman.scorer.train")
    @patch("batman.scorer.peptide2index")
    @patch("batman.tcrdist.TCRrep")
    def test_score_includes_tcrdist_score_field(self, mock_tcrep, mock_p2i, mock_train):
        """ScoreResponse includes tcrdist_score field."""
        import numpy as np
        from batman.cache import CachedModel
        from batman_server import _cache, app
        from fastapi.testclient import TestClient

        # Pre-load a model in the cache
        _cache.put(CachedModel("TCR_TD", "GILGFVFTL", np.ones((1, 9)), np.eye(20)))
        try:
            mock_p2i.return_value = np.array([3.0])

            mock_tr = MagicMock()
            mock_tr.pw_beta = np.array([[0, 12], [12, 0]])
            mock_tcrep.return_value = mock_tr

            client = TestClient(app)
            payload = {
                "tcr_id": "TCR_TD",
                "index_peptide": "GILGFVFTL",
                "candidate_peptide": "AILGFVFTL",
                "hla_allele": "HLA-A*02:01",
                "cdr3_b": "CASSIRSSYEQYF",
                "v_gene": "TRBV19",
                "j_gene": "TRBJ2-1",
            }
            r = client.post("/score", json=payload)
            assert r.status_code == 200
            data = r.json()
            assert "tcrdist_score" in data
        finally:
            _cache.clear()

    def test_score_without_cdr3_tcrdist_is_none(self):
        """When cdr3_b is omitted, tcrdist_score is null."""
        from batman_server import app
        from fastapi.testclient import TestClient
        client = TestClient(app)
        payload = {
            "tcr_id": "UNKNOWN_TCR",
            "index_peptide": "GILGFVFTL",
            "candidate_peptide": "AILGFVFTL",
            # No cdr3_b, v_gene, j_gene
        }
        r = client.post("/score", json=payload)
        assert r.status_code == 200
        data = r.json()
        assert "tcrdist_score" in data
        assert data["tcrdist_score"] is None
