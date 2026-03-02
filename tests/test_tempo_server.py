"""Integration tests for tempo_server FastAPI endpoints."""
from __future__ import annotations
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).parent.parent / "servers"))

import pytest
import numpy as np
from unittest.mock import patch, MagicMock
from fastapi.testclient import TestClient


@pytest.fixture
def client():
    from tempo_server import app
    return TestClient(app)


class TestHealthEndpoint:
    def test_health_returns_ok(self, client):
        r = client.get("/health")
        assert r.status_code == 200
        assert r.json()["status"] == "ok"
        assert r.json()["server"] == "tempo"


class TestTempoScoreEndpoint:
    def test_score_single_chain(self, client):
        payload = {
            "v_gene": "TRAV12-2",
            "j_gene": "TRAJ30",
            "cdr3": "CAVGDDKIIF",
            "chain": "alpha",
            "species": "human",
        }
        r = client.post("/tempo/score", json=payload)
        assert r.status_code == 200
        data = r.json()
        assert "log_likelihood" in data

    def test_score_paired(self, client):
        payload = {
            "v_a": "TRAV12-2", "j_a": "TRAJ30", "cdr3_a": "CAVGDDKIIF",
            "v_b": "TRBV28", "j_b": "TRBJ2-7", "cdr3_b": "CASTPQTAYEQYF",
            "species": "human",
        }
        r = client.post("/tempo/score_paired", json=payload)
        assert r.status_code == 200
        data = r.json()
        assert "log_likelihood" in data

    def test_score_missing_required_field(self, client):
        r = client.post("/tempo/score", json={"v_gene": "TRAV12-2"})
        assert r.status_code == 422


class TestTSPEndpoint:
    def test_tsp_profile(self, client):
        payload = {
            "tcrs": [
                {"v_gene": "TRBV28", "j_gene": "TRBJ2-7", "cdr3": "CASTPQTAYEQYF"},
                {"v_gene": "TRBV28", "j_gene": "TRBJ2-7", "cdr3": "CASSIRSSYEQYF"},
            ],
            "chain": "beta",
            "species": "human",
        }
        r = client.post("/tsp/profile", json=payload)
        assert r.status_code == 200
        data = r.json()
        assert "v_enrichment" in data
        assert "j_enrichment" in data
        assert "cdr3_length_dist" in data

    def test_tsp_empty_tcrs(self, client):
        payload = {"tcrs": [], "chain": "beta", "species": "human"}
        r = client.post("/tsp/profile", json=payload)
        assert r.status_code == 200
        assert r.json()["n_tcrs"] == 0


class TestCrossReactEndpoint:
    def test_crossreact_predict(self, client):
        payload = {
            "reference_peptide": "LLWNGPMAV",
            "variant_peptides": ["ALWNGPMAV", "LLWAGPMAV"],
            "mhc_allele": "HLA-A*02:01",
            "mhc_class": "I",
        }
        r = client.post("/crossreact/predict", json=payload)
        assert r.status_code == 200
        data = r.json()
        assert "results" in data
        assert len(data["results"]) == 2
        for result in data["results"]:
            assert "score" in result
            assert 0.0 <= result["score"] <= 3.0
            assert "category" in result
            assert result["category"] in ("High", "Medium", "Low")
            assert "normalized_score" in result


class TestBatcaveEndpoint:
    def test_batcave_variants(self, client):
        r = client.get("/batcave/variants", params={"reference_peptide": "GILGFVFTL"})
        assert r.status_code == 200
        data = r.json()
        assert "variants" in data
