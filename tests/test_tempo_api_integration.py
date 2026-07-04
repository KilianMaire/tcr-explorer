"""Integration tests for TEMPO features in main API."""
from __future__ import annotations
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).parent.parent / "src"))

import pytest
from unittest.mock import patch, AsyncMock, MagicMock
from fastapi.testclient import TestClient


@pytest.fixture
def client():
    with patch.dict("os.environ", {"TEMPO_ENABLE": "true", "BATMAN_ENABLE": "false"}):
        from tcr_explorer.api import app
        yield TestClient(app)


class TestConfig:
    def test_tempo_config_defaults(self):
        from tcr_explorer.config import Settings
        s = Settings()
        assert hasattr(s, "tempo_server_url")
        assert hasattr(s, "tempo_enable")
        assert hasattr(s, "tempo_timeout")


class TestModels:
    def test_gene_record_has_tempo_fields(self):
        from tcr_explorer.models import GeneRecord
        rec = GeneRecord(source="vdjdb", gene_name="test", sequence="CASS")
        assert hasattr(rec, "tempo_score")
        assert hasattr(rec, "tempo_rank")


class TestCrossReactivityEndpoint:
    @patch("tcr_explorer.api.tempo_client")
    def test_predict_crossreactivity(self, mock_tempo, client):
        mock_tempo.post_crossreact = AsyncMock(return_value={
            "results": [
                {
                    "variant_peptide": "ALWNGPMAV",
                    "score": 2.5,
                    "category": "High",
                    "normalized_score": 0.833,
                    "variant_positions": [],
                    "mhc_interface_conserved": True,
                    "length_match": True,
                },
            ],
        })
        payload = {
            "reference_peptide": "LLWNGPMAV",
            "variant_peptides": ["ALWNGPMAV"],
            "mhc_allele": "HLA-A*02:01",
            "mhc_class": "I",
        }
        r = client.post("/predict/crossreactivity", json=payload)
        assert r.status_code == 200
        data = r.json()
        assert "results" in data

    def test_predict_crossreactivity_with_tempo_disabled(self):
        with patch.dict("os.environ", {"TEMPO_ENABLE": "false", "BATMAN_ENABLE": "false"}):
            from importlib import reload
            import tcr_explorer.config
            reload(tcr_explorer.config)
            from tcr_explorer.api import app
            c = TestClient(app)
            payload = {
                "reference_peptide": "LLWNGPMAV",
                "variant_peptides": ["ALWNGPMAV"],
                "mhc_allele": "HLA-A*02:01",
            }
            r = c.post("/predict/crossreactivity", json=payload)
            # Should still work (returns empty or fallback)
            assert r.status_code == 200
