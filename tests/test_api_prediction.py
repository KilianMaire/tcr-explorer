"""Tests for TCRpredictor prediction API endpoints."""
import pytest
from fastapi.testclient import TestClient

from src.imgt_app.models import (
    BindingPredictionRequest,
    BindingPredictionResponse,
    BatchPredictionRequest,
    ModelStatusResponse,
)


def test_binding_prediction_request():
    req = BindingPredictionRequest(
        cdr3_beta="CASSLAPGATNEKLFF",
        epitope="PKYVKQNTLKLAT",
        mhc_allele="HLA-DRB1*04:01",
    )
    assert req.cdr3_beta == "CASSLAPGATNEKLFF"
    assert req.cdr3_alpha is None


def test_binding_prediction_response():
    resp = BindingPredictionResponse(
        tier1_score=0.85,
        composite_score=0.82,
        mhc_class="II",
        confidence="high",
    )
    assert resp.tier1_score == 0.85
    assert resp.confidence == "high"


def test_batch_prediction_request():
    req = BatchPredictionRequest(
        predictions=[
            BindingPredictionRequest(
                cdr3_beta="CASSLAPGATNEKLFF",
                epitope="PKYVKQNTLKLAT",
                mhc_allele="HLA-DRB1*04:01",
            ),
            BindingPredictionRequest(
                cdr3_beta="CASSLGQAYEQYF",
                epitope="GILGFVFTL",
                mhc_allele="HLA-A*02:01",
            ),
        ],
        tier1_threshold=0.5,
    )
    assert len(req.predictions) == 2


def test_model_status_response():
    resp = ModelStatusResponse(
        tier1_loaded=False,
        tier1_model="BindingPredictor",
        tier2_components=["batman", "tempo"],
        tier3_tools=["tcrmodel2", "af3"],
    )
    assert resp.tier1_loaded is False
    assert len(resp.tier2_components) == 2
