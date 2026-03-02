"""BATMAN Scoring Microservice — port 8105.

Endpoints:
    GET  /health          → status check
    POST /score           → get composite BATMAN score for TCR+epitope pair
    POST /train           → train a BATMAN model for a given TCR
"""
from __future__ import annotations

import logging
import os
import threading
from typing import Optional

import pandas as pd
from fastapi import FastAPI, HTTPException
from pydantic import BaseModel, Field

from batman.cache import CachedModel, ModelCache
from batman.scorer import BATMANScorer
from batman.tcrdist import TCRdistScorer

logger = logging.getLogger(__name__)

app = FastAPI(title="BATMAN Scoring Server", version="0.1.0")

_cache = ModelCache(
    cache_dir=os.getenv("BATMAN_CACHE_DIR", "./batman_cache"),
    max_ram=int(os.getenv("BATMAN_RAM_MODELS", "50")),
)
_scorer = BATMANScorer()
_tcrdist_scorer = TCRdistScorer()
_train_lock = threading.Lock()  # Global lock: pybatman uses shared C state; training is not thread-safe
_BATMAN_STEPS = int(os.getenv("BATMAN_STEPS", "20000"))


# ---------------------------------------------------------------------------
# Request / Response models
# ---------------------------------------------------------------------------
class ScoreRequest(BaseModel):
    tcr_id: str
    index_peptide: str
    candidate_peptide: str
    hla_allele: Optional[str] = None
    # Optional TCR chain identity for tcrdist3 scoring
    cdr3_b: Optional[str] = Field(default=None, min_length=1)
    v_gene: Optional[str] = Field(default=None, min_length=1)
    j_gene: Optional[str] = Field(default=None, min_length=1)


class ScoreResponse(BaseModel):
    tcr_id: str
    candidate_peptide: str
    batman_score: Optional[float] = None  # None = model not yet trained
    tcrdist_score: Optional[float] = None  # None if cdr3_b not provided or no reference
    cache_status: str  # "hit" | "miss"


class ActivationRow(BaseModel):
    peptide: str
    activation: int


class TrainRequest(BaseModel):
    tcr_id: str
    index_peptide: str
    activation_data: list[ActivationRow] = Field(..., min_length=1)


class TrainResponse(BaseModel):
    tcr_id: str
    status: str   # "trained" | "queued"
    n_samples: int


# ---------------------------------------------------------------------------
# Endpoints
# ---------------------------------------------------------------------------
@app.get("/health")
def health() -> dict[str, str]:
    return {"status": "ok", "server": "batman"}


@app.post("/score", response_model=ScoreResponse)
def score(req: ScoreRequest) -> ScoreResponse:
    """Return BATMAN score for a TCR-epitope pair.

    If no trained model is cached for this TCR, returns batman_score=None
    and cache_status="miss". The caller should POST /train to build the model.
    """
    model = _cache.get(req.tcr_id)

    batman_score = None
    cache_status = "miss"
    if model is not None:
        cache_status = "hit"
        try:
            batman_score = _scorer.score(
                index_peptide=model.index_peptide,
                candidate_peptide=req.candidate_peptide,
                weights=model.weights,
                aa_matrix=model.aa_matrix,
            )
        except Exception as exc:
            logger.warning("Scoring failed for TCR %s: %s", req.tcr_id, exc)

    tcrdist_score = None
    if req.cdr3_b is not None and req.v_gene is not None and req.j_gene is not None:
        try:
            tcrdist_score = _tcrdist_scorer.score(
                cdr3_b=req.cdr3_b,
                v_gene=req.v_gene,
                j_gene=req.j_gene,
                epitope=req.candidate_peptide,
            )
        except Exception as exc:
            logger.warning("TCRdist scoring failed for %s: %s", req.tcr_id, exc)

    return ScoreResponse(
        tcr_id=req.tcr_id,
        candidate_peptide=req.candidate_peptide,
        batman_score=batman_score,
        tcrdist_score=tcrdist_score,
        cache_status=cache_status,
    )


@app.post("/train", response_model=TrainResponse)
def train_model(req: TrainRequest) -> TrainResponse:
    """Train a BATMAN model for a TCR from provided activation data.

    Training is synchronous (blocks until complete). For async background
    training from IEDB data, use the batman_prefetch background task.
    """
    df = pd.DataFrame([row.model_dump() for row in req.activation_data])
    df["tcr"] = req.tcr_id
    df["index"] = req.index_peptide

    with _train_lock:
        try:
            result = _scorer.train(df, steps=_BATMAN_STEPS)
        except ValueError as exc:
            raise HTTPException(status_code=422, detail=str(exc)) from exc
        except Exception as exc:
            logger.exception("Training failed for TCR %s", req.tcr_id)
            raise HTTPException(status_code=500, detail="Training failed") from exc

        model = CachedModel(
            tcr_id=req.tcr_id,
            index_peptide=req.index_peptide,
            weights=result.weights,
            aa_matrix=result.aa_matrix,
        )
        _cache.put(model)

    return TrainResponse(
        tcr_id=req.tcr_id,
        status="trained",
        n_samples=len(req.activation_data),
    )
