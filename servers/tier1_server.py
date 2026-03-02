"""Tier 1 ML Screening FastAPI microservice.

Standalone server for fast TCR-epitope-HLA binding screening
using the custom Transformer + ESM-2 model.
Port: 8110
"""
from __future__ import annotations

import logging
from typing import Optional

from fastapi import FastAPI
from pydantic import BaseModel

logger = logging.getLogger(__name__)

app = FastAPI(title="TCRpredictor Tier 1 Screening", version="0.1.0")


class ScreenRequest(BaseModel):
    cdr3_beta: str
    epitope: str
    mhc_allele: str
    cdr3_alpha: Optional[str] = None


class ScreenResponse(BaseModel):
    score: float
    mhc_class: str


@app.post("/screen", response_model=ScreenResponse)
async def screen(req: ScreenRequest):
    """Screen a single TCR-epitope-HLA triad."""
    # Inline MHC class detection (was in pipeline.orchestrator, now removed)
    allele_upper = (req.mhc_allele or "").upper()
    if any(allele_upper.startswith(p) for p in ("HLA-A", "HLA-B", "HLA-C")):
        mhc_class = "I"
    elif any(allele_upper.startswith(p) for p in ("HLA-D", "HLA-DR", "HLA-DQ", "HLA-DP")):
        mhc_class = "II"
    else:
        mhc_class = "I"  # default

    # Placeholder: model loading and inference
    logger.warning("Tier 1 model not yet trained, returning default score")
    return ScreenResponse(score=0.5, mhc_class=mhc_class)


@app.get("/health")
async def health():
    return {"status": "ok", "model_loaded": False}


if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8110)
