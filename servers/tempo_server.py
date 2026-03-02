"""TEMPO Scoring Microservice — port 8106.

Endpoints:
    GET  /health              → status check
    POST /tempo/score         → TEMPO single-chain log-likelihood score
    POST /tempo/score_paired  → TEMPO paired TCRαβ score
    POST /tsp/profile         → Compute TSP for epitope-specific TCRs
    POST /crossreact/predict  → MixTCRcross cross-reactivity prediction
    GET  /batcave/variants    → Query BATCAVE for epitope variants
"""
from __future__ import annotations

import logging
from typing import Optional

from fastapi import FastAPI
from pydantic import BaseModel, Field

from tempo.baseline import BaselineModel, load_baseline
from tempo.scorer import TempoScorer
from tempo.tsp import compute_tsp
from tempo.crossreact import predict_cross_reactivity
from tempo.batcave import BatcaveClient
from tempo.tempo_binary import is_available as tempo_binary_available, predict as tempo_binary_predict, perc_rank_to_score, list_epitopes

logger = logging.getLogger(__name__)

app = FastAPI(title="TEMPO Scoring Server", version="0.1.0")

# Lazy-loaded baselines
_baselines: dict[str, BaselineModel] = {}
_batcave = BatcaveClient()


def _get_baseline(species: str, chain: str) -> Optional[BaselineModel]:
    """Load baseline lazily, return None if unavailable."""
    key = f"{species}_{chain}"
    if key not in _baselines:
        try:
            _baselines[key] = load_baseline(species, chain)
        except FileNotFoundError:
            logger.warning("Baseline not found for %s_%s, using empty baseline", species, chain)
            _baselines[key] = BaselineModel(species=species, chain=chain)
    return _baselines[key]


# ---------------------------------------------------------------------------
# Request / Response models
# ---------------------------------------------------------------------------
class SingleChainScoreRequest(BaseModel):
    v_gene: str
    j_gene: str
    cdr3: str
    chain: str = "beta"
    species: str = "human"


class PairedScoreRequest(BaseModel):
    v_a: str
    j_a: str
    cdr3_a: str
    v_b: str
    j_b: str
    cdr3_b: str
    species: str = "human"


class ScoreResponse(BaseModel):
    log_likelihood: float
    v_contribution: float = 0.0
    j_contribution: float = 0.0
    length_contribution: float = 0.0
    cdr3_contribution: float = 0.0


class TSPRequest(BaseModel):
    tcrs: list[dict[str, str]]
    chain: str = "beta"
    species: str = "human"


class CrossReactRequest(BaseModel):
    reference_peptide: str
    variant_peptides: list[str]
    mhc_allele: str
    mhc_class: str = "I"


class CrossReactVariantResult(BaseModel):
    variant_peptide: str
    score: float  # 0-3 scale
    category: str = "Low"  # High/Medium/Low
    normalized_score: float = 0.0  # 0-1 for backward compatibility
    variant_positions: list[dict] = []
    mhc_interface_conserved: bool = True
    length_match: bool = True


class CrossReactResponse(BaseModel):
    results: list[CrossReactVariantResult]


# ---------------------------------------------------------------------------
# Endpoints
# ---------------------------------------------------------------------------
@app.get("/health")
def health() -> dict:
    return {"status": "ok", "server": "tempo", "binary_available": tempo_binary_available()}


@app.post("/tempo/score", response_model=ScoreResponse)
def tempo_score(req: SingleChainScoreRequest) -> ScoreResponse:
    bl = _get_baseline(req.species, req.chain)
    scorer = TempoScorer(
        alpha_baseline=bl if req.chain == "alpha" else None,
        beta_baseline=bl if req.chain == "beta" else None,
    )
    result = scorer.score_single_chain(req.v_gene, req.j_gene, req.cdr3, req.chain)
    return ScoreResponse(
        log_likelihood=result.log_likelihood,
        v_contribution=result.v_contribution,
        j_contribution=result.j_contribution,
        length_contribution=result.length_contribution,
        cdr3_contribution=result.cdr3_contribution,
    )


@app.post("/tempo/score_paired", response_model=ScoreResponse)
def tempo_score_paired(req: PairedScoreRequest) -> ScoreResponse:
    bl_a = _get_baseline(req.species, "alpha")
    bl_b = _get_baseline(req.species, "beta")
    scorer = TempoScorer(alpha_baseline=bl_a, beta_baseline=bl_b)
    result = scorer.score_paired(
        req.v_a, req.j_a, req.cdr3_a,
        req.v_b, req.j_b, req.cdr3_b,
    )
    return ScoreResponse(
        log_likelihood=result.log_likelihood,
        v_contribution=result.v_contribution,
        j_contribution=result.j_contribution,
        length_contribution=result.length_contribution,
        cdr3_contribution=result.cdr3_contribution,
    )


@app.post("/tsp/profile")
def tsp_profile(req: TSPRequest) -> dict:
    bl = _get_baseline(req.species, req.chain)
    profile = compute_tsp(req.tcrs, bl)
    return profile.to_dict()


@app.post("/crossreact/predict", response_model=CrossReactResponse)
def crossreact_predict(req: CrossReactRequest) -> CrossReactResponse:
    results = []
    for variant in req.variant_peptides:
        r = predict_cross_reactivity(
            reference_peptide=req.reference_peptide,
            variant_peptide=variant,
            mhc_allele=req.mhc_allele,
            mhc_class=req.mhc_class,
        )
        results.append(CrossReactVariantResult(
            variant_peptide=variant,
            score=r.score,
            category=r.category,
            normalized_score=r.normalized_score,
            variant_positions=r.variant_positions,
            mhc_interface_conserved=r.mhc_interface_conserved,
            length_match=r.length_match,
        ))
    return CrossReactResponse(results=results)


@app.get("/batcave/variants")
def batcave_variants(
    reference_peptide: Optional[str] = None,
    mhc_class: Optional[str] = None,
    mhc_allele: Optional[str] = None,
) -> dict:
    variants = _batcave.lookup(
        reference_peptide=reference_peptide,
        mhc_class=mhc_class,
        mhc_allele=mhc_allele,
    )
    return {
        "variants": [
            {
                "reference_peptide": v.reference_peptide,
                "variant_peptide": v.variant_peptide,
                "activation_score": v.activation_score,
                "mhc_allele": v.mhc_allele,
                "mhc_class": v.mhc_class,
                "mutation_position": v.mutation_position,
                "original_aa": v.original_aa,
                "mutant_aa": v.mutant_aa,
            }
            for v in variants
        ],
        "total": len(variants),
    }


# ---------------------------------------------------------------------------
# TEMPO binary endpoints
# ---------------------------------------------------------------------------
class BinaryPredictRequest(BaseModel):
    tcrs: list[dict[str, str]]
    epitope_id: str
    chain: str = "AB"
    species: str = "HomoSapiens"


class BinaryPredictResult(BaseModel):
    perc_rank: Optional[float] = None
    score: Optional[float] = None
    problem: str = ""


class BinaryPredictResponse(BaseModel):
    results: list[BinaryPredictResult]
    binary_available: bool
    epitope_id: str


@app.get("/tempo/epitopes")
def get_epitopes() -> dict:
    """List available TEMPO epitope models."""
    return {"epitopes": list_epitopes(), "binary_available": tempo_binary_available()}


@app.post("/tempo/predict", response_model=BinaryPredictResponse)
def tempo_predict_binary(req: BinaryPredictRequest) -> BinaryPredictResponse:
    """Run TEMPO binary prediction (real motif models)."""
    if not tempo_binary_available():
        return BinaryPredictResponse(
            results=[BinaryPredictResult(problem="TEMPO binary not available")
                     for _ in req.tcrs],
            binary_available=False,
            epitope_id=req.epitope_id,
        )
    try:
        preds = tempo_binary_predict(
            req.tcrs, req.epitope_id, req.chain, req.species,
        )
        results = []
        for p in preds:
            score = perc_rank_to_score(p.perc_rank) if p.perc_rank is not None else None
            results.append(BinaryPredictResult(
                perc_rank=p.perc_rank,
                score=score,
                problem=p.problem,
            ))
        return BinaryPredictResponse(
            results=results,
            binary_available=True,
            epitope_id=req.epitope_id,
        )
    except Exception as e:
        return BinaryPredictResponse(
            results=[BinaryPredictResult(problem=str(e)) for _ in req.tcrs],
            binary_available=True,
            epitope_id=req.epitope_id,
        )
