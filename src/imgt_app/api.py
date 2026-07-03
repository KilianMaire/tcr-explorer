from __future__ import annotations

import asyncio
import csv
import io
import logging
import math as _math
from contextlib import asynccontextmanager
from typing import Optional

import httpx
from pydantic import BaseModel

from fastapi import FastAPI, File, Form, Request, UploadFile
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import HTMLResponse, PlainTextResponse

from .cdr_enricher import get_cdr1_cdr2
from .config import settings
from .dossier_models import (
    AlignRequest,
    AskRequest,
    AskResponse,
    DossierRequest,
    MSAResult,
    RecordsRequest,
    RecordsResponse,
    SimilarRequest,
    SimilarResponse,
    TCRDossier,
)
from .fasta_parser import parse_cdr3_fasta, parse_fasta_bytes
from .file_ingest import parse_file, parse_vdjdb_tsv
from .mcp_clients import ToolServerClient
from .models import (
    CDRPredictResponse, GeneRecord, GeneSource, IEDBHit, IngestResponse, NLQueryRequest,
    ReconstructRequest, ReconstructResponse,
    SearchRequest, SearchResponse, Species,
)
from .nl_query import lmstudio_parse
from .reconstructor import reconstruct_tcr
from .search_index import SearchIndex

_logger = logging.getLogger(__name__)


@asynccontextmanager
async def lifespan(app: FastAPI):
    """Startup/shutdown lifespan: preload models."""
    # Load Tier 1 (backward compatibility)
    try:
        from pipeline.tier1_screening import load_model, is_loaded
        if not is_loaded():
            ok = load_model()
            if ok:
                _logger.info("Tier 1 model preloaded at startup")
            else:
                _logger.warning("Tier 1 model not loaded (no checkpoint or missing deps)")
    except Exception as exc:
        _logger.warning("Tier 1 preload skipped: %s", exc)

    # Load ensemble scorer (all 3 models + weights)
    try:
        from pipeline.ensemble_scorer import get_scorer
        scorer = get_scorer()
        ok = scorer.load()
        if ok:
            _logger.info("Ensemble scorer loaded: %s", scorer.available_models())
        else:
            _logger.warning("No ensemble model checkpoints found")
    except Exception as exc:
        _logger.warning("Ensemble preload skipped: %s", exc)

    yield


app = FastAPI(title="IMGT Search Engine", version="0.1.0", lifespan=lifespan)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)
index = SearchIndex(settings.database_path)
hla_client = ToolServerClient(settings.hla_server_url)
tcr_client = ToolServerClient(settings.tcr_server_url)
vdjdb_client = ToolServerClient(settings.vdjdb_server_url)
iedb_client = ToolServerClient(settings.iedb_server_url)
mhc_client = ToolServerClient(settings.mhc_server_url)


class _BatmanClient:
    """Thin HTTP client for batman_server scoring/training endpoints."""

    def __init__(self) -> None:
        self.base_url = settings.batman_server_url.rstrip("/")
        self.timeout = settings.batman_timeout

    async def post_score(self, payload: dict) -> dict:
        async with httpx.AsyncClient(timeout=self.timeout) as client:
            r = await client.post(f"{self.base_url}/score", json=payload)
            r.raise_for_status()
            return r.json()

    async def post_train(self, payload: dict) -> dict:
        async with httpx.AsyncClient(timeout=self.timeout) as client:
            r = await client.post(f"{self.base_url}/train", json=payload)
            r.raise_for_status()
            return r.json()


batman_client = _BatmanClient()


class _TempoClient:
    """Thin HTTP client for tempo_server endpoints."""

    def __init__(self) -> None:
        self.base_url = settings.tempo_server_url.rstrip("/")
        self.timeout = settings.tempo_timeout

    async def post_crossreact(self, payload: dict) -> dict:
        async with httpx.AsyncClient(timeout=self.timeout) as client:
            r = await client.post(f"{self.base_url}/crossreact/predict", json=payload)
            r.raise_for_status()
            return r.json()

    async def post_tempo_score(self, payload: dict) -> dict:
        async with httpx.AsyncClient(timeout=self.timeout) as client:
            r = await client.post(f"{self.base_url}/tempo/score", json=payload)
            r.raise_for_status()
            return r.json()

    async def post_tsp_profile(self, payload: dict) -> dict:
        async with httpx.AsyncClient(timeout=self.timeout) as client:
            r = await client.post(f"{self.base_url}/tsp/profile", json=payload)
            r.raise_for_status()
            return r.json()


tempo_client = _TempoClient()


def _composite_score(
    batman: float | None,
    pmhc: float | None,
    tcrdist: float | None,
) -> float | None:
    """Geometric mean of available scores. Returns None if all are None."""
    scores = [s for s in (batman, pmhc, tcrdist) if s is not None]
    if not scores:
        return None
    product = _math.prod(scores)
    return product ** (1.0 / len(scores))


async def _enrich_with_batman(
    records: list[GeneRecord],
    req: SearchRequest,
) -> list[GeneRecord]:
    """Call batman_server /score for each VDJdb record with an epitope.

    Attaches batman_score, pmhc_score, tcrdist_score, composite_score.
    Skips records without antigen_epitope. Silently ignores server errors.
    """
    if not settings.batman_enable:
        return records

    enriched = []
    for rec in records:
        epitope = rec.antigen_epitope
        if not epitope:
            enriched.append(rec)
            continue

        meta = rec.metadata or {}
        hla = meta.get("mhc_a") or getattr(req, "hla_allele", None) or ""

        payload = {
            "tcr_id": rec.gene_name,
            "index_peptide": epitope,
            "candidate_peptide": epitope,
            "hla_allele": hla or None,
            "cdr3_b": rec.sequence or None,
            "v_gene": meta.get("v_segm") or None,
            "j_gene": meta.get("j_segm") or None,
        }

        try:
            resp = await batman_client.post_score(payload)
            b_score = resp.get("batman_score")
            p_score = resp.get("pmhc_score")
            t_score = resp.get("tcrdist_score")
            comp = _composite_score(b_score, p_score, t_score)
            enriched.append(rec.model_copy(update={
                "batman_score": b_score,
                "pmhc_score": p_score,
                "tcrdist_score": t_score,
                "composite_score": comp,
            }))
        except Exception:
            enriched.append(rec)

    return enriched


def _coerce_search_response(value: "SearchResponse | dict") -> SearchResponse:
    """Accept either a SearchResponse or a raw dict (as returned by mocks in tests)."""
    if isinstance(value, dict):
        return SearchResponse(**value)
    return value


def _merge_results(local: SearchResponse, remote: SearchResponse, req: SearchRequest) -> SearchResponse:
    """Merge local index and remote tool-server results, then apply offset + limit."""
    combined = local.records + remote.records
    page = combined[req.offset: req.offset + req.limit]
    return SearchResponse(
        total=local.total + remote.total,
        records=page,
        limit=req.limit,
        offset=req.offset,
    )


_IEDB_HITS_CAP = 5


def _enrich_with_iedb(
    vdjdb_records: list[GeneRecord],
    iedb_response: SearchResponse,
) -> list[GeneRecord]:
    """Attach IEDB assay hits to VDJdb records sharing the same antigen epitope.
    IEDB hits with no matching VDJdb record become phantom GeneRecords
    (sequence="", gene_name=epitope) acting as extensible slots."""
    iedb_by_epitope: dict[str, list[IEDBHit]] = {}
    for rec in iedb_response.records:
        ep = rec.sequence.upper()
        if not ep:
            continue
        hit = IEDBHit(
            epitope_sequence=rec.sequence,
            mhc_allele=rec.gene_name if rec.gene_name and rec.gene_name != "unknown" else None,
            mhc_class=rec.metadata.get("mhc_class"),
            source_organism=rec.metadata.get("source_organism"),
            antigen_name=rec.metadata.get("antigen_name"),
            assay_type=rec.metadata.get("assay_type"),
            effector_cell_type=rec.metadata.get("effector_cell_type"),
            qualitative_measure=rec.metadata.get("qualitative_measure"),
        )
        iedb_by_epitope.setdefault(ep, [])
        if len(iedb_by_epitope[ep]) < _IEDB_HITS_CAP:
            iedb_by_epitope[ep].append(hit)

    matched_epitopes: set[str] = set()
    enriched: list[GeneRecord] = []
    for rec in vdjdb_records:
        ep_key = (rec.antigen_epitope or "").upper()
        if ep_key and ep_key in iedb_by_epitope:
            matched_epitopes.add(ep_key)
            enriched.append(rec.model_copy(update={"iedb_hits": iedb_by_epitope[ep_key]}))
        else:
            enriched.append(rec)

    for ep_key, hits in iedb_by_epitope.items():
        if ep_key not in matched_epitopes:
            phantom = GeneRecord(
                source="vdjdb",
                gene_name=hits[0].epitope_sequence,
                sequence="",
                iedb_hits=hits,
            )
            enriched.append(phantom)

    return enriched


class UserActivationRow(BaseModel):
    peptide: str
    activation: int  # 0=none, 1=weak, 2=strong


class PredictActivationRequest(BaseModel):
    tcr_cdr3: str
    tcr_v_gene: Optional[str] = None
    tcr_j_gene: Optional[str] = None
    hla_allele: Optional[str] = None
    candidate_peptides: list[str]
    user_activation_data: list[UserActivationRow] = []


class PeptideScore(BaseModel):
    peptide: str
    batman_score: Optional[float] = None
    pmhc_score: Optional[float] = None
    tcrdist_score: Optional[float] = None
    composite_score: Optional[float] = None
    mhc_class: Optional[str] = None   # "I" or "II"


class PredictActivationResponse(BaseModel):
    results: list[PeptideScore]
    tcr_id: str
    hla_allele: Optional[str] = None


class CrossReactivityRequest(BaseModel):
    reference_peptide: str
    variant_peptides: list[str]
    mhc_allele: str
    mhc_class: str = "I"


class CrossReactivityVariant(BaseModel):
    variant_peptide: str
    score: float  # 0-3 scale
    category: str = "Low"
    normalized_score: float = 0.0
    variant_positions: list[dict] = []
    mhc_interface_conserved: bool = True
    length_match: bool = True


class CrossReactivityResponse(BaseModel):
    results: list[CrossReactivityVariant]


class TSPProfileRequest(BaseModel):
    tcrs: list[dict[str, str]]
    chain: str = "beta"
    species: str = "human"


@app.get("/health")
def health() -> dict[str, str]:
    return {"status": "ok"}


@app.get("/predict/cdr", response_model=CDRPredictResponse)
def predict_cdr(v_gene: str, species: Species = "human") -> CDRPredictResponse:
    """Return CDR1 and CDR2 amino acid sequences for a TCR V gene using stitchr IMGT data."""
    result = get_cdr1_cdr2(v_gene, species)
    return CDRPredictResponse(v_gene=v_gene, species=species, **result)


@app.post("/ingest/fasta", response_model=IngestResponse)
async def ingest_fasta(
    source: GeneSource = Form(...),
    species: Species = Form("other"),
    file: UploadFile = File(...),
) -> IngestResponse:
    raw = await file.read()
    parsed = list(parse_fasta_bytes(raw, source=source, default_species=species))
    inserted = index.upsert_many(parsed)
    return IngestResponse(inserted=inserted, source=source, species=species)


@app.post("/ingest/file", response_model=IngestResponse)
async def ingest_file(
    source: GeneSource = Form(...),
    species: Species = Form("other"),
    file: UploadFile = File(...),
) -> IngestResponse:
    raw = await file.read()
    parsed = parse_file(raw, file.filename or "unknown", source=source, species=species)
    inserted = index.upsert_many(parsed)
    return IngestResponse(inserted=inserted, source=source, species=species)


@app.post("/ingest/vdjdb", response_model=IngestResponse)
async def ingest_vdjdb(file: UploadFile = File(...)) -> IngestResponse:
    """Bulk-load a VDJdb TSV (or CSV) export into the local search index."""
    raw = await file.read()
    records = parse_vdjdb_tsv(raw)
    inserted = index.upsert_many(records)
    return IngestResponse(inserted=inserted, source="vdjdb", species="other")


@app.post("/search", response_model=SearchResponse)
async def search(req: SearchRequest) -> SearchResponse:
    local = index.search(req)

    # Tool servers don't support offset; inflate their limit so we fetch enough
    # records to correctly apply offset after merging with the local results.
    tool_req = req.model_copy(update={"limit": req.limit + req.offset, "offset": 0}) if req.offset else req

    remote = SearchResponse(total=0, records=[], limit=req.limit, offset=0)
    if req.source == "hla":
        remote = _coerce_search_response(await hla_client.search(tool_req))
    elif req.source == "tcr":
        remote = _coerce_search_response(await tcr_client.search(tool_req))
    elif req.source == "vdjdb":
        iedb_req = tool_req.model_copy(update={"source": None})
        vdjdb_raw, iedb_raw = await asyncio.gather(
            vdjdb_client.search(tool_req),
            iedb_client.search(iedb_req),
            return_exceptions=True,
        )
        remote = _coerce_search_response(
            vdjdb_raw if not isinstance(vdjdb_raw, Exception)
            else {"total": 0, "records": [], "limit": req.limit, "offset": 0}
        )
        iedb_remote = _coerce_search_response(
            iedb_raw if not isinstance(iedb_raw, Exception)
            else {"total": 0, "records": [], "limit": req.limit, "offset": 0}
        )
        merged_local = _merge_results(local, remote, req)
        enriched_records = _enrich_with_iedb(merged_local.records, iedb_remote)
        # Batman composite scoring
        enriched_records = await _enrich_with_batman(enriched_records, req)
        return SearchResponse(
            total=merged_local.total,
            records=enriched_records,
            limit=req.limit,
            offset=req.offset,
        )
    elif req.source == "iedb":
        remote = _coerce_search_response(await iedb_client.search(tool_req))
    elif req.source == "mhc":
        remote = _coerce_search_response(await mhc_client.search(tool_req))

    return _merge_results(local, remote, req)


@app.post("/query/nl", response_model=SearchResponse)
async def search_nl(req: NLQueryRequest) -> SearchResponse:
    parsed = await lmstudio_parse(req.query)
    search_req = SearchRequest(**parsed.model_dump(), limit=req.limit)
    return await search(search_req)


@app.post("/search/fasta", response_model=SearchResponse)
async def search_fasta(
    file: UploadFile = File(...),
    source: Optional[GeneSource] = Form(None),
    species: Optional[Species] = Form(None),
    limit: int = Form(50),
) -> SearchResponse:
    """
    Accept a FASTA file whose sequences are CDR3 amino-acid strings.
    Searches each CDR3 against VDJdb (or the specified source) and returns
    merged, deduplicated results.  Limited to 20 sequences per request.
    """
    raw = await file.read()
    entries = parse_cdr3_fasta(raw)[:20]   # cap at 20 CDR3s per call

    if not entries:
        return SearchResponse(total=0, records=[], limit=limit, offset=0)

    tasks = [
        search(SearchRequest(
            source=source or "vdjdb",
            species=species,
            sequence_contains=cdr3,
            limit=limit,
        ))
        for _, cdr3 in entries
    ]
    results = await asyncio.gather(*tasks)

    seen: set[str] = set()
    merged = []
    for r in results:
        for rec in r.records:
            key = rec.sequence + "|" + rec.gene_name
            if key not in seen:
                seen.add(key)
                merged.append(rec)

    return SearchResponse(total=len(merged), records=merged[:limit], limit=limit, offset=0)


@app.post("/reconstruct", response_model=ReconstructResponse)
def reconstruct(req: ReconstructRequest) -> ReconstructResponse:
    """
    Reconstruct a full TCR coding sequence from V gene + CDR3 + J gene.

    Uses stitchr IMGT germline data for V/J regions; CDR3 is back-translated
    with human-optimised codons.  Assembly follows IMGT/VDJdb CDR3 boundaries:
    Cys104 … Phe/Trp118 (inclusive).
    """
    result = reconstruct_tcr(req.v_gene, req.j_gene, req.cdr3_aa, req.species)
    return ReconstructResponse(**result)


@app.post("/search/table")
async def search_table(req: SearchRequest, fmt: str = "csv"):
    result = await search(req)
    headers = ["source", "species", "gene_name", "allele_name", "region", "sequence"]

    rows = [
        [r.source, r.species, r.gene_name, r.allele_name or "", r.region or "", r.sequence]
        for r in result.records
    ]

    if fmt == "md":
        md = ["| " + " | ".join(headers) + " |", "|" + "|".join(["---"] * len(headers)) + "|"]
        for row in rows:
            md.append("| " + " | ".join(row) + " |")
        return PlainTextResponse("\n".join(md), media_type="text/markdown")

    out = io.StringIO()
    writer = csv.writer(out)
    writer.writerow(headers)
    writer.writerows(rows)
    return PlainTextResponse(out.getvalue(), media_type="text/csv")


@app.post("/predict/activation", response_model=PredictActivationResponse)
async def predict_activation(req: PredictActivationRequest) -> PredictActivationResponse:
    """Score a list of candidate peptides against a given TCR + HLA.

    If user_activation_data is provided, sends a /train request to batman_server
    first to fit a custom BATMAN model. Otherwise uses cached model.

    Results are sorted by composite_score descending.
    """
    import hashlib
    tcr_id = f"query_{hashlib.md5(req.tcr_cdr3.encode()).hexdigest()[:8]}"

    # Train if user supplied activation data
    if req.user_activation_data:
        if settings.batman_enable:
            train_payload = {
                "tcr_id": tcr_id,
                "index_peptide": req.tcr_cdr3,
                "activation_data": [row.model_dump() for row in req.user_activation_data],
            }
            try:
                await batman_client.post_train(train_payload)
            except Exception:
                pass  # continue without custom model
        else:
            import logging
            logging.getLogger(__name__).warning(
                "user_activation_data provided but BATMAN_ENABLE=false; training skipped"
            )

    # Empty peptide list → early return
    if not req.candidate_peptides:
        return PredictActivationResponse(
            results=[], tcr_id=tcr_id, hla_allele=req.hla_allele
        )

    if not settings.batman_enable:
        return PredictActivationResponse(
            results=[PeptideScore(peptide=p) for p in req.candidate_peptides],
            tcr_id=tcr_id,
            hla_allele=req.hla_allele,
        )

    # Determine MHC class from allele (inline to avoid path manipulation)
    mhc_class: Optional[str] = None
    if req.hla_allele:
        mhc_class = (
            "II"
            if req.hla_allele.upper().startswith(
                ("HLA-D", "DRA", "DRB", "DQA", "DQB", "DPA", "DPB")
            )
            else "I"
        )

    # Score each candidate peptide
    results: list[PeptideScore] = []
    for peptide in req.candidate_peptides:
        score_payload = {
            "tcr_id": tcr_id,
            "index_peptide": req.tcr_cdr3,
            "candidate_peptide": peptide,
            "hla_allele": req.hla_allele,
            "cdr3_b": req.tcr_cdr3,
            "v_gene": req.tcr_v_gene,
            "j_gene": req.tcr_j_gene,
        }
        try:
            resp = await batman_client.post_score(score_payload)
            b = resp.get("batman_score")
            p = resp.get("pmhc_score")
            t = resp.get("tcrdist_score")
            comp = _composite_score(b, p, t)
        except Exception:
            b = p = t = comp = None

        results.append(PeptideScore(
            peptide=peptide,
            batman_score=b,
            pmhc_score=p,
            tcrdist_score=t,
            composite_score=comp,
            mhc_class=mhc_class,
        ))

    # Sort by composite_score descending (None last)
    results.sort(
        key=lambda x: x.composite_score if x.composite_score is not None else -float("inf"),
        reverse=True,
    )

    return PredictActivationResponse(
        results=results, tcr_id=tcr_id, hla_allele=req.hla_allele
    )


@app.post("/predict/crossreactivity", response_model=CrossReactivityResponse)
async def predict_crossreactivity(req: CrossReactivityRequest) -> CrossReactivityResponse:
    """Predict cross-reactivity between a reference epitope and peptide variants."""
    if not settings.tempo_enable:
        return CrossReactivityResponse(results=[
            CrossReactivityVariant(variant_peptide=p, score=0.0)
            for p in req.variant_peptides
        ])
    try:
        resp = await tempo_client.post_crossreact(req.model_dump())
        return CrossReactivityResponse(
            results=[CrossReactivityVariant(**r) for r in resp["results"]]
        )
    except Exception:
        return CrossReactivityResponse(results=[
            CrossReactivityVariant(variant_peptide=p, score=0.0)
            for p in req.variant_peptides
        ])


@app.post("/tsp/profile")
async def tsp_profile(req: TSPProfileRequest) -> dict:
    """Generate a TCR Specificity Profile for epitope-specific TCRs."""
    if not settings.tempo_enable:
        return {"error": "TEMPO is disabled", "v_enrichment": {}, "j_enrichment": {},
                "cdr3_length_dist": {}, "cdr3_motifs": {}, "n_tcrs": 0}
    try:
        return await tempo_client.post_tsp_profile(req.model_dump())
    except Exception:
        return {"error": "TEMPO server unavailable", "v_enrichment": {}, "j_enrichment": {},
                "cdr3_length_dist": {}, "cdr3_motifs": {}, "n_tcrs": 0}


# ──────────────────────────────────────────────────────────
# TCRpredictor Prediction Endpoints
# ──────────────────────────────────────────────────────────

from .models import (
    BindingPredictionRequest, BindingPredictionResponse,
    BatchPredictionRequest, BatchPredictionResponse,
    ModelStatusResponse,
    PeptidePredictionRequest, PeptidePredictionResponse, PeptideCandidate,
)


@app.post("/predict/binding", response_model=BindingPredictionResponse)
async def predict_binding(req: BindingPredictionRequest):
    """Predict TCR-Epitope-HLA binding through full tiered pipeline."""
    from pipeline.orchestrator import PredictionRequest, run_prediction

    prediction_req = PredictionRequest(
        cdr3_beta=req.cdr3_beta,
        epitope=req.epitope,
        mhc_allele=req.mhc_allele,
        cdr3_alpha=req.cdr3_alpha or "",
        v_beta=req.v_beta or "",
        j_beta=req.j_beta or "",
        v_alpha=req.v_alpha or "",
        j_alpha=req.j_alpha or "",
    )

    result = await run_prediction(prediction_req)

    return BindingPredictionResponse(
        tier1_score=result.tier1_score,
        tier2_scores=result.tier2_scores,
        tier3_structural=result.tier3_structural,
        composite_score=result.composite_score,
        mhc_class=result.mhc_class,
        binding_core=result.binding_core,
        confidence=result.confidence,
    )


@app.post("/predict/screen", response_model=BindingPredictionResponse)
async def predict_screen(req: BindingPredictionRequest):
    """Fast Tier 1 screening only (no Tier 2/3)."""
    from pipeline.orchestrator import PredictionRequest, TierConfig, run_prediction

    prediction_req = PredictionRequest(
        cdr3_beta=req.cdr3_beta,
        epitope=req.epitope,
        mhc_allele=req.mhc_allele,
        cdr3_alpha=req.cdr3_alpha or "",
    )

    config = TierConfig(tier1_enabled=True, tier2_enabled=False, tier3_enabled=False)
    result = await run_prediction(prediction_req, config)

    return BindingPredictionResponse(
        tier1_score=result.tier1_score,
        composite_score=result.tier1_score,
        mhc_class=result.mhc_class,
        confidence=result.confidence,
    )


@app.post("/predict/batch", response_model=BatchPredictionResponse)
async def predict_batch(req: BatchPredictionRequest):
    """Batch prediction with Tier 1 screening and selective Tier 2."""
    from pipeline.orchestrator import PredictionRequest, TierConfig, run_batch_prediction

    requests = [
        PredictionRequest(
            cdr3_beta=p.cdr3_beta,
            epitope=p.epitope,
            mhc_allele=p.mhc_allele,
            cdr3_alpha=p.cdr3_alpha or "",
        )
        for p in req.predictions
    ]

    config = TierConfig(
        tier1_threshold=req.tier1_threshold,
        tier2_top_n=req.max_tier2,
        tier3_enabled=False,
    )
    results = await run_batch_prediction(requests, config)

    responses = [
        BindingPredictionResponse(
            tier1_score=r.tier1_score,
            tier2_scores=r.tier2_scores,
            composite_score=r.composite_score,
            mhc_class=r.mhc_class,
            binding_core=r.binding_core,
            confidence=r.confidence,
        )
        for r in results
    ]

    tier2_count = sum(1 for r in results if r.tier2_scores)

    return BatchPredictionResponse(
        results=responses,
        total=len(responses),
        tier2_processed=tier2_count,
    )


@app.get("/models/status", response_model=ModelStatusResponse)
async def model_status():
    """Get status of loaded ML models."""
    tier1_loaded = False
    tier1_model = "BindingPredictor (not loaded)"
    try:
        from pipeline.tier1_screening import is_loaded
        tier1_loaded = is_loaded()
        if tier1_loaded:
            tier1_model = "BindingPredictor (loaded)"
    except ImportError:
        pass

    tier2 = []
    if settings.batman_enable:
        tier2.append("batman")
    if settings.tempo_enable:
        tier2.append("tempo")
    tier2.extend(["tcrdist", "pmtnet_omni", "mixtcrpred"])

    import os
    tier3 = []
    if os.environ.get("TCRMODEL2_SERVER_URL"):
        tier3.append("tcrmodel2")
    if os.environ.get("TFOLD_SERVER_URL"):
        tier3.append("tfold")
    if os.environ.get("AF3_SERVER_URL"):
        tier3.append("af3")

    # Ensemble models
    ensemble_models = []
    try:
        from pipeline.ensemble_scorer import get_scorer
        ensemble_models = get_scorer().available_models()
    except Exception:
        pass

    return ModelStatusResponse(
        tier1_loaded=tier1_loaded,
        tier1_model=tier1_model,
        tier2_components=tier2,
        tier3_tools=tier3 or ["none configured"],
        ensemble_models=ensemble_models,
    )


@app.post("/predict/peptides", response_model=PeptidePredictionResponse)
async def predict_peptides(req: PeptidePredictionRequest):
    """Predict candidate peptides for a given TCR alpha/beta query.

    Pipeline:
    1. Search VDJdb + IEDB for known epitopes matching the species
    2. Collect unique candidate peptides
    3. Score each candidate with ensemble scorer (Transformer + TEPCam + PanPep)
    4. Return ranked peptide candidates sorted by ensemble score
    """
    # Step 1: Gather candidate epitopes from VDJdb + IEDB
    candidate_epitopes: dict[str, dict] = {}  # peptide -> metadata

    # Search VDJdb by CDR3 beta sequence
    try:
        vdjdb_req = SearchRequest(
            source="vdjdb",
            species=req.species,
            sequence_contains=req.cdr3_beta,
            limit=200,
        )
        vdjdb_resp = await search(vdjdb_req)
        for rec in vdjdb_resp.records:
            ep = rec.antigen_epitope
            if ep and ep not in candidate_epitopes:
                meta = rec.metadata or {}
                candidate_epitopes[ep] = {
                    "mhc_allele": meta.get("mhc_a", req.mhc_allele or ""),
                    "source": "vdjdb",
                }
    except Exception as exc:
        _logger.warning("VDJdb search failed: %s", exc)

    # If few candidates from CDR3 match, also search broadly by species
    if len(candidate_epitopes) < 10:
        try:
            broad_req = SearchRequest(
                source="vdjdb",
                species=req.species,
                limit=200,
            )
            broad_resp = await search(broad_req)
            for rec in broad_resp.records:
                ep = rec.antigen_epitope
                if ep and ep not in candidate_epitopes:
                    meta = rec.metadata or {}
                    candidate_epitopes[ep] = {
                        "mhc_allele": meta.get("mhc_a", req.mhc_allele or ""),
                        "source": "vdjdb",
                    }
        except Exception as exc:
            _logger.warning("Broad VDJdb search failed: %s", exc)

    if not candidate_epitopes:
        return PeptidePredictionResponse(
            tcr_cdr3_beta=req.cdr3_beta,
            tcr_cdr3_alpha=req.cdr3_alpha,
            mhc_allele=req.mhc_allele,
            total_screened=0,
            total_passing=0,
        )

    # Step 2: Score with ensemble scorer (preferred) or tiered pipeline (fallback)
    from pipeline.ensemble_scorer import get_scorer
    scorer = get_scorer()
    use_ensemble = len(scorer.available_models()) > 0

    candidates: list[PeptideCandidate] = []
    epitope_list = list(candidate_epitopes.items())[:req.max_candidates]

    if use_ensemble:
        # ── Ensemble path ──
        for epitope, meta in epitope_list:
            allele = meta.get("mhc_allele") or req.mhc_allele or "HLA-A*02:01"
            try:
                result = scorer.score(
                    cdr3_beta=req.cdr3_beta,
                    peptide=epitope,
                    mhc_allele=allele,
                )
                ensemble_score = result.get("ensemble", 0.0)
                confidence = (
                    "high" if ensemble_score >= 0.7
                    else "medium" if ensemble_score >= 0.5
                    else "low"
                )
                candidates.append(PeptideCandidate(
                    peptide=epitope,
                    mhc_allele=allele,
                    transformer_score=result.get("transformer", 0.0),
                    tepcam_score=result.get("tepcam", 0.0),
                    panpep_score=result.get("panpep", 0.0),
                    ensemble_score=ensemble_score,
                    composite_score=ensemble_score,
                    confidence=confidence,
                    source=meta.get("source", "unknown"),
                ))
            except Exception as exc:
                _logger.warning("Ensemble scoring failed for epitope %s: %s", epitope, exc)
    else:
        # ── Tiered pipeline fallback ──
        from pipeline.orchestrator import PredictionRequest as PipelineRequest, TierConfig, run_prediction
        config = TierConfig(
            tier1_enabled=True,
            tier2_enabled=True,
            tier3_enabled=False,
            tier1_threshold=req.tier1_threshold,
        )
        for epitope, meta in epitope_list:
            allele = meta.get("mhc_allele") or req.mhc_allele or "HLA-A*02:01"
            prediction_req = PipelineRequest(
                cdr3_beta=req.cdr3_beta,
                epitope=epitope,
                mhc_allele=allele,
                cdr3_alpha=req.cdr3_alpha or "",
                v_beta=req.v_beta or "",
                j_beta=req.j_beta or "",
                v_alpha=req.v_alpha or "",
                j_alpha=req.j_alpha or "",
            )
            try:
                result = await run_prediction(prediction_req, config)
                candidates.append(PeptideCandidate(
                    peptide=epitope,
                    mhc_allele=allele,
                    mhc_class=result.mhc_class,
                    tier1_score=result.tier1_score,
                    tier2_scores=result.tier2_scores,
                    composite_score=result.composite_score,
                    confidence=result.confidence,
                    binding_core=result.binding_core,
                    source=meta.get("source", "unknown"),
                ))
            except Exception as exc:
                _logger.warning("Prediction failed for epitope %s: %s", epitope, exc)

    # Step 3: Sort by ensemble_score (or composite_score) descending
    candidates.sort(
        key=lambda c: c.ensemble_score if c.ensemble_score > 0 else c.composite_score,
        reverse=True,
    )

    total_passing = sum(
        1 for c in candidates
        if (c.ensemble_score or c.composite_score) >= req.tier1_threshold
    )

    return PeptidePredictionResponse(
        tcr_cdr3_beta=req.cdr3_beta,
        tcr_cdr3_alpha=req.cdr3_alpha,
        mhc_allele=req.mhc_allele,
        candidates=candidates,
        total_screened=len(epitope_list),
        total_passing=total_passing,
    )


def _dossier_markdown(d: TCRDossier) -> str:
    lines = [f"# TCR Dossier ({d.status})", "", d.summary, "",
             f"- chain: {d.chain}", f"- species: {d.species}"]
    if d.genes.get("v"):
        lines.append(f"- V: {d.genes['v'].call} ({d.genes['v'].score_method})")
    if d.known_epitopes:
        lines.append(f"- known epitopes: {d.known_epitopes_total}")
    if d.warnings:
        lines.append("- warnings: " + ", ".join(w.code for w in d.warnings))
    return "\n".join(lines)


@app.post("/v1/tcr/dossier", response_model=TCRDossier)
def tcr_dossier(req: DossierRequest, request: Request):
    # Synchronous by design: build_dossier is fully synchronous and its epitope
    # lookup (dossier_epitopes._run_search) must create its own event loop, which
    # it can only do off the request's running loop. FastAPI runs a sync route in
    # a threadpool, so known_epitopes actually surface (an async route left them
    # permanently empty).
    from .dossier import build_dossier  # local import: avoids circular import
    # (dossier -> dossier_epitopes -> api.search/_IEDB_HITS_CAP)

    d = build_dossier(req)
    if "text/markdown" in request.headers.get("accept", ""):
        return PlainTextResponse(_dossier_markdown(d), media_type="text/markdown")
    return d


@app.post("/v1/tcr/similar", response_model=SimilarResponse)
def tcr_similar(req: SimilarRequest):
    # Synchronous by design, same rationale as /v1/tcr/dossier above: the engine
    # reads a parquet index off disk, so it runs in FastAPI's threadpool rather
    # than blocking the event loop.
    from .similarity import find_similar_tcrs  # local import: avoids import cycles

    neigh, engine, total, warnings = find_similar_tcrs(
        req.cdr3,
        req.v_gene,
        req.j_gene,
        species=req.species,
        top_k=req.top_k,
        min_similarity=req.min_similarity,
    )
    return SimilarResponse(neighbours=neigh, engine=engine, total_candidates=total, warnings=warnings)


@app.post("/v1/tcr/records", response_model=RecordsResponse)
def tcr_records(req: RecordsRequest):
    # Synchronous by design, same rationale as /v1/tcr/dossier above: retrieval
    # reads a parquet index off disk via sync helpers that must not run inside
    # a running event loop, so it runs in FastAPI's threadpool.
    from .records import retrieve_records  # local import: avoids import cycles

    return retrieve_records(req)


@app.post("/v1/tcr/align", response_model=MSAResult)
def tcr_align(req: AlignRequest):
    # Synchronous by design, same rationale as /v1/tcr/dossier above: align may
    # resolve germline sets from disk or shell out to clustalo, so it runs in
    # FastAPI's threadpool rather than blocking the event loop.
    from .msa import align  # local import: avoids import cycles

    return align(req)


@app.post("/v1/tcr/ask", response_model=AskResponse)
def tcr_ask(req: AskRequest):
    # Synchronous by design, same rationale as /v1/tcr/dossier and /v1/tcr/similar
    # above: the routed intent may hit build_dossier/find_similar_tcrs, both of
    # which are themselves synchronous for the same reasons.
    from .ask import answer  # local import: avoids circular import

    return answer(req)


_UI_HTML = """<!doctype html><html lang="en"><head><meta charset="utf-8">
<title>TCR Explorer</title><meta name="viewport" content="width=device-width, initial-scale=1">
<style>
 body{font-family:system-ui,sans-serif;max-width:900px;margin:2rem auto;padding:0 1rem;color:#111}
 h1{font-size:1.3rem} input,select,button{font-size:1rem;padding:.4rem}
 #q{width:60%} .card{border:1px solid #ddd;border-radius:8px;padding:1rem;margin:1rem 0}
 .warn{color:#a15c00} .syn{color:#7a3e00;font-style:italic} table{border-collapse:collapse;width:100%}
 .loading{color:#0b5;font-weight:600} .loading::after{content:'';display:inline-block;width:.7em;height:.7em;margin-left:.4em;border:2px solid #0b5;border-top-color:transparent;border-radius:50%;animation:spin .7s linear infinite;vertical-align:middle} @keyframes spin{to{transform:rotate(360deg)}} button:disabled{opacity:.6}
 td,th{border:1px solid #eee;padding:.3rem;text-align:left;font-size:.9rem} .muted{color:#666}
 h3{margin:.3rem 0}
 .rec{border:1px solid #e2e2e2;border-radius:6px;padding:.6rem .8rem;margin:.5rem 0;background:#fafafa}
 .rec.neigh{background:#f5f8ff;border-style:dashed}
 .badge{display:inline-block;background:#333;color:#fff;font-size:.75rem;font-weight:700;padding:.1rem .4rem;border-radius:4px;letter-spacing:.03em}
 .kind{color:#666;font-size:.8rem}
 .comp{font-family:monospace;font-size:.85rem;background:#eee;padding:.2rem .4rem;border-radius:4px;display:inline-block}
 .partner{margin-top:.4rem;padding-left:.6rem;border-left:3px solid #ccc}
</style></head><body>
<h1>TCR Explorer</h1>
<p class="muted">Ask about a TCR: a gene (TRBV20-1), a CDR3 (CASSLGTEAFF), a sequence, or a V+J+CDR3. Known epitopes are retrieved; similar-TCR epitopes are inferred.</p>
<form id="f"><input id="q" placeholder="e.g. TRBV20-1" value="TRBV20-1">
<select id="sp"><option>human</option><option>mouse</option></select>
<button type="submit">Ask</button></form>
<div id="out"></div>
<h1>Align a gene set</h1>
<p class="muted">Align a germline set (species + chain + segment) or a gene list. V/J/C come from the germline source; D is not available there.</p>
<form id="af">
<select id="a_sp"><option>human</option><option>mouse</option></select>
<input id="a_chain" placeholder="chain e.g. TRB" value="TRB">
<select id="a_seg"><option>V</option><option>D</option><option>J</option><option>C</option></select>
<label><input type="checkbox" id="a_translate"> translate</label>
<button type="submit">Align</button></form>
<div id="a_out"></div>
<script>
const f=document.getElementById('f'),out=document.getElementById('out');
f.addEventListener('submit',async e=>{e.preventDefault();const btn=f.querySelector('button');const t0=btn.textContent;btn.disabled=true;btn.textContent='Searching...';out.innerHTML='<p class="loading">Searching...</p>';
 try{
  const q=document.getElementById('q').value,sp=document.getElementById('sp').value;
  const r=await fetch('/v1/tcr/ask',{method:'POST',headers:{'Content-Type':'application/json'},
   body:JSON.stringify({query:q,species:sp})});
  if(!r.ok){out.innerHTML='<p class="warn">Request failed ('+r.status+')</p>';return;}
  const b=await r.json();out.innerHTML=render(b);
  try{
   const rr=await fetch('/v1/tcr/records',{method:'POST',headers:{'Content-Type':'application/json'},
    body:JSON.stringify({query:q,species:sp})});
   if(rr.ok){const rb=await rr.json();out.innerHTML+=renderRecords(rb);}
  }catch(recErr){/* records lookup is additive, never fatal to the ask flow */}
 }catch(err){out.innerHTML='<p class="warn">Error: '+esc(String(err))+'</p>';}
 finally{btn.disabled=false;btn.textContent=t0;}});
function esc(s){return (s==null?'':String(s)).replace(/[&<>]/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;'}[c]));}
function render(b){let h=`<div class="card"><h3>intent: ${esc(b.intent)} <span class="muted">(source ${esc(b.plan_source)}, llm ${b.llm_used})</span></h3>`;
 if(b.dossier){const d=b.dossier;h+=`<p><b>${esc(d.summary)}</b></p><p>chain: ${esc(d.chain)} · species: ${esc(d.species)} · status: ${esc(d.status)}</p>`;
  if(d.status==='partial' && !(d.genes&&d.genes.v) && !(d.known_epitopes&&d.known_epitopes.length)){h+='<p class="muted">A CDR3 on its own cannot identify V/D/J. Provide the V and J genes (V+J+CDR3), a gene name, or a full V(D)J sequence to get an annotation.</p>';}
  if(d.genes&&d.genes.v){h+=`<p>V: ${esc(d.genes.v.call)} (${esc(d.genes.v.score_method)})</p>`;}
  if(d.regions){h+='<p>';for(const k of ['cdr1','cdr2','cdr3']){if(d.regions[k]&&d.regions[k].aa)h+=`${k}: <code>${esc(d.regions[k].aa)}</code> `;}h+='</p>';}
  if(d.junction&&d.junction.cdr3_nt_is_synthetic){h+=`<p class="syn">cdr3_nt is synthetic (back-translated)</p>`;}
  if(d.known_epitopes&&d.known_epitopes.length){h+='<h3>Known epitopes (retrieved)</h3><table><tr><th>epitope</th><th>MHC</th><th>antigen</th></tr>';
   for(const e of d.known_epitopes)h+=`<tr><td>${esc(e.epitope_sequence)}</td><td>${esc(e.mhc_allele)}</td><td>${esc(e.antigen_name)}</td></tr>`;h+='</table>';}
  if(d.neighbours&&d.neighbours.length){h+=neighTable(d.neighbours);}
  if(d.warnings&&d.warnings.length){h+='<p class="warn">warnings: '+d.warnings.map(w=>esc(w.code)).join(', ')+'</p>';}}
 if(b.neighbours_result){const nr=b.neighbours_result;h+=neighTable(nr.neighbours);if(nr.warnings&&nr.warnings.length){h+='<p class="warn">'+nr.warnings.map(w=>esc(w.message||w.code)).join('; ')+'</p>';}}
 if(b.search_result){h+=`<p class="muted">search returned ${b.search_result.total} records</p>`;}
 return h+'</div>';}
function neighTable(ns){if(!ns||!ns.length)return '<p class="muted">No matching known TCRs found in the reference database for this CDR3.</p>';
 let h='<h3>Similar TCRs (inferred, not confirmed specificity)</h3><table><tr><th>CDR3</th><th>V</th><th>sim</th><th>epitope</th></tr>';
 for(const n of ns)h+=`<tr><td>${esc(n.cdr3_b_aa)}</td><td>${esc(n.v_b_gene)}</td><td>${esc(n.similarity)}</td><td>${esc(n.epitope_aa)}</td></tr>`;return h+'</table>';}
function sourceBadge(s){return '<span class="badge">'+esc((s||'').toUpperCase())+'</span>';}
function compStrip(c){if(!c)return '';
 return '<div class="comp">'+esc(c.v_germline_aa||'-')+' | '+esc(c.cdr3_aa||'-')+' | '+esc(c.j_germline_aa||'-')+'</div>';}
function partnerLine(rec,pairs){if(!pairs||!pairs.length)return '';
 for(const p of pairs){
  if(p.pairing_key!==rec.pairing_key)continue;
  const other=(rec.chain==='alpha')?p.beta:p.alpha;
  if(other&&other.cdr3_aa!==rec.cdr3_aa){
   return '<div class="partner">paired '+esc(other.chain)+': <code>'+esc(other.cdr3_aa)+'</code> ('+esc(other.source)+')</div>';
  }
 }
 return '';}
function recordCard(rec,pairs,cls){let h='<div class="rec '+cls+'">';
 h+=sourceBadge(rec.source)+' <span class="muted">'+esc(rec.chain)+' / '+esc(rec.species)+'</span>';
 if(rec.external_url)h+=' <a href="'+esc(rec.external_url)+'" target="_blank" rel="noopener">source</a>';
 h+='<br>aa: <code>'+esc(rec.full_aa||rec.cdr3_aa)+'</code>';
 if(rec.full_aa_kind)h+=' <span class="kind">('+esc(rec.full_aa_kind)+')</span>';
 const nt=rec.full_nt||rec.cdr3_nt;
 const ntKind=rec.full_nt?rec.full_nt_kind:rec.cdr3_nt_kind;
 if(nt&&ntKind){h+='<br>nt: <code>'+esc(nt)+'</code> <span class="kind">('+esc(ntKind)+')</span>';}
 const genes=[];
 if(rec.v_gene)genes.push('V '+esc(rec.v_gene));
 if(rec.d_gene)genes.push('D '+esc(rec.d_gene));
 if(rec.j_gene)genes.push('J '+esc(rec.j_gene));
 if(genes.length)h+='<br>'+genes.join(' &middot; ');
 const cdrs=[];
 if(rec.cdr1_aa)cdrs.push('CDR1 '+esc(rec.cdr1_aa));
 if(rec.cdr2_aa)cdrs.push('CDR2 '+esc(rec.cdr2_aa));
 if(rec.cdr3_aa)cdrs.push('CDR3 '+esc(rec.cdr3_aa));
 if(cdrs.length)h+='<br>'+cdrs.join(' &middot; ');
 if(rec.epitope_aa||rec.mhc_a){
  h+='<br>epitope: '+esc(rec.epitope_aa||'unknown');
  if(rec.mhc_class||rec.mhc_a)h+=' &middot; MHC '+esc(rec.mhc_class||'')+' '+esc(rec.mhc_a||'');
 }
 if(rec.composition)h+='<br>'+compStrip(rec.composition);
 h+=partnerLine(rec,pairs);
 return h+'</div>';}
function renderRecords(data){if(!data)return '';
 let h='<h3>Exact records</h3>';
 if(data.exact&&data.exact.length){for(const rec of data.exact)h+=recordCard(rec,data.pairs,'exact');}
 else h+='<p class="muted">No exact records found.</p>';
 h+='<h3>Near neighbours</h3>';
 if(data.neighbours&&data.neighbours.length){for(const rec of data.neighbours)h+=recordCard(rec,data.pairs,'neigh');}
 else h+='<p class="muted">No near neighbours found.</p>';
 if(data.warnings&&data.warnings.length){h+='<p class="warn">'+esc(data.warnings.map(w=>w.message||w.code).join('; '))+'</p>';}
 return h;}
const af=document.getElementById('af'),a_out=document.getElementById('a_out');
af.addEventListener('submit',async e=>{e.preventDefault();const btn=af.querySelector('button');const t0=btn.textContent;btn.disabled=true;btn.textContent='Aligning...';a_out.innerHTML='<p class="loading">Aligning...</p>';
 try{
  const r=await fetch('/v1/tcr/align',{method:'POST',headers:{'Content-Type':'application/json'},
   body:JSON.stringify({species:document.getElementById('a_sp').value,
    chain:document.getElementById('a_chain').value,
    segment:document.getElementById('a_seg').value,
    translate:document.getElementById('a_translate').checked})});
  if(!r.ok){a_out.innerHTML='<p class="warn">Request failed ('+r.status+')</p>';return;}
  const b=await r.json();a_out.innerHTML=renderAlign(b);
 }catch(err){a_out.innerHTML='<p class="warn">Error: '+esc(String(err))+'</p>';}
 finally{btn.disabled=false;btn.textContent=t0;}});
function shade(c){c=c||0;if(c>=0.9)return 'background:#08519c;color:#fff';if(c>=0.7)return 'background:#3182bd;color:#fff';if(c>=0.5)return 'background:#6baed6';if(c>=0.3)return 'background:#bdd7e7';return '';}
function cell(t,s){return '<span style="'+s+'">'+esc(t)+'</span>';}
function renderAlign(b){let h='<div class="card"><h3>engine: '+esc(b.engine)+' <span class="muted">('+esc(b.n_sequences)+' sequences, '+esc(b.mean_pct_identity)+'% identity, view '+esc(b.view)+')</span></h3>';
 h+='<p class="muted">Shading shows per-column conservation (darker is more conserved).</p><pre>';
 const cons=b.conservation||[];
 for(const rec of (b.records||[])){
  const name=esc(rec.name.padEnd(12));
  if(rec.aligned_aa && rec.aligned_nt){
   let aa=name+' aa  ',nt=name+' nt  ';
   for(let i=0;i<rec.aligned_aa.length;i++){aa+=cell(' '+rec.aligned_aa[i]+' ',shade(cons[i]));nt+=cell(rec.aligned_nt.substr(3*i,3),shade(cons[i]));}
   h+=aa+'<br>'+nt+'<br><br>';
  }else{
   const s=rec.aligned_aa||rec.aligned||rec.aligned_nt||'';let row=name+'  ';
   for(let i=0;i<s.length;i++){row+=cell(s[i],shade(cons[i]));}
   h+=row+'<br>';
  }
 }
 h+='</pre>';
 if(b.warnings&&b.warnings.length){h+='<p class="warn">warnings: '+b.warnings.map(w=>esc(w.code)).join(', ')+'</p>';}
 return h+'</div>';}
</script></body></html>"""


@app.get("/ui", response_class=HTMLResponse)
def ui() -> HTMLResponse:
    return HTMLResponse(_UI_HTML)
