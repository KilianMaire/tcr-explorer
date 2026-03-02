from __future__ import annotations

from typing import Any, Literal, Optional

from pydantic import BaseModel, Field


GeneSource = Literal["hla", "tcr", "vdjdb", "iedb", "mhc"]
Species = Literal["human", "mouse", "other"]


class IEDBHit(BaseModel):
    epitope_sequence: str
    mhc_allele: Optional[str] = None
    mhc_class: Optional[str] = None
    source_organism: Optional[str] = None
    antigen_name: Optional[str] = None
    assay_type: Optional[str] = None
    effector_cell_type: Optional[str] = None
    qualitative_measure: Optional[str] = None


class GeneRecord(BaseModel):
    source: GeneSource
    species: Species = "other"
    gene_name: str
    allele_name: Optional[str] = None
    region: Optional[str] = None
    sequence: str
    antigen_epitope: Optional[str] = None  # VDJdb: antigen peptide sequence
    metadata: dict[str, Any] = Field(default_factory=dict)
    iedb_hits: Optional[list[IEDBHit]] = None
    # Scoring fields (populated by batman_server enrichment)
    batman_score: Optional[float] = Field(default=None, ge=0.0, le=1.0)
    pmhc_score: Optional[float] = Field(default=None, ge=0.0, le=1.0)
    tcrdist_score: Optional[float] = Field(default=None, ge=0.0, le=1.0)
    composite_score: Optional[float] = Field(default=None, ge=0.0, le=1.0)
    # Scoring fields (tempo_server)
    tempo_score: Optional[float] = None
    tempo_rank: Optional[float] = Field(default=None, ge=0.0, le=100.0)


class SearchRequest(BaseModel):
    source: Optional[GeneSource] = None
    species: Optional[Species] = None
    gene_name: Optional[str] = None
    region: Optional[str] = None
    sequence_contains: Optional[str] = None
    antigen_epitope: Optional[str] = None  # VDJdb: filter by antigen peptide sequence
    limit: int = Field(default=50, ge=1, le=500)
    offset: int = Field(default=0, ge=0, le=10000)


class SearchResponse(BaseModel):
    total: int
    records: list[GeneRecord]
    limit: int = 50    # echoes the requested page size
    offset: int = 0    # echoes the requested page offset


class NLQueryRequest(BaseModel):
    query: str
    limit: int = Field(default=50, ge=1, le=500)


class ParseQueryResult(BaseModel):
    source: Optional[GeneSource] = None
    species: Optional[Species] = None
    gene_name: Optional[str] = None
    region: Optional[str] = None
    sequence_contains: Optional[str] = None
    antigen_epitope: Optional[str] = None


class IngestResponse(BaseModel):
    inserted: int
    source: GeneSource
    species: Species


class CDRPredictResponse(BaseModel):
    v_gene: str
    species: Species
    allele: Optional[str]
    cdr1_aa: Optional[str]
    cdr2_aa: Optional[str]
    cdr1_nt: Optional[str]
    cdr2_nt: Optional[str]


class ReconstructRequest(BaseModel):
    v_gene: str
    j_gene: str
    cdr3_aa: str
    species: Species = "human"


class ReconstructResponse(BaseModel):
    v_gene: str
    j_gene: str
    cdr3_aa: str
    species: Species
    full_nt: Optional[str]       # full TCR coding sequence (V + CDR3 + J)
    full_aa: Optional[str]       # translated protein sequence
    v_region_nt: Optional[str]   # raw V-REGION from stitchr
    cdr3_nt: str                 # back-translated CDR3
    j_region_nt: Optional[str]   # raw J-REGION from stitchr
    v_found: bool
    j_found: bool
    note: str


# ──────────────────────────────────────────────────────────
# TCRpredictor Prediction Models
# ──────────────────────────────────────────────────────────

class BindingPredictionRequest(BaseModel):
    """Request for TCR-Epitope-HLA binding prediction."""
    cdr3_beta: str
    epitope: str
    mhc_allele: str
    cdr3_alpha: Optional[str] = None
    v_beta: Optional[str] = None
    j_beta: Optional[str] = None
    v_alpha: Optional[str] = None
    j_alpha: Optional[str] = None


class BindingPredictionResponse(BaseModel):
    """Response from binding prediction."""
    tier1_score: float = 0.0
    tier2_scores: dict[str, float] = Field(default_factory=dict)
    tier3_structural: Optional[dict[str, Any]] = None
    composite_score: float = 0.0
    mhc_class: str = ""
    binding_core: str = ""
    confidence: str = "low"


class BatchPredictionRequest(BaseModel):
    """Request for batch TCR-Epitope-HLA predictions."""
    predictions: list[BindingPredictionRequest]
    tier1_threshold: float = Field(default=0.5, ge=0.0, le=1.0)
    max_tier2: int = Field(default=100, ge=1, le=10000)


class BatchPredictionResponse(BaseModel):
    """Response from batch prediction."""
    results: list[BindingPredictionResponse]
    total: int
    tier2_processed: int = 0


class ModelStatusResponse(BaseModel):
    """Status of loaded ML models."""
    tier1_loaded: bool = False
    tier1_model: str = ""
    tier2_components: list[str] = Field(default_factory=list)
    tier3_tools: list[str] = Field(default_factory=list)
    ensemble_models: list[str] = Field(default_factory=list)
    last_trained: Optional[str] = None


class PeptidePredictionRequest(BaseModel):
    """Request to predict candidate peptides for a given TCR alpha/beta pair.

    Given TCR CDR3 sequences (alpha and/or beta), searches known epitope
    databases (VDJdb, IEDB) for candidate peptides and scores each one
    through the tiered prediction pipeline.
    """
    cdr3_beta: str
    cdr3_alpha: Optional[str] = None
    v_beta: Optional[str] = None
    j_beta: Optional[str] = None
    v_alpha: Optional[str] = None
    j_alpha: Optional[str] = None
    mhc_allele: Optional[str] = None
    species: Optional[str] = "human"
    max_candidates: int = Field(default=50, ge=1, le=500)
    tier1_threshold: float = Field(default=0.3, ge=0.0, le=1.0)


class PeptideCandidate(BaseModel):
    """A scored peptide candidate predicted for a TCR query."""
    peptide: str
    mhc_allele: str = ""
    mhc_class: str = ""
    # Per-model scores (ensemble integration)
    transformer_score: float = 0.0
    tepcam_score: float = 0.0
    panpep_score: float = 0.0
    ensemble_score: float = 0.0
    # Legacy tiered pipeline scores (backward compatibility)
    tier1_score: float = 0.0
    tier2_scores: dict[str, float] = Field(default_factory=dict)
    composite_score: float = 0.0
    confidence: str = "low"  # "low" / "medium" / "high"
    binding_core: str = ""
    source: str = ""  # "vdjdb", "iedb", or "combined"


class PeptidePredictionResponse(BaseModel):
    """Response containing ranked peptide candidates for the queried TCR."""
    tcr_cdr3_beta: str
    tcr_cdr3_alpha: Optional[str] = None
    mhc_allele: Optional[str] = None
    candidates: list[PeptideCandidate] = Field(default_factory=list)
    total_screened: int = 0
    total_passing: int = 0
