from __future__ import annotations
from typing import Any, Literal, Optional
from pydantic import BaseModel, Field
from .models import IEDBHit

InputType = Literal["auto", "raw_nt", "raw_aa", "gene_name", "allele", "id"]
ModeType = Literal["fast", "full"]
IncludeFlag = Literal["sequences", "germline"]
ChainType = Literal["alpha", "beta", "gamma", "delta", "unknown"]
SpeciesType = Literal["human", "mouse"]
DossierStatus = Literal["complete", "partial"]
WarningCode = Literal[
    "igblast_unavailable", "source_unavailable", "ambiguous_gene",
    "ambiguous_alphabet", "unresolved_input_type", "d_segment_unresolved",
    "aa_annotation_limited", "back_translated_nt", "partial_annotation", "timeout",
]
ProvBlock = Literal["annotation", "germline", "regions", "junction", "full_sequence", "known_epitopes"]
ProvSource = Literal["igblast", "kmer_align", "cdr_enricher", "reconstructor", "vdjdb", "iedb"]
ProvConfidence = Literal["high", "medium", "low"]
ProvKind = Literal["observed", "germline_lookup", "back_translated"]

class DossierRequest(BaseModel):
    query: str
    input_type: InputType = "auto"
    species: SpeciesType = "human"
    mode: ModeType = "fast"
    include: list[IncludeFlag] = Field(default_factory=list)
    v_gene: Optional[str] = None
    j_gene: Optional[str] = None
    cdr3_aa: Optional[str] = None

class GeneCall(BaseModel):
    call: Optional[str] = None
    allele: Optional[str] = None
    score: Optional[float] = None
    score_method: Optional[Literal["igblast_pident", "kmer"]] = None
    germline_nt: Optional[str] = None
    germline_aa: Optional[str] = None

class RegionSeq(BaseModel):
    nt: Optional[str] = None
    aa: Optional[str] = None

class Junction(BaseModel):
    cdr3_nt: Optional[str] = None
    cdr3_aa: Optional[str] = None
    cdr3_nt_is_synthetic: bool = False
    cys104_index: Optional[int] = None
    phe_trp118_index: Optional[int] = None

class Provenance(BaseModel):
    block: ProvBlock
    source: ProvSource
    confidence: ProvConfidence
    kind: Optional[ProvKind] = None

class DossierWarning(BaseModel):
    code: WarningCode
    block: Optional[str] = None
    message: str

class TCRDossier(BaseModel):
    schema_version: str = "1.0"
    status: DossierStatus
    summary: str
    query_echo: dict[str, Any]
    chain: ChainType
    species: SpeciesType
    genes: dict[str, Optional[GeneCall]]
    regions: dict[str, Optional[RegionSeq]]
    junction: Optional[Junction] = None
    full_sequence: Optional[RegionSeq] = None
    known_epitopes: list[IEDBHit] = Field(default_factory=list)
    known_epitopes_total: int = 0
    provenance: list[Provenance] = Field(default_factory=list)
    warnings: list[DossierWarning] = Field(default_factory=list)
