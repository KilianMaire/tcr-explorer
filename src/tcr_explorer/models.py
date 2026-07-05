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


class SearchRequest(BaseModel):
    source: Optional[GeneSource] = Field(None, description="Data source to search: 'hla', 'tcr', 'vdjdb', 'iedb', or 'mhc'.", examples=["tcr", "vdjdb", "iedb"])
    species: Optional[Species] = Field(None, description="Filter by organism: 'human', 'mouse', or 'other'.", examples=["human", "mouse"])
    gene_name: Optional[str] = Field(None, description="IMGT gene name to search, e.g. TRBV20-1.", examples=["TRBV20-1", "TRBJ2-7"])
    region: Optional[str] = Field(None, description="Gene region, e.g. 'V-REGION' or 'J-REGION'.", examples=["V-REGION", "J-REGION"])
    sequence_contains: Optional[str] = Field(None, description="Nucleotide or amino acid substring to search within sequences.", examples=["CASS", "TGTGCG"])
    antigen_epitope: Optional[str] = Field(None, description="VDJdb/IEDB: filter by antigen peptide sequence.", examples=["GILGFVFTL", "MSKGEDFCIKVQCAKTGLSLV"])
    limit: int = Field(default=50, ge=1, le=500, description="Maximum number of records to return (1-500).", examples=[50, 100, 200])
    offset: int = Field(default=0, ge=0, le=10000, description="Pagination offset, e.g. skip the first N records (0-10000).", examples=[0, 50, 100])


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
    v_gene: Optional[str] = None   # optional: inferred from the CDR3 if omitted
    j_gene: Optional[str] = None   # optional: inferred from the CDR3 if omitted
    cdr3_aa: str
    species: Species = "human"


class ReconstructResponse(BaseModel):
    v_gene: str
    j_gene: str
    v_allele_used: Optional[str] = None   # germline V allele used (defaults to *01)
    j_allele_used: Optional[str] = None   # germline J allele used (defaults to *01)
    cdr3_aa: str
    species: Species
    full_nt: Optional[str]       # variable-domain coding sequence (V + CDR3 + J)
    full_aa: Optional[str]       # translated variable-domain protein
    full_chain_aa: Optional[str] = None   # variable domain + membrane-bound constant
    constant_source: Optional[str] = None  # provenance of the appended constant
    v_region_nt: Optional[str]   # raw V-REGION from stitchr
    cdr3_nt: str                 # back-translated CDR3
    j_region_nt: Optional[str]   # raw J-REGION from stitchr
    v_prefix_nt: Optional[str] = None   # in-frame V up to Cys104
    j_suffix_nt: Optional[str] = None   # in-frame J FR4 after Phe/Trp118
    v_found: bool
    j_found: bool
    genes_inferred: bool = False           # True when V/J were inferred from the CDR3
    inference_support: Optional[int] = None  # records backing the chosen V/J pairing
    inference_alternatives: Optional[list[str]] = None  # other pairings, e.g. "TRBV19/TRBJ2-7 (n=3)"
    note: str


# ──────────────────────────────────────────────────────────
# TCR Allele Assignment Models
# ──────────────────────────────────────────────────────────

class AlleleCallModel(BaseModel):
    """A germline allele call: candidate alleles tied at the best identity."""
    alleles: list[str]
    identity: float
    aligned_span: list[int]


class DCallModel(AlleleCallModel):
    """D-gene call. D assignment is inherently noisy (short segment, heavy
    exonuclease trimming), so callers get an explicit low_confidence flag."""
    low_confidence: bool


class AssignRequest(BaseModel):
    sequence: str
    species: Optional[Species] = None
    chain: Optional[str] = None
    want_d: bool = False


class AssignResponse(BaseModel):
    input_kind: str
    species: str
    chain: Optional[str] = None
    v_call: Optional[AlleleCallModel] = None
    j_call: Optional[AlleleCallModel] = None
    d_call: Optional[DCallModel] = None
    constant_call: Optional[AlleleCallModel] = None
    regions: dict[str, float] = Field(default_factory=dict)
    cdr3_aa: Optional[str] = None
    v_determinable: bool = True
    v_reason: Optional[str] = None
    v_db_inference: Optional[list] = None
    reconstruction: Optional[dict] = None
    warnings: list[str] = Field(default_factory=list)


# ──────────────────────────────────────────────────────────
# Query Router Models
# ──────────────────────────────────────────────────────────

class QueryRequest(BaseModel):
    query: str
    species: Optional[Species] = None
    force: Optional[str] = None


class QueryBlock(BaseModel):
    tool: str
    title: str
    data: dict


class QueryUnderstanding(BaseModel):
    input: str
    detected_type: str
    species: Optional[str] = None
    tools: list[str]
    note: str


class QueryResponse(BaseModel):
    understood: QueryUnderstanding
    blocks: list[QueryBlock]
    warnings: list[str] = Field(default_factory=list)
