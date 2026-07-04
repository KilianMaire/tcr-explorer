from __future__ import annotations
from typing import Any, Literal, Optional
from pydantic import BaseModel, Field
from .models import IEDBHit

InputType = Literal["auto", "raw_nt", "raw_aa", "gene_name", "allele", "id"]
ModeType = Literal["fast", "full"]
IncludeFlag = Literal["sequences", "germline", "neighbours"]
ChainType = Literal["alpha", "beta", "gamma", "delta", "unknown"]
SpeciesType = Literal["human", "mouse"]
DossierStatus = Literal["complete", "partial"]
WarningCode = Literal[
    "igblast_unavailable", "source_unavailable", "ambiguous_gene",
    "ambiguous_alphabet", "unresolved_input_type", "d_segment_unresolved",
    "aa_annotation_limited", "back_translated_nt", "partial_annotation", "timeout",
    "tcrdist_unavailable", "similarity_index_unavailable", "no_reference_candidates",
    "species_unsupported", "segment_unavailable", "too_few_sequences", "alignment_failed",
    "no_pairing_found", "records_index_unavailable", "records_index_stale",
    "nonstandard_residues",
]
ProvBlock = Literal["annotation", "germline", "regions", "junction", "full_sequence", "known_epitopes", "neighbours", "alignment"]
ProvSource = Literal["igblast", "kmer_align", "cdr_enricher", "reconstructor", "vdjdb", "iedb", "unitcr", "tcrdist", "blosum_cdr3"]
ProvConfidence = Literal["high", "medium", "low"]
ProvKind = Literal["observed", "germline_lookup", "back_translated", "neighbor_inferred"]

class DossierRequest(BaseModel):
    query: str = Field(..., description="A TCR query: a raw nucleotide or amino acid sequence, an IMGT gene name or allele, or a namespaced database id (vdjdb:/iedb:).", examples=["TRBV20-1", "CASSLGTEAFF", "vdjdb:12345"])
    input_type: InputType = Field("auto", description="Force the interpretation of query; 'auto' runs the router.", examples=["auto", "gene_name", "raw_aa"])
    species: SpeciesType = Field("human", description="Organism of the query.", examples=["human", "mouse"])
    mode: ModeType = Field("fast", description="'fast' uses gene/id/CDR3 paths only; 'full' enables raw-sequence annotation and neighbours.", examples=["fast", "full"])
    include: list[IncludeFlag] = Field(default_factory=list, description="Optional extras: 'sequences' (long nt), 'germline', 'neighbours'.", examples=[["sequences"], ["neighbours"]])
    v_gene: Optional[str] = Field(None, description="Optional V gene for reconstruction/neighbours.", examples=["TRBV20-1"])
    j_gene: Optional[str] = Field(None, description="Optional J gene for reconstruction/neighbours.", examples=["TRBJ2-7"])
    cdr3_aa: Optional[str] = Field(None, description="Optional CDR3 amino acids for reconstruction/neighbours.", examples=["CASSLGTEAFF"])

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

class Neighbour(BaseModel):
    cdr3_b_aa: str
    v_b_gene: Optional[str] = None
    j_b_gene: Optional[str] = None
    similarity: float
    """Within-query relative similarity: normalised to the per-query candidate
    maximum distance, so only similarity==1.0 (identical) is absolute. Callers
    comparing across queries should threshold on the absolute `distance` field,
    not on `similarity`."""
    distance: float
    epitope_aa: Optional[str] = None
    mhc_class: Optional[str] = None
    mhc_a: Optional[str] = None
    antigen: Optional[str] = None
    antigen_organism: Optional[str] = None
    cluster_id: Optional[int] = None

class SimilarRequest(BaseModel):
    cdr3: str = Field(..., description="CDR3 amino acid sequence for similarity search.", examples=["CASSLGTEAFF", "CASSPYLGTNQYF"])
    v_gene: str = Field(..., description="V gene name, e.g. TRBV20-1", examples=["TRBV20-1", "TRBV5-1"])
    j_gene: str = Field(..., description="J gene name, e.g. TRBJ2-7", examples=["TRBJ2-7", "TRBJ1-2"])
    species: SpeciesType = Field("human", description="Organism of the query.", examples=["human", "mouse"])
    top_k: int = Field(default=10, ge=1, le=200, description="Maximum number of similar neighbours to return (1-200).", examples=[10, 50, 100])
    min_similarity: float = Field(default=0.0, ge=0.0, le=1.0, description="Minimum relative similarity (0.0-1.0) for a match; only 1.0 is truly identical.", examples=[0.0, 0.5, 0.8])

class SimilarResponse(BaseModel):
    neighbours: list[Neighbour] = Field(default_factory=list)
    engine: str
    total_candidates: int = 0
    warnings: list["DossierWarning"] = Field(default_factory=list)

class AskRequest(BaseModel):
    query: str = Field(..., description="A free-text question about a TCR.", examples=["annotate CASSLGTEAFF", "what does TRBV20-1 recognise?"])
    species: SpeciesType = "human"

class AskResponse(BaseModel):
    intent: str
    plan_source: str          # "llm" | "heuristic"
    llm_used: bool
    dossier: Optional["TCRDossier"] = None
    neighbours_result: Optional["SimilarResponse"] = None
    search_result: Optional[Any] = None   # SearchResponse
    warnings: list["DossierWarning"] = Field(
        default_factory=list,
        description=(
            "Warnings from the dossier path only. Similar/search-path warnings "
            "are nested in neighbours_result.warnings / search_result respectively, "
            "not surfaced here."
        ),
    )

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
    neighbours: Optional[list[Neighbour]] = None
    records: list["TCRRecord"] = Field(default_factory=list)

class AlignedRecord(BaseModel):
    name: str
    aligned: str
    aligned_aa: Optional[str] = None   # gapped amino acids (codon-aware view)
    aligned_nt: Optional[str] = None   # gapped nucleotides, 3 per aa column, in register

class ProvidedSeq(BaseModel):
    name: str
    seq: str

class AlignRequest(BaseModel):
    species: SpeciesType = "human"
    chain: Optional[str] = Field(None, description="IMGT chain locus for a germline set: TRA/TRB/TRG/TRD.", examples=["TRB"])
    segment: Optional[str] = Field(None, description="Segment for a germline set: V/D/J/C. D is not in the germline source.", examples=["J"])
    genes: Optional[list[str]] = Field(None, description="Explicit gene/allele names to resolve and align.", examples=[["TRBJ1-1", "TRBJ2-7"]])
    sequences: Optional[list[ProvidedSeq]] = Field(None, description="Provided named sequences to align verbatim.")
    seq_type: str = Field("nt", description="'nt' or 'aa'.", examples=["nt", "aa"])
    translate: bool = Field(False, description="Translate germline nt to aa before aligning.")

class MSAResult(BaseModel):
    engine: str
    seq_type: str
    n_sequences: int
    alignment_length: int
    records: list[AlignedRecord] = Field(default_factory=list)
    consensus: str = ""
    mean_pct_identity: float = 0.0
    conservation: list[float] = Field(default_factory=list)  # per aa column, 0..1
    view: str = "nt"                                          # "aa_nt" | "aa" | "nt"
    provenance: list[Provenance] = Field(default_factory=list)
    warnings: list[DossierWarning] = Field(default_factory=list)


class Composition(BaseModel):
    v_germline_aa: Optional[str] = Field(None, description="V gene germline amino acid sequence used to reconstruct this record.")
    cdr3_aa: str = Field(..., description="CDR3 amino acid sequence carried through the reconstruction.")
    j_germline_aa: Optional[str] = Field(None, description="J gene germline amino acid sequence used to reconstruct this record.")
    v_region_nt: Optional[str] = Field(None, description="V region germline nucleotide sequence used to reconstruct this record.")
    cdr3_nt: Optional[str] = Field(None, description="CDR3 nucleotide sequence, back translated when not deposited.")
    j_region_nt: Optional[str] = Field(None, description="J region germline nucleotide sequence used to reconstruct this record.")
    note: Optional[str] = Field(None, description="Free text explanation of how this sequence was reconstructed.")


class TCRRecord(BaseModel):
    source: str
    source_record_id: str
    external_url: str = Field(..., description="Link to view this record on the source database's own website.")
    pairing_key: str = Field(..., description="Key shared by the alpha and beta chain of the same receptor, for pairing.")
    chain: str
    species: str
    cdr3_aa: str
    cdr3_nt: Optional[str] = None
    cdr3_nt_kind: Optional[Literal["deposited", "reconstructed"]] = Field(None, description="Whether cdr3_nt was deposited by the source database or back translated from the amino acids.")
    full_aa: Optional[str] = None
    full_aa_kind: Optional[Literal["deposited", "reconstructed"]] = Field(None, description="Whether full_aa was deposited by the source database or reconstructed from germline plus CDR3.")
    full_nt: Optional[str] = None
    full_nt_kind: Optional[Literal["deposited", "reconstructed"]] = Field(None, description="Whether full_nt was deposited by the source database or reconstructed from germline plus CDR3.")
    nt_is_synthetic: bool = Field(False, description="True when any nucleotide sequence on this record was back translated rather than deposited; a synthetic nucleotide sequence is not the receptor's real DNA.")
    v_gene: Optional[str] = None
    d_gene: Optional[str] = None
    j_gene: Optional[str] = None
    cdr1_aa: Optional[str] = Field(None, description="CDR1 amino acid sequence, germline encoded and largely V gene determined.")
    cdr2_aa: Optional[str] = Field(None, description="CDR2 amino acid sequence, germline encoded and largely V gene determined.")
    epitope_aa: Optional[str] = Field(None, description="Amino acid sequence of the epitope this receptor was observed to recognize.")
    antigen: Optional[str] = None
    antigen_organism: Optional[str] = None
    mhc_class: Optional[str] = None
    mhc_a: Optional[str] = None
    mhc_b: Optional[str] = None
    pdb_id: Optional[str] = Field(None, description="RCSB PDB structure id for a solved complex containing this receptor, if any.")
    reference_pmid: Optional[str] = None
    score: Optional[float] = Field(None, description="Source database's own confidence or quality score for this record; scale varies by source.")
    composition: Optional[Composition] = None
    match_kind: Literal["exact", "neighbour"] = Field("exact", description="'exact' means the CDR3 matched the query verbatim; 'neighbour' means it was found by similarity search instead.")
    similarity: Optional[float] = Field(None, description="Within-query relative similarity to the queried CDR3, set only for neighbour matches; only 1.0 is truly identical.")
    concordance: int = Field(1, description="Number of distinct source databases that independently report this exact CDR3.")
    mhc_organism: Optional[str] = Field(None, description="Organism inferred from the mhc_a allele name ('human' for HLA*, 'mouse' for H2*/H-2*), independent of the record's own species field.")
    mhc_is_cross_species: bool = Field(False, description="True when mhc_organism is known and differs from this record's species, e.g. an HLA-transgenic mouse.")


class PairedRecord(BaseModel):
    pairing_key: str
    source: str
    alpha: Optional[TCRRecord] = None
    beta: Optional[TCRRecord] = None


class RecordsRequest(BaseModel):
    query: Optional[str] = Field(None, description="Free query: a CDR3, a gene, a V+J+CDR3 phrase, or a database id.", examples=["CASSLGTEAFF", "TRBV20-1", "vdjdb:c123"])
    cdr3_aa: Optional[str] = Field(None, description="A CDR3 amino acid sequence to retrieve records for.", examples=["CASSLGTEAFF"])
    cdr3_aa_b: Optional[str] = Field(None, description="A second-chain CDR3 for a pair query (alpha+beta).", examples=["CAVRDSNYQLIW"])
    v_gene: Optional[str] = Field(None, examples=["TRBV20-1"])
    j_gene: Optional[str] = Field(None, examples=["TRBJ2-7"])
    species: Optional[SpeciesType] = Field(None, description="Filter by organism; omit for all species.", examples=["human", "mouse"])
    top_k: int = Field(50, ge=1, le=500, description="Max exact records to return.", examples=[50])
    include_neighbours: bool = Field(True, description="Also return BLOSUM near neighbours (kept separate from exact hits).")
    include_cross_species_mhc: bool = Field(False, description="Include records whose MHC allele belongs to a different species than the query (e.g. human HLA on a mouse record); off by default.")


class RecordsResponse(BaseModel):
    query_echo: dict[str, Any]
    exact: list[TCRRecord] = Field(default_factory=list)
    neighbours: list[TCRRecord] = Field(default_factory=list)
    pairs: list[PairedRecord] = Field(default_factory=list)
    total_exact: int = 0
    sources_searched: list[str] = Field(default_factory=list)
    warnings: list[DossierWarning] = Field(default_factory=list)
