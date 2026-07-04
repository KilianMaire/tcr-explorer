"""MCP server (stdio) exposing the TCR dossier, similarity, and ask tools.

Each tool is a thin wrapper over the same pure function the REST API uses, so
there is no logic divergence between REST and MCP. Run with:
    python -m tcr_explorer.mcp_server
"""
from __future__ import annotations
import dataclasses
from typing import Optional
from mcp.server.fastmcp import FastMCP
from mcp.types import ToolAnnotations

from .dossier_models import DossierRequest, AskRequest, RecordsRequest, SimilarResponse
from .dossier_models import AlignRequest as _AlignRequest
from .dossier import build_dossier, find_similar_tcrs as find_similar_tcrs_fn
from .ask import answer as answer_fn
from .records import retrieve_records as retrieve_records_fn
from .tcr_align import assign as assign_fn

mcp = FastMCP("tcr-explorer")

# All tools are read-only lookups: they mutate nothing, are safe to retry, and
# reach external sources (NCBI/IEDB) for the dossier, so openWorldHint is true.
_READ_ONLY = ToolAnnotations(readOnlyHint=True, destructiveHint=False,
                             idempotentHint=True, openWorldHint=True)

@mcp.tool(annotations=_READ_ONLY)
def get_tcr_dossier(query: str, species: str = "human", input_type: str = "auto",
                    mode: str = "fast", include: Optional[list[str]] = None,
                    v_gene: Optional[str] = None, j_gene: Optional[str] = None,
                    cdr3_aa: Optional[str] = None) -> dict:
    """Return a consolidated TCR dossier for a query (a raw nt/aa sequence, an IMGT
    gene name or allele, a namespaced id vdjdb:/iedb:, or a V+J+CDR3 triple).
    Includes chain, V/D/J/C annotation, CDR1/2/3, sequences (synthetic ones are
    tagged), and KNOWN retrieved epitopes only (never predicted), with provenance
    and warnings. Set mode='full' and include=['neighbours'] to add nearest known
    TCRs. Species is 'human' or 'mouse'."""
    req = DossierRequest(query=query, species=species, input_type=input_type,
                         mode=mode, include=include or [], v_gene=v_gene,
                         j_gene=j_gene, cdr3_aa=cdr3_aa)
    return build_dossier(req).model_dump()

@mcp.tool(annotations=_READ_ONLY)
def find_similar_tcrs(cdr3: str, v_gene: str, j_gene: str, species: str = "human",
                      top_k: int = 10, min_similarity: float = 0.0) -> dict:
    """Return the nearest KNOWN TCRs to a query CDR3 (with V and J) from the unitcr
    reference set, with their epitopes, MHC, and antigens. These are an INFERRED,
    weaker signal (epitopes of similar TCRs), NOT confirmed specificity of the
    query. Each neighbour carries similarity and distance. Human beta only."""
    neigh, engine, total, warnings = find_similar_tcrs_fn(
        cdr3, v_gene, j_gene, species=species, top_k=top_k, min_similarity=min_similarity)
    return SimilarResponse(neighbours=neigh, engine=engine,
                           total_candidates=total, warnings=warnings).model_dump()

@mcp.tool(annotations=_READ_ONLY)
def ask_tcr(query: str, species: str = "human") -> dict:
    """Answer a free-text TCR question by routing it to the dossier, similarity, or
    search capability and returning the structured result. Prefer the specific
    tools when you already know the intent; use this for free text. No prose is
    generated; the response carries intent, plan_source, and llm_used."""
    return answer_fn(AskRequest(query=query, species=species)).model_dump()

@mcp.tool(annotations=_READ_ONLY)
def align_tcr_genes(species: str = "human", chain: Optional[str] = None,
                    segment: Optional[str] = None, genes: Optional[list[str]] = None,
                    sequences: Optional[list[dict]] = None, seq_type: str = "nt",
                    translate: bool = False) -> dict:
    """Multiple-sequence-align a TCR gene set: a germline set (species+chain+segment,
    e.g. mouse TRB J), an explicit gene list, or provided sequences. Returns the
    alignment, consensus, and mean pairwise identity. V/J/C are available from the
    germline source; D is not (supply D sequences directly)."""
    from .msa import align
    return align(_AlignRequest(species=species, chain=chain, segment=segment, genes=genes,
                               sequences=sequences, seq_type=seq_type, translate=translate)).model_dump()

@mcp.tool(annotations=_READ_ONLY)
def retrieve_tcr_records(query: str = "", cdr3: str = "", cdr3_b: str = "",
                         v_gene: str = "", j_gene: str = "", species: str = "",
                         top_k: int = 50) -> dict:
    """Retrieve raw federated TCR records (VDJdb/IEDB/McPAS/TCR3d) for a query:
    a CDR3, a V/J gene pair, an alpha+beta pair (cdr3 + cdr3_b), or a namespaced
    database id (vdjdb:/iedb:/mcpas:/tcr3d:) passed as query. Returns exact hits
    and BLOSUM near neighbours kept strictly separate, plus any alpha/beta pairs
    found by shared pairing_key. Species is 'human' or 'mouse'; leave blank for
    all species."""
    req = RecordsRequest(
        query=query or None,
        cdr3_aa=cdr3 or None,
        cdr3_aa_b=cdr3_b or None,
        v_gene=v_gene or None,
        j_gene=j_gene or None,
        species=species or None,
        top_k=top_k,
    )
    return retrieve_records_fn(req).model_dump()

@mcp.tool(annotations=_READ_ONLY)
def assign_tcr_alleles(sequence: str, species: str = "", chain: str = "", want_d: bool = False) -> dict:
    """Identify the composing IMGT germline V and J alleles (and D, constant when
    present) for any TCR sequence, extract the CDR3, and reconstruct the canonical
    chain. A bare CDR3 yields a J allele only."""
    result = assign_fn(sequence, species=species or None, chain=chain or None, want_d=want_d)
    return dataclasses.asdict(result)


def main() -> None:
    mcp.run()

if __name__ == "__main__":
    main()
