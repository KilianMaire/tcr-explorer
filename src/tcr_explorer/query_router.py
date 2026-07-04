from __future__ import annotations

import dataclasses
from dataclasses import dataclass, field
from typing import Optional

from .input_router import route as classify
from .query_nl import parse_query, _is_cdr3_token
from .records import retrieve_records
from .tcr_align import assign
from .dossier import build_dossier, find_similar_tcrs
from .ask import answer
from .dossier_models import RecordsRequest, DossierRequest, AskRequest, SimilarResponse


@dataclass
class Block:
    tool: str
    title: str
    data: dict


@dataclass
class QueryResult:
    input: str
    detected_type: str
    species: Optional[str]
    tools: list[str]
    note: str
    blocks: list[Block] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)


_FORCE = {"records", "assign", "dossier", "similar"}


def _records(**kw) -> Block:
    return Block("records", "Known database records",
                 retrieve_records(RecordsRequest(**kw)).model_dump())


def _assign(seq: str, species: Optional[str]) -> Block:
    return Block("assign", "Germline assignment",
                 dataclasses.asdict(assign(seq, species=species)))


def _dossier(query: str, species: Optional[str]) -> Block:
    return Block("dossier", "Dossier",
                 build_dossier(DossierRequest(query=query, species=species or "human")).model_dump())


def _ask(query: str, species: Optional[str]) -> Block:
    return Block("ask", "Answer",
                 answer(AskRequest(query=query, species=species or "human")).model_dump())


def _similar(cdr3, v, j, species, warnings) -> Optional[Block]:
    if not (cdr3 and v and j):
        warnings.append("similarity needs a CDR3 with V and J genes")
        return None
    neigh, engine, total, w = find_similar_tcrs(cdr3, v, j, species=species or "human")
    data = SimilarResponse(neighbours=neigh, engine=engine,
                           total_candidates=total, warnings=w).model_dump()
    return Block("similar", "Similar known TCRs", data)


def route_query(query: str, species: Optional[str] = None, force: Optional[str] = None) -> QueryResult:
    q = (query or "").strip()
    parsed = parse_query(q)
    sp = species or parsed["species"]
    dtype = classify(q).detected_type
    blocks: list[Block] = []
    warnings: list[str] = []

    if force in _FORCE:
        if force == "records":
            blocks.append(_records(query=q, species=sp))
        elif force == "assign":
            blocks.append(_assign(q, sp))
        elif force == "dossier":
            blocks.append(_dossier(q, sp))
        elif force == "similar":
            b = _similar(parsed["cdr3_aa"], parsed["v_gene"], parsed["j_gene"], sp, warnings)
            if b:
                blocks.append(b)
        return QueryResult(q, dtype, sp, [force], f"forced {force}", blocks, warnings)

    if dtype == "id":
        blocks.append(_records(query=q, species=sp))
        note = "database id lookup"
        tools = ["records"]
    elif dtype in ("gene_name", "allele"):
        blocks.append(_records(query=q, species=sp))
        note = "gene lookup"
        tools = ["records"]
    elif dtype == "raw_nt":
        blocks.append(_assign(q, sp))
        note = "nucleotide sequence, germline assignment"
        tools = ["assign"]
    elif dtype == "raw_aa" and _is_cdr3_token(q):
        blocks.append(_records(cdr3_aa=q, species=sp))
        blocks.append(_assign(q, sp))
        note = "amino acid CDR3, known records and germline assignment (J only, V needs framework)"
        tools = ["records", "assign"]
    elif dtype == "raw_aa":
        blocks.append(_assign(q, sp))
        note = "amino acid sequence, germline assignment"
        tools = ["assign"]
    elif parsed["cdr3_aa"] or parsed["v_gene"] or parsed["j_gene"] or parsed["record_id"]:
        blocks.append(_records(query=q, species=sp))
        note = "parsed a phrase into records filters"
        tools = ["records"]
    else:
        blocks.append(_ask(q, sp))
        note = "free text, routed to ask"
        tools = ["ask"]

    return QueryResult(q, dtype, sp, tools, note, blocks, warnings)
