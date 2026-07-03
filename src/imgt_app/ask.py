"""Free-text `/v1/tcr/ask` shim: routes a natural-language TCR question to the
existing dossier/similar/search machinery, via an optional LLM plan with a
deterministic heuristic fallback.

Honesty: this module never generates prose. `AskResponse` always carries
`intent`, `plan_source` ("llm" | "heuristic"), and `llm_used` so a caller can
tell exactly how the routing decision was made; malformed or unusable LLM
output falls back to the heuristic rather than silently claiming "llm".
"""
from __future__ import annotations

from .dossier_models import AskRequest, AskResponse, DossierRequest
from .llm_client import llm_available, llm_json
from .input_router import route
from .dossier import build_dossier, find_similar_tcrs
from .models import SearchRequest

_SYS = (
    "You route a TCR question. Return JSON with key 'intent' one of "
    "'dossier','similar','search', plus any of: query, v_gene, j_gene, cdr3_aa, gene_name. "
    "Use 'similar' for nearest-neighbour requests, 'dossier' for one TCR/sequence/gene, "
    "'search' for database filters. JSON only."
)


def _run_search_sync(req: SearchRequest):
    from .dossier_epitopes import _run_search  # reuse the sync loop helper

    return _run_search(req)


def _heuristic_plan(query: str, species: str) -> dict:
    q = query.lower()
    if any(k in q for k in ("similar", "nearest", "neighbour", "neighbor", "close to")):
        return {"intent": "similar", "query": query}
    r = route(query.split()[0] if query.split() else query, "auto")
    if r.detected_type in ("raw_nt", "raw_aa", "gene_name", "allele", "id"):
        return {"intent": "dossier", "query": query.split()[0]}
    return {"intent": "search", "query": query}


def plan_intent(query: str, species: str):
    if llm_available():
        raw = llm_json(_SYS, query)
        if isinstance(raw, dict) and raw.get("intent") in ("dossier", "similar", "search"):
            return raw, "llm", True
    return _heuristic_plan(query, species), "heuristic", False


def _tokens(query):
    return query.replace(",", " ").split()


def answer(request: AskRequest) -> AskResponse:
    plan, source, llm_used = plan_intent(request.query, request.species)
    intent = plan.get("intent", "search")
    if intent == "dossier":
        q = plan.get("query") or plan.get("gene_name") or request.query
        d = build_dossier(
            DossierRequest(
                query=q,
                species=request.species,
                v_gene=plan.get("v_gene"),
                j_gene=plan.get("j_gene"),
                cdr3_aa=plan.get("cdr3_aa"),
            )
        )
        return AskResponse(
            intent="dossier", plan_source=source, llm_used=llm_used, dossier=d, warnings=d.warnings
        )
    if intent == "similar":
        toks = _tokens(request.query)
        cdr3 = plan.get("cdr3_aa") or next((t for t in toks if t.upper().startswith("CAS")), "")
        vg = plan.get("v_gene") or next((t for t in toks if t.upper().startswith("TRBV")), "")
        jg = plan.get("j_gene") or next((t for t in toks if t.upper().startswith("TRBJ")), "")
        neigh, engine, total, w = find_similar_tcrs(cdr3, vg, jg, species=request.species)
        from .dossier_models import SimilarResponse

        return AskResponse(
            intent="similar",
            plan_source=source,
            llm_used=llm_used,
            neighbours_result=SimilarResponse(
                neighbours=neigh, engine=engine, total_candidates=total, warnings=w
            ),
        )
    # search
    sr = SearchRequest(gene_name=plan.get("gene_name"), species=request.species)
    result = _run_search_sync(sr)
    return AskResponse(intent="search", plan_source=source, llm_used=llm_used, search_result=result)
