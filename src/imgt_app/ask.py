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
from .nl_query import heuristic_parse

_SYS = (
    "You route a TCR question. Return JSON with key 'intent' one of "
    "'dossier','similar','search', plus any of: query, v_gene, j_gene, cdr3_aa, gene_name. "
    "Use 'similar' for nearest-neighbour requests, 'dossier' for one TCR/sequence/gene, "
    "'search' for database filters. JSON only."
)


def _run_search_sync(req: SearchRequest):
    from .dossier_epitopes import _run_search  # reuse the sync loop helper

    return _run_search(req)


_DOSSIER_TYPES = ("raw_nt", "raw_aa", "gene_name", "allele", "id")


def _looks_like_cdr3(tok: str) -> bool:
    """A protein token of CDR3-ish length: the thing a user wants looked up."""
    r = route(tok, "auto")
    return r.detected_type == "raw_aa" and 8 <= len(tok) <= 22


def _find_gene(tokens: list[str], seg: str) -> str:
    """Return the first token that looks like a TCR gene of segment seg (V/J)."""
    for t in tokens:
        u = t.upper()
        if len(u) >= 5 and u[:2] == "TR" and u[2] in "ABGD" and u[3:4] == seg:
            return t
    return ""


def _heuristic_plan(query: str, species: str) -> dict:
    q = query.lower()
    toks = query.split()
    cdr3 = next((t for t in toks if _looks_like_cdr3(t)), "")
    vg, jg = _find_gene(toks, "V"), _find_gene(toks, "J")
    forced_similar = any(k in q for k in ("similar", "nearest", "neighbour", "neighbor", "close to"))
    # A CDR3 is present: look it up against the known-TCR reference (nearest known
    # TCRs and their epitopes). This is the core purpose of the tool for a CDR3, so
    # it takes priority over germline annotation (which needs V and J anyway).
    if cdr3 or forced_similar:
        return {"intent": "similar", "cdr3_aa": cdr3, "v_gene": vg, "j_gene": jg}
    # Otherwise a gene name or a raw sequence: annotate it in a dossier. Only accept
    # an unambiguous classification (no router warnings) so English filler words
    # spellable in the amino-acid alphabet do not get misrouted ahead of the real
    # gene/sequence token.
    for tok in toks:
        r = route(tok, "auto")
        if r.detected_type in _DOSSIER_TYPES and not r.warnings:
            return {"intent": "dossier", "query": tok, "v_gene": vg, "j_gene": jg}
    return {"intent": "search", "query": query}


def plan_intent(query: str, species: str):
    if llm_available():
        raw = llm_json(_SYS, query)
        if isinstance(raw, dict) and raw.get("intent") in ("dossier", "similar", "search"):
            return raw, "llm", True
    return _heuristic_plan(query, species), "heuristic", False


def _tokens(query):
    return query.replace(",", " ").split()


def _execute(request: AskRequest, plan: dict, source: str, llm_used: bool) -> AskResponse:
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
    gene_name = plan.get("gene_name")
    if gene_name:
        sr = SearchRequest(gene_name=gene_name, species=request.species)
    else:
        parsed = heuristic_parse(request.query)
        sr = SearchRequest(
            gene_name=parsed.gene_name,
            region=parsed.region,
            sequence_contains=parsed.sequence_contains,
            antigen_epitope=parsed.antigen_epitope,
            source=parsed.source,
            species=request.species,
        )
    result = _run_search_sync(sr)
    return AskResponse(intent="search", plan_source=source, llm_used=llm_used, search_result=result)


def answer(request: AskRequest) -> AskResponse:
    plan, source, llm_used = plan_intent(request.query, request.species)
    if source == "llm":
        try:
            return _execute(request, plan, source, llm_used)
        except Exception:
            plan = _heuristic_plan(request.query, request.species)
            return _execute(request, plan, "heuristic", False)
    return _execute(request, plan, source, llm_used)
