from imgt_app.dossier_models import AskRequest
from imgt_app import ask


def test_heuristic_gene_routes_to_dossier(monkeypatch):
    monkeypatch.setattr(ask, "llm_available", lambda: False)
    resp = ask.answer(AskRequest(query="TRBV20-1"))
    assert resp.intent == "dossier"
    assert resp.llm_used is False and resp.plan_source == "heuristic"
    assert resp.dossier is not None


def test_heuristic_similar_keyword_routes_to_similar(monkeypatch):
    monkeypatch.setattr(ask, "llm_available", lambda: False)
    monkeypatch.setattr(
        ask, "find_similar_tcrs", lambda *a, **k: ([], "blosum_cdr3", 0, [])
    )
    resp = ask.answer(
        AskRequest(query="find TCRs similar to CASSLGTEAFF TRBV20-1 TRBJ1-1")
    )
    assert resp.intent == "similar"
    assert resp.neighbours_result is not None


def test_llm_plan_used_when_available(monkeypatch):
    monkeypatch.setattr(ask, "llm_available", lambda: True)
    monkeypatch.setattr(
        ask, "llm_json", lambda s, u: {"intent": "search", "gene_name": "TRBV20-1"}
    )
    monkeypatch.setattr(
        ask,
        "_run_search_sync",
        lambda req: __import__(
            "imgt_app.models", fromlist=["SearchResponse"]
        ).SearchResponse(total=0, records=[]),
    )
    resp = ask.answer(AskRequest(query="anything"))
    assert resp.intent == "search" and resp.llm_used is True and resp.plan_source == "llm"
