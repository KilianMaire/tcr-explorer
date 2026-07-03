from imgt_app.dossier_models import AskRequest
from imgt_app import ask


def test_llm_malformed_dossier_params_fall_back_to_heuristic(monkeypatch):
    monkeypatch.setattr(ask, "llm_available", lambda: True)
    monkeypatch.setattr(
        ask, "llm_json", lambda s, u: {"intent": "dossier", "v_gene": 123}
    )
    resp = ask.answer(AskRequest(query="TRBV20-1"))
    assert resp.plan_source == "heuristic"
    assert resp.llm_used is False
    assert resp.dossier is not None


def test_llm_malformed_similar_params_fall_back_to_heuristic(monkeypatch):
    monkeypatch.setattr(ask, "llm_available", lambda: True)
    monkeypatch.setattr(
        ask,
        "llm_json",
        lambda s, u: {"intent": "similar", "v_gene": 123, "cdr3_aa": "CASSLGTEAFF"},
    )
    resp = ask.answer(
        AskRequest(query="find TCRs similar to CASSLGTEAFF TRBV20-1 TRBJ1-1")
    )
    assert resp.plan_source == "heuristic"
    assert resp.llm_used is False
    assert resp.neighbours_result is not None


def test_bare_cdr3_routes_to_similarity_lookup(monkeypatch):
    # A bare CDR3 (no genes) is the tool's core lookup case: find known TCRs (and
    # their epitopes) matching it, rather than a dead-end germline annotation.
    monkeypatch.setattr(ask, "llm_available", lambda: False)
    captured = {}

    def fake_similar(cdr3, v, j, species="human", **k):
        captured.update(cdr3=cdr3)
        return ([], "blosum_cdr3", 0, [])

    monkeypatch.setattr(ask, "find_similar_tcrs", fake_similar)
    resp = ask.answer(AskRequest(query="find this CDR3 CASSLGTEAFF"))
    assert resp.intent == "similar"
    assert resp.neighbours_result is not None
    assert captured["cdr3"] == "CASSLGTEAFF"


def test_cdr3_with_genes_routes_to_dossier(monkeypatch):
    # A CDR3 supplied with V (and/or J) genes is characterized in a dossier
    # (annotation / reconstruction), which is species-agnostic (covers mouse).
    monkeypatch.setattr(ask, "llm_available", lambda: False)
    resp = ask.answer(AskRequest(
        query="CASSLGTEAFF TRBV20-1 TRBJ2-7", species="human"))
    assert resp.intent == "dossier"
    assert resp.dossier is not None
    assert resp.dossier.genes["v"] and resp.dossier.genes["v"].call == "TRBV20-1"


def test_heuristic_scans_all_tokens_for_gene_question(monkeypatch):
    monkeypatch.setattr(ask, "llm_available", lambda: False)
    resp = ask.answer(AskRequest(query="what does TRBV20-1 recognise?"))
    assert resp.intent == "dossier"
    assert resp.dossier is not None
    assert resp.dossier.chain == "beta"


def test_heuristic_search_uses_query_text(monkeypatch):
    monkeypatch.setattr(ask, "llm_available", lambda: False)
    captured = {}

    def fake_run_search(req):
        captured["req"] = req
        from imgt_app.models import SearchResponse

        return SearchResponse(total=0, records=[])

    monkeypatch.setattr(ask, "_run_search_sync", fake_run_search)
    resp = ask.answer(AskRequest(query="database records for HLA-A"))
    assert resp.intent == "search"
    assert resp.search_result is not None
    req = captured["req"]
    assert req is not None
    # heuristic_parse must have derived something from the query text itself
    # (species is always populated from the request regardless, so exclude it).
    assert any(
        getattr(req, f) is not None
        for f in ("gene_name", "sequence_contains", "antigen_epitope", "source")
    )


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
