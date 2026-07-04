from fastapi.testclient import TestClient

from tcr_explorer.api import app

client = TestClient(app)


def test_ask_endpoint_routes_gene_to_dossier():
    r = client.post("/v1/tcr/ask", json={"query": "TRBV20-1"})
    assert r.status_code == 200
    body = r.json()
    assert body["intent"] in ("dossier", "search")
    assert "llm_used" in body and "plan_source" in body
