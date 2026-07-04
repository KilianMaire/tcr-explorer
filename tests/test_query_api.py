from fastapi.testclient import TestClient
from tcr_explorer.api import app

client = TestClient(app)


def test_query_cdr3_returns_two_blocks():
    r = client.post("/v1/tcr/query", json={"query": "CASSLGTEAFF", "species": "human"})
    assert r.status_code == 200
    b = r.json()
    assert b["understood"]["tools"] == ["records", "assign"]
    assert [x["tool"] for x in b["blocks"]] == ["records", "assign"]


def test_query_gene_returns_records():
    r = client.post("/v1/tcr/query", json={"query": "TRBV20-1"})
    assert r.status_code == 200
    assert r.json()["understood"]["tools"] == ["records"]


def test_query_force_single_tool():
    r = client.post("/v1/tcr/query", json={"query": "CASSLGTEAFF", "force": "assign"})
    assert r.status_code == 200
    body = r.json()
    assert body["understood"]["tools"] == ["assign"] and len(body["blocks"]) == 1
