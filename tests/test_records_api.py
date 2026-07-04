from fastapi.testclient import TestClient
from tcr_explorer.api import app

client = TestClient(app)


def test_records_route_returns_response_shape():
    r = client.post("/v1/tcr/records", json={"cdr3_aa": "CASSLGTEAFF"})
    assert r.status_code == 200
    body = r.json()
    assert "exact" in body and "neighbours" in body and "sources_searched" in body
    # never mix exact and neighbours
    exact_ids = {x["source_record_id"] for x in body["exact"]}
    neigh_ids = {x["source_record_id"] for x in body["neighbours"]}
    assert exact_ids.isdisjoint(neigh_ids)


def test_records_route_id_lookup():
    r = client.post("/v1/tcr/records", json={"query": "vdjdb:c1"})
    assert r.status_code == 200
