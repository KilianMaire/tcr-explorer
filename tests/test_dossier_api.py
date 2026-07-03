from fastapi.testclient import TestClient
from imgt_app.api import app
from imgt_app import dossier_epitopes
from imgt_app.models import GeneRecord, IEDBHit, SearchResponse

client = TestClient(app)

def test_dossier_endpoint_json():
    r = client.post("/v1/tcr/dossier", json={"query": "TRBV20-1", "species": "human"})
    assert r.status_code == 200
    body = r.json()
    assert body["schema_version"] == "1.0"
    assert "summary" in body and "warnings" in body and "provenance" in body
    assert set(body["genes"].keys()) == {"v", "d", "j", "c"}

def test_dossier_endpoint_markdown():
    r = client.post("/v1/tcr/dossier", json={"query": "TRBV20-1"},
                    headers={"Accept": "text/markdown"})
    assert r.status_code == 200
    assert r.headers["content-type"].startswith("text/markdown")


def test_known_epitopes_surface_over_http(monkeypatch):
    # Regression for the async-route bug: build_dossier ran inside the request's
    # event loop, so dossier_epitopes._run_search bailed on the running loop and
    # known_epitopes was ALWAYS empty over HTTP. With a synchronous route FastAPI
    # runs it in a threadpool where _run_search can create its own loop.
    async def _fake_search(req):
        return SearchResponse(
            total=1,
            records=[GeneRecord(
                source="vdjdb", species="human", gene_name="TRBV20-1", sequence="",
                iedb_hits=[IEDBHit(epitope_sequence="NLVPMVATV")],
            )],
        )

    # Patch the search seam used by the real lookup path (default-bound in build_dossier).
    monkeypatch.setattr(dossier_epitopes, "search", _fake_search)

    r = client.post("/v1/tcr/dossier", json={
        "query": "CASSFGTEAFF", "input_type": "raw_aa",
        "v_gene": "TRBV20-1", "j_gene": "TRBJ2-7", "cdr3_aa": "CASSFGTEAFF",
    })
    assert r.status_code == 200
    body = r.json()
    assert body["known_epitopes"], "known epitopes must surface in the HTTP JSON"
    assert body["known_epitopes"][0]["epitope_sequence"] == "NLVPMVATV"
