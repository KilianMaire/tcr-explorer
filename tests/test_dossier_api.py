from fastapi.testclient import TestClient
from imgt_app.api import app

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
