from fastapi.testclient import TestClient
from tcr_explorer.api import app
client = TestClient(app)

def test_align_endpoint_provided_sequences():
    r = client.post("/v1/tcr/align", json={"sequences":[{"name":"a","seq":"CASSLGTEAFF"},{"name":"b","seq":"CASSLGTEAF"}], "seq_type":"aa"})
    assert r.status_code == 200
    b = r.json()
    assert b["n_sequences"] == 2 and b["engine"] in ("center_star","clustalo")
    assert len(b["records"]) == 2

def test_align_endpoint_too_few():
    r = client.post("/v1/tcr/align", json={"sequences":[{"name":"a","seq":"CASS"}], "seq_type":"aa"})
    assert r.status_code == 200
    assert any(w["code"] == "too_few_sequences" for w in r.json()["warnings"])
