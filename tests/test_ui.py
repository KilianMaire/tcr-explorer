from fastapi.testclient import TestClient
from imgt_app.api import app
client = TestClient(app)

def test_ui_served():
    r = client.get("/ui")
    assert r.status_code == 200
    assert r.headers["content-type"].startswith("text/html")
    body = r.text
    assert "<form" in body and "/v1/tcr/ask" in body
    assert "TCR" in body
