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


def test_ui_escapes_rendered_data():
    """Guard against a silent XSS regression: the render path must define esc()
    and route data through it before writing to innerHTML."""
    body = client.get("/ui").text
    assert "function esc(" in body
    # every user/data-derived interpolation in the render helpers goes through esc()
    assert "esc(b.intent)" in body
    assert "esc(e.epitope_sequence)" in body
    assert "esc(n.cdr3_b_aa)" in body
    # innerHTML is only assigned from strings built by render()/neighTable()
    assert ".innerHTML=render(" in body or "innerHTML=render(" in body
