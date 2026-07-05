from fastapi.testclient import TestClient
from tcr_explorer.api import app
client = TestClient(app)

def test_ui_served():
    r = client.get("/ui")
    assert r.status_code == 200
    assert r.headers["content-type"].startswith("text/html")
    body = r.text
    assert "<form" in body and "/v1/tcr/query" in body
    assert "TCR" in body


def test_ui_has_paired_similarity_form():
    body = client.get("/ui").text
    # the paired chip, its form inputs, the endpoint call, and the render helper.
    assert 'data-tool="paired"' in body
    assert 'id="pairform"' in body
    assert 'id="p_ca"' in body and 'id="p_vb"' in body
    assert "/v1/tcr/similar_paired" in body
    assert "function pairedTable(" in body
    # paired render must escape data-derived cells too.
    assert "esc(n.cdr3_a_aa)" in body


def test_ui_escapes_rendered_data():
    """Guard against a silent XSS regression: the render path must define esc()
    and route data through it before writing to innerHTML."""
    body = client.get("/ui").text
    assert "function esc(" in body
    # every user/data-derived interpolation in the render helpers goes through esc()
    assert "esc(b.intent)" in body
    assert "esc(e.epitope_sequence)" in body
    assert "esc(n.cdr3_b_aa)" in body
    # the adaptive result area is only assigned from strings built by
    # renderBlock() (which itself only calls renderRecords/renderAssign/render(),
    # all esc()-based), never from a raw template interpolation
    assert "out.innerHTML=h" in body
    assert "for(const block of (b.blocks||[]))h+=renderBlock(block)" in body
