from fastapi.testclient import TestClient

from tcr_explorer.api import app

client = TestClient(app)


def test_similar_endpoint_uses_tiny_index(monkeypatch):
    import tcr_explorer.api as api  # noqa: F401
    import tcr_explorer.similarity as similarity

    # force the tiny fixture + BLOSUM fallback via the engine's default; patch index path
    monkeypatch.setenv("UNITCR_INDEX_PATH", "tests/fixtures/unitcr_tiny.parquet")
    similarity._load_index.cache_clear()
    r = client.post(
        "/v1/tcr/similar",
        json={"cdr3": "CASSLGTEAFF", "v_gene": "TRBV20-1", "j_gene": "TRBJ1-1", "top_k": 2},
    )
    assert r.status_code == 200
    body = r.json()
    assert "neighbours" in body and "engine" in body
