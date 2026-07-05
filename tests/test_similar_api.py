import pytest
from fastapi.testclient import TestClient

from tcr_explorer.api import app
from tcr_explorer import tcrdist_engine

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


@pytest.mark.skipif(not tcrdist_engine.tcrdist_available(),
                    reason="pwseqdist not installed (tcrdist extra)")
def test_similar_paired_endpoint(monkeypatch, tmp_path):
    import pandas as pd
    import tcr_explorer.similarity as similarity

    rows = []
    for pk, va, ca, vb, cb, epi in [
        ("pk1", "TRAV12-1", "CAVNFGGGKLIF", "TRBV19", "CASSIRSSYEQYF", "NLVPMVATV"),
        ("pk2", "TRAV1-2", "CAVRDSNYQLIW", "TRBV28", "CASSLGQAYEQYF", "GILGFVFTL"),
    ]:
        rows.append({"pairing_key": pk, "chain": "alpha", "cdr3_aa": ca, "v_gene": va,
                     "species": "human", "epitope_aa": epi})
        rows.append({"pairing_key": pk, "chain": "beta", "cdr3_aa": cb, "v_gene": vb,
                     "species": "human", "epitope_aa": epi})
    p = tmp_path / "paired.parquet"
    pd.DataFrame(rows).to_parquet(p, index=False)

    monkeypatch.setenv("UNITCR_INDEX_PATH", str(p))
    similarity._load_index.cache_clear()
    r = client.post(
        "/v1/tcr/similar_paired",
        json={"cdr3_a": "CAVNFGGGKLIF", "v_a": "TRAV12-1",
              "cdr3_b": "CASSIRSSYEQYF", "v_b": "TRBV19", "top_k": 2},
    )
    assert r.status_code == 200
    body = r.json()
    assert body["engine"] == "tcrdist"
    assert body["neighbours"][0]["cdr3_a_aa"] == "CAVNFGGGKLIF"
    assert body["neighbours"][0]["distance"] == 0.0
