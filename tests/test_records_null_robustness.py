"""Regression tests for null-robustness against the real vendored index.

The real data/records_index.parquet carries NaN in v_gene, j_gene, and
species for some rows (3,828 null v_gene, 21,303 null j_gene, 6 null
species). The test fixtures used elsewhere never exercised those nulls, so
they slipped through and 500'd the live endpoints. These tests build a tiny
index with the same null shapes and assert the retrieval/similarity/API
paths tolerate them.
"""
import pandas as pd
import pytest
from fastapi.testclient import TestClient

from imgt_app import records as R
from imgt_app import similarity as S
from imgt_app.records_build import SCHEMA_COLUMNS
from imgt_app.dossier_models import RecordsRequest


@pytest.fixture
def null_idx(tmp_path):
    rows = [
        # null v_gene AND null j_gene: must not crash Neighbour(v_b_gene=..., j_b_gene=...)
        dict(source="vdjdb", source_record_id="vdjdb:1", pairing_key="vdjdb:u1", species="human",
             chain="beta", cdr3_aa="CASSLGTEAFF", v_gene=None, j_gene=None,
             epitope_aa="NLVPMVATV", external_url="u"),
        # a near neighbour with real genes, so the null-gene row's similarity search has company
        dict(source="vdjdb", source_record_id="vdjdb:2", pairing_key="vdjdb:u2", species="human",
             chain="beta", cdr3_aa="CASSLGTEAYF", v_gene="TRBV19", j_gene="TRBJ2-7",
             epitope_aa="GILGFVFTL", external_url="u"),
        # null species: must not crash build_record(species=row["species"])
        dict(source="iedb", source_record_id="iedb:3", pairing_key="iedb:3", species=None,
             chain="beta", cdr3_aa="ASGDTGGYEQY", v_gene="TRBV7-2", j_gene="TRBJ2-1",
             epitope_aa=None, external_url="u"),
    ]
    df = pd.DataFrame(rows)
    for c in SCHEMA_COLUMNS:
        if c not in df.columns:
            df[c] = None
    df = df[SCHEMA_COLUMNS]
    p = tmp_path / "null_idx.parquet"
    df.to_parquet(p, index=False)
    return str(p)


def test_find_similar_tcrs_tolerates_null_genes(null_idx):
    neigh, engine, total, warnings = S.find_similar_tcrs(
        "CASSLGTEAFF", "TRBV19", "TRBJ2-7", species="human", index_path=null_idx
    )
    assert total > 0
    assert neigh  # did not raise, returned candidates
    null_gene_hits = [n for n in neigh if n.cdr3_b_aa == "CASSLGTEAFF"]
    assert null_gene_hits
    assert null_gene_hits[0].v_b_gene is None
    assert null_gene_hits[0].j_b_gene is None


def test_retrieve_records_tolerates_null_species(null_idx):
    resp = R.retrieve_records(RecordsRequest(cdr3_aa="ASGDTGGYEQY"), index_path=null_idx)
    assert resp.total_exact == 1
    assert resp.exact[0].species == ""


def test_retrieve_records_id_lookup_tolerates_null_species(null_idx):
    resp = R.retrieve_records(RecordsRequest(query="iedb:3"), index_path=null_idx)
    assert resp.total_exact == 1
    assert resp.exact[0].species == ""


def test_api_records_endpoint_returns_200_over_null_species(null_idx, monkeypatch):
    monkeypatch.setenv("RECORDS_INDEX_PATH", null_idx)
    from imgt_app.api import app

    client = TestClient(app)
    resp = client.post("/v1/tcr/records", json={"cdr3_aa": "ASGDTGGYEQY"})
    assert resp.status_code == 200
    body = resp.json()
    assert body["total_exact"] == 1
    assert body["exact"][0]["species"] == ""


def test_api_similar_endpoint_returns_200_over_null_genes(null_idx, monkeypatch):
    monkeypatch.setenv("UNITCR_INDEX_PATH", null_idx)
    from imgt_app.api import app

    client = TestClient(app)
    resp = client.post(
        "/v1/tcr/similar",
        json={"cdr3": "CASSLGTEAFF", "v_gene": "TRBV19", "j_gene": "TRBJ2-7", "species": "human"},
    )
    assert resp.status_code == 200
    body = resp.json()
    assert body["neighbours"]
