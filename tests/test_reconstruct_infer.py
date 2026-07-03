"""Reconstruct a full chain from a CDR3 alone by inferring V and J.

V/J are inferred by tallying deposited records that carry the exact same CDR3
(most common pairing wins). The result is labeled inferred, with the supporting
record count and the alternative pairings, so it is never confused with an
explicit germline assignment.
"""
import pandas as pd
import pytest
from fastapi.testclient import TestClient

from imgt_app import records as R
from imgt_app.api import app
from imgt_app.records_build import SCHEMA_COLUMNS

client = TestClient(app)


@pytest.fixture
def idx(tmp_path):
    rows = [
        dict(source="vdjdb", source_record_id="v1", pairing_key="p1", species="human",
             chain="beta", cdr3_aa="CASSLGTEAFF", v_gene="TRBV19", j_gene="TRBJ2-7"),
        dict(source="mcpas", source_record_id="m2", pairing_key="p2", species="human",
             chain="beta", cdr3_aa="CASSLGTEAFF", v_gene="TRBV19", j_gene="TRBJ2-7"),
        dict(source="vdjdb", source_record_id="v3", pairing_key="p3", species="human",
             chain="beta", cdr3_aa="CASSLGTEAFF", v_gene="TRBV5-1*01", j_gene="TRBJ2-7"),
        dict(source="vdjdb", source_record_id="v4", pairing_key="p4", species="mouse",
             chain="beta", cdr3_aa="CASSLGTEAFF", v_gene="TRBV13-1", j_gene="TRBJ1-1"),
    ]
    df = pd.DataFrame(rows)
    for c in SCHEMA_COLUMNS:
        if c not in df.columns:
            df[c] = None
    df = df[SCHEMA_COLUMNS]
    p = tmp_path / "idx.parquet"
    df.to_parquet(p, index=False)
    return str(p)


def test_infer_ranks_pairings_by_support(idx):
    out = R.infer_vj_from_cdr3("CASSLGTEAFF", "human", index_path=idx)
    assert out[0]["v_gene"] == "TRBV19" and out[0]["j_gene"] == "TRBJ2-7"
    assert out[0]["count"] == 2 and out[0]["chain"] == "beta"
    # the allele-suffixed second pairing is normalized to its base
    assert any(c["v_gene"] == "TRBV5-1" for c in out)


def test_infer_respects_species(idx):
    out = R.infer_vj_from_cdr3("CASSLGTEAFF", "mouse", index_path=idx)
    assert len(out) == 1 and out[0]["v_gene"] == "TRBV13-1"


def test_infer_empty_for_unknown_cdr3(idx):
    assert R.infer_vj_from_cdr3("CAAAAAAAAAF", "human", index_path=idx) == []
    assert R.infer_vj_from_cdr3("", "human", index_path=idx) == []


def test_reconstruct_endpoint_infers_from_cdr3_alone():
    # CASSLGTEAFF is in the vendored index; top human pairing is TRBV4-1/TRBJ1-1
    r = client.post("/reconstruct", json={"cdr3_aa": "CASSLGTEAFF", "species": "human"})
    assert r.status_code == 200
    b = r.json()
    assert b["genes_inferred"] is True
    assert b["v_gene"].startswith("TRBV4-1") and b["j_gene"].startswith("TRBJ1-1")
    assert b["v_found"] and b["j_found"] and b["full_chain_aa"]
    assert "CASSLGTEAFF" in b["full_aa"]
    assert b["inference_support"] >= 1
    assert b["inference_alternatives"] and "n=" in b["inference_alternatives"][0]


def test_reconstruct_endpoint_reports_when_cdr3_not_inferable():
    r = client.post("/reconstruct", json={"cdr3_aa": "CAAAWWWKKKF", "species": "human"})
    assert r.status_code == 200
    b = r.json()
    assert b["genes_inferred"] is False
    assert b["v_found"] is False and b["j_found"] is False
    assert "infer" in b["note"].lower()


def test_reconstruct_endpoint_explicit_vj_unchanged():
    r = client.post("/reconstruct", json={
        "v_gene": "TRBV19", "j_gene": "TRBJ1-4", "cdr3_aa": "CASSMADRKFF", "species": "mouse"})
    assert r.status_code == 200
    b = r.json()
    assert b["genes_inferred"] is False
    assert b["full_chain_aa"] and b["full_chain_aa"].startswith(b["full_aa"])
