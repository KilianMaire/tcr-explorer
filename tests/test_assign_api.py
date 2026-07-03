from fastapi.testclient import TestClient
from imgt_app.api import app

client = TestClient(app)


def test_assign_full_chain_returns_alleles_and_cdr3():
    from imgt_app.reconstructor import reconstruct_tcr
    chain_aa = reconstruct_tcr("TRBV19", "TRBJ1-4", "CASSMADRKFF", "mouse")["full_chain_aa"]
    r = client.post("/v1/tcr/assign", json={"sequence": chain_aa, "species": "mouse"})
    assert r.status_code == 200
    b = r.json()
    assert b["chain"] == "TRB" and b["v_determinable"] is True
    assert any(n.startswith("TRBV19") for n in b["v_call"]["alleles"])
    assert b["cdr3_aa"] == "CASSMADRKFF"
    assert b["reconstruction"]["full_chain_aa"]


def test_assign_bare_cdr3_refuses_v():
    r = client.post("/v1/tcr/assign", json={"sequence": "CASSLGTEAFF", "species": "human"})
    assert r.status_code == 200
    b = r.json()
    assert b["v_determinable"] is False and b["v_call"] is None
    assert b["j_call"]["alleles"]
    assert b["v_db_inference"] is not None
