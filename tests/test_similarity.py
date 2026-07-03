from imgt_app import similarity
from imgt_app.similarity import find_similar_tcrs, cdr3_distance

FIX = "tests/fixtures/unitcr_tiny.parquet"


def test_self_match_is_nearest(monkeypatch):
    monkeypatch.setattr(similarity, "tcrdist3_available", lambda: False)
    neigh, engine, ncand, warns = find_similar_tcrs(
        "CASSLGTEAFF", "TRBV20-1", "TRBJ1-1", top_k=3, index_path=FIX)
    assert engine == "blosum_cdr3"
    assert neigh[0].cdr3_b_aa == "CASSLGTEAFF"
    assert neigh[0].similarity >= 0.99
    assert any(w.code == "tcrdist_unavailable" for w in warns)


def test_distance_zero_for_identical():
    assert cdr3_distance("CASSLGTEAFF", "CASSLGTEAFF") == 0.0


def test_missing_index_is_graceful():
    neigh, engine, ncand, warns = find_similar_tcrs(
        "CASSLGTEAFF", "TRBV20-1", "TRBJ1-1", index_path="/no/such/index.parquet")
    assert neigh == []
    assert any(w.code == "similarity_index_unavailable" for w in warns)
