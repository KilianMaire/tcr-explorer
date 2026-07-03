from imgt_app.dossier_models import Neighbour, SimilarRequest, DossierWarning

def test_neighbour_minimal():
    n = Neighbour(cdr3_b_aa="CASSLGTEAFF", v_b_gene="TRBV20-1", j_b_gene="TRBJ1-1",
                  similarity=1.0, distance=0.0)
    assert n.epitope_aa is None

def test_similar_request_defaults():
    r = SimilarRequest(cdr3="CASSLGTEAFF", v_gene="TRBV20-1", j_gene="TRBJ1-1")
    assert r.species == "human" and r.top_k == 10 and r.min_similarity == 0.0

def test_new_warning_codes_valid():
    for c in ("tcrdist_unavailable", "similarity_index_unavailable", "no_reference_candidates"):
        assert DossierWarning(code=c, message="m").code == c
