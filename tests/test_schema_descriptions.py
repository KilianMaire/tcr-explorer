from imgt_app.dossier_models import DossierRequest, SimilarRequest
from imgt_app.models import SearchRequest


def test_dossier_request_query_is_described():
    s = DossierRequest.model_json_schema()
    q = s["properties"]["query"]
    assert q.get("description") and len(q["description"]) > 20
    assert q.get("examples")


def test_dossier_request_input_type_is_described():
    s = DossierRequest.model_json_schema()
    it = s["properties"]["input_type"]
    assert it.get("description") and len(it["description"]) > 5
    assert it.get("examples")


def test_dossier_request_species_is_described():
    s = DossierRequest.model_json_schema()
    sp = s["properties"]["species"]
    assert sp.get("description") and len(sp["description"]) > 5
    assert sp.get("examples")


def test_dossier_request_mode_is_described():
    s = DossierRequest.model_json_schema()
    m = s["properties"]["mode"]
    assert m.get("description") and len(m["description"]) > 5
    assert m.get("examples")


def test_dossier_request_include_is_described():
    s = DossierRequest.model_json_schema()
    inc = s["properties"]["include"]
    assert inc.get("description") and len(inc["description"]) > 5
    assert inc.get("examples")


def test_dossier_request_v_gene_is_described():
    s = DossierRequest.model_json_schema()
    vg = s["properties"]["v_gene"]
    assert vg.get("description") and len(vg["description"]) > 5
    assert vg.get("examples")


def test_dossier_request_j_gene_is_described():
    s = DossierRequest.model_json_schema()
    jg = s["properties"]["j_gene"]
    assert jg.get("description") and len(jg["description"]) > 5
    assert jg.get("examples")


def test_dossier_request_cdr3_aa_is_described():
    s = DossierRequest.model_json_schema()
    cdr3 = s["properties"]["cdr3_aa"]
    assert cdr3.get("description") and len(cdr3["description"]) > 5
    assert cdr3.get("examples")


def test_similar_request_cdr3_is_described():
    s = SimilarRequest.model_json_schema()
    cdr3 = s["properties"]["cdr3"]
    assert cdr3.get("description") and len(cdr3["description"]) > 5
    assert cdr3.get("examples")


def test_similar_request_v_gene_is_described():
    s = SimilarRequest.model_json_schema()
    vg = s["properties"]["v_gene"]
    assert vg.get("description") and len(vg["description"]) > 5
    assert vg.get("examples")


def test_similar_request_j_gene_is_described():
    s = SimilarRequest.model_json_schema()
    jg = s["properties"]["j_gene"]
    assert jg.get("description") and len(jg["description"]) > 5
    assert jg.get("examples")


def test_similar_request_species_is_described():
    s = SimilarRequest.model_json_schema()
    sp = s["properties"]["species"]
    assert sp.get("description") and len(sp["description"]) > 5
    assert sp.get("examples")


def test_similar_request_top_k_is_described():
    s = SimilarRequest.model_json_schema()
    tk = s["properties"]["top_k"]
    assert tk.get("description") and len(tk["description"]) > 5
    assert tk.get("examples")


def test_similar_request_min_similarity_is_described():
    s = SimilarRequest.model_json_schema()
    ms = s["properties"]["min_similarity"]
    assert ms.get("description") and len(ms["description"]) > 5
    assert ms.get("examples")


def test_search_request_source_is_described():
    s = SearchRequest.model_json_schema()
    src = s["properties"]["source"]
    assert src.get("description") and len(src["description"]) > 5
    assert src.get("examples")


def test_search_request_species_is_described():
    s = SearchRequest.model_json_schema()
    sp = s["properties"]["species"]
    assert sp.get("description") and len(sp["description"]) > 5
    assert sp.get("examples")


def test_search_request_gene_name_is_described():
    s = SearchRequest.model_json_schema()
    gn = s["properties"]["gene_name"]
    assert gn.get("description") and len(gn["description"]) > 5
    assert gn.get("examples")


def test_search_request_region_is_described():
    s = SearchRequest.model_json_schema()
    r = s["properties"]["region"]
    assert r.get("description") and len(r["description"]) > 5
    assert r.get("examples")


def test_search_request_sequence_contains_is_described():
    s = SearchRequest.model_json_schema()
    sc = s["properties"]["sequence_contains"]
    assert sc.get("description") and len(sc["description"]) > 5
    assert sc.get("examples")


def test_search_request_antigen_epitope_is_described():
    s = SearchRequest.model_json_schema()
    ae = s["properties"]["antigen_epitope"]
    assert ae.get("description") and len(ae["description"]) > 5
    assert ae.get("examples")


def test_search_request_limit_is_described():
    s = SearchRequest.model_json_schema()
    lim = s["properties"]["limit"]
    assert lim.get("description") and len(lim["description"]) > 5
    assert lim.get("examples")


def test_search_request_offset_is_described():
    s = SearchRequest.model_json_schema()
    off = s["properties"]["offset"]
    assert off.get("description") and len(off["description"]) > 5
    assert off.get("examples")
