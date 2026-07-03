from imgt_app.dossier_models import DossierRequest
from imgt_app.dossier import build_dossier


def test_neighbours_populated_and_separated(monkeypatch):
    import imgt_app.dossier as D
    from imgt_app.dossier_models import Neighbour

    def fake_similar(cdr3, v, j, species="human", top_k=10, min_similarity=0.0, index_path=None):
        return ([Neighbour(cdr3_b_aa="CASSLGTEAFF", v_b_gene="TRBV20-1", j_b_gene="TRBJ1-1",
                           similarity=0.9, distance=1.0, epitope_aa="NLVPMVATV")], "blosum_cdr3", 1, [])
    monkeypatch.setattr(D, "find_similar_tcrs", fake_similar)
    req = DossierRequest(query="CASSLGTEAFF", input_type="raw_aa",
                         v_gene="TRBV20-1", j_gene="TRBJ1-1", cdr3_aa="CASSLGTEAFF",
                         include=["neighbours"])
    d = build_dossier(req)
    assert d.neighbours and d.neighbours[0].epitope_aa == "NLVPMVATV"
    # honesty: the neighbour epitope must NOT be in known_epitopes
    assert all(h.epitope_sequence != "NLVPMVATV" for h in d.known_epitopes)
    assert any(p.block == "neighbours" and p.kind == "neighbor_inferred" for p in d.provenance)


def test_neighbours_absent_when_not_requested():
    d = build_dossier(DossierRequest(query="TRBV20-1"))
    assert d.neighbours is None
