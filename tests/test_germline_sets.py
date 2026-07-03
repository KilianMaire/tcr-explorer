import pytest
from imgt_app.cdr_enricher import _stitchr_data_dir
from imgt_app.dossier_models import AlignRequest
from imgt_app.germline_sets import load_segment_map, resolve_sequences

pytestmark = pytest.mark.skipif(_stitchr_data_dir() is None, reason="stitchr germline not installed")

def test_mouse_trbj_segment_nonempty():
    m = load_segment_map("TRB", "MOUSE", "J")
    assert len(m) >= 5
    assert all(isinstance(v, str) and v for v in m.values())

def test_germline_set_request_resolves_j():
    seqs, warns = resolve_sequences(AlignRequest(species="mouse", chain="TRB", segment="J"))
    assert len(seqs) >= 5 and all(len(s[1]) > 0 for s in seqs)

def test_d_segment_unavailable_warns():
    seqs, warns = resolve_sequences(AlignRequest(species="mouse", chain="TRB", segment="D"))
    assert seqs == []
    assert any(w.code == "segment_unavailable" for w in warns)

def test_gene_name_list_resolves_and_skips_unknown():
    seqs, warns = resolve_sequences(AlignRequest(genes=["TRBJ1-1", "TRBJ2-7", "TRBJ9-9"]))
    names = [n for n, _ in seqs]
    assert "TRBJ1-1" in names or len(names) >= 1
    assert any(w.code == "ambiguous_gene" for w in warns)

def test_provided_sequences_pass_through():
    req = AlignRequest(sequences=[{"name": "a", "seq": "CASS"}, {"name": "b", "seq": "CASF"}], seq_type="aa")
    seqs, warns = resolve_sequences(req)
    assert seqs == [("a", "CASS"), ("b", "CASF")]
